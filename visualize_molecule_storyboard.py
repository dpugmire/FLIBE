"""
Blender Python script: Molecular visualizer + storyboard animator
-----------------------------------------------------------------
Run inside Blender via the Scripting workspace:
  1. Open Blender → Scripting tab
  2. Open this file (or paste it in)
  3. Adjust the USER-CONFIGURABLE SETTINGS block below
  4. Click "Run Script"

Frame mapping
  Blender frame 1  =  DCD frame DCD_FRAME_OFFSET (default 700).
  The Blender timeline spans DCD frames DCD_FRAME_OFFSET through SEG5_DCD_END.

Storyboard (boundaries given as DCD frame numbers):
  Seg 1  DCD  700– 900 : all atoms, animated
  Seg 2  DCD  900–1000 : start zoom/pan toward H + 2F focus group
  Seg 3  DCD 1000–1150 : close-up on H + 2F
  Seg 4  DCD 1150–1350 : continue close-up
  Seg 5  DCD 1350–1550 : animate user-defined FOCUS_INDICES list
"""

import bpy
import os
import math
import struct
import numpy as np
from mathutils import Vector, Matrix


# ┌─────────────────────────────────────────────────────────────────────────────┐
# │                        USER-CONFIGURABLE SETTINGS                          │
# └─────────────────────────────────────────────────────────────────────────────┘

# ── File paths ────────────────────────────────────────────────────────────────
# Leave as "" to auto-detect files next to the saved .blend file.
PDB_FILE  = ""   # e.g. "/home/user/sim/structure.pdb"
TRAJ_FILE = ""   # e.g. "/home/user/sim/traj_nvt.dcd"

# ── Trajectory sampling ───────────────────────────────────────────────────────
# Use every Nth DCD frame.
#   FRAME_STEP = 1   → all frames          (smooth but slow to scrub)
#   FRAME_STEP = 10  → every 10th frame    (good balance)
#   FRAME_STEP = 40  → fast preview
FRAME_STEP = 1

# Interpolation mode for atom positions between sampled DCD frames:
#   "nearest" → step/hold (no interpolation; legacy behavior)
#   "linear"  → linear interpolation between sampled frames
#   "spline"  → Catmull-Rom cubic spline interpolation
INTERP_MODE = "linear"

# ── DCD ↔ Blender frame mapping ───────────────────────────────────────────────
# Blender frame 1 corresponds to this DCD frame number.
DCD_FRAME_OFFSET = 700

# ── Storyboard segment boundaries (DCD frame numbers) ────────────────────────
SEG1_DCD_START = 700    # Seg 1 begins: show all atoms
SEG2_DCD_START = 800    # Seg 2 begins: zoom → tritium + F_1 + F_2
SEG3_DCD_START = 1000   # Seg 3 begins: animate tritium + F_1 + F_2
SEG4_DCD_START = 1150   # Seg 4 begins: continue tritium + F_1 + F_2
SEG5_DCD_START = 1350   # Seg 5 begins: animate FOCUS_INDICES
SEG5_DCD_END   = 1550   # end of animation

# ── Zoom / visibility timing controls (DCD frame numbers) ─────────────────────
# Start camera zoom/pan at this frame (all atoms can still be visible here).
ZOOM_START_DCD = SEG2_DCD_START
# Hide non-focus atoms at this frame and keep only H + 2F atoms visible.
HIDE_OTHERS_DCD_START = SEG3_DCD_START

# ── Key atom indices (0-based, matching PDB ATOM/HETATM record order) ─────────
# NOTE: PDB HETATM serial numbers are 1-based; subtract 1 to get 0-based index.
#   e.g. PDB serial 1681 (TR1/Tritium)  → index 1680
#        PDB serial 1682 (F5/Fluoride)  → index 1681
#        PDB serial 86   (F/Fluoride)   → index 85
TRITIUM_IDX    = 1680
FLUORIDE_1_IDX = 86
FLUORIDE_2_IDX = 856

# ── Zoom-focus atom list (H + 2F by default) ──────────────────────────────────
# The script guarantees TRITIUM_IDX, FLUORIDE_1_IDX, FLUORIDE_2_IDX are included.
ZOOM_FOCUS_INDICES = [TRITIUM_IDX, FLUORIDE_1_IDX, FLUORIDE_2_IDX]

# ── Segment-5 focus atom list ─────────────────────────────────────────────────
# Replace these with your actual indices (0-based). Any number of indices is fine.
FOCUS_INDICES = [
    1259, 552, 1006, 1273, 89, 1434, 1159, 592, 1551, 588, 589, 590,
    591, 856, 1155, 1156, 1157, 1158, 1547, 1548, 1549, 1550, 1680,
]

# ── Camera rig ────────────────────────────────────────────────────────────────
# The camera is parented to the AtomCentre empty and sits at
# normalize(initial_camera_offset) * <distance> in the empty's local space.
# Changing a distance value zooms in or out for that segment.
# The camera always points back at AtomCentre via a TrackTo constraint.

# Optional explicit initial camera position at segment start (Blender units).
# If set, this XYZ defines the viewing direction (camera ray) from AtomCentre.
# Segment-1/2/3/4/5 distances are still controlled by CAM_DIST_SEG1..SEG5.
# In other words, CAM_INITIAL_XYZ chooses "where from", CAM_DIST_* chooses
# "how far".
# CAM_INITIAL_SPACE controls whether coordinates are read in AtomCentre local
# space ("LOCAL") or Blender world space ("WORLD").
#
# If left as None, a legacy fallback direction is used.
CAM_INITIAL_XYZ   = (2, -10, 4)       # e.g. (12.0, -10.0, 6.5)
CAM_INITIAL_SPACE = "WORLD"           # "WORLD" or "LOCAL"

# Distance (Blender units; 1 BU = 10 Å) from AtomCentre for each segment.
# Larger = further away / more zoomed out.
CAM_DIST_SEG1 = 15.0   # Seg 1: wide shot — all atoms
CAM_DIST_SEG2 = 15.0   # At ZOOM_START_DCD: start zoom/pan toward H + 2F focus
CAM_DIST_SEG3 =  2.0   # Seg 3: close-up — H + two fluorides
CAM_DIST_SEG4 =  2.0   # Seg 4: same close-up — H + two fluorides
CAM_DIST_SEG5 =  4.0 #8.0   # Seg 5: backed off — focus cluster

# ── Display ───────────────────────────────────────────────────────────────────
SCALE = 0.1          # Å → Blender units  (1 Å = 0.1 BU)

SPHERE_SEGMENTS = 32   # longitude divisions (higher = smoother spheres)
SPHERE_RINGS    = 16   # latitude divisions

# ── Atom radii (Å) — van der Waals defaults ───────────────────────────────────
ATOM_RADII = {
    "F"  : 1.47,   "Li" : 1.82,   "Be" : 1.53,   "H"  : 1.20,
    "B"  : 1.92,   "C"  : 1.70,   "N"  : 1.55,   "O"  : 1.52,
    "Na" : 2.27,   "Mg" : 1.73,   "Al" : 1.84,   "Si" : 2.10,
    "P"  : 1.80,   "S"  : 1.80,   "Cl" : 1.75,   "K"  : 2.75,
    "Ca" : 2.31,
}
DEFAULT_RADIUS = 1.50

ATOM_RADII['F'] *= 0.5
ATOM_RADII['H'] *= 0.25
ATOM_RADII['Be'] *= 0.5
ATOM_RADII['Li'] *= 0.5

# ── CPK colours (RGBA, linear sRGB, values 0.0–1.0) ──────────────────────────
ATOM_COLORS = {
    "F"  : (0.565, 0.878, 0.314, 1.0),   # green
    "Li" : (0.784, 0.502, 1.000, 1.0),   # violet
    "Be" : (0.765, 1.000, 0.000, 1.0),   # yellow-green
    "H"  : (1.000, 1.000, 1.000, 1.0),   # white
    "B"  : (1.000, 0.671, 0.475, 1.0),
    "C"  : (0.565, 0.565, 0.565, 1.0),
    "N"  : (0.188, 0.314, 0.973, 1.0),
    "O"  : (1.000, 0.051, 0.051, 1.0),
    "Na" : (0.671, 0.361, 0.949, 1.0),
    "Mg" : (0.541, 1.000, 0.000, 1.0),
    "Al" : (0.749, 0.651, 0.651, 1.0),
    "Si" : (0.941, 0.784, 0.627, 1.0),
    "P"  : (1.000, 0.502, 0.000, 1.0),
    "S"  : (1.000, 1.000, 0.188, 1.0),
    "Cl" : (0.122, 0.941, 0.122, 1.0),
    "K"  : (0.561, 0.251, 0.831, 1.0),
    "Ca" : (0.239, 1.000, 0.000, 1.0),
}
DEFAULT_COLOR = (0.800, 0.800, 0.800, 1.0)

# Optional custom colors for the two highlighted fluorides.
# Set to None to fall back to the default fluorine color from ATOM_COLORS["F"].
FLUORIDE_1_COLOR = (0.565, 0.878, 0.314, 1.0)
FLUORIDE_2_COLOR = (0.365, 0.820, 0.560, 1.0)


# ┌─────────────────────────────────────────────────────────────────────────────┐
# │                             INTERNAL FUNCTIONS                             │
# └─────────────────────────────────────────────────────────────────────────────┘

# ── Frame mapping ─────────────────────────────────────────────────────────────

def dcd_to_bl(dcd_frame):
    """Convert a DCD frame number to a Blender frame number."""
    return dcd_frame - DCD_FRAME_OFFSET + 1


def bl_to_sample_pos(bl_frame, n_traj_frames):
    """Convert a Blender frame number to a fractional coords-array index (clamped)."""
    dcd_frame = bl_frame + DCD_FRAME_OFFSET - 1
    idx_f = dcd_frame / FRAME_STEP
    return min(max(idx_f, 0.0), float(n_traj_frames - 1))


def get_interpolated_coords(coords, bl_frame):
    """
    Return interpolated atom coordinates for a Blender frame.

    coords shape: (n_sampled, natoms, 3), sampled every FRAME_STEP DCD frames.
    """
    n_traj = coords.shape[0]
    s_pos  = bl_to_sample_pos(bl_frame, n_traj)
    i0     = int(math.floor(s_pos))

    if i0 >= n_traj - 1:
        return coords[n_traj - 1]

    t = s_pos - i0
    if t <= 1.0e-12:
        return coords[i0]

    mode = INTERP_MODE.strip().lower()

    if mode == "nearest":
        return coords[i0]

    i1 = i0 + 1
    if mode == "linear":
        return (1.0 - t) * coords[i0] + t * coords[i1]

    if mode == "spline":
        im1 = max(i0 - 1, 0)
        i2  = min(i1 + 1, n_traj - 1)

        p0 = coords[im1]
        p1 = coords[i0]
        p2 = coords[i1]
        p3 = coords[i2]

        t2 = t * t
        t3 = t2 * t
        return 0.5 * (
            (2.0 * p1)
            + (-p0 + p2) * t
            + (2.0 * p0 - 5.0 * p1 + 4.0 * p2 - p3) * t2
            + (-p0 + 3.0 * p1 - 3.0 * p2 + p3) * t3
        )

    raise ValueError(f'INTERP_MODE must be "nearest", "linear", or "spline" (got {INTERP_MODE!r}).')


def get_zoom_focus_indices():
    """
    Zoom-focus atom list, guaranteed to include tritium + both fluorides.

    Preserves user-provided ZOOM_FOCUS_INDICES order and appends missing
    required indices at the end.
    """
    zoom_indices = list(ZOOM_FOCUS_INDICES)
    for required_idx in (TRITIUM_IDX, FLUORIDE_1_IDX, FLUORIDE_2_IDX):
        if required_idx not in zoom_indices:
            zoom_indices.append(required_idx)
    return zoom_indices


def get_seg5_focus_indices():
    """
    Segment-5 atom list, guaranteed to include tritium + both fluorides.

    Preserves user-provided FOCUS_INDICES order and appends missing required
    indices at the end.
    """
    seg5_indices = list(FOCUS_INDICES)
    for required_idx in (TRITIUM_IDX, FLUORIDE_1_IDX, FLUORIDE_2_IDX):
        if required_idx not in seg5_indices:
            seg5_indices.append(required_idx)
    return seg5_indices


# ── I/O ───────────────────────────────────────────────────────────────────────

def resolve_path(user_path, filename):
    """Return absolute path: use user_path if given, else look beside the .blend."""
    if user_path:
        return os.path.abspath(user_path)
    blend_dir = bpy.path.abspath("//")
    candidate = os.path.join(blend_dir, filename)
    return candidate if os.path.isfile(candidate) else None


def parse_pdb(filepath):
    """
    Parse ATOM/HETATM records.
    Returns:
      atoms_by_element : {element: [(x,y,z), …]}
      atom_order       : [element, …]  — one entry per atom, in PDB order
      flat_positions   : [(x,y,z), …] — positions in PDB order (Å)
    """
    atoms_by_element = {}
    atom_order       = []
    flat_positions   = []

    with open(filepath, "r") as fh:
        for line in fh:
            rec = line[:6].strip()
            if rec not in ("ATOM", "HETATM"):
                continue
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue

            element = line[76:78].strip() if len(line) > 76 else ""
            if not element:
                raw = line[12:16].strip()
                element = ''.join(c for c in raw if c.isalpha())
            element = element.capitalize()

            atoms_by_element.setdefault(element, []).append((x, y, z))
            atom_order.append(element)
            flat_positions.append((x, y, z))

    return atoms_by_element, atom_order, flat_positions


def parse_dcd(filepath, frame_step=1):
    """
    Parse a CHARMM/NAMD/OpenMM DCD file.
    Returns coords : np.ndarray  shape (n_sampled, natoms, 3)  float32  Å
    """
    def read_fortran_record(f):
        size = struct.unpack("i", f.read(4))[0]
        data = f.read(size)
        f.read(4)
        return data, size

    with open(filepath, "rb") as f:
        hdr, _         = read_fortran_record(f)
        n_frames       = struct.unpack_from("i", hdr, 4)[0]
        read_fortran_record(f)                           # title block
        natoms_data, _ = read_fortran_record(f)
        natoms         = struct.unpack("i", natoms_data)[0]
        print(f"  DCD: {n_frames} frames, {natoms} atoms")

        frame_start = f.tell()
        first_size  = struct.unpack("i", f.read(4))[0]
        f.seek(frame_start)
        has_cell = (first_size == 48)

        sampled_indices  = range(0, n_frames, frame_step)
        n_sampled        = len(sampled_indices)
        coords           = np.zeros((n_sampled, natoms, 3), dtype=np.float32)

        frame_size_bytes = (
            (8 + 48 if has_cell else 0) +
            3 * (8 + natoms * 4)
        )

        for out_i, traj_i in enumerate(sampled_indices):
            f.seek(frame_start + traj_i * frame_size_bytes)
            if has_cell:
                read_fortran_record(f)
            x_data, _ = read_fortran_record(f)
            y_data, _ = read_fortran_record(f)
            z_data, _ = read_fortran_record(f)
            coords[out_i, :, 0] = np.frombuffer(x_data, dtype=np.float32)
            coords[out_i, :, 1] = np.frombuffer(y_data, dtype=np.float32)
            coords[out_i, :, 2] = np.frombuffer(z_data, dtype=np.float32)
            if out_i % 100 == 0:
                print(f"  DCD: loaded frame {out_i}/{n_sampled}", end="\r")

        print(f"  DCD: {n_sampled} sampled frames  ({n_frames} total, step={frame_step})")
        return coords


# ── Blender scene helpers ─────────────────────────────────────────────────────

def make_material(element):
    """Create (or reuse) a Principled BSDF material for the given element."""
    mat_name = f"Atom_{element}"
    if mat_name in bpy.data.materials:
        return bpy.data.materials[mat_name]
    mat  = bpy.data.materials.new(name=mat_name)
    mat.use_nodes = True
    bsdf = mat.node_tree.nodes.get("Principled BSDF")
    bsdf.inputs["Base Color"].default_value = ATOM_COLORS.get(element, DEFAULT_COLOR)
    bsdf.inputs["Metallic"].default_value   = 0.15
    bsdf.inputs["Roughness"].default_value  = 0.25
    return mat


def make_material_with_color(mat_name, color_rgba):
    """Create (or reuse) a material with an explicit RGBA base color."""
    if mat_name in bpy.data.materials:
        return bpy.data.materials[mat_name]
    mat  = bpy.data.materials.new(name=mat_name)
    mat.use_nodes = True
    bsdf = mat.node_tree.nodes.get("Principled BSDF")
    bsdf.inputs["Base Color"].default_value = color_rgba
    bsdf.inputs["Metallic"].default_value   = 0.15
    bsdf.inputs["Roughness"].default_value  = 0.25
    return mat


def create_element_objects(element, positions, collection, object_name=None, material=None):
    """Build a vertex-instancer for one element. Returns (parent_obj, sphere_obj)."""
    import bmesh

    obj_name = object_name if object_name else element

    mesh  = bpy.data.meshes.new(obj_name)
    verts = [(x * SCALE, y * SCALE, z * SCALE) for x, y, z in positions]
    mesh.from_pydata(verts, [], [])
    mesh.update()

    parent_obj = bpy.data.objects.new(obj_name, mesh)
    collection.objects.link(parent_obj)
    parent_obj.instance_type = "VERTS"

    radius_bu   = ATOM_RADII.get(element, DEFAULT_RADIUS) * SCALE
    sphere_mesh = bpy.data.meshes.new(f"{obj_name}_sphere")
    bm = bmesh.new()
    bmesh.ops.create_uvsphere(bm, u_segments=SPHERE_SEGMENTS,
                              v_segments=SPHERE_RINGS, radius=radius_bu)
    bm.to_mesh(sphere_mesh)
    bm.free()
    sphere_mesh.polygons.foreach_set("use_smooth", [True] * len(sphere_mesh.polygons))
    if hasattr(sphere_mesh, "use_auto_smooth"):
        sphere_mesh.use_auto_smooth   = True
        sphere_mesh.auto_smooth_angle = math.radians(180)
    sphere_mesh.update()

    sphere_obj        = bpy.data.objects.new(f"{obj_name}_template", sphere_mesh)
    sphere_obj.parent = parent_obj
    collection.objects.link(sphere_obj)

    # Hide the template object itself in the viewport.
    # IMPORTANT: for child-vertex instancing, hide_render on the template can
    # also hide all of its instances in final renders.
    sphere_obj.hide_set(True)
    sphere_obj.hide_render = False

    mat = material if material is not None else make_material(element)
    if sphere_obj.data.materials:
        sphere_obj.data.materials[0] = mat
    else:
        sphere_obj.data.materials.append(mat)

    return parent_obj, sphere_obj


def set_bezier(obj):
    """Set all keyframe interpolations on obj to BEZIER for smooth motion."""
    if obj.animation_data and obj.animation_data.action:
        for fc in obj.animation_data.action.fcurves:
            for kp in fc.keyframe_points:
                kp.interpolation = 'BEZIER'
                kp.handle_left_type = 'AUTO_CLAMPED'
                kp.handle_right_type = 'AUTO_CLAMPED'


def get_initial_camera_local_offset(target_world):
    """
    Return (segment1_local_offset, camera_direction_unit_vector).

    The returned direction vector is used for distance-driven keyframes in
    later storyboard segments.
    """
    if CAM_INITIAL_XYZ is None:
        # Backward-compatible default view direction if no explicit XYZ is given.
        fallback_dir = Vector((1.0, -1.0, 0.7)).normalized()
        return fallback_dir * CAM_DIST_SEG1, fallback_dir

    if len(CAM_INITIAL_XYZ) != 3:
        raise ValueError("CAM_INITIAL_XYZ must be a 3-value tuple/list or None.")

    initial_vec = Vector(CAM_INITIAL_XYZ)
    space = CAM_INITIAL_SPACE.upper()
    if space == "WORLD":
        local_offset = initial_vec - target_world
    elif space == "LOCAL":
        local_offset = initial_vec
    else:
        raise ValueError('CAM_INITIAL_SPACE must be "WORLD" or "LOCAL".')

    if local_offset.length < 1.0e-9:
        raise ValueError("CAM_INITIAL_XYZ cannot be at AtomCentre (zero camera-target distance).")

    cam_dir = local_offset.normalized()
    return cam_dir * CAM_DIST_SEG1, cam_dir


def setup_camera_animation(cam, target, focus_pos, centroid, cluster_center):
    """
    Set up the camera rig and keyframe the storyboard animation.

    Camera rig
    ----------
    The camera is parented to AtomCentre (target). Its LOCAL location uses the
    direction implied by the initial camera offset, scaled by per-segment
    distances, so zooming does not depend on a separate direction setting.
    A TrackTo constraint keeps the lens pointed at the target regardless of
    where in the scene the target sits.

    Animation overview
    ------------------
    target.location keyframes move the look-at point:
        Seg 1–2  →  molecule centroid      (wide overview)
        Seg 3–4  →  zoom-focus centroid    (H + 2F close-up)
        Seg 5    →  cluster centroid    (focus group)

    cam.location keyframes (LOCAL space) change the orbital distance:
        Seg 1–2  →  CAM_DIST_SEG1 / CAM_DIST_SEG2   (wide to zoom start)
        Seg 3–4  →  CAM_DIST_SEG3 / CAM_DIST_SEG4   (close)
        Seg 5    →  CAM_DIST_SEG5                    (mid)

    Both sets of keyframes use Bezier interpolation for smooth motion.
    """
    seg1_local, cam_dir = get_initial_camera_local_offset(centroid)

    seg1_bl     = dcd_to_bl(SEG1_DCD_START)
    zoom_start_bl = dcd_to_bl(ZOOM_START_DCD)
    seg3_bl     = dcd_to_bl(SEG3_DCD_START)
    seg5_bl     = dcd_to_bl(SEG5_DCD_START)
    seg5_end_bl = dcd_to_bl(SEG5_DCD_END)

    # ── Camera LOCAL location keyframes (distance from target) ────────────────
    # Because the camera is parented to target, cam.location is in LOCAL space,
    # so these keyframes purely control distance — world position follows target.

    cam.location = seg1_local
    cam.keyframe_insert(data_path="location", frame=seg1_bl)

    cam.location = cam_dir * CAM_DIST_SEG2
    cam.keyframe_insert(data_path="location", frame=zoom_start_bl)   # transition starts

    cam.location = cam_dir * CAM_DIST_SEG3
    cam.keyframe_insert(data_path="location", frame=seg3_bl)   # zoomed in

    cam.location = cam_dir * CAM_DIST_SEG4
    cam.keyframe_insert(data_path="location", frame=seg5_bl)   # hold through seg 4

    cam.location = cam_dir * CAM_DIST_SEG5
    cam.keyframe_insert(data_path="location", frame=seg5_end_bl)  # backed off for cluster

    # ── Target WORLD location keyframes (pan / look-at point) ─────────────────
    target.location = centroid
    target.keyframe_insert(data_path="location", frame=seg1_bl)
    target.keyframe_insert(data_path="location", frame=zoom_start_bl)   # hold at centroid

    target.location = focus_pos
    target.keyframe_insert(data_path="location", frame=seg3_bl)   # pan to zoom-focus group
    target.keyframe_insert(data_path="location", frame=seg5_bl)   # hold through seg 4

    target.location = cluster_center
    target.keyframe_insert(data_path="location", frame=seg5_end_bl)  # pan to cluster

    # ── Smooth Bezier interpolation on all keyframes ──────────────────────────
    set_bezier(cam)
    set_bezier(target)


# ── Storyboard frame handler ──────────────────────────────────────────────────

def register_storyboard_handler(coords, instancer_objects, instancer_index_map):
    """
    Register a frame_change_pre handler that drives the storyboard each frame.

    Atom visibility by segment:
      Seg 1 (DCD 700– 900) : all atoms
      Seg 2 (DCD 900–1000) : all atoms (zoom can already be in progress)
      Seg 3 (DCD 1000–1150): TRITIUM_IDX, FLUORIDE_1_IDX, FLUORIDE_2_IDX
      Seg 4 (DCD 1150–1350): TRITIUM_IDX, FLUORIDE_1_IDX, FLUORIDE_2_IDX
      Seg 5 (DCD 1350–1550): FOCUS_INDICES (+ TRITIUM_IDX/F1/F2 guaranteed)
    """
    n_traj = coords.shape[0]

    # Pre-built visibility arrays for each segment
    seg23_vis = np.array(get_zoom_focus_indices(), dtype=np.int32)
    seg4_vis  = np.array(get_zoom_focus_indices(), dtype=np.int32)
    seg5_vis  = np.array(get_seg5_focus_indices(), dtype=np.int32)

    # Precompute Blender boundary frames
    hide_others_bl = dcd_to_bl(HIDE_OTHERS_DCD_START)
    seg4_bl = dcd_to_bl(SEG4_DCD_START)
    seg5_bl = dcd_to_bl(SEG5_DCD_START)

    def update_frame(scene, depsgraph=None):
        bf        = scene.frame_current
        frame_xyz = get_interpolated_coords(coords, bf)   # (natoms, 3) in Å

        # Determine visibility mode for this frame
        if bf < hide_others_bl:
            show_all = True
            vis_idx  = None
        else:
            show_all = False
            if bf < seg4_bl:       # zoom-focus phase before seg 4 boundary
                vis_idx = seg23_vis
            elif bf < seg5_bl:     # seg 4
                vis_idx = seg4_vis
            else:                  # seg 5
                vis_idx = seg5_vis

        # Update each instancer mesh from its own global atom index list
        for obj_name, obj in instancer_objects.items():
            g_idxs = instancer_index_map.get(obj_name)
            if g_idxs is None:
                continue

            if show_all:
                positions = frame_xyz[g_idxs] * SCALE
            else:
                mask      = np.isin(g_idxs, vis_idx)
                positions = frame_xyz[g_idxs[mask]] * SCALE

            mesh = obj.data
            mesh.clear_geometry()
            if len(positions):
                mesh.vertices.add(len(positions))
                mesh.vertices.foreach_set("co", positions.ravel())
            mesh.update()

    # Remove any previously registered handler from this script
    hlist  = bpy.app.handlers.frame_change_pre
    hlist[:] = [h for h in hlist if getattr(h, "__name__", "") != "mol_storyboard"]
    update_frame.__name__ = "mol_storyboard"
    hlist.append(update_frame)

    # Trigger immediately so the starting frame renders correctly
    update_frame(bpy.context.scene)
    print(f"Storyboard handler registered. {n_traj} sampled frames available.")


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    if FRAME_STEP < 1:
        raise ValueError("FRAME_STEP must be >= 1.")
    if INTERP_MODE.strip().lower() not in {"nearest", "linear", "spline"}:
        raise ValueError('INTERP_MODE must be "nearest", "linear", or "spline".')
    if not (SEG1_DCD_START <= ZOOM_START_DCD <= SEG5_DCD_END):
        raise ValueError("ZOOM_START_DCD must be between SEG1_DCD_START and SEG5_DCD_END.")
    if not (SEG1_DCD_START <= HIDE_OTHERS_DCD_START <= SEG5_DCD_END):
        raise ValueError("HIDE_OTHERS_DCD_START must be between SEG1_DCD_START and SEG5_DCD_END.")
    if ZOOM_START_DCD > HIDE_OTHERS_DCD_START:
        raise ValueError("ZOOM_START_DCD must be <= HIDE_OTHERS_DCD_START.")

    # 1. Clear scene
    bpy.ops.object.select_all(action="SELECT")
    bpy.ops.object.delete(use_global=False)
    for blk in list(bpy.data.meshes):    bpy.data.meshes.remove(blk)
    for blk in list(bpy.data.materials): bpy.data.materials.remove(blk)

    # 2. Parse PDB
    pdb_path = resolve_path(PDB_FILE, "structure.pdb")
    if not pdb_path:
        raise FileNotFoundError(
            "structure.pdb not found. Save your .blend next to it, or set PDB_FILE.")
    print(f"PDB: {pdb_path}")
    atoms_by_element, atom_order, flat_positions = parse_pdb(pdb_path)
    for elem, pos in atoms_by_element.items():
        print(f"  {elem:4s}: {len(pos):5d} atoms")

    # 3. Build Blender collection
    col_name = "Molecule"
    if col_name in bpy.data.collections:
        bpy.data.collections.remove(bpy.data.collections[col_name])
    mol_col = bpy.data.collections.new(col_name)
    bpy.context.scene.collection.children.link(mol_col)

    # 4. Build instancer objects (one per element), with dedicated objects for
    #    FLUORIDE_1_IDX and FLUORIDE_2_IDX so they can have distinct colors.
    element_global_indices = {}
    for g_idx, elem in enumerate(atom_order):
        element_global_indices.setdefault(elem, []).append(g_idx)

    instancer_objects = {}
    instancer_index_map = {}

    for element in sorted(element_global_indices.keys()):
        g_idxs = element_global_indices[element]

        if element != "F":
            positions = [flat_positions[i] for i in g_idxs]
            parent_obj, _ = create_element_objects(element, positions, mol_col, object_name=element)
            instancer_objects[element] = parent_obj
            instancer_index_map[element] = np.array(g_idxs, dtype=np.int32)
            print(f"  Instancer: {element} ({len(g_idxs)} atoms)")
            continue

        # Split fluorine into bulk + two singled-out fluorides.
        special_f = {FLUORIDE_1_IDX, FLUORIDE_2_IDX}
        f_bulk_idxs = [i for i in g_idxs if i not in special_f]

        if f_bulk_idxs:
            positions = [flat_positions[i] for i in f_bulk_idxs]
            parent_obj, _ = create_element_objects(element, positions, mol_col, object_name="F")
            instancer_objects["F"] = parent_obj
            instancer_index_map["F"] = np.array(f_bulk_idxs, dtype=np.int32)
            print(f"  Instancer: F ({len(f_bulk_idxs)} atoms)")

        if FLUORIDE_1_IDX in g_idxs:
            f1_color = FLUORIDE_1_COLOR if FLUORIDE_1_COLOR is not None else ATOM_COLORS.get("F", DEFAULT_COLOR)
            f1_mat   = make_material_with_color("Atom_F_1", f1_color)
            positions = [flat_positions[FLUORIDE_1_IDX]]
            parent_obj, _ = create_element_objects(
                "F", positions, mol_col, object_name="F_1", material=f1_mat
            )
            instancer_objects["F_1"] = parent_obj
            instancer_index_map["F_1"] = np.array([FLUORIDE_1_IDX], dtype=np.int32)
            print(f"  Instancer: F_1 (index {FLUORIDE_1_IDX})")
        else:
            print(f"WARNING: FLUORIDE_1_IDX {FLUORIDE_1_IDX} is not a fluorine atom in this PDB.")

        if FLUORIDE_2_IDX in g_idxs:
            f2_color = FLUORIDE_2_COLOR if FLUORIDE_2_COLOR is not None else ATOM_COLORS.get("F", DEFAULT_COLOR)
            f2_mat   = make_material_with_color("Atom_F_2", f2_color)
            positions = [flat_positions[FLUORIDE_2_IDX]]
            parent_obj, _ = create_element_objects(
                "F", positions, mol_col, object_name="F_2", material=f2_mat
            )
            instancer_objects["F_2"] = parent_obj
            instancer_index_map["F_2"] = np.array([FLUORIDE_2_IDX], dtype=np.int32)
            print(f"  Instancer: F_2 (index {FLUORIDE_2_IDX})")
        else:
            print(f"WARNING: FLUORIDE_2_IDX {FLUORIDE_2_IDX} is not a fluorine atom in this PDB.")

    # 5. Compute molecule centroid and zoom-focus position (from static PDB)
    n_atoms  = len(flat_positions)
    cx = sum(p[0] for p in flat_positions) / n_atoms * SCALE
    cy = sum(p[1] for p in flat_positions) / n_atoms * SCALE
    cz = sum(p[2] for p in flat_positions) / n_atoms * SCALE
    centroid = Vector((cx, cy, cz))

    valid_zoom_focus = [i for i in get_zoom_focus_indices() if i < n_atoms]
    if valid_zoom_focus:
        zfx = sum(flat_positions[i][0] for i in valid_zoom_focus) / len(valid_zoom_focus) * SCALE
        zfy = sum(flat_positions[i][1] for i in valid_zoom_focus) / len(valid_zoom_focus) * SCALE
        zfz = sum(flat_positions[i][2] for i in valid_zoom_focus) / len(valid_zoom_focus) * SCALE
        focus_pos = Vector((zfx, zfy, zfz))
    else:
        print("WARNING: zoom-focus indices out of range — using centroid for zoom.")
        focus_pos = centroid

    # Centroid of the seg-5 focus cluster (from static PDB positions)
    valid_focus = [i for i in get_seg5_focus_indices() if i < n_atoms]
    if valid_focus:
        fcx = sum(flat_positions[i][0] for i in valid_focus) / len(valid_focus) * SCALE
        fcy = sum(flat_positions[i][1] for i in valid_focus) / len(valid_focus) * SCALE
        fcz = sum(flat_positions[i][2] for i in valid_focus) / len(valid_focus) * SCALE
        cluster_center = Vector((fcx, fcy, fcz))
    else:
        print("WARNING: FOCUS_INDICES out of range — using molecule centroid for seg 5.")
        cluster_center = centroid

    # 6. AtomCentre empty (look-at target) + camera rig + sun light
    target = bpy.data.objects.new("AtomCentre", None)
    target.empty_display_type = "SPHERE"
    target.empty_display_size  = 0.3
    target.location = centroid
    mol_col.objects.link(target)

    cam_data = bpy.data.cameras.new("MolCamera")
    cam      = bpy.data.objects.new("MolCamera", cam_data)
    mol_col.objects.link(cam)
    bpy.context.scene.camera = cam

    # Parent camera to AtomCentre so that cam.location is a LOCAL offset from the
    # target.  Clearing the parent-inverse matrix keeps that offset clean.
    cam.parent                 = target
    cam.matrix_parent_inverse  = Matrix.Identity(4)

    # Initial local position (matches seg-1 camera keyframe setup)
    initial_cam_local, _ = get_initial_camera_local_offset(centroid)
    cam.location = initial_cam_local

    # TrackTo keeps the lens pointed at AtomCentre regardless of orbital distance
    c = cam.constraints.new(type="TRACK_TO")
    c.target     = target
    c.track_axis = "TRACK_NEGATIVE_Z"
    c.up_axis    = "UP_Y"

    extent          = max(abs(cx), abs(cy), abs(cz))
    sun_data        = bpy.data.lights.new("MolSun", type="SUN")
    sun_data.energy = 3.0
    sun             = bpy.data.objects.new("MolSun", sun_data)
    sun.location    = (cx, cy, cz + extent * 4.5)
    mol_col.objects.link(sun)

    # 7. Parse DCD trajectory
    traj_path = resolve_path(TRAJ_FILE, "traj_nvt.dcd")
    if not traj_path:
        raise FileNotFoundError(
            "traj_nvt.dcd not found. Save your .blend next to it, or set TRAJ_FILE.")
    print(f"DCD: {traj_path}")
    coords = parse_dcd(traj_path, frame_step=FRAME_STEP)

    # 8. Set up camera animation keyframes
    setup_camera_animation(cam, target, focus_pos, centroid, cluster_center)

    # 9. Set Blender timeline
    bl_start = dcd_to_bl(SEG1_DCD_START)
    bl_end   = dcd_to_bl(SEG5_DCD_END)
    bpy.context.scene.frame_start   = bl_start
    bpy.context.scene.frame_end     = bl_end
    bpy.context.scene.frame_current = bl_start
    print(f"Timeline: Blender frames {bl_start}–{bl_end}  "
          f"(DCD frames {SEG1_DCD_START}–{SEG5_DCD_END})")

    # 10. Register storyboard frame handler
    register_storyboard_handler(coords, instancer_objects, instancer_index_map)

    # 11. Switch viewport to Material Preview (GUI only)
    if bpy.context.screen is not None:
        for area in bpy.context.screen.areas:
            if area.type == "VIEW_3D":
                for space in area.spaces:
                    if space.type == "VIEW_3D":
                        space.shading.type = "MATERIAL"

    print("Done!  Press Space to play.  Z → Rendered mode for full quality.")


if __name__ == "__main__":
    main()
