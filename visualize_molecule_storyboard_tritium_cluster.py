"""
Blender Python script: Molecular visualizer + storyboard animator
-----------------------------------------------------------------
Run inside Blender via the Scripting workspace:
  1. Open Blender → Scripting tab
  2. Open this file (or paste it in)
  3. Adjust the USER-CONFIGURABLE SETTINGS block below
  4. Click "Run Script"

Storyboard (all timings are DCD frame numbers):
  Phase 1  all atoms + zoom to tritium
  Phase 2  tritium-focused hold
  Phase 3  fade out non-cluster atoms while fading sphere in
  Phase 4  fade sphere out immediately after sphere fade-in completes
  Phase 5  cluster+tritium tail
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
TRAJ_FILE = ""
DEFAULT_TRAJ_FILENAME = "traj_nvt_interp.dcd"

# Optional setup .blend containing camera, empty, track-to constraint, and lights.
# If left empty, the script creates a default camera/empty/light rig from scratch.
SETUP_BLEND_FILE  = ""   # e.g. "//setup_base.blend" or "/abs/path/setup_base.blend"
SETUP_CAMERA_NAME = ""   # optional exact camera object name in setup file
SETUP_EMPTY_NAME  = ""   # optional exact empty target object name in setup file

# Auto-light rig used when no setup file is provided (or setup has no lights).
# Lights are positioned relative to the camera view direction and TrackTo AtomCentre.
AUTO_LIGHT_DISTANCE_FACTOR = 1.6   # multiplied by current camera-target distance
AUTO_LIGHT_MIN_DISTANCE    = 8.0
AUTO_LIGHT_SIZE            = 2.4
AUTO_LIGHT_POWER_RIGHT     = 650.0
AUTO_LIGHT_POWER_LEFT      = 60.0
AUTO_LIGHT_POWER_TOP       = 220.0
AUTO_LIGHT_SPREAD_DEGREES  = 85.0

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

# ── Story timing (DCD frames) ─────────────────────────────────────────────────
# Storyboard timeline window.
STORY_DCD_START = 900
STORY_DCD_STEPS = 660

# Storyboard phase timing (DCD frames). Remaining frames become phase 6 tail.
PHASE1_ZOOM_TO_T_DCD_FRAMES = 200
PHASE2_T_FOCUS_HOLD_DCD_FRAMES = 0

# Fade knobs (DCD frames).
FADE_NON_CLUSTER_OUT_DCD_FRAMES = 60
FADE_SPHERE_IN_DCD_FRAMES = 110
FADE_SPHERE_OUT_DCD_FRAMES = 40

# ── DCD ↔ Blender frame mapping ───────────────────────────────────────────────
# Blender frame 1 corresponds to this DCD frame number.
DCD_FRAME_OFFSET = STORY_DCD_START

# ── Key atom indices (0-based, matching PDB ATOM/HETATM record order) ─────────
# NOTE: PDB HETATM serial numbers are 1-based; subtract 1 to get 0-based index.
#   e.g. PDB serial 1681 (TR1/Tritium)  → index 1680
TRITIUM_IDX = 1680

# ── Cluster atom list ─────────────────────────────────────────────────────────
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
# Storyboard distances are controlled by CAM_DIST_ALL_ATOMS/TRITIUM/SPHERE/CLUSTER.
# In other words, CAM_INITIAL_XYZ chooses "where from", CAM_DIST_* chooses
# "how far".
# CAM_INITIAL_SPACE controls whether coordinates are read in AtomCentre local
# space ("LOCAL") or Blender world space ("WORLD").
#
# If left as None, a legacy fallback direction is used.
CAM_INITIAL_XYZ   = (-3.7, -9.4, 2.3) # tuned left/up view angle
CAM_INITIAL_SPACE = "WORLD"           # "WORLD" or "LOCAL"

# Distance (Blender units; 1 BU = 10 Å) from AtomCentre.
# Larger = further away / more zoomed out.
CAM_DIST_ALL_ATOMS = 10.0  # phase 1 start (wide shot)
CAM_DIST_TRITIUM   = 4.0   # zoomed-in tritium view
CAM_DIST_SPHERE    = 4.0   # zoomed-out distance before sphere fade-in
CAM_DIST_CLUSTER   = 4.0   # cluster view after blend-in

# ── Display ───────────────────────────────────────────────────────────────────
SCALE = 0.1          # Å → Blender units  (1 Å = 0.1 BU)

SPHERE_SEGMENTS = 64   # longitude divisions (higher = smoother spheres)
SPHERE_RINGS    = 32   # latitude divisions

# Boundary style shown during sphere phases (3-5).
# Options: "sphere", "none"
BOUNDARY_MODE = "sphere"

# Static reference sphere (centered at the target empty).
SPHERE_RADIUS_ANGSTROM = 7.7
SPHERE_COLOR = (0.82, 0.82, 0.82, 1.0)  # light gray
SPHERE_FINAL_ALPHA = 0.5                # 50% transparent

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
ATOM_RADII['Li'] *= 0.4

# ── CPK colours (RGBA, linear sRGB, values 0.0–1.0) ──────────────────────────
ATOM_COLORS = {
    "F"  : (0.000, 1.000, 0.000, 1.0),   # green
    "Li" : (0.000, 0.000, 1.000, 1.0),   # blue
    "Be" : (1.000, 1.000, 0.000, 1.0),   # yellow
    "H"  : (1.000, 0.000, 0.000, 1.0),   # red (tritium)
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


def get_cluster_focus_indices():
    """
    Cluster atom list, guaranteed to include TRITIUM_IDX.

    Preserves user-provided FOCUS_INDICES order and appends TRITIUM_IDX if
    missing.
    """
    cluster_indices = list(FOCUS_INDICES)
    if TRITIUM_IDX not in cluster_indices:
        cluster_indices.append(TRITIUM_IDX)
    return cluster_indices


def get_story_dcd_boundaries():
    """
    Return storyboard boundary frames in DCD units.

    Keys:
      start
      zoom_to_t_end
      t_focus_hold_end
      transition_start
      transition_end
      with_sphere_hold_end
      sphere_fade_out_end
      end
    """
    start = int(STORY_DCD_START)
    end = start + int(STORY_DCD_STEPS)
    zoom_to_t_end = start + int(PHASE1_ZOOM_TO_T_DCD_FRAMES)
    t_focus_hold_end = zoom_to_t_end + int(PHASE2_T_FOCUS_HOLD_DCD_FRAMES)
    transition_start = t_focus_hold_end
    transition_end = transition_start + max(
        int(FADE_NON_CLUSTER_OUT_DCD_FRAMES),
        int(FADE_SPHERE_IN_DCD_FRAMES),
    )
    # No plateau: as soon as sphere fade-in reaches SPHERE_FINAL_ALPHA, fade-out starts.
    with_sphere_hold_end = transition_start + int(FADE_SPHERE_IN_DCD_FRAMES)
    sphere_fade_out_end = with_sphere_hold_end + int(FADE_SPHERE_OUT_DCD_FRAMES)
    return {
        "start": start,
        "zoom_to_t_end": zoom_to_t_end,
        "t_focus_hold_end": t_focus_hold_end,
        "transition_start": transition_start,
        "transition_end": transition_end,
        "with_sphere_hold_end": with_sphere_hold_end,
        "sphere_fade_out_end": sphere_fade_out_end,
        "end": end,
    }


# ── I/O ───────────────────────────────────────────────────────────────────────

def resolve_path(user_path, filename):
    """Return absolute path: use user_path if given, else look beside the .blend."""
    if user_path:
        return os.path.abspath(user_path)
    blend_dir = bpy.path.abspath("//")
    candidate = os.path.join(blend_dir, filename)
    return candidate if os.path.isfile(candidate) else None


def resolve_setup_blend_path(user_path):
    """Resolve optional setup .blend path. Returns absolute path or None."""
    if not user_path:
        return None

    # Supports Blender-style // relative paths and absolute paths.
    candidate = os.path.abspath(bpy.path.abspath(user_path))
    if os.path.isfile(candidate):
        return candidate

    candidate = os.path.abspath(user_path)
    if os.path.isfile(candidate):
        return candidate

    return None


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
    bsdf.inputs["Metallic"].default_value   = 0.0
    bsdf.inputs["Roughness"].default_value  = 0.90
    return mat


def make_material_with_color(mat_name, color_rgba):
    """Create (or reuse) a material with an explicit RGBA base color."""
    if mat_name in bpy.data.materials:
        return bpy.data.materials[mat_name]
    mat  = bpy.data.materials.new(name=mat_name)
    mat.use_nodes = True
    bsdf = mat.node_tree.nodes.get("Principled BSDF")
    bsdf.inputs["Base Color"].default_value = color_rgba
    bsdf.inputs["Metallic"].default_value   = 0.0
    bsdf.inputs["Roughness"].default_value  = 0.90
    return mat


def set_material_alpha(mat, alpha):
    """Set material alpha (0..1) for fade transitions."""
    if mat is None:
        return

    a = max(0.0, min(1.0, float(alpha)))
    if not mat.use_nodes:
        mat.use_nodes = True
    bsdf = mat.node_tree.nodes.get("Principled BSDF")
    if bsdf is not None and "Alpha" in bsdf.inputs:
        bsdf.inputs["Alpha"].default_value = a

    # Eevee: avoid BLEND sorting artifacts on many overlapping transparent spheres.
    # Prefer dithered/hashed-style transparency during fades.
    fading = (a < 0.999)

    if hasattr(mat, "blend_method"):
        try:
            mat.blend_method = "HASHED" if fading else "OPAQUE"
        except Exception:
            # Fallback for Blender builds that do not expose HASHED here.
            try:
                mat.blend_method = "BLEND" if fading else "OPAQUE"
            except Exception:
                pass

    # Blender 4.x Eevee Next path (when available).
    if hasattr(mat, "surface_render_method"):
        try:
            mat.surface_render_method = "DITHERED" if fading else "OPAQUE"
        except Exception:
            pass

    if hasattr(mat, "shadow_method"):
        try:
            mat.shadow_method = "HASHED" if fading else "OPAQUE"
        except Exception:
            pass

    if hasattr(mat, "use_transparent_shadow"):
        try:
            mat.use_transparent_shadow = fading
        except Exception:
            pass

    if hasattr(mat, "show_transparent_back"):
        try:
            mat.show_transparent_back = not fading
        except Exception:
            pass


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


def create_reference_sphere(collection, center_empty):
    """Create a semi-transparent reference sphere centered on the given empty."""
    import bmesh

    mesh = bpy.data.meshes.new("TritiumReferenceSphere")
    bm = bmesh.new()
    bmesh.ops.create_uvsphere(
        bm,
        u_segments=SPHERE_SEGMENTS,
        v_segments=SPHERE_RINGS,
        radius=SPHERE_RADIUS_ANGSTROM * SCALE,
    )
    bm.to_mesh(mesh)
    bm.free()
    mesh.polygons.foreach_set("use_smooth", [True] * len(mesh.polygons))
    if hasattr(mesh, "use_auto_smooth"):
        mesh.use_auto_smooth = True
        mesh.auto_smooth_angle = math.radians(180)
    mesh.update()

    sphere_obj = bpy.data.objects.new("TritiumReferenceSphere", mesh)
    sphere_obj.parent = center_empty
    sphere_obj.matrix_parent_inverse = Matrix.Identity(4)
    sphere_obj.location = Vector((0.0, 0.0, 0.0))
    collection.objects.link(sphere_obj)

    sphere_mat = make_material_with_color("Atom_TritiumReferenceSphere", SPHERE_COLOR)
    if sphere_obj.data.materials:
        sphere_obj.data.materials[0] = sphere_mat
    else:
        sphere_obj.data.materials.append(sphere_mat)
    set_material_alpha(sphere_mat, 0.0)
    return sphere_obj, sphere_mat


def load_setup_rig(setup_blend_path, collection):
    """
    Append camera/empty/light objects from a setup .blend and return (camera, target, has_lights).
    """
    with bpy.data.libraries.load(setup_blend_path, link=False) as (data_from, data_to):
        data_to.objects = data_from.objects

    appended = [obj for obj in data_to.objects if obj is not None]
    allowed_types = {"CAMERA", "EMPTY", "LIGHT"}
    rig_objects = []

    for obj in appended:
        if obj.type in allowed_types:
            if collection not in obj.users_collection:
                collection.objects.link(obj)
            rig_objects.append(obj)
        else:
            bpy.data.objects.remove(obj, do_unlink=True)

    cameras = [o for o in rig_objects if o.type == "CAMERA"]
    empties = [o for o in rig_objects if o.type == "EMPTY"]
    lights  = [o for o in rig_objects if o.type == "LIGHT"]

    if not cameras:
        raise ValueError(f"Setup file '{setup_blend_path}' has no camera objects.")
    if not empties:
        raise ValueError(f"Setup file '{setup_blend_path}' has no empty objects for TrackTo target.")

    if SETUP_CAMERA_NAME:
        cam = next((o for o in cameras if o.name == SETUP_CAMERA_NAME), None)
        if cam is None:
            raise ValueError(
                f"SETUP_CAMERA_NAME '{SETUP_CAMERA_NAME}' not found in setup file '{setup_blend_path}'."
            )
    else:
        cam = next((o for o in cameras if o.name == "Camera"), cameras[0])

    if SETUP_EMPTY_NAME:
        target = next((o for o in empties if o.name == SETUP_EMPTY_NAME), None)
        if target is None:
            raise ValueError(
                f"SETUP_EMPTY_NAME '{SETUP_EMPTY_NAME}' not found in setup file '{setup_blend_path}'."
            )
    else:
        track = next((c for c in cam.constraints if c.type == "TRACK_TO" and c.target is not None), None)
        if track and track.target.type == "EMPTY":
            target = track.target
        else:
            target = next((o for o in empties if o.name == "Empty"), empties[0])

    # Ensure camera tracks the chosen target.
    ensure_track_to_constraint(cam, target)

    return cam, target, bool(lights)


def clear_object_animation(obj):
    """Remove existing animation data from an object, if any."""
    if obj is not None and obj.animation_data:
        obj.animation_data_clear()


def ensure_track_to_constraint(obj, target):
    """Ensure obj has a TrackTo constraint aimed at target."""
    track = next((c for c in obj.constraints if c.type == "TRACK_TO"), None)
    if track is None:
        track = obj.constraints.new(type="TRACK_TO")
    track.target     = target
    track.track_axis = "TRACK_NEGATIVE_Z"
    track.up_axis    = "UP_Y"
    return track


def create_area_light(name, location, power, size, collection, target):
    """Create one area light at location, linked to collection, tracking target."""
    light_data = bpy.data.lights.new(name, type="AREA")
    light_data.energy = power
    light_data.shape  = "SQUARE"
    light_data.size   = size
    if hasattr(light_data, "spread"):
        light_data.spread = math.radians(float(AUTO_LIGHT_SPREAD_DEGREES))
    if hasattr(light_data, "use_shadow"):
        light_data.use_shadow = True

    light_obj = bpy.data.objects.new(name, light_data)
    light_obj.location = location
    collection.objects.link(light_obj)
    ensure_track_to_constraint(light_obj, target)
    return light_obj


def create_default_three_point_lights(collection, target, cam):
    """
    Create right/left/top area lights based on camera view direction.
    """
    target_pos = target.matrix_world.translation.copy()
    cam_pos    = cam.matrix_world.translation.copy()
    cam_vec    = cam_pos - target_pos
    if cam_vec.length < 1.0e-9:
        cam_vec = Vector((0.0, -1.0, 0.0))
    view_dir = cam_vec.normalized()  # direction from target to camera

    world_up = Vector((0.0, 0.0, 1.0))
    side_dir = view_dir.cross(world_up)
    if side_dir.length < 1.0e-9:
        side_dir = view_dir.cross(Vector((1.0, 0.0, 0.0)))
    side_dir.normalize()

    top_dir = world_up - view_dir * world_up.dot(view_dir)
    if top_dir.length < 1.0e-9:
        top_dir = Vector((0.0, 1.0, 0.0))
    top_dir.normalize()

    dist = max(cam_vec.length * AUTO_LIGHT_DISTANCE_FACTOR, AUTO_LIGHT_MIN_DISTANCE)

    # Asymmetric placement: strong camera-side key, weak opposite fill, and top rim.
    right_pos = target_pos + ( side_dir * 1.05 + top_dir * 0.25 + view_dir * 0.45) * dist
    left_pos  = target_pos + (-side_dir * 1.10 + top_dir * 0.10 - view_dir * 0.05) * dist
    top_pos   = target_pos + ( top_dir * 1.20 - view_dir * 0.10) * dist

    create_area_light("MolAreaRight", right_pos, AUTO_LIGHT_POWER_RIGHT, AUTO_LIGHT_SIZE, collection, target)
    create_area_light("MolAreaLeft",  left_pos,  AUTO_LIGHT_POWER_LEFT,  AUTO_LIGHT_SIZE, collection, target)
    create_area_light("MolAreaTop",   top_pos,   AUTO_LIGHT_POWER_TOP,   AUTO_LIGHT_SIZE, collection, target)


def force_black_world_background(scene):
    """Force world background to black and remove world light contribution."""
    world = scene.world
    if world is None:
        world = bpy.data.worlds.new("World")
        scene.world = world

    world.use_nodes = True
    bg = world.node_tree.nodes.get("Background")
    if bg is not None:
        bg.inputs["Color"].default_value = (0.0, 0.0, 0.0, 1.0)
        bg.inputs["Strength"].default_value = 0.0
    else:
        world.color = (0.0, 0.0, 0.0)


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
        return fallback_dir * CAM_DIST_ALL_ATOMS, fallback_dir

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
    return cam_dir * CAM_DIST_ALL_ATOMS, cam_dir


def setup_camera_animation(cam, target, centroid, tritium_anchor, cluster_center, story_dcd):
    """
    Set up the camera rig and keyframe the storyboard animation.
    """
    cam_start_local, cam_dir = get_initial_camera_local_offset(centroid)

    start_bl = dcd_to_bl(story_dcd["start"])
    zoom_to_t_end_bl = dcd_to_bl(story_dcd["zoom_to_t_end"])
    t_focus_hold_end_bl = dcd_to_bl(story_dcd["t_focus_hold_end"])
    transition_end_bl = dcd_to_bl(story_dcd["transition_end"])
    with_sphere_hold_end_bl = dcd_to_bl(story_dcd["with_sphere_hold_end"])
    sphere_fade_out_end_bl = dcd_to_bl(story_dcd["sphere_fade_out_end"])
    end_bl = dcd_to_bl(story_dcd["end"])

    # Camera LOCAL location keyframes (distance from target)
    cam.location = cam_start_local
    cam.keyframe_insert(data_path="location", frame=start_bl)
    cam.location = cam_dir * CAM_DIST_TRITIUM
    cam.keyframe_insert(data_path="location", frame=zoom_to_t_end_bl)
    cam.location = cam_dir * CAM_DIST_TRITIUM
    cam.keyframe_insert(data_path="location", frame=t_focus_hold_end_bl)
    cam.location = cam_dir * CAM_DIST_SPHERE
    cam.keyframe_insert(data_path="location", frame=transition_end_bl)
    cam.location = cam_dir * CAM_DIST_SPHERE
    cam.keyframe_insert(data_path="location", frame=with_sphere_hold_end_bl)
    cam.location = cam_dir * CAM_DIST_CLUSTER
    cam.keyframe_insert(data_path="location", frame=sphere_fade_out_end_bl)
    cam.location = cam_dir * CAM_DIST_CLUSTER
    cam.keyframe_insert(data_path="location", frame=end_bl)

    # Target WORLD location keyframes (look-at point)
    target.location = centroid
    target.keyframe_insert(data_path="location", frame=start_bl)
    target.location = tritium_anchor
    target.keyframe_insert(data_path="location", frame=zoom_to_t_end_bl)
    target.location = tritium_anchor
    target.keyframe_insert(data_path="location", frame=t_focus_hold_end_bl)
    target.location = cluster_center
    target.keyframe_insert(data_path="location", frame=transition_end_bl)
    target.location = cluster_center
    target.keyframe_insert(data_path="location", frame=with_sphere_hold_end_bl)
    target.location = cluster_center
    target.keyframe_insert(data_path="location", frame=sphere_fade_out_end_bl)
    target.location = cluster_center
    target.keyframe_insert(data_path="location", frame=end_bl)

    # Smooth Bezier interpolation on all keyframes
    set_bezier(cam)
    set_bezier(target)


# ── Storyboard frame handler ──────────────────────────────────────────────────

def register_storyboard_handler(
    coords,
    instancer_objects,
    instancer_index_map,
    noncluster_object_names,
    boundary_materials,
    target,
    dynamic_center_indices,
    follow_cluster_center,
    story_dcd,
):
    """
    Register a frame_change_pre handler that drives the storyboard each frame.
    """
    n_traj = coords.shape[0]

    transition_start_bl = dcd_to_bl(story_dcd["transition_start"])
    transition_end_bl = dcd_to_bl(story_dcd["transition_end"])
    with_sphere_hold_end_bl = dcd_to_bl(story_dcd["with_sphere_hold_end"])
    sphere_fade_out_end_bl = dcd_to_bl(story_dcd["sphere_fade_out_end"])
    noncluster_fade_out_end_bl = transition_start_bl + int(FADE_NON_CLUSTER_OUT_DCD_FRAMES)
    sphere_fade_in_end_bl = transition_start_bl + int(FADE_SPHERE_IN_DCD_FRAMES)

    # Cache per-instancer material for fast alpha updates.
    obj_materials = {}
    for obj_name, obj in instancer_objects.items():
        mat = None
        for child in obj.children:
            if child.type == "MESH" and child.data and child.data.materials:
                mat = child.data.materials[0]
                break
        obj_materials[obj_name] = mat

    def phase_progress(frame, start_frame, end_frame):
        if end_frame <= start_frame:
            return 1.0 if frame >= end_frame else 0.0
        t = (frame - start_frame) / float(end_frame - start_frame)
        return max(0.0, min(1.0, t))

    def update_frame(scene, depsgraph=None):
        bf = scene.frame_current
        frame_xyz = get_interpolated_coords(coords, bf)   # (natoms, 3) in Å

        # During transition and sphere phases, keep AtomCentre on the live cluster centroid so
        # the reference sphere tracks the moving cluster.
        if follow_cluster_center and bf >= transition_start_bl and len(dynamic_center_indices):
            center_xyz = frame_xyz[dynamic_center_indices].mean(axis=0) * SCALE
            target.location = Vector(
                (float(center_xyz[0]), float(center_xyz[1]), float(center_xyz[2]))
            )

        # Non-cluster atoms: stay visible through zoom-to-T/focus, then fade out.
        if bf < transition_start_bl:
            noncluster_alpha = 1.0
        elif bf < noncluster_fade_out_end_bl:
            noncluster_alpha = 1.0 - phase_progress(
                bf, transition_start_bl, noncluster_fade_out_end_bl
            )
        else:
            noncluster_alpha = 0.0

        # Sphere boundary: fade in during transition, then immediately fade out.
        if bf < transition_start_bl:
            boundary_alpha = 0.0
        elif bf < sphere_fade_in_end_bl:
            boundary_alpha = SPHERE_FINAL_ALPHA * phase_progress(
                bf, transition_start_bl, sphere_fade_in_end_bl
            )
        elif bf < sphere_fade_out_end_bl:
            boundary_alpha = SPHERE_FINAL_ALPHA * (
                1.0 - phase_progress(bf, with_sphere_hold_end_bl, sphere_fade_out_end_bl)
            )
        else:
            boundary_alpha = 0.0

        for obj_name, mat in obj_materials.items():
            if obj_name in noncluster_object_names:
                set_material_alpha(mat, noncluster_alpha)
            else:
                set_material_alpha(mat, 1.0)

        for mat in boundary_materials:
            set_material_alpha(mat, boundary_alpha)

        # Update each instancer mesh from its own global atom index list
        for obj_name, obj in instancer_objects.items():
            g_idxs = instancer_index_map.get(obj_name)
            if g_idxs is None:
                continue

            if obj_name in noncluster_object_names and bf >= transition_end_bl:
                positions = np.empty((0, 3), dtype=np.float32)
            else:
                positions = frame_xyz[g_idxs] * SCALE

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
    if STORY_DCD_STEPS < 1:
        raise ValueError("STORY_DCD_STEPS must be >= 1.")
    if DCD_FRAME_OFFSET != STORY_DCD_START:
        raise ValueError("DCD_FRAME_OFFSET must match STORY_DCD_START in this storyboard script.")
    boundary_mode = BOUNDARY_MODE.strip().lower()
    if boundary_mode not in {"sphere", "none"}:
        raise ValueError('BOUNDARY_MODE must be one of: "sphere", "none".')
    timing_values = [
        PHASE1_ZOOM_TO_T_DCD_FRAMES,
        PHASE2_T_FOCUS_HOLD_DCD_FRAMES,
        FADE_NON_CLUSTER_OUT_DCD_FRAMES,
        FADE_SPHERE_IN_DCD_FRAMES,
        FADE_SPHERE_OUT_DCD_FRAMES,
    ]
    if any(p < 0 for p in timing_values):
        raise ValueError("All PHASE*/FADE*_DCD_FRAMES values must be >= 0.")
    transition_frames = max(int(FADE_NON_CLUSTER_OUT_DCD_FRAMES), int(FADE_SPHERE_IN_DCD_FRAMES))
    phase_total = (
        int(PHASE1_ZOOM_TO_T_DCD_FRAMES)
        + int(PHASE2_T_FOCUS_HOLD_DCD_FRAMES)
        + transition_frames
        + int(FADE_SPHERE_OUT_DCD_FRAMES)
    )
    if not (0.0 <= SPHERE_FINAL_ALPHA <= 1.0):
        raise ValueError("SPHERE_FINAL_ALPHA must be in [0, 1].")
    if SPHERE_RADIUS_ANGSTROM <= 0.0 and boundary_mode == "sphere":
        raise ValueError("SPHERE_RADIUS_ANGSTROM must be > 0.")

    story_dcd = get_story_dcd_boundaries()
    if story_dcd["sphere_fade_out_end"] > story_dcd["end"]:
        raise ValueError(
            "Phase durations exceed STORY_DCD_STEPS: "
            f"total={phase_total}, STORY_DCD_STEPS={STORY_DCD_STEPS}. "
            "Reduce PHASE*/FADE*_DCD_FRAMES or increase STORY_DCD_STEPS."
        )

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

    # 4. Build instancer objects.
    #    - One dedicated tritium object (always visible).
    #    - Per-element cluster objects (always visible after transition).
    #    - Per-element non-cluster objects (faded out during transition).
    element_global_indices = {}
    for g_idx, elem in enumerate(atom_order):
        element_global_indices.setdefault(elem, []).append(g_idx)

    if not (0 <= TRITIUM_IDX < len(atom_order)):
        raise ValueError(f"TRITIUM_IDX {TRITIUM_IDX} is out of range for this PDB ({len(atom_order)} atoms).")

    n_atoms = len(flat_positions)
    valid_cluster = [i for i in get_cluster_focus_indices() if 0 <= i < n_atoms]
    invalid_cluster = [i for i in get_cluster_focus_indices() if not (0 <= i < n_atoms)]
    if invalid_cluster:
        print(f"WARNING: Ignoring out-of-range FOCUS_INDICES: {invalid_cluster}")
    cluster_focus_set = set(valid_cluster)
    cluster_focus_set.discard(TRITIUM_IDX)

    instancer_objects = {}
    instancer_index_map = {}
    noncluster_object_names = set()

    tritium_element = atom_order[TRITIUM_IDX]
    tritium_obj_name = f"Tritium_{TRITIUM_IDX}_{tritium_element}"
    tritium_color = ATOM_COLORS.get(tritium_element, DEFAULT_COLOR)
    tritium_mat = make_material_with_color(f"Atom_{tritium_obj_name}", tritium_color)
    tritium_parent, _ = create_element_objects(
        tritium_element,
        [flat_positions[TRITIUM_IDX]],
        mol_col,
        object_name=tritium_obj_name,
        material=tritium_mat,
    )
    instancer_objects[tritium_obj_name] = tritium_parent
    instancer_index_map[tritium_obj_name] = np.array([TRITIUM_IDX], dtype=np.int32)
    print(f"  Instancer: {tritium_obj_name} (index {TRITIUM_IDX})")

    for element in sorted(element_global_indices.keys()):
        g_idxs = element_global_indices[element]
        cluster_idxs = [i for i in g_idxs if i in cluster_focus_set]
        noncluster_idxs = [i for i in g_idxs if i != TRITIUM_IDX and i not in cluster_focus_set]

        if cluster_idxs:
            cluster_name = f"Cluster_{element}"
            cluster_positions = [flat_positions[i] for i in cluster_idxs]
            cluster_color = ATOM_COLORS.get(element, DEFAULT_COLOR)
            cluster_mat = make_material_with_color(f"Atom_{cluster_name}", cluster_color)
            parent_obj, _ = create_element_objects(
                element,
                cluster_positions,
                mol_col,
                object_name=cluster_name,
                material=cluster_mat,
            )
            instancer_objects[cluster_name] = parent_obj
            instancer_index_map[cluster_name] = np.array(cluster_idxs, dtype=np.int32)
            print(f"  Instancer: {cluster_name} ({len(cluster_idxs)} atoms)")

        if noncluster_idxs:
            bulk_name = f"Bulk_{element}"
            bulk_positions = [flat_positions[i] for i in noncluster_idxs]
            bulk_color = ATOM_COLORS.get(element, DEFAULT_COLOR)
            bulk_mat = make_material_with_color(f"Atom_{bulk_name}", bulk_color)
            parent_obj, _ = create_element_objects(
                element,
                bulk_positions,
                mol_col,
                object_name=bulk_name,
                material=bulk_mat,
            )
            instancer_objects[bulk_name] = parent_obj
            instancer_index_map[bulk_name] = np.array(noncluster_idxs, dtype=np.int32)
            noncluster_object_names.add(bulk_name)
            print(f"  Instancer: {bulk_name} ({len(noncluster_idxs)} atoms)")

    # 5. Compute static centroids from PDB
    cx = sum(p[0] for p in flat_positions) / n_atoms * SCALE
    cy = sum(p[1] for p in flat_positions) / n_atoms * SCALE
    cz = sum(p[2] for p in flat_positions) / n_atoms * SCALE
    centroid = Vector((cx, cy, cz))
    if valid_cluster:
        fcx = sum(flat_positions[i][0] for i in valid_cluster) / len(valid_cluster) * SCALE
        fcy = sum(flat_positions[i][1] for i in valid_cluster) / len(valid_cluster) * SCALE
        fcz = sum(flat_positions[i][2] for i in valid_cluster) / len(valid_cluster) * SCALE
        cluster_center = Vector((fcx, fcy, fcz))
    else:
        print("WARNING: FOCUS_INDICES out of range — using molecule centroid for cluster center.")
        cluster_center = centroid

    # 6. Camera/target/light rig from optional setup file, else create from scratch
    setup_blend_path = resolve_setup_blend_path(SETUP_BLEND_FILE)
    if SETUP_BLEND_FILE and not setup_blend_path:
        raise FileNotFoundError(
            f"SETUP_BLEND_FILE not found: {SETUP_BLEND_FILE!r}. "
            "Use absolute path or Blender-relative //path."
        )

    if setup_blend_path:
        print(f"Setup rig: {setup_blend_path}")
        cam, target, has_lights = load_setup_rig(setup_blend_path, mol_col)
        bpy.context.scene.camera = cam

        # Place look-at target on molecule centroid and drive camera in local space.
        target.location = centroid
        cam.parent = target
        cam.matrix_parent_inverse = Matrix.Identity(4)

        # Keep segment-1 start offset logic consistent with non-setup mode.
        initial_cam_local, _ = get_initial_camera_local_offset(centroid)
        cam.location = initial_cam_local

        # Avoid stale keyframes from the setup file.
        clear_object_animation(cam)
        clear_object_animation(target)

        if not has_lights:
            print("WARNING: Setup file has no lights — creating default 3 area lights.")
            create_default_three_point_lights(mol_col, target, cam)
    else:
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
        # target. Clearing the parent-inverse matrix keeps that offset clean.
        cam.parent                 = target
        cam.matrix_parent_inverse  = Matrix.Identity(4)

        # Initial local position (matches seg-1 camera keyframe setup)
        initial_cam_local, _ = get_initial_camera_local_offset(centroid)
        cam.location = initial_cam_local

        # TrackTo keeps the lens pointed at AtomCentre regardless of orbital distance
        ensure_track_to_constraint(cam, target)
        create_default_three_point_lights(mol_col, target, cam)

    # 7. Parse DCD trajectory
    traj_path = resolve_path(TRAJ_FILE, DEFAULT_TRAJ_FILENAME)
    if not traj_path:
        raise FileNotFoundError(
            f"{DEFAULT_TRAJ_FILENAME} not found. Save your .blend next to it, or set TRAJ_FILE."
        )
    print(f"DCD: {traj_path}")
    coords = parse_dcd(traj_path, frame_step=FRAME_STEP)
    if TRITIUM_IDX >= coords.shape[1]:
        raise ValueError(
            f"TRITIUM_IDX {TRITIUM_IDX} is out of range for DCD with {coords.shape[1]} atoms."
        )

    # 8. Dynamic anchor points from trajectory
    tritium_zoom_end_bl = dcd_to_bl(story_dcd["zoom_to_t_end"])
    tritium_zoom_end_pos = get_interpolated_coords(coords, tritium_zoom_end_bl)[TRITIUM_IDX] * SCALE
    tritium_anchor = Vector((float(tritium_zoom_end_pos[0]), float(tritium_zoom_end_pos[1]), float(tritium_zoom_end_pos[2])))

    # Optional boundary visual (sphere).
    boundary_materials = []

    if boundary_mode == "sphere":
        _, sphere_material = create_reference_sphere(mol_col, target)
        boundary_materials.append(sphere_material)

    # 9. Set up camera animation keyframes
    setup_camera_animation(cam, target, centroid, tritium_anchor, cluster_center, story_dcd)

    # 10. Set Blender timeline
    bl_start = dcd_to_bl(story_dcd["start"])
    bl_end   = dcd_to_bl(story_dcd["end"])
    bpy.context.scene.frame_start   = bl_start
    bpy.context.scene.frame_end     = bl_end
    bpy.context.scene.frame_current = bl_start
    force_black_world_background(bpy.context.scene)
    print(f"Timeline: Blender frames {bl_start}–{bl_end}  "
          f"(DCD frames {story_dcd['start']}–{story_dcd['end']})")

    # 11. Register storyboard frame handler
    register_storyboard_handler(
        coords,
        instancer_objects,
        instancer_index_map,
        noncluster_object_names,
        boundary_materials,
        target,
        np.array(valid_cluster, dtype=np.int32),
        boundary_mode == "sphere",
        story_dcd,
    )

    # 12. Switch viewport to Material Preview (GUI only)
    if bpy.context.screen is not None:
        for area in bpy.context.screen.areas:
            if area.type == "VIEW_3D":
                for space in area.spaces:
                    if space.type == "VIEW_3D":
                        space.shading.type = "MATERIAL"

    print("Done!  Press Space to play.  Z → Rendered mode for full quality.")


if __name__ == "__main__":
    main()
