"""
Blender Python script: Molecular visualizer + trajectory animator
-----------------------------------------------------------------
Run inside Blender via the Scripting workspace:
  1. Open Blender → Scripting tab
  2. Open this file (or paste it in)
  3. Adjust settings in the USER-CONFIGURABLE SETTINGS block below
  4. Click "Run Script"

Static mode  — only structure.pdb is needed.
Animated mode — also requires traj_nvt.dcd (or any CHARMM/NAMD/OpenMM DCD).

Animation works via a frame_change_pre handler: Blender calls it every time
the timeline moves, and the handler swaps the instancer-mesh vertices to the
coordinates for that trajectory frame.  No keyframes are stored — scrubbing
and rendering both work normally.
"""

import bpy
import os
import math
import struct
import numpy as np


# ┌─────────────────────────────────────────────────────────────────────────────┐
# │                        USER-CONFIGURABLE SETTINGS                          │
# └─────────────────────────────────────────────────────────────────────────────┘

# ── File paths ────────────────────────────────────────────────────────────────
# Leave as "" to auto-detect files next to the saved .blend file.
PDB_FILE  = ""   # e.g. "/home/user/sim/structure.pdb"
TRAJ_FILE = ""   # e.g. "/home/user/sim/traj_nvt.dcd"  — leave "" to skip animation

# ── Trajectory playback ───────────────────────────────────────────────────────
# Use every Nth frame of the DCD.
#   FRAME_STEP = 1   → all 8160 frames  (smooth but slow to scrub)
#   FRAME_STEP = 10  → 816 frames       (good balance)
#   FRAME_STEP = 40  → 204 frames       (fast preview)
FRAME_STEP = 10

# ── Neighbour filter ──────────────────────────────────────────────────────────
# Show only atoms within NEIGHBOR_CUTOFF Ångströms of NEIGHBOR_ELEMENT.
# Set NEIGHBOR_CUTOFF = 0 to disable the filter and show all atoms.
#
# The H atom is the main feature of this simulation — a good starting cutoff
# is 5–8 Å, which captures the first coordination shell.
#
NEIGHBOR_ELEMENT = "H"    # centre atom type
NEIGHBOR_CUTOFF  = 4.0    # Å  — set to 0 to show everything

# ── Display ───────────────────────────────────────────────────────────────────
# Scale: Ångströms → Blender units  (0.1 → 1 Å = 0.1 BU)
SCALE = 0.1

# Sphere tessellation — increase for smoother spheres (costs render time)
SPHERE_SEGMENTS = 32   # longitude divisions
SPHERE_RINGS    = 16   # latitude divisions

# ── Atom radii (Ångströms) — EDIT THESE ──────────────────────────────────────
#
# These are the atoms present in this simulation.  Change the values below
# to rescale any atom type.  Larger number = bigger sphere.
#
#   Radius styles for reference:
#     Van der Waals  — "space-filling" look  (default below)
#     Covalent       — smaller, shows bonds more clearly
#     Ionic          — reflects charge state
#
ATOM_RADII = {
    # Element : radius in Å
    "F"  : 1.47,   # Fluorine   — van der Waals
    "Li" : 1.82,   # Lithium    — van der Waals
    "Be" : 1.53,   # Beryllium  — van der Waals
    "H"  : 1.20,   # Hydrogen   — van der Waals

    # ── Additional elements (used if present in the PDB) ─────────────────────
    "B"  : 1.92,
    "C"  : 1.70,
    "N"  : 1.55,
    "O"  : 1.52,
    "Na" : 2.27,
    "Mg" : 1.73,
    "Al" : 1.84,
    "Si" : 2.10,
    "P"  : 1.80,
    "S"  : 1.80,
    "Cl" : 1.75,
    "K"  : 2.75,
    "Ca" : 2.31,
}

# Fallback radius for any element not listed above
DEFAULT_RADIUS = 1.50  # Å

# ── CPK colours (RGBA, linear sRGB) — EDIT THESE ─────────────────────────────
#
# Tuple format: (Red, Green, Blue, Alpha)  — values from 0.0 to 1.0
#
ATOM_COLORS = {
    # Element : (R,    G,    B,    A)
    "F"  : (0.565, 0.878, 0.314, 1.0),   # green
    "Li" : (0.784, 0.502, 1.000, 1.0),   # violet
    "Be" : (0.765, 1.000, 0.000, 1.0),   # yellow-green
    "H"  : (1.000, 1.000, 1.000, 1.0),   # white

    # ── Additional elements ───────────────────────────────────────────────────
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

# Fallback colour for any element not listed above
DEFAULT_COLOR = (0.800, 0.800, 0.800, 1.0)   # light grey


# ┌─────────────────────────────────────────────────────────────────────────────┐
# │                             INTERNAL FUNCTIONS                             │
# └─────────────────────────────────────────────────────────────────────────────┘

def resolve_path(user_path, filename):
    """Return absolute path: use user_path if given, else look beside .blend."""
    if user_path:
        return os.path.abspath(user_path)
    blend_dir = bpy.path.abspath("//")
    candidate = os.path.join(blend_dir, filename)
    if os.path.isfile(candidate):
        return candidate
    return None


def parse_pdb(filepath):
    """
    Parse ATOM / HETATM records.
    Returns:
      atoms_by_element : {element: [(x,y,z), ...]}
      atom_order       : [element, ...]  — one entry per atom, in PDB order
    """
    atoms_by_element = {}
    atom_order = []

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

    return atoms_by_element, atom_order


def parse_dcd(filepath, frame_step=1):
    """
    Parse a CHARMM/NAMD/OpenMM DCD file.

    Returns:
      coords : np.ndarray, shape (n_sampled_frames, natoms, 3), float32, in Å
    """
    def read_fortran_record(f):
        size = struct.unpack("i", f.read(4))[0]
        data = f.read(size)
        f.read(4)
        return data, size

    with open(filepath, "rb") as f:
        # ── Header block 1 ────────────────────────────────────────────────────
        hdr, _ = read_fortran_record(f)
        n_frames = struct.unpack_from("i", hdr, 4)[0]

        # ── Header block 2: title ─────────────────────────────────────────────
        read_fortran_record(f)

        # ── Header block 3: natoms ────────────────────────────────────────────
        natoms_data, _ = read_fortran_record(f)
        natoms = struct.unpack("i", natoms_data)[0]

        print(f"  DCD: {n_frames} frames, {natoms} atoms")

        # ── Detect whether a unit-cell block precedes X/Y/Z each frame ────────
        # Peek at the first record size: 48 bytes = 6 doubles = cell block
        frame_start = f.tell()
        first_size = struct.unpack("i", f.read(4))[0]
        f.seek(frame_start)
        has_cell = (first_size == 48)

        # ── Allocate output array ─────────────────────────────────────────────
        sampled_indices = range(0, n_frames, frame_step)
        n_sampled = len(sampled_indices)
        coords = np.zeros((n_sampled, natoms, 3), dtype=np.float32)

        # ── Read frames ───────────────────────────────────────────────────────
        frame_size_bytes = (
            (8 + 48 if has_cell else 0) +   # cell block
            3 * (8 + natoms * 4)             # X, Y, Z blocks
        )

        for out_i, traj_i in enumerate(sampled_indices):
            # Seek directly to the frame start (fast random access)
            f.seek(frame_start + traj_i * frame_size_bytes)

            if has_cell:
                read_fortran_record(f)   # skip unit cell

            x_data, _ = read_fortran_record(f)
            y_data, _ = read_fortran_record(f)
            z_data, _ = read_fortran_record(f)

            coords[out_i, :, 0] = np.frombuffer(x_data, dtype=np.float32)
            coords[out_i, :, 1] = np.frombuffer(y_data, dtype=np.float32)
            coords[out_i, :, 2] = np.frombuffer(z_data, dtype=np.float32)

            if out_i % 100 == 0:
                print(f"  DCD: loaded frame {out_i}/{n_sampled}", end="\r")

        print(f"  DCD: loaded {n_sampled} frames ({n_frames} total, step={frame_step})")
        return coords


def make_material(element):
    """Create (or reuse) a Principled BSDF material for the given element."""
    mat_name = f"Atom_{element}"
    if mat_name in bpy.data.materials:
        return bpy.data.materials[mat_name]

    mat = bpy.data.materials.new(name=mat_name)
    mat.use_nodes = True
    bsdf = mat.node_tree.nodes.get("Principled BSDF")
    color = ATOM_COLORS.get(element, DEFAULT_COLOR)
    bsdf.inputs["Base Color"].default_value = color
    bsdf.inputs["Metallic"].default_value   = 0.15
    bsdf.inputs["Roughness"].default_value  = 0.25
    return mat


def create_element_objects(element, positions, collection):
    """
    Build a vertex-instancer for one element.
    Returns (parent_obj, sphere_obj).
    """
    # Parent mesh: vertices = atom positions
    mesh = bpy.data.meshes.new(f"{element}")
    verts = [(x * SCALE, y * SCALE, z * SCALE) for x, y, z in positions]
    mesh.from_pydata(verts, [], [])
    mesh.update()

    parent_obj = bpy.data.objects.new(f"{element}", mesh)
    collection.objects.link(parent_obj)
    parent_obj.instance_type = "VERTS"

    # Child sphere template
    import bmesh
    radius_bu = ATOM_RADII.get(element, DEFAULT_RADIUS) * SCALE
    sphere_mesh = bpy.data.meshes.new(f"{element}_sphere")
    bm = bmesh.new()
    bmesh.ops.create_uvsphere(
        bm,
        u_segments=SPHERE_SEGMENTS,
        v_segments=SPHERE_RINGS,
        radius=radius_bu,
    )
    bm.to_mesh(sphere_mesh)
    bm.free()

    # Smooth shading
    sphere_mesh.polygons.foreach_set(
        "use_smooth", [True] * len(sphere_mesh.polygons)
    )
    if hasattr(sphere_mesh, "use_auto_smooth"):
        sphere_mesh.use_auto_smooth = True
        sphere_mesh.auto_smooth_angle = math.radians(180)
    sphere_mesh.update()

    sphere_obj = bpy.data.objects.new(f"{element}_template", sphere_mesh)
    collection.objects.link(sphere_obj)
    sphere_obj.parent = parent_obj

    mat = make_material(element)
    if sphere_obj.data.materials:
        sphere_obj.data.materials[0] = mat
    else:
        sphere_obj.data.materials.append(mat)

    return parent_obj, sphere_obj


def register_trajectory_handler(coords, atom_order, instancer_objects):
    """
    Register a frame_change_pre handler that updates instancer mesh vertices
    from the preloaded trajectory array on every frame change.

    If NEIGHBOR_CUTOFF > 0, only atoms within that distance of any atom of
    type NEIGHBOR_ELEMENT are shown.  The vertex count of each instancer mesh
    is rebuilt each frame to reflect the filtered set.

    coords           : np.ndarray (n_frames, natoms, 3) in Å
    atom_order       : list of element symbols, length natoms, in PDB atom order
    instancer_objects: {element: parent_obj}
    """
    n_traj_frames = coords.shape[0]

    # Build element → numpy index array into the global atom list
    element_indices = {}
    for global_idx, elem in enumerate(atom_order):
        element_indices.setdefault(elem, []).append(global_idx)
    element_indices_np = {
        elem: np.array(idxs) for elem, idxs in element_indices.items()
    }

    # Indices of the centre atom(s) used for neighbour filtering
    centre_indices = element_indices_np.get(NEIGHBOR_ELEMENT, None)
    use_filter = (NEIGHBOR_CUTOFF > 0) and (centre_indices is not None)

    if use_filter:
        print(f"  Neighbour filter: showing atoms within {NEIGHBOR_CUTOFF} Å "
              f"of {NEIGHBOR_ELEMENT} ({len(centre_indices)} centre atom(s))")
    else:
        print("  Neighbour filter: disabled — showing all atoms")

    def update_frame(scene, depsgraph=None):
        blender_frame = scene.frame_current
        traj_frame    = min(blender_frame, n_traj_frames - 1)
        frame_coords  = coords[traj_frame]          # (natoms, 3) in Å

        # ── Compute neighbour mask ────────────────────────────────────────────
        if use_filter:
            # Positions of all centre atoms this frame
            centre_pos = frame_coords[centre_indices]   # (n_centre, 3)

            # Distance from every atom to every centre atom; keep atom if any
            # centre is within cutoff.  Broadcasting: (natoms,3) vs (n_centre,3)
            diff      = frame_coords[:, None, :] - centre_pos[None, :, :]  # (N,C,3)
            dist_sq   = np.einsum("nci,nci->nc", diff, diff)               # (N,C)
            within    = np.any(dist_sq <= NEIGHBOR_CUTOFF ** 2, axis=1)    # (N,)
        else:
            within = np.ones(len(atom_order), dtype=bool)

        # ── Update each element's instancer mesh ──────────────────────────────
        for elem, obj in instancer_objects.items():
            if elem not in element_indices_np:
                continue

            idxs         = element_indices_np[elem]
            visible_mask = within[idxs]
            visible_pos  = frame_coords[idxs[visible_mask]] * SCALE  # (n_vis, 3)

            mesh = obj.data
            mesh.clear_geometry()
            if len(visible_pos) > 0:
                mesh.vertices.add(len(visible_pos))
                mesh.vertices.foreach_set("co", visible_pos.ravel())
            mesh.update()

    # Remove any previously registered handler from this script
    handlers = bpy.app.handlers.frame_change_pre
    handlers[:] = [h for h in handlers
                   if getattr(h, "__name__", "") != "mol_traj_update"]
    update_frame.__name__ = "mol_traj_update"
    handlers.append(update_frame)

    # Trigger immediately so frame 0 shows correctly
    update_frame(bpy.context.scene)
    print(f"  Animation handler registered ({n_traj_frames} frames available).")


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    # 1. Clear the scene
    bpy.ops.object.select_all(action="SELECT")
    bpy.ops.object.delete(use_global=False)
    for block in list(bpy.data.meshes):
        bpy.data.meshes.remove(block)
    for block in list(bpy.data.materials):
        bpy.data.materials.remove(block)

    # 2. Parse PDB
    pdb_path = resolve_path(PDB_FILE, "structure.pdb")
    if not pdb_path:
        raise FileNotFoundError(
            "Could not locate structure.pdb. "
            "Save your .blend next to it, or set PDB_FILE above."
        )
    print(f"PDB: {pdb_path}")
    atoms_by_element, atom_order = parse_pdb(pdb_path)
    for elem, pos in atoms_by_element.items():
        print(f"  {elem:4s}: {len(pos):5d} atoms")

    # 3. Create collection
    col_name = "Molecule"
    if col_name in bpy.data.collections:
        bpy.data.collections.remove(bpy.data.collections[col_name])
    mol_col = bpy.data.collections.new(col_name)
    bpy.context.scene.collection.children.link(mol_col)

    # 4. Build instancer objects for each element
    instancer_objects = {}
    for element, positions in sorted(atoms_by_element.items()):
        parent_obj, _ = create_element_objects(element, positions, mol_col)
        instancer_objects[element] = parent_obj
        print(f"  Created instancer: {element} ({len(positions)} atoms)")

    # 5. Camera + empty setup
    all_positions = [p for positions in atoms_by_element.values() for p in positions]
    cx = sum(p[0] for p in all_positions) / len(all_positions) * SCALE
    cy = sum(p[1] for p in all_positions) / len(all_positions) * SCALE
    cz = sum(p[2] for p in all_positions) / len(all_positions) * SCALE

    target = bpy.data.objects.new("AtomCentre", None)
    target.empty_display_type = "SPHERE"
    target.empty_display_size = 0.3
    target.location = (cx, cy, cz)
    mol_col.objects.link(target)

    offset = max(cx, cy, cz) * 4.5
    cam_data = bpy.data.cameras.new("MolCamera")
    cam = bpy.data.objects.new("MolCamera", cam_data)
    cam.location = (cx + offset, cy - offset, cz + offset * 0.7)
    mol_col.objects.link(cam)
    bpy.context.scene.camera = cam

    constraint = cam.constraints.new(type="TRACK_TO")
    constraint.target     = target
    constraint.track_axis = "TRACK_NEGATIVE_Z"
    constraint.up_axis    = "UP_Y"

    sun_data = bpy.data.lights.new("MolSun", type="SUN")
    sun_data.energy = 3.0
    sun = bpy.data.objects.new("MolSun", sun_data)
    sun.location = (cx, cy, cz + offset)
    mol_col.objects.link(sun)

    # 6. Load trajectory and register animation handler (if DCD is available)
    traj_path = resolve_path(TRAJ_FILE, "traj_nvt.dcd")
    if traj_path:
        print(f"DCD:  {traj_path}")
        coords = parse_dcd(traj_path, frame_step=FRAME_STEP)
        n_frames = coords.shape[0]

        # Set Blender timeline to match the number of sampled frames
        bpy.context.scene.frame_start = 0
        bpy.context.scene.frame_end   = n_frames - 1
        bpy.context.scene.frame_current = 0

        register_trajectory_handler(coords, atom_order, instancer_objects)
        print(f"  Timeline set to 0 – {n_frames - 1} frames.")
    else:
        print("  No DCD file found — showing static structure only.")

    # 7. Viewport shading
    for area in bpy.context.screen.areas:
        if area.type == "VIEW_3D":
            for space in area.spaces:
                if space.type == "VIEW_3D":
                    space.shading.type = "MATERIAL"

    print("Done!  Press Space to play.  Z → Rendered for full quality.")


if __name__ == "__main__":
    main()
