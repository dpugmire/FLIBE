"""
register_handler.py
-------------------
Pass this to Blender with -P when batch-rendering to re-register the
trajectory animation handler.  The .blend file must have been saved after
running visualize_molecule.py (so the scene objects already exist).

Usage:
  blender -b setup.blend -P register_handler.py -o ./frames/frame.#### -s 0 -e 815 -j 10 -a

Edit the settings below to match what you used in visualize_molecule.py.
"""

import bpy
import os
import sys
import struct
import numpy as np

# ── Settings — keep these in sync with visualize_molecule.py ─────────────────

TRAJ_FILE        = ""       # leave "" to auto-detect traj_nvt.dcd beside the .blend
FRAME_STEP       = 10       # must match what was used when the .blend was built
NEIGHBOR_ELEMENT = "H"
NEIGHBOR_CUTOFF  = 6.0      # Å  — set to 0 to show all atoms
SCALE            = 0.1      # Ångströms → Blender units

# ─────────────────────────────────────────────────────────────────────────────

def resolve_path(user_path, filename):
    if user_path:
        return os.path.abspath(user_path)
    blend_dir = bpy.path.abspath("//")
    candidate = os.path.join(blend_dir, filename)
    return candidate if os.path.isfile(candidate) else None


def parse_dcd(filepath, frame_step=1):
    def read_record(f):
        size = struct.unpack("i", f.read(4))[0]
        data = f.read(size); f.read(4)
        return data, size

    with open(filepath, "rb") as f:
        hdr, _        = read_record(f)
        n_frames      = struct.unpack_from("i", hdr, 4)[0]
        read_record(f)                                      # title
        natoms_data,_ = read_record(f)
        natoms        = struct.unpack("i", natoms_data)[0]
        print(f"  DCD: {n_frames} frames, {natoms} atoms")

        frame_start   = f.tell()
        first_size    = struct.unpack("i", f.read(4))[0]
        f.seek(frame_start)
        has_cell      = (first_size == 48)

        sampled       = range(0, n_frames, frame_step)
        coords        = np.zeros((len(sampled), natoms, 3), dtype=np.float32)
        frame_bytes   = (8 + 48 if has_cell else 0) + 3 * (8 + natoms * 4)

        for out_i, traj_i in enumerate(sampled):
            f.seek(frame_start + traj_i * frame_bytes)
            if has_cell:
                read_record(f)
            xd,_ = read_record(f)
            yd,_ = read_record(f)
            zd,_ = read_record(f)
            coords[out_i,:,0] = np.frombuffer(xd, dtype=np.float32)
            coords[out_i,:,1] = np.frombuffer(yd, dtype=np.float32)
            coords[out_i,:,2] = np.frombuffer(zd, dtype=np.float32)
            if out_i % 200 == 0:
                print(f"  DCD: loaded {out_i}/{len(sampled)}", end="\r")

    print(f"  DCD: {len(sampled)} frames loaded (step={frame_step})")
    return coords, natoms


def get_atom_order_from_scene():
    """
    Reconstruct atom_order list from the instancer objects already in the scene.
    Each instancer mesh named 'F', 'Li', 'Be', etc. has N vertices = N atoms.
    We sort instancer objects by their first vertex index stored as a custom
    attribute — but since we don't have that, we fall back to reading the PDB.
    """
    # Simpler: read the PDB again (it's tiny)
    pdb_path = resolve_path("", "structure.pdb")
    if not pdb_path:
        raise FileNotFoundError("Cannot find structure.pdb beside the .blend file.")
    atom_order = []
    with open(pdb_path) as fh:
        for line in fh:
            rec = line[:6].strip()
            if rec not in ("ATOM", "HETATM"):
                continue
            element = line[76:78].strip() if len(line) > 76 else ""
            if not element:
                element = ''.join(c for c in line[12:16].strip() if c.isalpha())
            atom_order.append(element.capitalize())
    return atom_order


def register(coords, atom_order):
    n_traj_frames = coords.shape[0]

    element_indices = {}
    for i, elem in enumerate(atom_order):
        element_indices.setdefault(elem, []).append(i)
    elem_idx_np = {e: np.array(v) for e, v in element_indices.items()}

    centre_indices = elem_idx_np.get(NEIGHBOR_ELEMENT)
    use_filter = (NEIGHBOR_CUTOFF > 0) and (centre_indices is not None)

    # Collect instancer objects from the scene (named by element symbol)
    instancer_objects = {
        obj.name: obj
        for obj in bpy.data.objects
        if obj.type == "MESH" and obj.instance_type == "VERTS"
    }
    print(f"  Found instancer objects: {list(instancer_objects.keys())}")

    def update_frame(scene, depsgraph=None):
        traj_frame   = min(scene.frame_current, n_traj_frames - 1)
        frame_coords = coords[traj_frame]

        if use_filter:
            centre_pos = frame_coords[centre_indices]
            diff       = frame_coords[:, None, :] - centre_pos[None, :, :]
            dist_sq    = np.einsum("nci,nci->nc", diff, diff)
            within     = np.any(dist_sq <= NEIGHBOR_CUTOFF ** 2, axis=1)
        else:
            within = np.ones(len(atom_order), dtype=bool)

        for elem, obj in instancer_objects.items():
            if elem not in elem_idx_np:
                continue
            idxs        = elem_idx_np[elem]
            visible_pos = frame_coords[idxs[within[idxs]]] * SCALE
            mesh        = obj.data
            mesh.clear_geometry()
            if len(visible_pos) > 0:
                mesh.vertices.add(len(visible_pos))
                mesh.vertices.foreach_set("co", visible_pos.ravel())
            mesh.update()

    handlers = bpy.app.handlers.frame_change_pre
    handlers[:] = [h for h in handlers
                   if getattr(h, "__name__", "") != "mol_traj_update"]
    update_frame.__name__ = "mol_traj_update"
    handlers.append(update_frame)

    # Also hook render_pre so the handler fires during rendering
    render_handlers = bpy.app.handlers.render_pre
    render_handlers[:] = [h for h in render_handlers
                          if getattr(h, "__name__", "") != "mol_traj_update"]
    render_handlers.append(update_frame)

    update_frame(bpy.context.scene)
    print(f"  Handler registered. {n_traj_frames} frames available.")


def main():
    traj_path = resolve_path(TRAJ_FILE, "traj_nvt.dcd")
    if not traj_path:
        print("ERROR: cannot find traj_nvt.dcd — set TRAJ_FILE above.")
        sys.exit(1)

    print(f"DCD: {traj_path}")
    coords, natoms = parse_dcd(traj_path, frame_step=FRAME_STEP)
    atom_order     = get_atom_order_from_scene()

    if len(atom_order) != natoms:
        print(f"WARNING: PDB has {len(atom_order)} atoms but DCD has {natoms}.")

    register(coords, atom_order)
    print("Done — starting render.")


main()
