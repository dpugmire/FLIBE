# flibe

Blender Python scripts for molecular visualization and storyboard rendering.

## Scripts

- `visualize_molecule.py`
- `visualize_molecule_storyboard.py`
- `register_handler.py`

## Notes

Large trajectory/data/render artifacts are excluded from version control via `.gitignore`.

## Running `visualize_molecule_storyboard.py`

The storyboard script supports two setup modes:

1. Script-only mode (no setup `.blend`)
2. Setup-file mode (reuse camera/empty/lights from a `.blend`)

### 1) Script-only mode (auto rig)

In `visualize_molecule_storyboard.py` set:

```python
SETUP_BLEND_FILE = ""
```

Then run the script in Blender. It will:

- Create camera + `AtomCentre` empty + camera TrackTo
- Create default 3-point area lights (all TrackTo `AtomCentre`)
- Build atom geometry and animation

### 2) Setup-file mode (reuse your rig)

In `visualize_molecule_storyboard.py` set:

```python
SETUP_BLEND_FILE = "//setup.blend"   # or absolute path
SETUP_CAMERA_NAME = ""               # optional explicit camera object name
SETUP_EMPTY_NAME  = ""               # optional explicit empty object name
```

Then run the script. It will append camera/empty/lights from that setup file and reuse them.
If no lights are found in the setup file, the script falls back to auto 3-point lights.

### Command-line example

```bash
/Applications/Blender.app/Contents/MacOS/Blender \
  -b --factory-startup \
  -P /Users/dpn/proj/flibe/OneDrive_1_4-9-2026/FLIBE/visualize_molecule_storyboard.py \
  -a
```
