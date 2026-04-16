# flibe

Blender Python scripts for molecular visualization and storyboard rendering.

## Scripts

- `visualize_molecule.py`
- `visualize_molecule_storyboard.py`
- `visualize_molecule_storyboard_tritium_cluster.py`
- `register_handler.py`

## Notes

Large trajectory/data/render artifacts are excluded from version control via `.gitignore`.

## Running `visualize_molecule_storyboard_tritium_cluster.py`

`visualize_molecule_storyboard_tritium_cluster.py` is a dedicated storyboard for:

1. all atoms + zoom toward tritium
2. tritium-focused hold
3. fade out non-cluster atoms while sphere fades in
4. cluster + tritium hold with semi-transparent sphere
5. sphere fade-out
6. cluster + tritium tail

Key defaults and knobs are in the script header:

- Trajectory window: `STORY_DCD_START`, `STORY_DCD_STEPS`
- Phase timing: `PHASE1_ZOOM_TO_T_DCD_FRAMES`, `PHASE2_T_FOCUS_HOLD_DCD_FRAMES`,
  `PHASE4_WITH_SPHERE_HOLD_DCD_FRAMES`
- Fade timing: `FADE_NON_CLUSTER_OUT_DCD_FRAMES`, `FADE_SPHERE_IN_DCD_FRAMES`,
  `FADE_SPHERE_OUT_DCD_FRAMES`
- Camera distances: `CAM_DIST_ALL_ATOMS`, `CAM_DIST_TRITIUM`, `CAM_DIST_SPHERE`, `CAM_DIST_CLUSTER`
- Sphere style: `SPHERE_RADIUS_ANGSTROM`, `SPHERE_COLOR`, `SPHERE_FINAL_ALPHA`
- Boundary mode: `BOUNDARY_MODE = "sphere" | "none"`

Notes:

- `TRITIUM_IDX = 1680`
- Cluster atoms come from `FOCUS_INDICES`
- From phase 3 onward, the reference sphere center follows the live cluster centroid
  (mean of `FOCUS_INDICES` atom positions) via `AtomCentre`
- If `TRAJ_FILE` is blank, the script looks for `traj_nvt_interp.dcd`

Minimal run:

```bash
blender --python /Users/dpn/proj/flibe/OneDrive_1_4-9-2026/FLIBE/visualize_molecule_storyboard_tritium_cluster.py
```

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

### Transition fade controls (zoom to F,F,H and cluster)

In `visualize_molecule_storyboard.py`, these knobs control the zoom/visibility transition:

```python
ZOOM_START_DCD = SEG2_DCD_START
HIDE_OTHERS_DCD_START = SEG3_DCD_START
FADE_OTHERS_FRAMES = 30
FADE_CLUSTER_IN_FRAMES = 30
```

- `HIDE_OTHERS_DCD_START`: frame where non-focus atoms begin to fade out.
- `FADE_OTHERS_FRAMES`: number of Blender frames for non-focus fade-out (`0` = instant hide).
- `SEG5_DCD_START`: cluster phase begins.
- `FADE_CLUSTER_IN_FRAMES`: number of Blender frames for cluster fade-in during seg 5 (`0` = instant pop-in).

The focus trio (tritium + `F_1` + `F_2`) stays visible while the box fades out and while the cluster blends in.

### Command-line render (full movie)

```bash
cd /Users/dpn/proj/flibe/OneDrive_1_4-9-2026/FLIBE
mkdir -p output

/Applications/Blender.app/Contents/MacOS/Blender \
  -b --factory-startup \
  -P /Users/dpn/proj/flibe/OneDrive_1_4-9-2026/FLIBE/visualize_molecule_storyboard.py \
  --python-expr "import bpy; bpy.context.scene.render.resolution_percentage=50" \
  -o /Users/dpn/proj/flibe/OneDrive_1_4-9-2026/FLIBE/output/frame_##### \
  -F PNG -x 1 \
  -a
```

This renders the full storyboard range set in the script (`SEG1_DCD_START` to `SEG5_DCD_END`).

### Command-line render (subrange / frame skip)

Add `-s`, `-e`, and `-j` before `-a`:

```bash
blender \
  -b --factory-startup \
  -P /Users/dpn/proj/flibe/OneDrive_1_4-9-2026/FLIBE/visualize_molecule_storyboard.py \
  --python-expr "import bpy; bpy.context.scene.render.resolution_percentage=50" \
  -o /Users/dpn/proj/flibe/OneDrive_1_4-9-2026/nvt_630_TF_large/output/frame_##### \
  -F PNG -x 1 \
  -a
```

`-j 2` renders every other frame.

### Minimal command-line example

```bash
blender \
  -b --factory-startup \
  -P /Users/dpn/proj/flibe/OneDrive_1_4-9-2026/FLIBE/visualize_molecule_storyboard.py \
  -a
```

## Transparent label overlays

Use `make_label_overlays.py` to generate transparent PNG overlays with:

- Right-side labels in order: `F`, `T`, `Be`, `Li`
- Label colors matched to atom colors
- `Be`/`Li` fade out during zoom focus transition
- `Be`/`Li` fade back in during cluster fade-in
- Bottom-center white text: `Single Cluster` (fades in with cluster)

### Generate overlay PNG sequence

Install Pillow once:

```bash
python3 -m pip install pillow
```

Generate overlays for rendered frames:

```bash
python3 /Users/dpn/proj/flibe/OneDrive_1_4-9-2026/FLIBE/make_label_overlays.py \
  --frames-dir /Users/dpn/proj/flibe/OneDrive_1_4-9-2026/nvt_630_TF_large/output \
  --output-dir /Users/dpn/proj/flibe/OneDrive_1_4-9-2026/nvt_630_TF_large/output/labels \
  --storyboard-script /Users/dpn/proj/flibe/OneDrive_1_4-9-2026/FLIBE/visualize_molecule_storyboard.py
```

By default, the overlay script auto-loads timing constants from
`visualize_molecule_storyboard.py` and prints the resolved fade frame ranges.
You can still override any timing directly from CLI.

Timing constants used:

- `HIDE_OTHERS_DCD_START = 1000`
- `FADE_OTHERS_FRAMES = 30`
- `SEG5_DCD_START = 1350`
- `FADE_CLUSTER_IN_FRAMES = 30`
- `DCD_FRAME_OFFSET = 700`

You can override timings with:
`--dcd-offset`, `--hide-others-dcd-start`, `--fade-others-frames`,
`--seg5-dcd-start`, and `--fade-cluster-in-frames`.

If your frame filenames were renumbered and no longer match Blender frame numbers,
use `--frame-index-offset` so label fades line up with atom fades.

### Composite overlays over rendered frames

```bash
ffmpeg \
  -framerate 30 -start_number 1 \
  -i /Users/dpn/proj/flibe/OneDrive_1_4-9-2026/nvt_630_TF_large/output/frame_%05d.png \
  -framerate 30 -start_number 1 \
  -i /Users/dpn/proj/flibe/OneDrive_1_4-9-2026/nvt_630_TF_large/output/labels/labels_%05d.png \
  -filter_complex "[0:v][1:v]overlay=0:0:format=auto:shortest=1" \
  -c:v libx264 -pix_fmt yuv420p \
  /Users/dpn/proj/flibe/OneDrive_1_4-9-2026/nvt_630_TF_large/output/movie_with_labels.mp4
```
