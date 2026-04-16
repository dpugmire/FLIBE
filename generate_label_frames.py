#!/usr/bin/env python3
"""
Generate transparent label overlay frames for storyboard compositing.

Creates:
  - Always-visible right-side labels: T, Li, Be, F (atom colors)
  - Bottom text "Single Cluster" from a chosen storyboard boundary onward

The script reads timing knobs and ATOM_COLORS from the storyboard script so
the overlay sequence stays synchronized with your render frames.
"""

from __future__ import annotations

import argparse
import ast
from pathlib import Path
from typing import Dict, Tuple

try:
    from PIL import Image, ImageDraw, ImageFont
except Exception as exc:  # pragma: no cover
    raise SystemExit(
        "Pillow is required. Install with: python3 -m pip install pillow"
    ) from exc


def _load_top_level_literals(script_path: Path) -> Dict[str, object]:
    src = script_path.read_text(encoding="utf-8")
    tree = ast.parse(src, filename=str(script_path))
    out: Dict[str, object] = {}

    for node in tree.body:
        target_name = None
        if isinstance(node, ast.Assign) and len(node.targets) == 1 and isinstance(node.targets[0], ast.Name):
            target_name = node.targets[0].id
            value_node = node.value
        elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name) and node.value is not None:
            target_name = node.target.id
            value_node = node.value
        else:
            continue

        try:
            out[target_name] = ast.literal_eval(value_node)
        except Exception:
            # Skip non-literal expressions.
            pass

    return out


def _to_rgba255(color) -> Tuple[int, int, int, int]:
    if not isinstance(color, (tuple, list)) or len(color) < 3:
        return (255, 255, 255, 255)
    r = int(round(max(0.0, min(1.0, float(color[0]))) * 255))
    g = int(round(max(0.0, min(1.0, float(color[1]))) * 255))
    b = int(round(max(0.0, min(1.0, float(color[2]))) * 255))
    if len(color) >= 4:
        a = int(round(max(0.0, min(1.0, float(color[3]))) * 255))
    else:
        a = 255
    return (r, g, b, a)


def _compute_story_boundaries_dcd(vals: Dict[str, object]) -> Dict[str, int]:
    start = int(vals["STORY_DCD_START"])
    end = start + int(vals["STORY_DCD_STEPS"])
    pre = int(vals["PRE_ZOOM_ALIGN_DCD_FRAMES"])
    phase1 = int(vals["PHASE1_ZOOM_TO_T_DCD_FRAMES"])
    phase2 = int(vals["PHASE2_T_FOCUS_HOLD_DCD_FRAMES"])
    fade_non = int(vals["FADE_NON_CLUSTER_OUT_DCD_FRAMES"])
    fade_in = int(vals["FADE_SPHERE_IN_DCD_FRAMES"])
    fade_out = int(vals["FADE_SPHERE_OUT_DCD_FRAMES"])

    pre_zoom_align_end = start + pre
    zoom_to_t_end = pre_zoom_align_end + phase1
    t_focus_hold_end = zoom_to_t_end + phase2
    transition_start = t_focus_hold_end
    transition_end = transition_start + max(fade_non, fade_in)
    with_sphere_hold_end = transition_start + fade_in
    sphere_fade_out_end = with_sphere_hold_end + fade_out

    return {
        "start": start,
        "pre_zoom_align_end": pre_zoom_align_end,
        "zoom_to_t_end": zoom_to_t_end,
        "t_focus_hold_end": t_focus_hold_end,
        "transition_start": transition_start,
        "transition_end": transition_end,
        "with_sphere_hold_end": with_sphere_hold_end,
        "sphere_fade_out_end": sphere_fade_out_end,
        "end": end,
    }


def _dcd_to_bl(dcd_frame: int, dcd_offset: int) -> int:
    return int(dcd_frame) - int(dcd_offset) + 1


def _load_font(font_size: int, font_path: str | None):
    if font_path:
        return ImageFont.truetype(font_path, font_size)
    for candidate in (
        "/System/Library/Fonts/Supplemental/Arial Bold.ttf",
        "/Library/Fonts/Arial Bold.ttf",
        "/System/Library/Fonts/Supplemental/Helvetica.ttc",
        "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
    ):
        p = Path(candidate)
        if p.exists():
            try:
                return ImageFont.truetype(str(p), font_size)
            except Exception:
                pass
    return ImageFont.load_default()


def _draw_overlay(
    width: int,
    height: int,
    labels: Tuple[Tuple[str, Tuple[int, int, int, int]], ...],
    show_cluster_text: bool,
    cluster_text: str,
    font_label,
    font_cluster,
) -> Image.Image:
    img = Image.new("RGBA", (width, height), (0, 0, 0, 0))
    draw = ImageDraw.Draw(img)

    scale = min(width, height) / 1080.0
    right_margin = int(round(70 * scale))
    top_margin = int(round(80 * scale))
    line_spacing = int(round(62 * scale))
    stroke_w = max(1, int(round(3 * scale)))

    for i, (text, rgba) in enumerate(labels):
        y = top_margin + i * line_spacing
        bbox = draw.textbbox((0, 0), text, font=font_label, stroke_width=stroke_w)
        tw = bbox[2] - bbox[0]
        x = width - right_margin - tw
        draw.text(
            (x, y),
            text,
            fill=rgba,
            font=font_label,
            stroke_width=stroke_w,
            stroke_fill=(0, 0, 0, 230),
        )

    if show_cluster_text:
        stroke_w2 = max(1, int(round(4 * scale)))
        bottom_margin = int(round(60 * scale))
        bbox = draw.textbbox((0, 0), cluster_text, font=font_cluster, stroke_width=stroke_w2)
        tw = bbox[2] - bbox[0]
        th = bbox[3] - bbox[1]
        x = (width - tw) // 2
        y = height - bottom_margin - th
        draw.text(
            (x, y),
            cluster_text,
            fill=(255, 255, 255, 255),
            font=font_cluster,
            stroke_width=stroke_w2,
            stroke_fill=(0, 0, 0, 230),
        )

    return img


def main() -> None:
    repo = Path(__file__).resolve().parent
    default_story = repo / "visualize_molecule_storyboard_tritium_cluster.py"

    ap = argparse.ArgumentParser(description="Generate transparent label overlay frames.")
    ap.add_argument("--story-script", default=str(default_story), help="Path to storyboard .py script")
    ap.add_argument("--width", type=int, default=1080, help="Overlay width (px)")
    ap.add_argument("--height", type=int, default=1080, help="Overlay height (px)")
    ap.add_argument("--output-dir", default=str(repo / "output_tritium_labels"), help="Output folder")
    ap.add_argument("--prefix", default="label_", help="Filename prefix")
    ap.add_argument("--digits", type=int, default=5, help="Frame number zero padding")
    ap.add_argument(
        "--cluster-label-start",
        choices=["transition_start", "transition_end", "with_sphere_hold_end"],
        default="transition_start",
        help="Storyboard boundary where 'Single Cluster' begins",
    )
    ap.add_argument(
        "--cluster-text",
        default="Single Cluster",
        help="Bottom label text",
    )
    ap.add_argument("--font-path", default="", help="Optional TTF/OTF font path")
    args = ap.parse_args()

    story_script = Path(args.story_script).resolve()
    if not story_script.exists():
        raise SystemExit(f"Storyboard script not found: {story_script}")

    vals = _load_top_level_literals(story_script)
    required = [
        "STORY_DCD_START",
        "STORY_DCD_STEPS",
        "PRE_ZOOM_ALIGN_DCD_FRAMES",
        "PHASE1_ZOOM_TO_T_DCD_FRAMES",
        "PHASE2_T_FOCUS_HOLD_DCD_FRAMES",
        "FADE_NON_CLUSTER_OUT_DCD_FRAMES",
        "FADE_SPHERE_IN_DCD_FRAMES",
        "FADE_SPHERE_OUT_DCD_FRAMES",
        "ATOM_COLORS",
    ]
    missing = [k for k in required if k not in vals]
    if missing:
        raise SystemExit(f"Missing constants in storyboard script: {', '.join(missing)}")

    # In the storyboard script this is commonly defined as:
    #   DCD_FRAME_OFFSET = STORY_DCD_START
    # which is not a literal and is skipped by ast.literal_eval. Fall back to
    # STORY_DCD_START in that case.
    vals.setdefault("DCD_FRAME_OFFSET", vals["STORY_DCD_START"])

    boundaries_dcd = _compute_story_boundaries_dcd(vals)
    dcd_offset = int(vals["DCD_FRAME_OFFSET"])
    boundaries_bl = {k: _dcd_to_bl(v, dcd_offset) for k, v in boundaries_dcd.items()}

    start_bl = boundaries_bl["start"]
    end_bl = boundaries_bl["end"]
    cluster_start_bl = boundaries_bl[args.cluster_label_start]

    atom_colors = vals["ATOM_COLORS"]
    labels = (
        ("T", _to_rgba255(atom_colors.get("H", (1.0, 0.0, 0.0, 1.0)))),
        ("Li", _to_rgba255(atom_colors.get("Li", (0.0, 0.0, 1.0, 1.0)))),
        ("Be", _to_rgba255(atom_colors.get("Be", (1.0, 1.0, 0.0, 1.0)))),
        ("F", _to_rgba255(atom_colors.get("F", (0.0, 1.0, 0.0, 1.0)))),
    )

    scale = min(args.width, args.height) / 1080.0
    font_label = _load_font(max(18, int(round(52 * scale))), args.font_path or None)
    font_cluster = _load_font(max(18, int(round(64 * scale))), args.font_path or None)

    out_dir = Path(args.output_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    for bf in range(start_bl, end_bl + 1):
        img = _draw_overlay(
            width=int(args.width),
            height=int(args.height),
            labels=labels,
            show_cluster_text=(bf >= cluster_start_bl),
            cluster_text=args.cluster_text,
            font_label=font_label,
            font_cluster=font_cluster,
        )
        filename = f"{args.prefix}{bf:0{int(args.digits)}d}.png"
        img.save(out_dir / filename)

    print(f"Storyboard script: {story_script}")
    print(f"Frames generated: {start_bl}..{end_bl} ({end_bl - start_bl + 1} frames)")
    print(f"'Single Cluster' starts at frame: {cluster_start_bl} ({args.cluster_label_start})")
    print(f"Output: {out_dir}")


if __name__ == "__main__":
    main()
