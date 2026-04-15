#!/usr/bin/env python3
"""
Create transparent PNG overlay frames with atom labels for compositing.

Behavior:
- Right-side labels in order: F, T, Be, Li
- F/T stay fully visible
- Be/Li fade with explicit output-frame ranges:
  --fade-out f0 f1
  --fade-in  f0 f1
- "Single Cluster" text fades in using --fade-in
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

try:
    from PIL import Image, ImageDraw, ImageFont
except Exception:
    print("ERROR: Pillow is required. Install with: python3 -m pip install pillow")
    raise


LABELS = [
    ("F", (0.565, 0.878, 0.314)),
    ("T", (1.000, 0.000, 0.000)),
    ("Be", (0.765, 1.000, 0.000)),
    ("Li", (0.784, 0.502, 1.000)),
]

FONT_CANDIDATES = [
    "/System/Library/Fonts/Supplemental/Arial Bold.ttf",
    "/System/Library/Fonts/SFNS.ttf",
    "/Library/Fonts/Arial Bold.ttf",
]


def clamp01(x: float) -> float:
    return max(0.0, min(1.0, x))


def to_rgba(rgbf: tuple[float, float, float], alpha: float) -> tuple[int, int, int, int]:
    a = int(round(clamp01(alpha) * 255.0))
    return (
        int(round(clamp01(rgbf[0]) * 255.0)),
        int(round(clamp01(rgbf[1]) * 255.0)),
        int(round(clamp01(rgbf[2]) * 255.0)),
        a,
    )


def load_font(size: int, font_path: str | None) -> ImageFont.FreeTypeFont | ImageFont.ImageFont:
    if font_path:
        return ImageFont.truetype(font_path, size=size)
    for candidate in FONT_CANDIDATES:
        p = Path(candidate)
        if p.exists():
            try:
                return ImageFont.truetype(str(p), size=size)
            except Exception:
                pass
    return ImageFont.load_default()


def parse_frame_index(path: Path) -> int | None:
    m = re.search(r"(\d+)(?=\.[^.]+$)", path.name)
    if not m:
        return None
    return int(m.group(1))


def discover_frames(frames_dir: Path, pattern: str) -> list[tuple[int, Path]]:
    out: list[tuple[int, Path]] = []
    for p in frames_dir.glob(pattern):
        idx = parse_frame_index(p)
        if idx is not None:
            out.append((idx, p))
    out.sort(key=lambda t: t[0])
    return out


def ramp_alpha(frame: int, f0: int | None, f1: int | None, start_alpha: float, end_alpha: float) -> float | None:
    if f0 is None or f1 is None:
        return None
    if frame < f0:
        return start_alpha
    if f1 <= f0:
        return end_alpha if frame >= f0 else start_alpha
    if frame >= f1:
        return end_alpha
    t = (frame - f0) / float(f1 - f0)
    return (1.0 - t) * start_alpha + t * end_alpha


def compute_be_li_alpha(frame: int, fade_out: tuple[int, int] | None, fade_in: tuple[int, int] | None) -> float:
    # Default fully visible.
    alpha = 1.0

    # Apply fade-out first.
    if fade_out is not None:
        out_alpha = ramp_alpha(frame, fade_out[0], fade_out[1], 1.0, 0.0)
        if out_alpha is not None:
            alpha = out_alpha

    # Fade-in overrides from its start onward.
    if fade_in is not None and frame >= fade_in[0]:
        in_alpha = ramp_alpha(frame, fade_in[0], fade_in[1], 0.0, 1.0)
        if in_alpha is not None:
            alpha = in_alpha

    return clamp01(alpha)


def compute_cluster_alpha(frame: int, fade_in: tuple[int, int] | None) -> float:
    if fade_in is None:
        return 0.0
    a = ramp_alpha(frame, fade_in[0], fade_in[1], 0.0, 1.0)
    return 0.0 if a is None else clamp01(a)


def draw_right_label(
    draw: ImageDraw.ImageDraw,
    text: str,
    color_rgb: tuple[float, float, float],
    alpha: float,
    x_right: int,
    y_top: int,
    font: ImageFont.FreeTypeFont | ImageFont.ImageFont,
    stroke: int,
) -> None:
    fill = to_rgba(color_rgb, alpha)
    stroke_fill = (0, 0, 0, int(round(alpha * 220.0)))
    bbox = draw.textbbox((0, 0), text, font=font, stroke_width=stroke)
    w = bbox[2] - bbox[0]
    x = x_right - w
    draw.text(
        (x, y_top),
        text,
        font=font,
        fill=fill,
        stroke_width=stroke,
        stroke_fill=stroke_fill,
    )


def main() -> int:
    p = argparse.ArgumentParser(description="Generate transparent label overlay PNGs.")
    p.add_argument("--frames-dir", required=True, help="Directory containing source frames (e.g., frame_#####.png).")
    p.add_argument("--output-dir", default="", help="Output directory for transparent label frames.")
    p.add_argument("--pattern", default="frame_*.png", help="Glob pattern in --frames-dir.")
    p.add_argument("--start", type=int, default=0, help="First frame index (0 = auto).")
    p.add_argument("--end", type=int, default=0, help="Last frame index (0 = auto).")
    p.add_argument("--step", type=int, default=1, help="Frame step.")
    p.add_argument("--font", default="", help="Optional .ttf/.otf font path.")
    p.add_argument("--cluster-text", default="Single Cluster", help="Bottom-center text.")
    p.add_argument("--fade-out", nargs=2, type=int, metavar=("F0", "F1"),
                   help="Be/Li fade-out frame range (inclusive start, exclusive end).")
    p.add_argument("--fade-in", nargs=2, type=int, metavar=("F0", "F1"),
                   help="Be/Li + cluster text fade-in frame range (inclusive start, exclusive end).")
    args = p.parse_args()

    if args.step < 1:
        print("ERROR: --step must be >= 1")
        return 2

    fade_out = tuple(args.fade_out) if args.fade_out else None
    fade_in = tuple(args.fade_in) if args.fade_in else None

    frames_dir = Path(args.frames_dir).expanduser().resolve()
    if not frames_dir.exists():
        print(f"ERROR: frames dir does not exist: {frames_dir}")
        return 2

    found = discover_frames(frames_dir, args.pattern)
    if not found:
        print(f"ERROR: no frames found in {frames_dir} matching {args.pattern!r}")
        return 2

    first_idx = found[0][0]
    last_idx = found[-1][0]
    start = args.start if args.start > 0 else first_idx
    end = args.end if args.end > 0 else last_idx
    selected = [t for t in found if start <= t[0] <= end and ((t[0] - start) % args.step == 0)]
    if not selected:
        print("ERROR: no selected frames after start/end/step filtering")
        return 2

    output_dir = Path(args.output_dir).expanduser().resolve() if args.output_dir else (frames_dir / "labels")
    output_dir.mkdir(parents=True, exist_ok=True)

    with Image.open(selected[0][1]) as im0:
        width, height = im0.size

    label_font_size = max(20, int(round(height * 0.070)))
    cluster_font_size = max(18, int(round(height * 0.060)))
    label_font = load_font(label_font_size, args.font or None)
    cluster_font = load_font(cluster_font_size, args.font or None)

    x_right = width - int(round(width * 0.045))
    y_start = int(round(height * 0.18))
    line_h = int(round(label_font_size * 1.35))
    stroke = max(1, label_font_size // 18)
    cluster_stroke = max(1, cluster_font_size // 20)

    print("Using frame-index timing:")
    print(f"  fade-out = {fade_out}")
    print(f"  fade-in  = {fade_in}")

    total = len(selected)
    for i, (frame_idx, _) in enumerate(selected, start=1):
        overlay = Image.new("RGBA", (width, height), (0, 0, 0, 0))
        draw = ImageDraw.Draw(overlay)

        be_li_alpha = compute_be_li_alpha(frame_idx, fade_out, fade_in)
        cluster_alpha = compute_cluster_alpha(frame_idx, fade_in)

        # Required order: F, T, Be, Li
        for row, (name, rgb) in enumerate(LABELS):
            alpha = 1.0 if name in {"F", "T"} else be_li_alpha
            y = y_start + row * line_h
            draw_right_label(draw, name, rgb, alpha, x_right, y, label_font, stroke)

        if cluster_alpha > 0.0:
            text = args.cluster_text
            bbox = draw.textbbox((0, 0), text, font=cluster_font, stroke_width=cluster_stroke)
            tw = bbox[2] - bbox[0]
            tx = (width - tw) // 2
            ty = int(round(height * 0.90 - cluster_font_size))
            draw.text(
                (tx, ty),
                text,
                font=cluster_font,
                fill=(255, 255, 255, int(round(cluster_alpha * 255.0))),
                stroke_width=cluster_stroke,
                stroke_fill=(0, 0, 0, int(round(cluster_alpha * 220.0))),
            )

        out_name = f"labels_{frame_idx:05d}.png"
        overlay.save(output_dir / out_name)

        if i == 1 or i == total or i % 100 == 0:
            print(f"[{i}/{total}] wrote {out_name}")

    print(f"Done. Wrote {total} overlay frames to: {output_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
