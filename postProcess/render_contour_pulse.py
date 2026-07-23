#!/usr/bin/env python3
"""
# Render Drop-Injection Contour Pulse Files

Render the lightweight `interface-latest.dat` files written by contour-mode
simulations. The renderer needs only Matplotlib; it does not restore Basilisk
snapshots or compile helpers, so scheduled monitors can run safely away from
the cluster.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def parse_meta(path: Path) -> dict[str, str]:
    """Read the simple key-value pulse metadata file."""
    values: dict[str, str] = {}
    if not path.exists():
        return values
    for line in path.read_text().splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            values[key.strip()] = value.strip()
    return values


def parse_classification(path: Path) -> dict[str, str]:
    """Read a case's latest classification row when available."""
    if not path.exists():
        return {}
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    return rows[-1] if rows else {}


def read_segments(path: Path) -> list[tuple[tuple[float, float], tuple[float, float]]]:
    """Parse Basilisk `output_facets()` pairs separated by blank lines."""
    segments: list[tuple[tuple[float, float], tuple[float, float]]] = []
    points: list[tuple[float, float]] = []
    for raw in path.read_text().splitlines() + [""]:
        line = raw.strip()
        if line:
            fields = line.split()
            if len(fields) >= 2:
                points.append((float(fields[0]), float(fields[1])))
        elif len(points) >= 2:
            segments.append((points[0], points[1]))
            points = []
        else:
            points = []
    return segments


def render_case(case_dir: Path, output_dir: Path) -> dict[str, str]:
    """Render one staged case and return a manifest row."""
    facets = case_dir / "interface-latest.dat"
    meta = parse_meta(case_dir / "interface-latest.meta")
    classification = parse_classification(case_dir / "classification.status")
    segments = read_segments(facets)
    if not segments:
        raise ValueError(f"No interface segments in {facets}")

    output_dir.mkdir(parents=True, exist_ok=True)
    figure, axis = plt.subplots(figsize=(7.2, 5.4), constrained_layout=True)
    for (x0, y0), (x1, y1) in segments:
        axis.plot((x0, x1), (y0, y1), color="#8b1e3f", linewidth=0.8)
        if max(y0, y1) > 1e-12:
            axis.plot((x0, x1), (-y0, -y1), color="#8b1e3f", linewidth=0.8)

    xs = [x for segment in segments for point in segment for x in (point[0],)]
    ys = [abs(y) for segment in segments for point in segment for y in (point[1],)]
    xmin, xmax = min(xs), max(xs)
    ymax = max(ys)
    xpad = max(0.25, 0.05 * (xmax - xmin))
    ypad = max(0.25, 0.08 * ymax)
    axis.set_xlim(xmin - xpad, xmax + xpad)
    axis.set_ylim(-ymax - ypad, ymax + ypad)
    axis.set_aspect("equal", adjustable="box")
    axis.set_xlabel("axial coordinate")
    axis.set_ylabel("radial coordinate")
    axis.grid(alpha=0.18, linewidth=0.5)

    case_name = case_dir.name
    status_id = classification.get("id", "-1")
    status_text = {"1": "drop", "0": "no drop", "-1": "running"}.get(
        status_id, "unknown"
    )
    title = (
        f"{case_name}: Rr={meta.get('Rr', '?')}, Oh={meta.get('Oh', '?')}, "
        f"t={meta.get('t', '?')} [{status_text}]"
    )
    axis.set_title(title)

    destination = output_dir / f"{case_name}.png"
    temporary = destination.with_suffix(".png.tmp")
    figure.savefig(temporary, dpi=180, format="png")
    plt.close(figure)
    temporary.replace(destination)
    return {
        "case": case_name,
        "Rr": meta.get("Rr", ""),
        "Oh": meta.get("Oh", ""),
        "t": meta.get("t", ""),
        "id": status_id,
        "reason": classification.get("reason", ""),
        "image": destination.name,
    }


def main() -> int:
    """Render every staged case with an available pulse file."""
    parser = argparse.ArgumentParser()
    parser.add_argument("staged_root", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()

    case_dirs = sorted(
        path.parent for path in args.staged_root.glob("*/interface-latest.dat")
    )
    if not case_dirs:
        parser.error(f"no pulse files below {args.staged_root}")

    rows = []
    failures = []
    for case_dir in case_dirs:
        try:
            rows.append(render_case(case_dir, args.output_dir))
        except Exception as error:  # keep the other 15 renders useful
            failures.append((case_dir.name, str(error)))

    manifest = args.output_dir / "render-manifest.csv"
    with manifest.open("w", newline="") as stream:
        writer = csv.DictWriter(
            stream, fieldnames=("case", "Rr", "Oh", "t", "id", "reason", "image")
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"rendered={len(rows)} failed={len(failures)} manifest={manifest}")
    for case, error in failures:
        print(f"failed {case}: {error}")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
