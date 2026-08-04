#!/usr/bin/env python3
"""Collect retained case evidence into one contour-attempt result table."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

try:
    from result_quality import case_quality
except ModuleNotFoundError:
    from contourWorkflow.result_quality import case_quality


def read_status(path: Path) -> dict[str, str]:
    """Return the last CSV status row, or an empty mapping when absent."""
    if not path.exists():
        return {}
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    return rows[-1] if rows else {}


def read_key_values(path: Path) -> dict[str, str]:
    """Read the atomic key=value runner status."""
    if not path.exists():
        return {}
    return dict(
        line.split("=", 1)
        for line in path.read_text().splitlines()
        if "=" in line
    )


def collect(attempt_root: Path) -> Path:
    """Write results.csv from the attempt's immutable case manifest."""
    manifest_path = attempt_root / "case-manifest.json"
    if not manifest_path.exists():
        raise FileNotFoundError(manifest_path)
    manifest = json.loads(manifest_path.read_text())
    if not manifest:
        raise ValueError(f"empty case manifest: {manifest_path}")

    rows = []
    for item in manifest:
        case_dir = Path(item["case_dir"])
        status = read_status(case_dir / "classification.status")
        runner = read_key_values(case_dir / "runner.status")
        rows.append(
            {
                "caseId": item["CaseNo"],
                "x": item["Rr"],
                "y": item["OhOut"],
                "id": status.get("id", "-1"),
                "runner_state": runner.get("state", "missing"),
                "exit_code": runner.get("exit_code", ""),
                "reason": status.get("reason", "missing_classification"),
                "t": status.get("t", ""),
                "drop_volume": status.get("drop_volume", ""),
                "drop_radius": status.get("drop_radius", ""),
                "drop_axial_position": status.get("drop_axial_position", ""),
                "drop_axial_velocity": status.get("drop_axial_velocity", ""),
                **case_quality(case_dir),
            }
        )

    destination = attempt_root / "results.csv"
    temporary = destination.with_suffix(".csv.tmp")
    with temporary.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(destination)
    resolved = sum(
        row["id"] in {"0", "1"} and row["runner_state"] == "complete"
        for row in rows
    )
    print(f"results={destination} resolved={resolved}/{len(rows)}")
    return destination


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("attempt_root", type=Path)
    args = parser.parse_args()
    collect(args.attempt_root.resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
