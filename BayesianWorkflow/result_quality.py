#!/usr/bin/env python3
"""Compute conservative quality evidence for one contour simulation case."""

from __future__ import annotations

import math
from pathlib import Path


MAX_KE_LIMIT = 1000.0
FACET_LINE_LIMIT = 8000


def kinetic_energy_evidence(path: Path) -> tuple[float | None, list[str]]:
    """Return maximum finite kinetic energy and any corruption reasons."""
    if not path.is_file():
        return None, ["missing_log"]
    values: list[float] = []
    reasons: list[str] = []
    for line in path.read_text(errors="replace").splitlines():
        fields = line.split()
        if len(fields) < 4:
            continue
        try:
            float(fields[0])
            float(fields[1])
            float(fields[2])
            value = float(fields[3])
        except ValueError:
            continue
        if not math.isfinite(value):
            reasons.append("nonfinite_ke")
        else:
            values.append(value)
    if not values:
        reasons.append("missing_ke")
        return None, reasons
    maximum = max(values)
    if maximum > MAX_KE_LIMIT:
        reasons.append("max_ke_exceeds_1000")
    return maximum, reasons


def facet_evidence(path: Path) -> tuple[int | None, list[str]]:
    """Return nonblank facet-line count and any corruption reasons."""
    if not path.is_file():
        return None, ["missing_facets"]
    lines = [line for line in path.read_text(errors="replace").splitlines() if line.strip()]
    reasons: list[str] = []
    for line in lines:
        try:
            values = [float(value) for value in line.split()]
        except ValueError:
            reasons.append("nonfinite_facets")
            break
        if not values or any(not math.isfinite(value) for value in values):
            reasons.append("nonfinite_facets")
            break
    if len(lines) > FACET_LINE_LIMIT:
        reasons.append("facet_lines_exceeds_8000")
    return len(lines), reasons


def case_quality(case_dir: Path) -> dict[str, str]:
    """Return CSV-ready quality columns for one completed case directory."""
    max_ke, ke_reasons = kinetic_energy_evidence(case_dir / "log")
    facet_lines, facet_reasons = facet_evidence(case_dir / "interface-latest.dat")
    reasons = list(dict.fromkeys((*ke_reasons, *facet_reasons)))
    return {
        "max_ke": "" if max_ke is None else f"{max_ke:.12g}",
        "facet_lines": "" if facet_lines is None else str(facet_lines),
        "quality_state": "fail" if reasons else "pass",
        "quality_reason": ";".join(reasons) if reasons else "within_gates",
    }
