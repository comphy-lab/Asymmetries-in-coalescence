#!/usr/bin/env python3
"""
# Materialise a Drop-Injection Contour Batch

Validate one Bayesian proposal or selective-retry subset and create
collision-safe Hamilton case directories. Radius ratios are restricted to
initial-shape files that actually exist in `simulationCases/DataFiles`.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import shutil
from pathlib import Path


DEFAULTS = {
    "RhoIn": "1e-3",
    "MAXlevel": "12",
    "tmax": "1.0",
    "zWall": "0.05",
    # One physical threshold across the map: two finest cells in the largest
    # (Ldomain=16, MAXlevel=12) campaign domain.
    "dropRadiusMin": "0.0078125",
    "dropPersistence": "3",
    "snapshotInterval": "0.05",
    "drillAMR": "0",
    "drillMaxlevelStart": "9",
    "drillMaxlevelFocus": "10",
    "drillNcells": "5",
    "drillRegionMinX": "-2.1",
    "drillArmSteps": "5",
    "drillArmTime": "0",
    "drillCoarsenTime": "0",
    "drillRegionMaxX": "3",
    "drillRegionRadius": "1.5",
    "drillFireX": "0.25",
    "drillTipRadius": "0.25",
    "drillRegionalOnly": "0",
}
PARAMETER_KEYS = (
    "CaseNo",
    "OhOut",
    "RhoIn",
    "Rr",
    "MAXlevel",
    "tmax",
    "zWall",
    "dropRadiusMin",
    "dropPersistence",
    "snapshotInterval",
    "drillAMR",
    "drillMaxlevelStart",
    "drillMaxlevelFocus",
    "drillNcells",
    "drillRegionMinX",
    "drillArmSteps",
    "drillArmTime",
    "drillCoarsenTime",
    "drillRegionMaxX",
    "drillRegionRadius",
    "drillFireX",
    "drillTipRadius",
    "drillRegionalOnly",
)

OH_MIN = 0.01
OH_MAX = 0.075
RADIUS_MATCH_TOLERANCE = 1e-12
SAFE_CASE_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")


def available_radius_ratios(data_dir: Path) -> dict[float, Path]:
    """Return exact numeric radius ratios backed by initial-shape files."""
    available: dict[float, Path] = {}
    for path in sorted(data_dir.glob("InitialConditionRr-*.dat")):
        value = float(path.stem.removeprefix("InitialConditionRr-"))
        if not math.isfinite(value):
            raise ValueError(f"non-finite radius ratio in {path.name}")
        if value in available:
            raise ValueError(
                f"ambiguous initial-shape files for Rr={value:g}: "
                f"{available[value].name}, {path.name}"
            )
        available[value] = path
    return available


def canonicalise(row: dict[str, str]) -> tuple[str, float, float]:
    """Read either physical (`Rr`, `Oh`) or predictor (`x`, `y`) columns."""
    case_id = row.get("caseId", "").strip()
    rr_text = row.get("Rr", row.get("x", "")).strip()
    oh_text = row.get("Oh", row.get("y", "")).strip()
    if not case_id or not rr_text or not oh_text:
        raise ValueError("each row needs caseId and either Rr/Oh or x/y")
    if not SAFE_CASE_ID.fullmatch(case_id):
        raise ValueError(f"unsafe caseId {case_id!r}")
    return case_id, float(rr_text), float(oh_text)


def params_text(values: dict[str, str]) -> str:
    """Return the canonical contents of one auditable parameter file."""
    lines = ["# Drop-injection Bayesian contour case"]
    lines.extend(f"{key}={value}" for key, value in values.items())
    return "\n".join(lines) + "\n"


def item_params(item: dict[str, str]) -> dict[str, str]:
    """Select case parameters in their stable on-disk order."""
    return {key: item[key] for key in PARAMETER_KEYS}


def matched_radius_ratio(rr: float, available: dict[float, Path]) -> float:
    """Match a proposal to a real initial shape without post-hoc rounding."""
    if not math.isfinite(rr):
        raise ValueError(f"Rr must be finite, got {rr}")
    matches = [
        value
        for value in available
        if math.isclose(rr, value, rel_tol=0.0, abs_tol=RADIUS_MATCH_TOLERANCE)
    ]
    if len(matches) != 1:
        raise ValueError(
            f"Rr={rr:.12g} does not exactly match one available initial shape; "
            "constrain Bayesian acquisition instead of rounding after scoring"
        )
    return matches[0]


def validate_existing_case(
    case_dir: Path,
    *,
    source_bytes: bytes,
    expected_params: bytes,
    data_dir: Path,
) -> bool:
    """Return true for a byte-identical reusable case, false for an empty slot."""
    if not case_dir.exists():
        return False
    if not case_dir.is_dir():
        raise ValueError(f"case path is not a directory: {case_dir}")
    if not any(case_dir.iterdir()):
        return False

    source = case_dir / "coalescenceBubble.c"
    params = case_dir / "case.params"
    data_link = case_dir / "DataFiles"
    intermediate = case_dir / "intermediate"
    if not source.is_file() or source.read_bytes() != source_bytes:
        raise ValueError(f"existing case has different source: {case_dir}")
    if not params.is_file() or params.read_bytes() != expected_params:
        raise ValueError(f"existing case has different parameters: {case_dir}")
    if not data_link.is_symlink() or data_link.resolve() != data_dir:
        raise ValueError(f"existing case has different DataFiles link: {case_dir}")
    if not intermediate.is_dir():
        raise ValueError(f"existing case lacks intermediate directory: {case_dir}")
    return True


def main() -> int:
    """Validate a bounded batch and create its case directories."""
    parser = argparse.ArgumentParser()
    parser.add_argument("cases_csv", type=Path)
    parser.add_argument("case_root", type=Path)
    parser.add_argument("--data-dir", required=True, type=Path)
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--expected", type=int, default=16)
    for key, value in DEFAULTS.items():
        parser.add_argument(f"--{key}", default=value)
    args = parser.parse_args()

    if not 1 <= args.expected <= 16:
        parser.error(f"expected row count must be in [1, 16], got {args.expected}")

    with args.cases_csv.open(newline="") as stream:
        raw_rows = list(csv.DictReader(stream))
    if len(raw_rows) != args.expected:
        parser.error(f"expected {args.expected} proposal rows, found {len(raw_rows)}")

    available = available_radius_ratios(args.data_dir)
    if not available:
        parser.error(f"no initial-shape files in {args.data_dir}")

    source_bytes = args.source.read_bytes()
    data_dir = args.data_dir.resolve()
    seen_ids: set[str] = set()
    prepared = []

    # Validate the complete proposal before creating any case directories. This
    # prevents one bad late row from leaving a half-materialised, unrestartable
    # batch behind.
    for raw in raw_rows:
        try:
            case_id, rr, oh = canonicalise(raw)
            matched_rr = matched_radius_ratio(rr, available)
        except ValueError as error:
            parser.error(str(error))
        rr_key = f"{matched_rr:.2f}"
        if not math.isfinite(oh) or not (OH_MIN <= oh <= OH_MAX):
            parser.error(f"Oh must be finite and in [{OH_MIN}, {OH_MAX}], got {oh}")
        oh_key = f"{oh:.12g}"
        if case_id in seen_ids:
            parser.error(f"duplicate caseId {case_id}")
        seen_ids.add(case_id)
        values = {
            "CaseNo": case_id,
            "OhOut": oh_key,
            "RhoIn": str(args.RhoIn),
            "Rr": rr_key,
            "MAXlevel": str(args.MAXlevel),
            "tmax": str(args.tmax),
            "zWall": str(args.zWall),
            "dropRadiusMin": str(args.dropRadiusMin),
            "dropPersistence": str(args.dropPersistence),
            "snapshotInterval": str(args.snapshotInterval),
            "drillAMR": str(args.drillAMR),
            "drillMaxlevelStart": str(args.drillMaxlevelStart),
            "drillMaxlevelFocus": str(args.drillMaxlevelFocus),
            "drillNcells": str(args.drillNcells),
            "drillRegionMinX": str(args.drillRegionMinX),
            "drillArmSteps": str(args.drillArmSteps),
            "drillArmTime": str(args.drillArmTime),
            "drillCoarsenTime": str(args.drillCoarsenTime),
            "drillRegionMaxX": str(args.drillRegionMaxX),
            "drillRegionRadius": str(args.drillRegionRadius),
            "drillFireX": str(args.drillFireX),
            "drillTipRadius": str(args.drillTipRadius),
            "drillRegionalOnly": str(args.drillRegionalOnly),
        }
        prepared.append(
            {
                "caseId": case_id,
                "case_dir": str((args.case_root / f"case-{case_id}").resolve()),
                **values,
            }
        )

    # Check every pre-existing target before writing any new one. Output files
    # from an interrupted allocation are allowed, but the immutable source,
    # parameters and initial-shape link must be byte-for-byte identical.
    reusable: set[Path] = set()
    for item in prepared:
        case_dir = Path(item["case_dir"])
        expected_params = params_text(item_params(item)).encode()
        try:
            if validate_existing_case(
                case_dir,
                source_bytes=source_bytes,
                expected_params=expected_params,
                data_dir=data_dir,
            ):
                reusable.add(case_dir)
        except ValueError as error:
            parser.error(str(error))

    args.case_root.mkdir(parents=True, exist_ok=True)
    for item in prepared:
        case_dir = Path(item["case_dir"])
        if case_dir in reusable:
            continue
        case_dir.mkdir(parents=True, exist_ok=True)
        (case_dir / "intermediate").mkdir(exist_ok=True)
        shutil.copy2(args.source, case_dir / "coalescenceBubble.c")
        (case_dir / "DataFiles").symlink_to(data_dir, target_is_directory=True)
        (case_dir / "case.params").write_text(params_text(item_params(item)))

    manifest = args.case_root.parent / "case-manifest.json"
    temporary = manifest.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(prepared, indent=2) + "\n")
    temporary.replace(manifest)
    print(
        f"materialised={len(prepared) - len(reusable)} "
        f"reused={len(reusable)} manifest={manifest}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
