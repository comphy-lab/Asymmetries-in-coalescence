#!/usr/bin/env python3
"""
# Rearmable Drop-Injection Contour Campaign

Drive a 16-case-per-iteration Bayesian contour campaign without keeping an
interactive agent alive. The state machine is idempotent: scheduled monitors
may call `advance` repeatedly, and a new batch is submitted only after the
previous batch has produced 16 resolved classifications.
"""

from __future__ import annotations

import argparse
import csv
import fcntl
import json
import math
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from contextlib import contextmanager
from typing import Iterator, Sequence

try:
    from result_quality import case_quality
except ModuleNotFoundError:
    from contourWorkflow.result_quality import case_quality


BATCH_SIZE = 16
DEFAULT_PHASE_ONE_END = 8
DEFAULT_FINAL_ITERATION = 16
DEFAULT_CASE_ID_START = 5000
Y_MIN = 0.01
Y_MAX = 0.075
DROP_RADIUS_MIN = 0.0078125
TERMINAL_SLURM_STATES = {
    "BOOT_FAIL",
    "CANCELLED",
    "COMPLETED",
    "DEADLINE",
    "FAILED",
    "NODE_FAIL",
    "OUT_OF_MEMORY",
    "PREEMPTED",
    "TIMEOUT",
}
TERMINAL_EXECUTION_STATES = TERMINAL_SLURM_STATES | {"SUCCESS"}


@dataclass(frozen=True)
class Campaign:
    """Resolved campaign paths and configuration."""

    root: Path
    project_root: Path
    predictor_root: Path
    backend: str = "slurm"

    @property
    def completed(self) -> Path:
        return self.root / "completed"

    @property
    def proposals(self) -> Path:
        return self.root / "proposals"

    @property
    def contours(self) -> Path:
        return self.root / "contours"

    @property
    def iterations(self) -> Path:
        return self.root / "iterations"

    @property
    def measurements(self) -> Path:
        return self.root / "measurements"


def campaign_config(campaign: Campaign) -> dict[str, object]:
    """Read the immutable campaign configuration."""
    return json.loads((campaign.root / "campaign-config.json").read_text())


def atomic_write(path: Path, text: str) -> None:
    """Replace a small campaign state file atomically."""
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(text)
    temporary.replace(path)


@contextmanager
def campaign_lock(root: Path) -> Iterator[None]:
    """Serialise scheduled ticks so proposal generation cannot race."""
    root.mkdir(parents=True, exist_ok=True)
    with (root / "campaign.lock").open("a+") as stream:
        fcntl.flock(stream.fileno(), fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(stream.fileno(), fcntl.LOCK_UN)


def read_rows(path: Path) -> list[dict[str, str]]:
    """Read a CSV into dictionaries."""
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def write_rows(path: Path, rows: Sequence[dict[str, str]], fields: Sequence[str]) -> None:
    """Write a CSV atomically with an explicit schema."""
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def run(command: Sequence[str], *, cwd: Path | None = None) -> subprocess.CompletedProcess[str]:
    """Run one checked subprocess and retain text output for audit logs."""
    return subprocess.run(
        list(command), cwd=cwd, text=True, capture_output=True, check=True
    )


def available_x_values(project_root: Path) -> list[float]:
    """Return simulatable radius ratios within the campaign domain."""
    values = []
    data_dir = project_root / "simulationCases" / "DataFiles"
    for path in data_dir.glob("InitialConditionRr-*.dat"):
        value = float(path.stem.removeprefix("InitialConditionRr-"))
        if 1.0 <= value <= 16.0:
            values.append(value)
    return sorted(set(values))


def initialise(
    campaign: Campaign,
    seed: Path,
    exclude_x: Sequence[float],
    *,
    final_iteration: int = DEFAULT_FINAL_ITERATION,
    manual_checkpoint_after: int | None = DEFAULT_PHASE_ONE_END,
    case_id_start: int = DEFAULT_CASE_ID_START,
    n_new: int = BATCH_SIZE,
    n_repeats: int = 0,
    max_level: int = 12,
    drop_radius_min: float = DROP_RADIUS_MIN,
    workers: int = 3,
    threads_per_case: int = 8,
    max_threads: int = 48,
    unit_prefix: str = "dropinj",
    allow_unbracketed_edges: bool = False,
    posterior_samples: int = 8,
) -> None:
    """Convert the physical seed table to canonical predictor columns."""
    predictor_help = run(
        [sys.executable, str(campaign.predictor_root / "propose_next_sweep.py"), "--help"]
    ).stdout
    if "--x-candidates" not in predictor_help:
        raise RuntimeError(
            "predictor lacks --x-candidates; use Bayesian-Contour-Predictor PR #3 "
            "or a later release before initialising this campaign"
        )
    seed_rows = read_rows(seed)
    canonical = []
    excluded = []
    for row in seed_rows:
        case_id = row.get("caseId", "").strip()
        x = row.get("x", row.get("Rr", "")).strip()
        y = row.get("y", row.get("Oh", "")).strip()
        label = row.get("id", "").strip()
        if not case_id or not x or not y or label not in {"0", "1"}:
            raise ValueError(f"invalid seed row: {row}")
        canonical_row = {"caseId": case_id, "x": x, "y": y, "id": label}
        if any(abs(float(x) - value) <= 1e-12 for value in exclude_x):
            excluded.append(canonical_row)
        else:
            canonical.append(canonical_row)
    if not canonical:
        raise ValueError("seed table is empty")

    allowed = available_x_values(campaign.project_root)
    if not allowed:
        raise ValueError("no simulatable radius ratios in [1, 16]")
    predictor_commit = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=campaign.predictor_root,
        text=True, capture_output=True, check=False
    ).stdout.strip()
    if final_iteration < 1:
        raise ValueError("final_iteration must be positive")
    if (
        manual_checkpoint_after is not None
        and not 1 <= manual_checkpoint_after < final_iteration
    ):
        raise ValueError("manual checkpoint must lie before final_iteration")
    if n_new < 0 or n_repeats < 0 or n_new + n_repeats != BATCH_SIZE:
        raise ValueError(f"n_new + n_repeats must equal {BATCH_SIZE}")
    if max_level < 1 or workers < 1 or threads_per_case < 1:
        raise ValueError("max_level, workers and threads_per_case must be positive")
    if posterior_samples < 8:
        raise ValueError("posterior_samples must be at least 8")
    if max_threads > 48 or workers * threads_per_case > max_threads:
        raise ValueError("local concurrency exceeds the 48-thread workstation ceiling")
    if case_id_start <= max(int(row["caseId"]) for row in canonical):
        raise ValueError("case_id_start must exceed every numeric seed caseId")
    if not unit_prefix or any(
        character not in "abcdefghijklmnopqrstuvwxyz0123456789-"
        for character in unit_prefix
    ):
        raise ValueError("unit_prefix must contain only lowercase letters, digits and hyphens")

    config = {
        "batch_size": BATCH_SIZE,
        "manual_checkpoint_after": manual_checkpoint_after,
        "final_iteration": final_iteration,
        "case_id_start": case_id_start,
        "n_new": n_new,
        "n_repeats": n_repeats,
        "max_level": max_level,
        "workers": workers,
        "threads_per_case": threads_per_case,
        "max_threads": max_threads,
        "unit_prefix": unit_prefix,
        "allow_unbracketed_edges": allow_unbracketed_edges,
        "posterior_samples": posterior_samples,
        "x_bounds": [1.0, 16.0],
        "y_bounds": [Y_MIN, Y_MAX],
        "allowed_x_values": allowed,
        "excluded_seed_x_values": list(exclude_x),
        "excluded_seed_rows": len(excluded),
        "drop_radius_min": drop_radius_min,
        "drop_radius_resolution": "2*max(Ldomain)/2^MAXlevel",
        "drop_persistence": 3,
        "snapshot_interval": 0.05,
        "project_root": str(campaign.project_root),
        "predictor_root": str(campaign.predictor_root),
        "predictor_commit": predictor_commit or "unavailable",
    }

    destination = campaign.completed / "Sweep-0_completed.csv"
    config_path = campaign.root / "campaign-config.json"
    state_path = campaign.root / "campaign-state.json"
    existing_payload = [
        path for path in campaign.root.iterdir()
        if path.name != "campaign.lock"
    ] if campaign.root.exists() else []
    if existing_payload:
        if not config_path.exists() or not destination.exists() or not state_path.exists():
            raise FileExistsError(
                f"refusing to initialise non-empty or partial campaign {campaign.root}"
            )
        if json.loads(config_path.read_text()) != config:
            raise FileExistsError(f"existing campaign config differs: {config_path}")
        if read_rows(destination) != canonical:
            raise FileExistsError(f"existing campaign seed differs: {destination}")
        print(f"already initialised campaign={campaign.root}")
        return

    campaign.root.mkdir(parents=True, exist_ok=True)
    for directory in (
        campaign.completed,
        campaign.proposals,
        campaign.contours,
        campaign.iterations,
        campaign.measurements,
    ):
        directory.mkdir(exist_ok=True)
    write_rows(destination, canonical, ("caseId", "x", "y", "id"))
    atomic_write(config_path, json.dumps(config, indent=2) + "\n")
    atomic_write(
        state_path,
        json.dumps({"state": "ready", "last_completed_iteration": 0}, indent=2) + "\n",
    )
    print(
        f"initialised seed_rows={len(canonical)} excluded={len(excluded)} "
        f"campaign={campaign.root}"
    )


def completed_iterations(campaign: Campaign) -> list[int]:
    """List completed Bayesian batch numbers."""
    values = []
    for path in campaign.completed.glob("Sweep-*_completed.csv"):
        try:
            values.append(int(path.stem.split("-")[1].split("_")[0]))
        except (IndexError, ValueError):
            continue
    return sorted(set(values))


def last_completed_iteration(campaign: Campaign) -> int:
    """Return the last contiguous promoted iteration, rejecting gaps."""
    completed = completed_iterations(campaign)
    if not completed:
        return 0
    expected = list(range(0, max(completed) + 1))
    if completed != expected:
        raise ValueError(f"non-contiguous completed iterations: {completed}")
    return completed[-1]


def model_args(campaign: Campaign) -> list[str]:
    """Return the pinned production model configuration."""
    config = campaign_config(campaign)
    return [
        "--mode", "monotone-y",
        "--monotone-direction", "decreasing",
        "--contour-fit", "local-linear",
        "--x-scale", "log10",
        "--y-scale", "log10",
        "--x-min", "1",
        "--x-max", "16",
        "--y-min", "0.01",
        "--y-max", "0.075",
        "--grid-size", "21",
        "--posterior-samples", str(config.get("posterior_samples", 8)),
        "--transition-width", "0.04",
        "--label-noise", "0.005",
        "--length-scale-x", "0.18",
        "--scarcity-fraction", "0.125",
        "--scarcity-candidate-bins", "12",
        "--scarcity-fan-width", "0.08",
    ]


def completed_files(campaign: Campaign) -> list[Path]:
    """Return completed tables in numerical iteration order."""
    return [
        campaign.completed / f"Sweep-{iteration}_completed.csv"
        for iteration in completed_iterations(campaign)
    ]


def assess(campaign: Campaign, iteration: int, candidate: Path) -> None:
    """Refresh contour state transactionally before promoting a batch."""
    contour = campaign.contours / f"Sweep-{iteration}_contour.csv"
    contour_pending = contour.with_name(f".{contour.name}.pending")
    state = campaign.root / "predictor-state.json"
    state_pending = state.with_name(f".{state.name}.pending")
    if state.exists():
        shutil.copy2(state, state_pending)
    else:
        state_pending.unlink(missing_ok=True)
    command = [
        sys.executable,
        str(campaign.predictor_root / "assess_contour.py"),
        *(str(path) for path in completed_files(campaign)),
        str(candidate),
        "--outfile",
        str(contour_pending),
        "--state",
        str(state_pending),
        "--iteration",
        str(iteration),
        *model_args(campaign),
    ]
    if campaign_config(campaign).get("allow_unbracketed_edges", False):
        command.append("--allow-unbracketed-edges")
    try:
        result = run(command)
        if not contour_pending.exists() or not state_pending.exists():
            raise RuntimeError("contour assessment did not produce both outputs")
        contour_pending.replace(contour)
        state_pending.replace(state)
        if result.stdout:
            print(result.stdout.strip())
    finally:
        contour_pending.unlink(missing_ok=True)
        state_pending.unlink(missing_ok=True)


def validate_proposal(campaign: Campaign, rows: Sequence[dict[str, str]]) -> None:
    """Enforce the immutable 16-point simulation-domain contract."""
    if len(rows) != BATCH_SIZE:
        raise ValueError(f"proposal has {len(rows)} rows, expected {BATCH_SIZE}")
    allowed = set(float(value) for value in campaign_config(campaign)["allowed_x_values"])
    seen_ids: set[str] = set()
    seen_points: set[tuple[float, float]] = set()
    repeated_points = 0
    for row in rows:
        case_id = row.get("caseId", "").strip()
        try:
            x = float(row.get("x", row.get("Rr", "nan")))
            y = float(row.get("y", row.get("Oh", "nan")))
        except ValueError as error:
            raise ValueError(f"non-numeric proposal row: {row}") from error
        if not case_id or case_id in seen_ids:
            raise ValueError(f"missing or duplicate caseId: {case_id!r}")
        if not math.isfinite(x) or x not in allowed:
            raise ValueError(f"proposal contains unavailable Rr={x}")
        if not math.isfinite(y) or not Y_MIN <= y <= Y_MAX:
            raise ValueError(f"proposal Oh={y} lies outside [{Y_MIN}, {Y_MAX}]")
        if row.get("id", "").strip() != "-1":
            raise ValueError("unrun proposal rows must have id=-1")
        point = (x, y)
        if point in seen_points:
            repeated_points += 1
        seen_ids.add(case_id)
        seen_points.add(point)
    repeat_quota = int(campaign_config(campaign).get("n_repeats", 0))
    if repeated_points > repeat_quota:
        raise ValueError(
            f"proposal contains {repeated_points} repeated points, "
            f"exceeding configured repeat quota {repeat_quota}"
        )


def proposal_for(
    campaign: Campaign, iteration: int, *, output: Path | None = None
) -> Path:
    """Generate a constrained 16-point proposal for one iteration."""
    config = campaign_config(campaign)
    allowed = ",".join(f"{value:g}" for value in config["allowed_x_values"])
    output = output or campaign.proposals / f"Sweep-{iteration}_proposed.csv"
    if output.exists():
        validate_proposal(campaign, read_rows(output))
        return output
    building = output.with_name(f".{output.name}.building")
    building.unlink(missing_ok=True)
    command = [
        sys.executable,
        str(campaign.predictor_root / "propose_next_sweep.py"),
        *(str(path) for path in completed_files(campaign)),
        "--outfile",
        str(building),
        "--n-simulations",
        str(BATCH_SIZE),
        "--n-new",
        str(config.get("n_new", BATCH_SIZE)),
        "--n-repeats",
        str(config.get("n_repeats", 0)),
        "--seed",
        str(12 + iteration),
        "--x-candidates",
        allowed,
        *model_args(campaign),
    ]
    try:
        result = run(command)
        rows = read_rows(building)
        for index, row in enumerate(rows):
            row["caseId"] = str(
                int(config.get("case_id_start", DEFAULT_CASE_ID_START))
                + (iteration - 1) * BATCH_SIZE
                + index
            )
        validate_proposal(campaign, rows)
        fields = tuple(rows[0].keys())
        write_rows(building, rows, fields)
        building.replace(output)
        if result.stdout:
            print(result.stdout.strip())
    finally:
        building.unlink(missing_ok=True)
    return output


def iteration_dir(campaign: Campaign, iteration: int) -> Path:
    """Return the stable directory for one simulation batch."""
    return campaign.iterations / f"iteration-{iteration:02d}"


def case_key(row: dict[str, str]) -> tuple[str, float, float]:
    """Return the immutable identity of one proposed or result row."""
    try:
        return (row["caseId"], float(row["x"]), float(row["y"]))
    except (KeyError, ValueError) as error:
        raise ValueError(f"invalid case row: {row}") from error


def quarantined_evidence(campaign: Campaign) -> set[tuple[int, int, str]]:
    """Return exact iteration/attempt/case records excluded by manual review."""
    path = campaign.root / "quality-quarantine.csv"
    if not path.exists():
        return set()
    quarantine: set[tuple[int, int, str]] = set()
    for row in read_rows(path):
        try:
            iteration = int(row["iteration"])
            attempt = int(row["attempt"])
            case_id = row["caseId"].strip()
        except (AttributeError, KeyError, TypeError, ValueError) as error:
            raise ValueError(f"invalid quality quarantine row: {row}") from error
        if iteration < 1 or attempt < 1 or not case_id:
            raise ValueError(f"invalid quality quarantine row: {row}")
        quarantine.add((iteration, attempt, case_id))
    return quarantine


def effective_quality_state(
    row: dict[str, str], attempt_root: Path
) -> str:
    """Return explicit quality or backfill it from retained case evidence."""
    explicit = row.get("quality_state", "").strip().lower()
    if explicit:
        return explicit
    case_root = attempt_root / "cases"
    if not case_root.exists():
        # Grandfather archived result tables whose case evidence was not retained.
        return "legacy"
    return case_quality(case_root / f"case-{row['caseId']}")["quality_state"]


def merged_resolved_results(
    campaign: Campaign, iteration: int, through_attempt: int
) -> dict[str, dict[str, str]]:
    """Merge valid resolved labels across attempts, rejecting stale evidence."""
    root = iteration_dir(campaign, iteration)
    canonical = read_rows(root / "cases.csv")
    canonical_by_id = {row["caseId"]: case_key(row) for row in canonical}
    if len(canonical_by_id) != len(canonical):
        raise ValueError(f"iteration {iteration} canonical cases contain duplicate caseIds")

    quarantine = quarantined_evidence(campaign)
    resolved: dict[str, dict[str, str]] = {}
    for attempt in range(1, through_attempt + 1):
        attempt_root = root / f"attempt-{attempt:02d}"
        attempt_cases_path = attempt_root / "cases.csv"
        if not attempt_cases_path.exists():
            raise FileNotFoundError(attempt_cases_path)
        attempt_cases = read_rows(attempt_cases_path)
        attempt_keys = [case_key(row) for row in attempt_cases]
        attempt_ids = [key[0] for key in attempt_keys]
        if not attempt_cases or len(set(attempt_ids)) != len(attempt_ids):
            raise ValueError(f"iteration {iteration} attempt {attempt} has invalid cases.csv")
        for key in attempt_keys:
            if canonical_by_id.get(key[0]) != key:
                raise ValueError(
                    f"iteration {iteration} attempt {attempt} contains non-canonical case {key[0]}"
                )

        results_path = attempt_root / "results.csv"
        if not results_path.exists():
            continue
        rows = read_rows(results_path)
        if [case_key(row) for row in rows] != attempt_keys:
            raise ValueError(
                f"iteration {iteration} attempt {attempt} results do not match its cases.csv"
            )
        for row in rows:
            if (iteration, attempt, row["caseId"]) in quarantine:
                continue
            if effective_quality_state(row, attempt_root) == "fail":
                continue
            if row.get("id") not in {"0", "1"} or row.get("runner_state") != "complete":
                continue
            candidate = dict(row)
            candidate["source_attempt"] = str(attempt)
            previous = resolved.get(row["caseId"])
            if previous and previous["id"] != row["id"]:
                raise ValueError(
                    f"conflicting resolved labels for case {row['caseId']}: "
                    f"{previous['id']} and {row['id']}"
                )
            resolved.setdefault(row["caseId"], candidate)
    return resolved


def unresolved_cases(
    campaign: Campaign, iteration: int, through_attempt: int
) -> list[dict[str, str]]:
    """Return canonical cases without a valid result in any completed attempt."""
    root = iteration_dir(campaign, iteration)
    canonical = read_rows(root / "cases.csv")
    resolved = merged_resolved_results(campaign, iteration, through_attempt)
    return [row for row in canonical if row["caseId"] not in resolved]


def stage_proposal(campaign: Campaign, iteration: int, proposal: Path) -> Path:
    """Copy an approved proposal into its immutable run directory."""
    validate_proposal(campaign, read_rows(proposal))
    target_dir = iteration_dir(campaign, iteration)
    target_dir.mkdir(parents=True, exist_ok=True)
    target = target_dir / "cases.csv"
    if target.exists() and target.read_text() != proposal.read_text():
        raise FileExistsError(f"iteration {iteration} already has a different cases.csv")
    if not target.exists():
        shutil.copy2(proposal, target)
    return target_dir


def submit(
    campaign: Campaign,
    iteration: int,
    proposal: Path,
    *,
    new_attempt: bool = False,
    attempt_rows: Sequence[dict[str, str]] | None = None,
) -> str:
    """Submit one isolated attempt, idempotently by default."""
    target_dir = stage_proposal(campaign, iteration, proposal)
    current_file = target_dir / "current-attempt.json"
    current = json.loads(current_file.read_text()) if current_file.exists() else None
    if current and not new_attempt:
        return str(current["job_id"])
    if new_attempt and not current:
        raise ValueError(f"iteration {iteration} has no prior attempt to retry")

    attempt = int(current["attempt"]) + 1 if current else 1
    if attempt_rows is None:
        if current:
            attempt_rows = unresolved_cases(campaign, iteration, int(current["attempt"]))
        else:
            attempt_rows = read_rows(target_dir / "cases.csv")
    if not attempt_rows:
        raise ValueError("all simulations are resolved; repair collection and run advance")

    attempt_root = target_dir / f"attempt-{attempt:02d}"
    attempt_root.mkdir(parents=True, exist_ok=True)
    attempt_cases = attempt_root / "cases.csv"
    fields = tuple(attempt_rows[0].keys())
    if attempt_cases.exists() and read_rows(attempt_cases) != list(attempt_rows):
        raise FileExistsError(f"attempt {attempt} contains a different cases.csv")
    if not attempt_cases.exists():
        write_rows(attempt_cases, attempt_rows, fields)
    job_file = attempt_root / f"{campaign.backend}-job.json"
    if job_file.exists():
        raise FileExistsError(f"refusing to replace attempt record {job_file}")

    execution_result: subprocess.CompletedProcess[str] | None = None
    if campaign.backend == "slurm":
        result = run(
            ["sbatch", "runContourHamilton.sbatch", str(attempt_root)],
            cwd=campaign.project_root,
        )
        words = result.stdout.strip().split()
        if len(words) < 4 or not words[-1].isdigit():
            raise RuntimeError(f"could not parse sbatch output: {result.stdout!r}")
        job_id = words[-1]
    elif campaign.backend == "local":
        config = campaign_config(campaign)
        unit = f"{config.get('unit_prefix', 'dropinj')}-i{iteration:02d}-a{attempt:02d}"
        run(
            [
                "systemd-run", "--user", "--unit", unit,
                "--property", "Nice=10",
                "--property", "IOSchedulingClass=best-effort",
                "--property", "IOSchedulingPriority=7",
                "--setenv", f"SYSTEMD_UNIT={unit}.service",
                "--setenv", f"CONTOUR_MAXLEVEL={int(config.get('max_level', 12))}",
                "--setenv", f"CONTOUR_DROP_RADIUS_MIN={float(config.get('drop_radius_min', DROP_RADIUS_MIN)):g}",
                "--setenv", f"CONTOUR_WORKERS={int(config.get('workers', 3))}",
                "--setenv", f"CONTOUR_THREADS_PER_CASE={int(config.get('threads_per_case', 8))}",
                "--setenv", f"CONTOUR_MAX_THREADS={int(config.get('max_threads', 48))}",
                str(campaign.project_root / "runContourLocal.sh"),
                str(attempt_root),
            ],
            cwd=campaign.project_root,
        )
        job_id = f"local:{unit}.service"
    elif campaign.backend == "inline":
        config = campaign_config(campaign)
        job_id = f"inline:{config.get('unit_prefix', 'dropinj')}-i{iteration:02d}-a{attempt:02d}"
        environment = dict(os.environ)
        environment.update(
            {
                "CONTOUR_MAXLEVEL": str(int(config.get("max_level", 12))),
                "CONTOUR_DROP_RADIUS_MIN": f"{float(config.get('drop_radius_min', DROP_RADIUS_MIN)):g}",
                "CONTOUR_WORKERS": str(int(config.get("workers", 3))),
                "CONTOUR_THREADS_PER_CASE": str(
                    int(config.get("threads_per_case", 8))
                ),
                "CONTOUR_MAX_THREADS": str(int(config.get("max_threads", 48))),
            }
        )
        execution_result = subprocess.run(
            [str(campaign.project_root / "runContourLocal.sh"), str(attempt_root)],
            cwd=campaign.project_root,
            env=environment,
            text=True,
            capture_output=True,
            check=False,
        )
        if execution_result.stdout:
            print(execution_result.stdout.strip())
        if execution_result.stderr:
            print(execution_result.stderr.strip(), file=sys.stderr)
    else:
        raise ValueError(f"unsupported execution backend: {campaign.backend}")
    atomic_write(
        job_file,
        json.dumps(
            {
                "job_id": job_id,
                "iteration": iteration,
                "attempt": attempt,
                "backend": campaign.backend,
            },
            indent=2,
        )
        + "\n",
    )
    atomic_write(current_file, job_file.read_text())
    atomic_write(
        campaign.root / "campaign-state.json",
        json.dumps(
            {
                "state": "running",
                "iteration": iteration,
                "attempt": attempt,
                "job_id": job_id,
            },
            indent=2,
        )
        + "\n",
    )
    print(f"submitted iteration={iteration} attempt={attempt} job={job_id}")
    return job_id


def slurm_state(job_id: str) -> str:
    """Read the current or accounting state of one Slurm allocation."""
    queued = subprocess.run(
        ["squeue", "-h", "-j", job_id, "-o", "%T"],
        text=True,
        capture_output=True,
        check=False,
    ).stdout.strip()
    if queued:
        return queued.splitlines()[0].split()[0].upper()
    accounted = subprocess.run(
        ["sacct", "-n", "-X", "-j", job_id, "--format=State"],
        text=True,
        capture_output=True,
        check=False,
    ).stdout.strip()
    return accounted.splitlines()[0].split()[0].split("+")[0].upper() if accounted else "UNKNOWN"


def local_state(job_id: str) -> str:
    """Read the state of one user-systemd workstation batch."""
    unit = job_id.removeprefix("local:")
    result = subprocess.run(
        [
            "systemctl", "--user", "show", unit,
            "--property=LoadState", "--property=ActiveState", "--property=Result",
        ],
        text=True,
        capture_output=True,
        check=False,
    )
    fields = dict(
        line.split("=", 1)
        for line in result.stdout.splitlines()
        if "=" in line
    )
    if fields.get("LoadState") == "not-found":
        return "UNKNOWN"
    active = fields.get("ActiveState", "unknown")
    if active in {"active", "activating", "deactivating"}:
        return "RUNNING"
    if active == "failed" or fields.get("Result") not in {"", "success"}:
        return "FAILED"
    if active == "inactive" and fields.get("Result") == "success":
        return "COMPLETED"
    return "UNKNOWN"


def execution_state(job_id: str) -> str:
    """Read a Slurm or local-systemd execution state from its stable job ID."""
    if job_id.startswith("inline:"):
        return "COMPLETED"
    return local_state(job_id) if job_id.startswith("local:") else slurm_state(job_id)


def collect(campaign: Campaign, iteration: int, attempt: int) -> Path:
    """Merge attempts and promote labels plus their measured drop sizes."""
    iteration_root = iteration_dir(campaign, iteration)
    expected = read_rows(iteration_root / "cases.csv")
    resolved = merged_resolved_results(campaign, iteration, attempt)
    if len(resolved) != BATCH_SIZE:
        raise ValueError(
            f"iteration {iteration} has {len(resolved)}/{BATCH_SIZE} resolved labels"
        )
    destination = campaign.completed / f"Sweep-{iteration}_completed.csv"
    measurements = campaign.measurements / f"Sweep-{iteration}_measurements.csv"
    canonical = [
        {
            "caseId": row["caseId"],
            "x": row["x"],
            "y": row["y"],
            "id": resolved[row["caseId"]]["id"],
        }
        for row in expected
    ]
    measurement_fields = (
        "iteration", "source_attempt", "caseId", "x", "y", "id",
        "drop_volume", "drop_radius", "drop_axial_position",
        "drop_axial_velocity", "reason", "t", "runner_state",
        "exit_code", "max_ke", "facet_lines", "quality_state", "quality_reason",
    )
    measurement_rows = [
        {**resolved[row["caseId"]], "iteration": str(iteration)}
        for row in expected
    ]
    canonical_measurements = [
        {field: row.get(field, "") for field in measurement_fields}
        for row in measurement_rows
    ]
    if destination.exists() and read_rows(destination) != canonical:
        raise FileExistsError(f"refusing to replace {destination}")
    if measurements.exists() and read_rows(measurements) != canonical_measurements:
        raise FileExistsError(f"refusing to replace {measurements}")
    if destination.exists():
        return destination
    pending = destination.with_name(f".{destination.name}.pending")
    measurements_pending = measurements.with_name(f".{measurements.name}.pending")
    try:
        write_rows(pending, canonical, ("caseId", "x", "y", "id"))
        write_rows(measurements_pending, measurement_rows, measurement_fields)
        assess(campaign, iteration, pending)
        pending.replace(destination)
        measurements_pending.replace(measurements)
    finally:
        pending.unlink(missing_ok=True)
        measurements_pending.unlink(missing_ok=True)
    return destination


def advance(campaign: Campaign, *, submit_jobs: bool) -> None:
    """Collect a finished batch and, when permitted, launch the next one."""
    config = campaign_config(campaign)
    final_iteration = int(config.get("final_iteration", DEFAULT_FINAL_ITERATION))
    manual_checkpoint_after = config.get(
        "manual_checkpoint_after", DEFAULT_PHASE_ONE_END
    )
    last_completed = last_completed_iteration(campaign)
    if last_completed >= final_iteration:
        atomic_write(
            campaign.root / "campaign-state.json",
            json.dumps({"state": "complete", "last_completed_iteration": last_completed}, indent=2) + "\n",
        )
        print(f"campaign complete at iteration {last_completed}")
        return

    current_iteration = last_completed + 1
    current_dir = iteration_dir(campaign, current_iteration)
    current_file = current_dir / "current-attempt.json"
    if current_file.exists():
        current = json.loads(current_file.read_text())
        job_id = str(current["job_id"])
        attempt = int(current["attempt"])
        state = execution_state(job_id)
        print(
            f"iteration={current_iteration} attempt={attempt} "
            f"job={job_id} execution_state={state}"
        )
        if state not in TERMINAL_EXECUTION_STATES:
            return
        try:
            collect(campaign, current_iteration, attempt)
        except (FileNotFoundError, ValueError, RuntimeError, subprocess.CalledProcessError) as error:
            atomic_write(
                campaign.root / "campaign-state.json",
                json.dumps(
                    {
                        "state": "needs_attention",
                        "iteration": current_iteration,
                        "attempt": attempt,
                        "job_id": job_id,
                        "execution_state": state,
                        "reason": str(error),
                        "recovery": "inspect evidence, then run retry --submit only if simulations are unresolved",
                    },
                    indent=2,
                )
                + "\n",
            )
            print(f"needs attention: {error}")
            return
        last_completed = current_iteration
        if last_completed >= final_iteration:
            atomic_write(
                campaign.root / "campaign-state.json",
                json.dumps(
                    {"state": "complete", "last_completed_iteration": last_completed},
                    indent=2,
                )
                + "\n",
            )
            print(f"campaign complete at iteration {last_completed}")
            return

    if (
        manual_checkpoint_after is not None
        and last_completed == int(manual_checkpoint_after)
    ):
        manual = campaign.proposals / "Sweep-9_manual-approved.csv"
        if not manual.exists():
            atomic_write(
                campaign.root / "campaign-state.json",
                json.dumps(
                    {
                        "state": "manual_checkpoint",
                        "last_completed_iteration": int(manual_checkpoint_after),
                        "next_required": "human reviews a proposed Sweep 9 and approves exactly 16 cases",
                    },
                    indent=2,
                )
                + "\n",
            )
            print("manual checkpoint after iteration 8; iteration 9 not submitted")
            return
        proposal = manual
        next_iteration = 9
    else:
        next_iteration = last_completed + 1
        proposal = campaign.proposals / f"Sweep-{next_iteration}_proposed.csv"
        proposal = proposal_for(campaign, next_iteration, output=proposal)

    if submit_jobs:
        submit(campaign, next_iteration, proposal)
    else:
        stage_proposal(campaign, next_iteration, proposal)
        print(f"staged iteration={next_iteration}; submission disabled")


def propose_manual_batch(campaign: Campaign) -> Path:
    """Generate, but never approve, the iteration-9 review candidate."""
    checkpoint = campaign_config(campaign).get(
        "manual_checkpoint_after", DEFAULT_PHASE_ONE_END
    )
    if checkpoint is None or last_completed_iteration(campaign) != int(checkpoint):
        raise ValueError("manual iteration-9 proposal is valid only after iteration 8")
    destination = campaign.proposals / "Sweep-9_manual-candidate.csv"
    proposal_for(campaign, 9, output=destination)
    print(f"manual candidate ready for review: {destination}")
    return destination


def retry_current_iteration(campaign: Campaign, *, submit_job: bool) -> str:
    """Create an isolated attempt for a terminal unresolved simulation batch."""
    iteration = last_completed_iteration(campaign) + 1
    root = iteration_dir(campaign, iteration)
    current_file = root / "current-attempt.json"
    if not current_file.exists():
        raise ValueError(f"iteration {iteration} has no submitted attempt")
    current = json.loads(current_file.read_text())
    state = execution_state(str(current["job_id"]))
    if state not in TERMINAL_EXECUTION_STATES:
        raise ValueError(f"job {current['job_id']} is still {state}; refusing retry")
    attempt = int(current["attempt"])
    remaining = unresolved_cases(campaign, iteration, attempt)
    if not remaining:
        raise ValueError("all simulations are resolved; repair collection and run advance")
    if not submit_job:
        raise ValueError("retry is explicit: pass --submit after inspecting the failed attempt")
    proposal = root / "cases.csv"
    return submit(
        campaign,
        iteration,
        proposal,
        new_attempt=True,
        attempt_rows=remaining,
    )


def approve_manual_batch(campaign: Campaign, source: Path) -> None:
    """Install Taylor's reviewed iteration-9 batch at the phase boundary."""
    config = campaign_config(campaign)
    checkpoint = config.get("manual_checkpoint_after", DEFAULT_PHASE_ONE_END)
    if checkpoint is None or last_completed_iteration(campaign) != int(checkpoint):
        raise ValueError("manual iteration-9 approval is valid only after iteration 8")
    rows = read_rows(source)
    if len(rows) != BATCH_SIZE:
        raise ValueError(f"manual batch needs {BATCH_SIZE} rows")
    for index, row in enumerate(rows):
        row["caseId"] = str(
            int(config.get("case_id_start", DEFAULT_CASE_ID_START))
            + 8 * BATCH_SIZE
            + index
        )
        row["x"] = row.get("x", row.get("Rr", ""))
        row["y"] = row.get("y", row.get("Oh", ""))
        row["id"] = row.get("id", "-1")
    validate_proposal(campaign, rows)
    destination = campaign.proposals / "Sweep-9_manual-approved.csv"
    fields = ("caseId", "x", "y", "id", *(
        key for key in rows[0].keys() if key not in {"caseId", "x", "y", "id"}
    ))
    if destination.exists() and read_rows(destination) != rows:
        raise FileExistsError(destination)
    if not destination.exists():
        write_rows(destination, rows, fields)
    print(f"approved manual batch {destination}")


def build_campaign(args: argparse.Namespace) -> Campaign:
    """Resolve user-supplied campaign paths."""
    return Campaign(
        root=args.campaign_root.resolve(),
        project_root=args.project_root.resolve(),
        predictor_root=args.predictor_root.resolve(),
        backend=args.backend,
    )


def main() -> int:
    """Dispatch campaign lifecycle commands."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--campaign-root", required=True, type=Path)
    parser.add_argument("--project-root", default=Path.cwd(), type=Path)
    parser.add_argument("--predictor-root", required=True, type=Path)
    parser.add_argument(
        "--backend", choices=("slurm", "local", "inline"), default="slurm"
    )
    subparsers = parser.add_subparsers(dest="command", required=True)
    init_parser = subparsers.add_parser("init")
    init_parser.add_argument("--seed", required=True, type=Path)
    init_parser.add_argument("--exclude-x", action="append", type=float, default=[])
    init_parser.add_argument("--final-iteration", type=int, default=DEFAULT_FINAL_ITERATION)
    init_parser.add_argument("--no-manual-checkpoint", action="store_true")
    init_parser.add_argument("--case-id-start", type=int, default=DEFAULT_CASE_ID_START)
    init_parser.add_argument("--n-new", type=int, default=BATCH_SIZE)
    init_parser.add_argument("--n-repeats", type=int, default=0)
    init_parser.add_argument("--max-level", type=int, default=12)
    init_parser.add_argument("--drop-radius-min", type=float, default=DROP_RADIUS_MIN)
    init_parser.add_argument("--workers", type=int, default=3)
    init_parser.add_argument("--threads-per-case", type=int, default=8)
    init_parser.add_argument("--max-threads", type=int, default=48)
    init_parser.add_argument("--unit-prefix", default="dropinj")
    init_parser.add_argument("--allow-unbracketed-edges", action="store_true")
    init_parser.add_argument("--posterior-samples", type=int, default=8)
    advance_parser = subparsers.add_parser("advance")
    advance_parser.add_argument("--submit", action="store_true")
    retry_parser = subparsers.add_parser("retry")
    retry_parser.add_argument("--submit", action="store_true")
    subparsers.add_parser("propose-manual-batch")
    approve_parser = subparsers.add_parser("approve-manual-batch")
    approve_parser.add_argument("source", type=Path)
    subparsers.add_parser("status")
    args = parser.parse_args()
    campaign = build_campaign(args)

    with campaign_lock(campaign.root):
        if args.command == "init":
            initialise(
                campaign,
                args.seed,
                args.exclude_x,
                final_iteration=args.final_iteration,
                manual_checkpoint_after=(
                    None if args.no_manual_checkpoint else DEFAULT_PHASE_ONE_END
                ),
                case_id_start=args.case_id_start,
                n_new=args.n_new,
                n_repeats=args.n_repeats,
                max_level=args.max_level,
                drop_radius_min=args.drop_radius_min,
                workers=args.workers,
                threads_per_case=args.threads_per_case,
                max_threads=args.max_threads,
                unit_prefix=args.unit_prefix,
                allow_unbracketed_edges=args.allow_unbracketed_edges,
                posterior_samples=args.posterior_samples,
            )
        elif args.command == "advance":
            advance(campaign, submit_jobs=args.submit)
        elif args.command == "retry":
            retry_current_iteration(campaign, submit_job=args.submit)
        elif args.command == "propose-manual-batch":
            propose_manual_batch(campaign)
        elif args.command == "approve-manual-batch":
            approve_manual_batch(campaign, args.source)
        elif args.command == "status":
            state = campaign.root / "campaign-state.json"
            print(state.read_text() if state.exists() else '{"state": "uninitialised"}')
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
