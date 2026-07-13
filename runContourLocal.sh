#!/usr/bin/env bash

# Run one contour attempt on a controlled Linux workstation. Cases are
# executed in bounded OpenMP batches; no MPI runtime is involved.

set -euo pipefail

project_root=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
iteration_dir=${1:?usage: runContourLocal.sh ITERATION_DIR}
if [[ "$iteration_dir" != /* ]]; then
  iteration_dir="${project_root}/${iteration_dir}"
fi
cases_csv="${iteration_dir}/cases.csv"
case_root="${iteration_dir}/cases"

threads=${CONTOUR_THREADS_PER_CASE:-8}
workers=${CONTOUR_WORKERS:-3}
max_threads=${CONTOUR_MAX_THREADS:-48}
max_level=${CONTOUR_MAXLEVEL:-}
drop_radius_min=${CONTOUR_DROP_RADIUS_MIN:-}

for value_name in threads workers max_threads; do
  value=${!value_name}
  if [[ ! "$value" =~ ^[1-9][0-9]*$ ]]; then
    echo "$value_name must be a positive integer, got $value" >&2
    exit 2
  fi
done
if (( max_threads > 48 )); then
  echo "CONTOUR_MAX_THREADS may not exceed the workstation ceiling of 48" >&2
  exit 2
fi
if (( workers * threads > max_threads )); then
  echo "requested $workers x $threads = $((workers * threads)) threads; ceiling is $max_threads" >&2
  exit 2
fi
if [[ ! -f "$cases_csv" ]]; then
  echo "Missing proposal: $cases_csv" >&2
  exit 2
fi

exec 9>"${iteration_dir}/local-run.lock"
if ! flock -n 9; then
  echo "Another local runner owns ${iteration_dir}/local-run.lock" >&2
  exit 2
fi

cd "$project_root"
case_count=$(python3 - "$cases_csv" <<'PY'
import csv
import sys

with open(sys.argv[1], newline="") as stream:
    count = sum(1 for _ in csv.DictReader(stream))
if not 1 <= count <= 16:
    raise SystemExit(f"cases.csv must contain 1..16 rows, found {count}")
print(count)
PY
)

# Compile against the project-local, ref-locked source tree. The qcc binary is
# host-built because a qcc copied from another Linux installation embeds its
# original include path.
export BASILISK="${project_root}/basilisk"
qcc_command=${CONTOUR_QCC:-}
if [[ -z "$qcc_command" ]]; then
  qcc_command=$(command -v qcc 2>/dev/null || true)
fi
if [[ -z "$qcc_command" && -x /home/vatsal/CMP-codes/basilisk/src/qcc ]]; then
  qcc_command=/home/vatsal/CMP-codes/basilisk/src/qcc
fi
if [[ -z "$qcc_command" || ! -x "$qcc_command" ]]; then
  echo "No host-built qcc found; set CONTOUR_QCC to its absolute path" >&2
  exit 2
fi
basilisk_ref=$(awk -F= '$1 == "ref" {print $2; exit}' "${project_root}/basilisk/.comphy-lock")

echo "Drop-injection contour workstation run"
echo "host=$(hostname) iteration_dir=${iteration_dir} cases=${case_count}"
echo "workers=${workers} threads_per_case=${threads} concurrent_threads=$((workers * threads)) ceiling=${max_threads}"
"$qcc_command" --version

python3 - "$iteration_dir" "$project_root" "$basilisk_ref" "$case_count" "$workers" "$threads" "$max_threads" "$max_level" "$drop_radius_min" "$qcc_command" <<'PY'
import json
import os
import socket
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

iteration_dir = Path(sys.argv[1])
project_root = Path(sys.argv[2])
commit = subprocess.run(
    ["git", "rev-parse", "HEAD"], cwd=project_root, text=True,
    capture_output=True, check=False,
).stdout.strip()
metadata = {
    "backend": "local-openmp",
    "started_at": datetime.now(timezone.utc).isoformat(),
    "source_commit": commit or "unavailable",
    "basilisk_ref": sys.argv[3],
    "node": socket.gethostname(),
    "cases": int(sys.argv[4]),
    "workers": int(sys.argv[5]),
    "threads_per_case": int(sys.argv[6]),
    "concurrent_threads": int(sys.argv[5]) * int(sys.argv[6]),
    "max_threads": int(sys.argv[7]),
    "max_level_override": sys.argv[8] or None,
    "drop_radius_min_override": sys.argv[9] or None,
    "qcc": sys.argv[10],
    "systemd_unit": os.environ.get("SYSTEMD_UNIT"),
}
path = iteration_dir / "run-metadata.json"
temporary = path.with_suffix(".json.tmp")
temporary.write_text(json.dumps(metadata, indent=2) + "\n")
temporary.replace(path)
PY

materialize=(
  python3 "${project_root}/contourWorkflow/materialize_cases.py"
  "$cases_csv" "$case_root"
  --data-dir "${project_root}/simulationCases/DataFiles"
  --source "${project_root}/simulationCases/coalescenceBubble.c"
  --expected "$case_count"
)
if [[ -n "$max_level" ]]; then
  materialize+=(--MAXlevel "$max_level")
fi
if [[ -n "$drop_radius_min" ]]; then
  materialize+=(--dropRadiusMin "$drop_radius_min")
fi
"${materialize[@]}"

build_dir="${iteration_dir}/build"
mkdir -p "$build_dir" "${iteration_dir}/logs"
cp "${project_root}/simulationCases/coalescenceBubble.c" "$build_dir/"
(
  cd "$build_dir"
  "$qcc_command" -I"${project_root}/src-local" -Wall -O2 -fopenmp -disable-dimensions \
    coalescenceBubble.c -o coalescenceBubbleContour -lm
)

mapfile -t case_dirs < <(
  python3 -c 'import json,sys; print("\n".join(row["case_dir"] for row in json.load(open(sys.argv[1]))))' \
    "${iteration_dir}/case-manifest.json"
)
if [[ ${#case_dirs[@]} -ne $case_count ]]; then
  echo "Expected ${case_count} materialised cases, got ${#case_dirs[@]}" >&2
  exit 2
fi

export OMP_NUM_THREADS="$threads"
export OMP_PLACES=cores
export OMP_PROC_BIND=close

rc=0
for ((batch_start = 0; batch_start < case_count; batch_start += workers)); do
  pids=()
  for ((slot = 0; slot < workers && batch_start + slot < case_count; slot++)); do
    case_dir=${case_dirs[$((batch_start + slot))]}
    case_name=$(basename "$case_dir")
    cpu_start=$((slot * threads))
    cpu_end=$((cpu_start + threads - 1))
    taskset -c "${cpu_start}-${cpu_end}" \
      env OMP_NUM_THREADS="$threads" OMP_PLACES=cores OMP_PROC_BIND=close \
      bash "${project_root}/contourWorkflow/run_one_contour_case.sh" \
        "$case_dir" "${build_dir}/coalescenceBubbleContour" \
        >"${iteration_dir}/logs/${case_name}.out" \
        2>"${iteration_dir}/logs/${case_name}.err" &
    pids+=("$!")
  done
  for pid_value in "${pids[@]}"; do
    wait "$pid_value" || rc=1
  done
done

python3 - "$iteration_dir" "$project_root" <<'PY'
import csv
import json
import sys
from pathlib import Path

root = Path(sys.argv[1])
sys.path.insert(0, str(Path(sys.argv[2]) / "contourWorkflow"))
from result_quality import case_quality

manifest = json.loads((root / "case-manifest.json").read_text())
rows = []
for item in manifest:
    case_dir = Path(item["case_dir"])
    status = {}
    runner = {}
    status_path = case_dir / "classification.status"
    if status_path.exists():
        with status_path.open(newline="") as stream:
            parsed = list(csv.DictReader(stream))
        if parsed:
            status = parsed[-1]
    runner_path = case_dir / "runner.status"
    if runner_path.exists():
        runner = dict(
            line.split("=", 1)
            for line in runner_path.read_text().splitlines()
            if "=" in line
        )
    rows.append({
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
        **case_quality(case_dir),
    })
destination = root / "results.csv"
temporary = destination.with_suffix(".csv.tmp")
with temporary.open("w", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
    writer.writeheader()
    writer.writerows(rows)
temporary.replace(destination)
print(
    f"results={destination} "
    f"resolved={sum(row['id'] in ('0', '1') and row['runner_state'] == 'complete' for row in rows)}/{len(rows)}"
)
PY

exit "$rc"
