#!/bin/bash
#SBATCH --job-name=dropContour
#SBATCH -N 1
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=8
#SBATCH --mem=0
#SBATCH --time=3-00:00:00
#SBATCH --gres=tmp:16G
#SBATCH -p multi
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=vatsal.sanjay@durham.ac.uk
#SBATCH --output=slurm-contour-%j.out
#SBATCH --error=slurm-contour-%j.err

# Hamilton's `tmp` resource is a consumable GRES. The job requests 16 GiB and
# each nested step requests 1 GiB; allowing every step to inherit the full GRES
# silently serialises the 16-way launch.

set -euo pipefail

project_root=${SLURM_SUBMIT_DIR}
iteration_dir=${1:?usage: sbatch runContourHamilton.sbatch ITERATION_DIR}
if [[ "$iteration_dir" != /* ]]; then
  iteration_dir="${project_root}/${iteration_dir}"
fi
cases_csv="${iteration_dir}/cases.csv"
case_root="${iteration_dir}/cases"

if [[ ! "$project_root" =~ ^/nobackup/ ]]; then
  echo "Hamilton campaign must run from /nobackup, got $project_root" >&2
  exit 2
fi
if [[ ! -f "$cases_csv" ]]; then
  echo "Missing proposal: $cases_csv" >&2
  exit 2
fi

cd "$project_root"

echo "Drop-injection contour batch"
echo "job=${SLURM_JOB_ID} node=$(hostname) iteration_dir=${iteration_dir}"
echo "allocated_tasks=${SLURM_NTASKS} cpus_per_task=${SLURM_CPUS_PER_TASK}"

module purge
module load gcc/11.2
module load python/3.10.8
python3 --version

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
echo "launched_cases=${case_count}"

source "${project_root}/src-local/basilisk_version.sh"
lock_file="${project_root}/basilisk/.comphy-lock"
lock_ref=""
if [[ -f "$lock_file" ]]; then
  lock_ref=$(awk -F= '$1 == "ref" {print $2; exit}' "$lock_file")
fi
if [[ ! -f "${project_root}/.project_config" || "$lock_ref" != "$BASILISK_REF" ]]; then
  install_basilisk --hard
else
  install_basilisk
fi
source "${project_root}/.project_config"
actual_ref=$(awk -F= '$1 == "ref" {print $2; exit}' "$lock_file")
if [[ "$actual_ref" != "$BASILISK_REF" ]]; then
  echo "Basilisk lock mismatch: expected $BASILISK_REF, got ${actual_ref:-missing}" >&2
  exit 2
fi
qcc --version

python3 - "$iteration_dir" "$project_root" "$BASILISK_REF" "$case_count" <<'PY'
import json
import subprocess
import sys
from pathlib import Path

iteration_dir = Path(sys.argv[1])
project_root = Path(sys.argv[2])
commit = subprocess.run(
    ["git", "rev-parse", "HEAD"], cwd=project_root, text=True,
    capture_output=True, check=False
).stdout.strip()
metadata = {
    "slurm_job_id": __import__("os").environ.get("SLURM_JOB_ID"),
    "source_commit": commit or "unavailable",
    "basilisk_ref": sys.argv[3],
    "node": __import__("socket").gethostname(),
    "cases": int(sys.argv[4]),
    "threads_per_case": int(__import__("os").environ["SLURM_CPUS_PER_TASK"]),
}
path = iteration_dir / "run-metadata.json"
temporary = path.with_suffix(".json.tmp")
temporary.write_text(json.dumps(metadata, indent=2) + "\n")
temporary.replace(path)
PY

python3 "${project_root}/contourWorkflow/materialize_cases.py" \
  "$cases_csv" "$case_root" \
  --data-dir "${project_root}/simulationCases/DataFiles" \
  --source "${project_root}/simulationCases/coalescenceBubble.c" \
  --expected "$case_count"

build_dir="${iteration_dir}/build"
mkdir -p "$build_dir" "${iteration_dir}/logs"
cp "${project_root}/simulationCases/coalescenceBubble.c" "$build_dir/"
(
  cd "$build_dir"
  qcc -I"${project_root}/src-local" -Wall -O2 -fopenmp -disable-dimensions \
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

threads=${SLURM_CPUS_PER_TASK:?}
export OMP_NUM_THREADS="$threads"
export OMP_PLACES=cores
export OMP_PROC_BIND=close

pids=()
for case_dir in "${case_dirs[@]}"; do
  case_name=$(basename "$case_dir")
  srun --exclusive --exact -N1 -n1 -c "$threads" --mem=15000M \
    --gres=tmp:1G --cpu-bind=cores \
    env OMP_NUM_THREADS="$threads" OMP_PLACES=cores OMP_PROC_BIND=close \
    bash "${project_root}/contourWorkflow/run_one_contour_case.sh" \
      "$case_dir" "${build_dir}/coalescenceBubbleContour" \
      >"${iteration_dir}/logs/${case_name}.out" \
      2>"${iteration_dir}/logs/${case_name}.err" &
  pids+=("$!")
done

rc=0
for pid_value in "${pids[@]}"; do
  wait "$pid_value" || rc=1
done

python3 - "$iteration_dir" <<'PY'
import csv
import json
import sys
from pathlib import Path

root = Path(sys.argv[1])
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
