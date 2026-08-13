#!/bin/bash
#SBATCH --job-name=dropContour
#SBATCH -N 1
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=12
#SBATCH --mem=0
#SBATCH --time=24:00:00
#SBATCH -p genoa
#SBATCH --output=slurm-contour-%j.out
#SBATCH --error=slurm-contour-%j.err

# Snellius (genoa) port of runContourHamilton.sbatch: pack one 16-case contour
# iteration onto a single node, one OpenMP case per task.
#
# Differences from the Hamilton original:
#   * `genoa` partition and no consumable `tmp` GRES (Snellius has none), so the
#     nested steps need no --gres request to run concurrently
#   * SURF module stack (2024 + OpenMPI/5.0.3-GCC-13.3.0) instead of gcc/python
#   * the project root lives on scratch; Snellius has no /nobackup
#   * the iteration directory and the initial-shape directory are supplied
#     explicitly, so run output never lands inside the source checkout
#   * per-case COALESCENCE_* solver knobs are exported to every step
#
# The accounting project is deliberately not hardcoded. Submit with an explicit
# account, or export SBATCH_ACCOUNT before submitting:
#
#   sbatch -A <account> runContourSnellius.sbatch ITERATION_DIR DATA_DIR \
#     [extra materialize_cases.py options...]
#
# Mail notifications are opt-in and left unset so this script is reusable and
# does not embed a personal address. Supply them via `sbatch --mail-type=... \
# --mail-user=...` if you want completion/failure emails.

set -euo pipefail

project_root=${SLURM_SUBMIT_DIR}
iteration_dir=${1:?usage: sbatch runContourSnellius.sbatch ITERATION_DIR DATA_DIR [materialize options...]}
data_dir=${2:?usage: sbatch runContourSnellius.sbatch ITERATION_DIR DATA_DIR [materialize options...]}
shift 2
# Anything left over is forwarded verbatim to materialize_cases.py, which is how
# a campaign pins MAXlevel, tmax, zWall, dropRadiusMin and the drill settings
# without editing this launcher.
materialize_options=("$@")

if [[ "$iteration_dir" != /* ]]; then
  iteration_dir="${project_root}/${iteration_dir}"
fi
if [[ "$data_dir" != /* ]]; then
  data_dir="${project_root}/${data_dir}"
fi
cases_csv="${iteration_dir}/cases.csv"
case_root="${iteration_dir}/cases"

# Scratch roots on Snellius. Override COALESCENCE_SCRATCH_PATTERN if the site
# exposes campaign scratch under a different prefix.
scratch_pattern=${COALESCENCE_SCRATCH_PATTERN:-'^(/scratch-shared|/scratch-local|/scratch-node|/gpfs/scratch1)/'}
if [[ ! "$project_root" =~ $scratch_pattern ]]; then
  echo "Snellius campaign must run from scratch, got $project_root" >&2
  exit 2
fi
if [[ ! -f "$cases_csv" ]]; then
  echo "Missing proposal: $cases_csv" >&2
  exit 2
fi
if [[ ! -d "$data_dir" ]]; then
  echo "Missing initial-shape directory: $data_dir" >&2
  exit 2
fi

cd "$project_root"

echo "Drop-injection contour batch (Snellius)"
echo "job=${SLURM_JOB_ID} node=$(hostname) iteration_dir=${iteration_dir}"
echo "account=${SLURM_JOB_ACCOUNT:-unset} data_dir=${data_dir}"
echo "allocated_tasks=${SLURM_NTASKS} cpus_per_task=${SLURM_CPUS_PER_TASK}"

module purge
module load 2024
module load OpenMPI/5.0.3-GCC-13.3.0
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

python3 "${project_root}/BayesianWorkflow/materialize_cases.py" \
  "$cases_csv" "$case_root" \
  --data-dir "$data_dir" \
  --source "${project_root}/simulationCases/coalescenceBubble.c" \
  --expected "$case_count" \
  ${materialize_options[@]+"${materialize_options[@]}"}

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

# Solver tolerances travel with each step rather than living in the binary, so a
# campaign can retune them at submission time.
coalescence_ferr=${COALESCENCE_FERR:-1e-3}
coalescence_interface_eps=${COALESCENCE_INTERFACE_EPS:-1e-6}
echo "COALESCENCE_FERR=${coalescence_ferr} COALESCENCE_INTERFACE_EPS=${coalescence_interface_eps}"

pids=()
for case_dir in "${case_dirs[@]}"; do
  case_name=$(basename "$case_dir")
  # No per-step --mem: the job requests the whole node (--mem=0) and per-case
  # RSS is small, so a large per-step request would exceed node RAM and
  # serialise the launches. Steps share the allocation's memory instead.
  srun --exclusive --exact -N1 -n1 -c "$threads" --cpu-bind=cores \
    env OMP_NUM_THREADS="$threads" OMP_PLACES=cores OMP_PROC_BIND=close \
        COALESCENCE_FERR="$coalescence_ferr" \
        COALESCENCE_INTERFACE_EPS="$coalescence_interface_eps" \
    bash "${project_root}/BayesianWorkflow/run_one_contour_case.sh" \
      "$case_dir" "${build_dir}/coalescenceBubbleContour" \
      >"${iteration_dir}/logs/${case_name}.out" \
      2>"${iteration_dir}/logs/${case_name}.err" &
  pids+=("$!")
done

rc=0
for pid_value in "${pids[@]}"; do
  wait "$pid_value" || rc=1
done

python3 - "$iteration_dir" "$project_root" <<'PY'
import csv
import json
import sys
from pathlib import Path

root = Path(sys.argv[1])
sys.path.insert(0, str(Path(sys.argv[2]) / "BayesianWorkflow"))
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
        "drop_axial_position": status.get("drop_axial_position", ""),
        "drop_axial_velocity": status.get("drop_axial_velocity", ""),
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
