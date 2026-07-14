#!/usr/bin/env bash

# Run one contour attempt on a controlled Linux workstation. Cases are
# executed through a machine-wide physical-CPU pool; no MPI runtime is involved.

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
cpu_start=${CONTOUR_CPU_START:-0}
cpu_count=${CONTOUR_CPU_COUNT:-24}
max_level=${CONTOUR_MAXLEVEL:-}
drop_radius_min=${CONTOUR_DROP_RADIUS_MIN:-}
drill_amr=${CONTOUR_DRILL_AMR:-0}
drill_start=${CONTOUR_DRILL_START:-9}
drill_focus=${CONTOUR_DRILL_FOCUS:-10}
drill_ncells=${CONTOUR_DRILL_NCELLS:-5}
drill_region_min_x=${CONTOUR_DRILL_REGION_MIN_X:--2.1}
drill_arm_steps=${CONTOUR_DRILL_ARM_STEPS:-5}
drill_arm_time=${CONTOUR_DRILL_ARM_TIME:-0}
drill_coarsen_time=${CONTOUR_DRILL_COARSEN_TIME:-0}
drill_region_max_x=${CONTOUR_DRILL_REGION_MAX_X:-3}
drill_region_radius=${CONTOUR_DRILL_REGION_RADIUS:-1.5}
drill_fire_x=${CONTOUR_DRILL_FIRE_X:-0.25}
drill_tip_radius=${CONTOUR_DRILL_TIP_RADIUS:-0.25}
drill_regional_only=${CONTOUR_DRILL_REGIONAL_ONLY:-0}
geometry_mode=${CONTOUR_GEOMETRY_MODE:-finite}
wall_clearance=${CONTOUR_WALL_CLEARANCE:--1}

for value_name in threads workers max_threads cpu_count; do
  value=${!value_name}
  if [[ ! "$value" =~ ^[1-9][0-9]*$ ]]; then
    echo "$value_name must be a positive integer, got $value" >&2
    exit 2
  fi
done
if [[ ! "$cpu_start" =~ ^[0-9]+$ ]]; then
  echo "CONTOUR_CPU_START must be a non-negative integer, got $cpu_start" >&2
  exit 2
fi
if (( max_threads > 48 )); then
  echo "CONTOUR_MAX_THREADS may not exceed the workstation ceiling of 48" >&2
  exit 2
fi
if (( workers * threads > max_threads )); then
  echo "requested $workers x $threads = $((workers * threads)) threads; ceiling is $max_threads" >&2
  exit 2
fi
if (( threads > cpu_count )); then
  echo "one case needs $threads CPUs but the shared pool contains only $cpu_count" >&2
  exit 2
fi
if (( cpu_count % threads != 0 )); then
  echo "CONTOUR_CPU_COUNT ($cpu_count) must be divisible by threads per case ($threads)" >&2
  exit 2
fi
if (( workers > cpu_count / threads )); then
  echo "requested $workers workers but the shared pool has only $((cpu_count / threads)) slots" >&2
  exit 2
fi
online_cpus=$(getconf _NPROCESSORS_ONLN)
if (( cpu_start + cpu_count > online_cpus )); then
  echo "CPU pool ${cpu_start}-$((cpu_start + cpu_count - 1)) exceeds ${online_cpus} online CPUs" >&2
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
echo "shared_cpu_pool=${cpu_start}-$((cpu_start + cpu_count - 1)) scheduler=rolling-global-locks"
"$qcc_command" --version

python3 - "$iteration_dir" "$project_root" "$basilisk_ref" "$case_count" "$workers" "$threads" "$max_threads" "$max_level" "$drop_radius_min" "$qcc_command" "$drill_amr" "$drill_start" "$drill_focus" "$drill_ncells" "$drill_region_min_x" "$drill_arm_steps" "$drill_arm_time" "$drill_coarsen_time" "$drill_region_max_x" "$drill_region_radius" "$drill_fire_x" "$drill_tip_radius" "$drill_regional_only" "$geometry_mode" "$wall_clearance" "$cpu_start" "$cpu_count" <<'PY'
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
    "drill_amr": int(sys.argv[11]),
    "drill_start": int(sys.argv[12]),
    "drill_focus": int(sys.argv[13]),
    "drill_ncells": float(sys.argv[14]),
    "drill_region_min_x": float(sys.argv[15]),
    "drill_arm_steps": int(sys.argv[16]),
    "drill_arm_time": float(sys.argv[17]),
    "drill_coarsen_time": float(sys.argv[18]),
    "drill_region_max_x": float(sys.argv[19]),
    "drill_region_radius": float(sys.argv[20]),
    "drill_fire_x": float(sys.argv[21]),
    "drill_tip_radius": float(sys.argv[22]),
    "drill_regional_only": int(sys.argv[23]),
    "geometry_mode": sys.argv[24],
    "wall_clearance": float(sys.argv[25]),
    "cpu_start": int(sys.argv[26]),
    "cpu_count": int(sys.argv[27]),
    "scheduler": "rolling-global-locks",
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
  --drillAMR "$drill_amr"
  --drillMaxlevelStart "$drill_start"
  --drillMaxlevelFocus "$drill_focus"
  --drillNcells "$drill_ncells"
  --drillRegionMinX "$drill_region_min_x"
  --drillArmSteps "$drill_arm_steps"
  --drillArmTime "$drill_arm_time"
  --drillCoarsenTime "$drill_coarsen_time"
  --drillRegionMaxX "$drill_region_max_x"
  --drillRegionRadius "$drill_region_radius"
  --drillFireX "$drill_fire_x"
  --drillTipRadius "$drill_tip_radius"
  --drillRegionalOnly "$drill_regional_only"
  --geometryMode "$geometry_mode"
  --wallClearance "$wall_clearance"
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

queue_file="${iteration_dir}/.local-run-next-case"
queue_lock="${iteration_dir}/.local-run-next-case.lock"
printf '0\n' >"$queue_file"
lock_root=${CONTOUR_CPU_LOCK_ROOT:-${XDG_RUNTIME_DIR:-/tmp}/drop-injection-contour-cpus}
mkdir -p "$lock_root"

claim_next_case() {
  (
    flock -x 200
    next=$(<"$queue_file")
    if (( next >= case_count )); then
      exit 1
    fi
    printf '%s\n' "$((next + 1))" >"$queue_file"
    printf '%s\n' "$next"
  ) 200>"$queue_lock"
}

acquire_cpu_block() {
  local candidate cpu fd acquired
  while true; do
    for ((candidate = cpu_start; candidate + threads <= cpu_start + cpu_count; candidate += threads)); do
      cpu_lock_fds=()
      acquired=1
      for ((cpu = candidate; cpu < candidate + threads; cpu++)); do
        exec {fd}>"${lock_root}/cpu-${cpu}.lock"
        if flock -n "$fd"; then
          cpu_lock_fds+=("$fd")
        else
          eval "exec ${fd}>&-"
          acquired=0
          break
        fi
      done
      if (( acquired )); then
        acquired_cpu_start=$candidate
        acquired_cpu_end=$((candidate + threads - 1))
        return 0
      fi
      for fd in "${cpu_lock_fds[@]}"; do
        eval "exec ${fd}>&-"
      done
    done
    sleep 1
  done
}

release_cpu_block() {
  local fd
  for fd in "${cpu_lock_fds[@]}"; do
    eval "exec ${fd}>&-"
  done
  cpu_lock_fds=()
}

worker_loop() {
  local index case_dir case_name worker_rc=0
  while index=$(claim_next_case); do
    case_dir=${case_dirs[$index]}
    case_name=$(basename "$case_dir")
    acquire_cpu_block
    echo "${case_name}: CPUs ${acquired_cpu_start}-${acquired_cpu_end}"
    if ! taskset -c "${acquired_cpu_start}-${acquired_cpu_end}" \
      env OMP_NUM_THREADS="$threads" OMP_PLACES=cores OMP_PROC_BIND=close \
      bash "${project_root}/contourWorkflow/run_one_contour_case.sh" \
        "$case_dir" "${build_dir}/coalescenceBubbleContour" \
        >"${iteration_dir}/logs/${case_name}.out" \
        2>"${iteration_dir}/logs/${case_name}.err"; then
      worker_rc=1
    fi
    release_cpu_block
  done
  return "$worker_rc"
}

rc=0
pids=()
for ((worker = 0; worker < workers; worker++)); do
  worker_loop &
  pids+=("$!")
done
for pid_value in "${pids[@]}"; do
  wait "$pid_value" || rc=1
done

python3 "${project_root}/contourWorkflow/collect_attempt_results.py" "$iteration_dir"

exit "$rc"
