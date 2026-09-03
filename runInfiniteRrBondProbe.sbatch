#!/bin/bash
#SBATCH --job-name=sbbInfRrBo
#SBATCH --account=nctth181
#SBATCH --partition=genoa
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=48
#SBATCH --time=72:00:00
#SBATCH --output=slurm-infrr-bond-%j.out
#SBATCH --error=slurm-infrr-bond-%j.err

# Stack-1 infinite-Rr Bond probe at the Oh=0.0460 rung (same physics as
# case 6520) with Bursting-Bubble Young-Laplace cavities at Bo=0.001 and
# Bo=0.01. Bond is argv 27: it selects DataFiles/Bo%5.4f.dat and applies
# G.x = -Bo via reduced.h. The Bo=0 control is the already-running case 6520.

set -euo pipefail

MODE="${MODE:-production}"
MAXLEVEL="${MAXLEVEL:-13}"
ZWALL="${ZWALL:-4.0}"
WALLCLEARANCE="${WALLCLEARANCE:--1}"
MURIN="${MURIN:-1e-2}"
OH="${OH:-0.0460}"
THREADS="${SLURM_CPUS_PER_TASK:-48}"

CASE_IDS=(6571 6572)
BOND_VALUES=(0.001 0.01)

case "$MODE" in
  preflight)
    TMAX=0.15; SNAPINT=0.05
    for i in 0 1; do CASE_IDS[$i]=$(( ${CASE_IDS[$i]} + 80 )); done
    ;;
  production)
    TMAX="${TMAX:-1.5}"; SNAPINT="${SNAPINT:-0.05}"
    ;;
  *) echo "ERROR: MODE must be preflight or production" >&2; exit 2 ;;
esac

ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
cd "$ROOT"

source "$ROOT/src-local/solver_args.sh"
solver_default_lookup() {
  printf '%s' "${COALESCENCE_SOLVER_ARG_DEFAULTS[$1]-}"
}

module purge
module load 2024
module load OpenMPI/5.0.3-GCC-13.3.0

BASILISK_ROOT="${BASILISK_ROOT:-${HOME}/basilisk}"
EXPECTED_BASILISK_REF="${EXPECTED_BASILISK_REF:-v2026-07-20}"
export BASILISK="${BASILISK_ROOT}/src"
export PATH="${BASILISK}:${PATH}"
if [[ ! -x "${BASILISK}/qcc" ]]; then
  echo "ERROR: no qcc at ${BASILISK}/qcc" >&2; exit 1
fi
if [[ ! -f "${BASILISK_ROOT}/.comphy-lock" ]]; then
  echo "ERROR: no ${BASILISK_ROOT}/.comphy-lock; Basilisk install is not ref-locked" >&2
  exit 1
fi
if ! grep -qx "ref=${EXPECTED_BASILISK_REF}" "${BASILISK_ROOT}/.comphy-lock"; then
  echo "ERROR: expected Basilisk ref=${EXPECTED_BASILISK_REF}, lock file has:" >&2
  grep '^ref=' "${BASILISK_ROOT}/.comphy-lock" >&2 || echo "  (no ref= line)" >&2
  exit 1
fi

echo "mode=$MODE tmax=$TMAX snapshotInterval=$SNAPINT MAXlevel=$MAXLEVEL zWall=$ZWALL wallClearance=$WALLCLEARANCE MuRin=$MURIN Oh=$OH threads_per_case=$THREADS"
echo "cases=${CASE_IDS[*]} bond=${BOND_VALUES[*]}"

for header in adapt_wavelet_limited.h jetFoot.h; do
  if [[ ! -f "$ROOT/src-local/$header" ]]; then
    echo "ERROR: staged tree is missing src-local/$header" >&2
    exit 1
  fi
done

cd simulationCases

BINARY="burstingBubbleInfiniteRr-Bond-${SLURM_JOB_ID:-$$}"
qcc -I. -I../src-local -O2 -Wall -disable-dimensions -fopenmp \
  burstingBubbleInfiniteRr.c -o "$BINARY" -lm

pids=()
for index in 0 1; do
  case_id="${CASE_IDS[$index]}"
  bond="${BOND_VALUES[$index]}"
  case_dir="$ROOT/simulationCases/$case_id"
  mkdir -p "$case_dir/intermediate"
  cp "$BINARY" "$case_dir/"
  ln -sfn "$ROOT/simulationCases/DataFiles" "$case_dir/DataFiles"

  coalescence_solver_args solver_default_lookup \
    "OhOut=$OH" "RhoIn=1e-3" "Rr=1000" "MAXlevel=$MAXLEVEL" \
    "tmax=$TMAX" "zWall=$ZWALL" \
    "dropRadiusMin=0.021005127" "dropPersistence=3" \
    "snapshotInterval=$SNAPINT" "drillAMR=0" \
    "geometryMode=halfspace" "wallClearance=$WALLCLEARANCE" \
    "interfaceFloor=1" "MuRin=$MURIN" "Bond=$bond"
  solver_args=("${COALESCENCE_SOLVER_ARGS[@]}")

  {
    printf 'CaseNo=%s\nOh=%s\nBond=%s\nMuRin=%s\nMAXlevel=%s\ntmax=%s\nzWall=%s\n' \
      "$case_id" "$OH" "$bond" "$MURIN" "$MAXLEVEL" "$TMAX" "$ZWALL"
    printf 'drillAMR=0\ngeometryMode=halfspace\nwallClearance=%s\n' "$WALLCLEARANCE"
    printf 'dropRadiusMin=0.021005127\ndropPersistence=3\nsnapshotInterval=%s\n' "$SNAPINT"
    printf 'initialShape=Bo%5.4f.dat\n' "$bond"
    printf 'solverStack=filtered+double-projection\nargv=%s\n' "${solver_args[*]}"
  } > "$case_dir/case.params"

  (
    cd "$case_dir"
    export OMP_NUM_THREADS="$THREADS"
    export OMP_PROC_BIND=close
    export OMP_PLACES=cores
    srun --exclusive --exact -N1 -n1 -c "$THREADS" --cpu-bind=cores \
      "./$BINARY" "${solver_args[@]}" \
        > runner.out 2> runner.err &
    echo $! > .srun.pid
    wait $!
    rm -f .srun.pid
  ) &
  pids+=("$!")
done

watchdog_nap() {
  local dir="$1" remaining="$2" chunk
  while (( remaining > 0 )); do
    [[ -e "$dir/.srun.pid" ]] || return 1
    chunk=$(( remaining < 15 ? remaining : 15 ))
    sleep "$chunk"
    remaining=$(( remaining - chunk ))
  done
  return 0
}

watchdog() {
  local dir="$1" last_t="" stalls=0 row t dt pid
  while :; do
    watchdog_nap "$dir" 600 || return 0
    pid=$(cat "$dir/.srun.pid" 2>/dev/null) || true
    [[ -n "${pid:-}" ]] && kill -0 "$pid" 2>/dev/null || return 0
    tr '\0' ' ' < "/proc/$pid/cmdline" 2>/dev/null \
      | grep -q -- "$BINARY" || return 0
    row=$(grep -E '^[0-9]+ ' "$dir/log" 2>/dev/null | tail -1) || continue
    [[ -n "$row" ]] || continue
    t=$(awk '{print $3}' <<< "$row")
    dt=$(awk '{print $2}' <<< "$row")
    if [[ "$t" == "$last_t" ]] && awk -v d="$dt" 'BEGIN{exit !(d+0 < 1e-9)}'; then
      stalls=$((stalls + 1))
      echo "WATCHDOG: $(basename "$dir") no progress at t=$t dt=$dt (${stalls}/3)" >&2
      if (( stalls >= 3 )); then
        echo "WATCHDOG: $(basename "$dir") stalled — terminating" >&2
        echo "watchdog=stalled t=$t dt=$dt" > "$dir/WATCHDOG"
        kill "$pid" 2>/dev/null || true
        sleep 30
        kill -9 "$pid" 2>/dev/null || true
        return 0
      fi
    else
      stalls=0
    fi
    last_t="$t"
  done
}

wpids=()
for case_id in "${CASE_IDS[@]}"; do
  watchdog "$ROOT/simulationCases/$case_id" &
  wpids+=("$!")
done

rc=0
for pid in "${pids[@]}"; do
  wait "$pid" || rc=1
done
for wpid in "${wpids[@]}"; do
  kill "$wpid" 2>/dev/null || true
done

for case_id in "${CASE_IDS[@]}"; do
  d="$ROOT/simulationCases/$case_id"
  echo "=== case=$case_id ==="
  test -s "$d/log"            || { echo "ERROR: no log for $case_id" >&2; rc=1; }
  test -s "$d/jet.log"        || { echo "ERROR: no jet.log for $case_id" >&2; rc=1; }
  test -s "$d/components.log" || { echo "ERROR: no components.log for $case_id" >&2; rc=1; }
  test -s "$d/foot.dat"       || { echo "ERROR: no foot.dat for $case_id" >&2; rc=1; }
  if [[ -s "$d/jet.log" ]]; then
    last_t=$(tail -1 "$d/jet.log" | cut -d, -f1)
    echo "  reached t=$last_t (horizon $TMAX)"
    awk -v t="$last_t" -v h="$TMAX" 'BEGIN {exit !(t < h - 0.02)}' \
      && { echo "  ERROR: $case_id stopped at t=$last_t before its $TMAX horizon" >&2
           [[ -s "$d/projection_guard.status" ]] && head -2 "$d/projection_guard.status" >&2
           rc=1; }
  fi
  [[ -s "$d/WATCHDOG" ]] && { cat "$d/WATCHDOG"; rc=1; }
  [[ -s "$d/classification.status" ]] && head -2 "$d/classification.status"
done

exit "$rc"
