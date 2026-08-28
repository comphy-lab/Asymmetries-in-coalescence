#!/bin/bash
#SBATCH --job-name=sbbDropMap
#SBATCH --account=nctth181
#SBATCH --partition=genoa
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=48
#SBATCH --time=72:00:00
#SBATCH --output=slurm-dropmap-%j.out
#SBATCH --error=slurm-dropmap-%j.err

# Matched-Delta r_d(Oh) ladder, bulk Bo=0, Rr=infinity (zWall=4 -> Ldomain=10).
# Every case runs to its horizon and logs EVERY liquid component plus a jet
# time series; the pinned classifier is recorded but never terminates a run.
#
# Four cases per node at 48 threads each fills the 192-CPU genoa node exactly.
# GROUP selects the Oh quartet; MODE=preflight runs the same quartet to a short
# horizon so configuration errors surface in minutes rather than days.

set -euo pipefail

MODE="${MODE:-production}"
GROUP="${GROUP:-A}"
MAXLEVEL="${MAXLEVEL:-13}"
THREADS="${SLURM_CPUS_PER_TASK:-48}"

# Matched Delta requires MAXlevel=13 with zWall=4: Delta = 10/2^13 = 1.2207e-3,
# identical to the completed cases 6211 (Oh=0.0345) and 6212 (Oh=0.0355).
case "$GROUP" in
  A) CASE_IDS=(6301 6302 6303 6304); OH_VALUES=(0.0100 0.0150 0.0200 0.0250) ;;
  B) CASE_IDS=(6305 6306 6307 6308); OH_VALUES=(0.0270 0.0280 0.0290 0.0300) ;;
  C) CASE_IDS=(6309 6310 6311 6312); OH_VALUES=(0.0310 0.0325 0.0333 0.0365) ;;
  *) echo "ERROR: GROUP must be A, B or C" >&2; exit 2 ;;
esac

case "$MODE" in
  preflight)
    TMAX=0.15; TSNAP=0.05
    for i in 0 1 2 3; do CASE_IDS[$i]=$(( ${CASE_IDS[$i]} + 80 )); done
    ;;
  production)
    TMAX="${TMAX:-1.5}"; TSNAP="${TSNAP:-0.05}"
    ;;
  *) echo "ERROR: MODE must be preflight or production" >&2; exit 2 ;;
esac

ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
cd "$ROOT"

module purge
module load 2024
module load OpenMPI/5.0.3-GCC-13.3.0

BASILISK_ROOT="${BASILISK_ROOT:-${HOME}/basilisk}"
EXPECTED_BASILISK_REF="${EXPECTED_BASILISK_REF:-v2026-07-20}"
export BASILISK="${BASILISK_ROOT}/src"
export PATH="${BASILISK}:${PATH}"
test -x "${BASILISK}/qcc"
test -f "${BASILISK_ROOT}/.comphy-lock"
grep -qx "ref=${EXPECTED_BASILISK_REF}" "${BASILISK_ROOT}/.comphy-lock"

# The base solver must be byte-identical to the one that produced 6211/6212,
# otherwise the ladder is not matched.
expected_solver='2ebb137102b81de6f2622489bc367c4056b728e43517f9f523fd7fb167b18daa'
actual_solver=$(sha256sum simulationCases/burstingBubble-drillResolution.c | awk '{print $1}')
test "$actual_solver" = "$expected_solver"

echo "mode=$MODE group=$GROUP tmax=$TMAX tsnap=$TSNAP MAXlevel=$MAXLEVEL threads_per_case=$THREADS"
echo "solver_sha256=$actual_solver basilisk_ref=$EXPECTED_BASILISK_REF"
echo "cases=${CASE_IDS[*]} oh=${OH_VALUES[*]}"

cd simulationCases
qcc -O2 -Wall -disable-dimensions -fopenmp -I../src-local \
  burstingBubbleDropMap.c -o bursting-bubble-dropmap -lm

pids=()
for index in 0 1 2 3; do
  case_id="${CASE_IDS[$index]}"
  oh="${OH_VALUES[$index]}"
  case_dir="$ROOT/simulationCases/$case_id"
  mkdir -p "$case_dir/intermediate"
  cp bursting-bubble-dropmap "$case_dir/"
  ln -sfn "$ROOT/simulationCases/DataFiles" "$case_dir/DataFiles"
  cat > "$case_dir/case.params" <<PARAMS
CaseNo=$case_id
Oh=$oh
Bond=0
OhRatio=2e-2
MAXlevel=$MAXLEVEL
MINlevel=4
init_grid_level=5
fErr=1e-3
VelErr=1e-3
KErr=1e-6
CFL=0.1
dtmax=1e-2
TOLERANCE=1e-4
tsnap=$TSNAP
tmax=$TMAX
zWall=4.0
drillAMR=1
drillMaxlevelStart=8
drillNcellsK=5
drillNcellsJet=5
drillRelaxLevel=-1
drillTsnapStages=1
drillTsnapMinFactor=0.03125
drillMaxlevelFocus=12
drillRemoveGasSize=3
keStopMax=1e4
keStopMin=1e-6
PARAMS

  (
    cd "$case_dir"
    export OMP_NUM_THREADS="$THREADS"
    export OMP_PROC_BIND=close
    export OMP_PLACES=cores
    srun --exclusive -N1 -n1 -c "$THREADS" --cpu-bind=cores \
      ./bursting-bubble-dropmap case.params > runner.out 2> runner.err
  ) &
  pids+=("$!")
done

rc=0
for pid in "${pids[@]}"; do
  wait "$pid" || rc=1
done

for case_id in "${CASE_IDS[@]}"; do
  d="$ROOT/simulationCases/$case_id"
  echo "=== case=$case_id ==="
  test -s "$d/log"            || { echo "ERROR: no log for $case_id" >&2; rc=1; }
  test -s "$d/jet.log"        || { echo "ERROR: no jet.log for $case_id" >&2; rc=1; }
  test -s "$d/components.log" || { echo "ERROR: no components.log for $case_id" >&2; rc=1; }
  if [[ -s "$d/jet.log" ]]; then
    last_t=$(tail -1 "$d/jet.log" | cut -d, -f1)
    ncomp_max=$(awk -F, 'NR>1 {if ($7+0 > m) m=$7+0} END {print m+0}' "$d/jet.log")
    echo "  reached t=$last_t (horizon $TMAX), max components seen=$ncomp_max"
    awk -v t="$last_t" -v h="$TMAX" 'BEGIN {exit !(t < h - 0.02)}' \
      && { echo "  WARNING: $case_id stopped short of its horizon" >&2; }
  fi
  [[ -s "$d/classification.status" ]] && cat "$d/classification.status"
done

exit "$rc"
