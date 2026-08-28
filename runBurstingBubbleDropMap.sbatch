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

# Matching 6211/6212 means matching the whole build, not just the solver .c.
# burstingBubble-drillResolution.c includes params.h and
# adapt_wavelet_limited.h, so a bundle staged with this repository's own
# src-local would compile against a different regional-AMR header and a
# missing params.h. Guarding all three turns that silent mismatch into an
# immediate failure.
guard_sha256() {
  local path="$1" expected="$2" actual
  actual=$(sha256sum "$path" | awk '{print $1}')
  if [[ "$actual" != "$expected" ]]; then
    echo "ERROR: $path is not the file 6211/6212 were built against." >&2
    echo "       expected $expected" >&2
    echo "       actual   $actual" >&2
    echo "       Stage the bursting-bubble bundle; do not run this from a" >&2
    echo "       bare checkout of this repository." >&2
    exit 3
  fi
  echo "guarded $path $actual"
}

guard_sha256 simulationCases/burstingBubble-drillResolution.c \
  2ebb137102b81de6f2622489bc367c4056b728e43517f9f523fd7fb167b18daa
guard_sha256 src-local/params.h \
  808dc389acb05fde6e7e0bba30be2b88b68181f65c621b6696b5de7bd97b7500
guard_sha256 src-local/adapt_wavelet_limited.h \
  16564742e060837637fbdc7eafb2e10ed94f33565d2c351a36d0ed65e5ccc7cc
actual_solver=2ebb137102b81de6f2622489bc367c4056b728e43517f9f523fd7fb167b18daa

echo "mode=$MODE group=$GROUP tmax=$TMAX tsnap=$TSNAP MAXlevel=$MAXLEVEL threads_per_case=$THREADS"
echo "solver_sha256=$actual_solver basilisk_ref=$EXPECTED_BASILISK_REF"
echo "cases=${CASE_IDS[*]} oh=${OH_VALUES[*]}"

cd simulationCases
# Deliberately no -DFILTERED and no $QCC_FLAGS. The comparator cases 6211/6212
# were built by this project's runBulkVisibleBracket.sbatch with exactly this
# bare line, and their case.params sets no FILTERED. Bursting-Bubble's
# runSnelliusDrillOhSweep.sbatch does set FILTERED=1 through
# append_solver_qcc_flags, but that is a different campaign; adding the flag
# here would smear the two-phase properties on one side of the comparison and
# break the match this ladder exists to make.
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
    srun --exclusive --exact -N1 -n1 -c "$THREADS" --cpu-bind=cores \
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
    # A short stop is a warning, not a failure: keStopMin=1e-6 in case.params
    # ends a run whose kinetic energy has decayed, which is a complete result.
    awk -v t="$last_t" -v h="$TMAX" 'BEGIN {exit !(t < h - 0.02)}' \
      && { echo "  NOTE: $case_id stopped at t=$last_t before its $TMAX horizon;" \
                "check the log for a keStop termination before treating it" \
                "as a failure" >&2; }
  fi
  [[ -s "$d/classification.status" ]] && cat "$d/classification.status"
done

exit "$rc"
