#!/bin/bash
#SBATCH --job-name=sbbDeltaAnchor
#SBATCH --account=nctth181
#SBATCH --partition=genoa
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=48
#SBATCH --time=72:00:00
#SBATCH --output=slurm-anchor-%j.out
#SBATCH --error=slurm-anchor-%j.err

# Delta-convergence anchor across solver stacks, at one Ohnesorge number.
#
# THE QUESTION
#
# At Oh = 0.0280, matched Delta = 1.2207e-3 (MAXlevel 13), the two campaigns
# disagree by 1.55-1.77x in peak jet speed, 2-3 vs 7-10 drops per burst, and
# a tip cap 7x fatter on the drill stack. The tip-cap diagnostic is
# geometrically identical between drivers, so the SOLUTIONS differ. Neither
# stack holds a Delta-convergence certificate at these settings, so neither
# number can yet be called the right one.
#
# THE DESIGN
#
# One Oh, one gas viscosity (the manuscript's 1e-2), three stacks at a coarse
# rung where all of them are healthy, plus the fine rung on the production
# stack:
#
#   6561  ML12  Stack 1 (FILTERED + double projection)
#   6562  ML12  Stack 2 (conserving, no FILTERED)      legal, dt-cheap at ML12
#   6563  ML12  FILTERED + conserving                  the drill pairing, via
#                                                      -DALLOW_FILTERED_CONSERVING
#   6564  ML14  Stack 1                                the convergence rung
#
# With the existing ML13 points (Stack 1: crossdriver 6506; drill pairing:
# ladder 6306/6326 died, 6307/6305 flank it; probe 6553/6554) this gives
# v_max(Delta) and t_pinch(Delta) per stack. Read it as:
#
#   - stacks converge to one v_max as Delta shrinks  -> stack choice is
#     economics; Stack 1 wins on cost (conserving dt is ~13x smaller at ML13).
#   - stacks converge to different v_max             -> one is wrong; the
#     (Oh, Bo) boundary of Walls et al. (2015) and the G&BR v_jet(Oh) form
#     arbitrate, fitted on our own prefactors.
#
# Horizon t = 0.60: past the focusing peak (t ~ 0.50) and the first pinch
# (t ~ 0.505), not the full burst. ML12 cases cost a few hours; the ML14 rung
# is the walltime driver (~5.7x the ML13 rate per unit t) and may stop on
# walltime near t ~ 0.55, which still contains everything this job needs, so
# a shortfall on 6564 alone is reported but does not fail the job.
#
# The pinned detector radius is passed as always, but at ML12 it spans only a
# few cells: drop COUNTS from 6561-6563 are not comparable to ML13 points and
# are not the point of this job. v_max, t_pinch and the neck history are.

set -euo pipefail

MODE="${MODE:-production}"
OH="${OH:-0.0280}"
MURIN="${MURIN:-1e-2}"
ZWALL="${ZWALL:-4.0}"          # Ldomain = min(zWall+6, 16) = 10
WALLCLEARANCE="${WALLCLEARANCE:--1}"
THREADS="${SLURM_CPUS_PER_TASK:-48}"

CASE_IDS=(6561 6562 6563 6564)
MAXLEVELS=(12   12   12   14)
STACK_FLAGS=(""
             "-DUSE_CONSERVING=1"
             "-DUSE_CONSERVING=1 -DFILTERED=1 -DALLOW_FILTERED_CONSERVING=1"
             "")
STACK_NAMES=("filtered+double-projection"
             "conserving"
             "filtered+conserving (ALLOW_FILTERED_CONSERVING)"
             "filtered+double-projection")

case "$MODE" in
  preflight)
    TMAX=0.15; SNAPINT=0.05
    for i in 0 1 2 3; do CASE_IDS[$i]=$(( ${CASE_IDS[$i]} + 80 )); done
    ;;
  production)
    TMAX="${TMAX:-0.60}"; SNAPINT="${SNAPINT:-0.05}"
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
# Fail loudly: `set -e` on a bare `grep -qx` aborts the job with no
# diagnostic, so a ref mismatch looks identical to a missing file.
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

echo "mode=$MODE Oh=$OH MuRin=$MURIN tmax=$TMAX zWall=$ZWALL threads_per_case=$THREADS"
echo "cases=${CASE_IDS[*]} maxlevels=${MAXLEVELS[*]}"

cd simulationCases

pids=()
for index in 0 1 2 3; do
  case_id="${CASE_IDS[$index]}"
  maxlevel="${MAXLEVELS[$index]}"
  flags="${STACK_FLAGS[$index]}"
  stack="${STACK_NAMES[$index]}"
  binary="burstingBubbleInfiniteRr-$case_id-${SLURM_JOB_ID:-$$}"
  case_dir="$ROOT/simulationCases/$case_id"

  # shellcheck disable=SC2086
  qcc -I. -I../src-local -O2 -Wall -disable-dimensions -fopenmp $flags \
    burstingBubbleInfiniteRr.c -o "$binary" -lm

  mkdir -p "$case_dir/intermediate"
  cp "$binary" "$case_dir/"
  ln -sfn "$ROOT/simulationCases/DataFiles" "$case_dir/DataFiles"

  coalescence_solver_args solver_default_lookup \
    "OhOut=$OH" "RhoIn=1e-3" "Rr=1000" "MAXlevel=$maxlevel" \
    "tmax=$TMAX" "zWall=$ZWALL" \
    "dropRadiusMin=0.021005127" "dropPersistence=3" \
    "snapshotInterval=$SNAPINT" "drillAMR=0" \
    "geometryMode=halfspace" "wallClearance=$WALLCLEARANCE" \
    "interfaceFloor=1" "MuRin=$MURIN"
  solver_args=("${COALESCENCE_SOLVER_ARGS[@]}")

  {
    printf 'CaseNo=%s\nOh=%s\nMuRin=%s\nMAXlevel=%s\ntmax=%s\nzWall=%s\n' \
      "$case_id" "$OH" "$MURIN" "$maxlevel" "$TMAX" "$ZWALL"
    printf 'drillAMR=0\ngeometryMode=halfspace\nwallClearance=%s\n' "$WALLCLEARANCE"
    printf 'dropRadiusMin=0.021005127\ndropPersistence=3\nsnapshotInterval=%s\n' "$SNAPINT"
    printf 'solverStack=%s\nqccFlags=%s\nargv=%s\n' "$stack" "$flags" "${solver_args[*]}"
  } > "$case_dir/case.params"

  (
    cd "$case_dir"
    export OMP_NUM_THREADS="$THREADS"
    export OMP_PROC_BIND=close
    export OMP_PLACES=cores
    srun --exclusive --exact -N1 -n1 -c "$THREADS" --cpu-bind=cores \
      "./$binary" "${solver_args[@]}" > runner.out 2> runner.err &
    echo $! > .srun.pid
    wait $!
    rm -f .srun.pid
  ) &
  pids+=("$!")
done

# Two outcomes here are RESULTS, not failures: a divergence through the ke
# gate (that is the question this job asks) and the ML14 rung running out of
# walltime after the focusing peak. Everything else -- a worker that could not
# start, a solver that crashed, a coarse rung stopping short -- has no benign
# explanation and must reach Slurm as a failure, or the job reports success
# while having measured nothing.
worker_rc=()
for index in 0 1 2 3; do
  if wait "${pids[$index]}"; then worker_rc[$index]=0; else worker_rc[$index]=$?; fi
done

echo
echo "=== delta-anchor summary ==="
printf '%-6s %-4s %-48s %-9s %s\n' case ML stack reached verdict
unexplained=0
for index in 0 1 2 3; do
  case_id="${CASE_IDS[$index]}"
  maxlevel="${MAXLEVELS[$index]}"
  d="$ROOT/simulationCases/$case_id"
  last_t="nodata"
  [[ -s "$d/jet.log" ]] && last_t=$(tail -1 "$d/jet.log" | cut -d, -f1)

  ke_gate=0
  grep -q 'kinetic energy blew up' "$d/runner.err" 2>/dev/null && ke_gate=1

  short=0
  if [[ "$last_t" == "nodata" ]] || \
     awk -v t="$last_t" -v h="$TMAX" 'BEGIN {exit !(t < h - 0.02)}'; then
    short=1
  fi

  if (( ke_gate )); then
    verdict="DIVERGED (ke gate) — a result"
  elif (( ! short )); then
    verdict="reached horizon"
  elif [[ "$maxlevel" == "14" ]]; then
    verdict="stopped short at t=$last_t (walltime-permitted on the fine rung)"
  else
    verdict="STOPPED SHORT with no benign explanation"
    unexplained=1
  fi

  # A worker that failed for a reason other than the ke gate is unexplained
  # regardless of how far it got.
  if (( worker_rc[index] != 0 && ! ke_gate )); then
    verdict="$verdict; worker exit ${worker_rc[index]}"
    unexplained=1
  fi

  printf '%-6s %-4s %-48s %-9s %s\n' \
    "$case_id" "$maxlevel" "${STACK_NAMES[$index]}" "$last_t" "$verdict"
done

exit "$unexplained"
