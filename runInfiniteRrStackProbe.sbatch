#!/bin/bash
#SBATCH --job-name=sbbStackProbe
#SBATCH --account=nctth181
#SBATCH --partition=genoa
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=48
#SBATCH --time=36:00:00
#SBATCH --output=slurm-stackprobe-%j.out
#SBATCH --error=slurm-stackprobe-%j.err

# Controlled discriminator for the drop-map ladder's two divergences.
#
# THE OBSERVATION
#
# On the drill solver (`burstingBubble-drillResolution.c`, stack
# `FILTERED` + `two-phase.h` + `navier-stokes/conserving.h`, single
# projection) the bulk half-space cases Oh = 0.0280 and Oh = 0.0325 both die
# shortly after the capillary-wave focusing event, at MAXlevel 13 and matched
# Delta = 1.2207e-3:
#
#   6326  Oh=0.0280  focusing at t~0.4865, ke 5.32 -> 1.09e6 over 38 steps
#                    ending t=0.50706, with r_base = 0.159 (~130 cells)
#   6330  Oh=0.0325  focusing at t~0.487,  ke -> 1.15e6 ending t=0.49394,
#                    with r_base = 0.0409 (~33 cells)
#
# Their immediate neighbours (0.0270, 0.0290, 0.0310, 0.0333) are healthy, so
# the failure is not monotone in Oh. Group D re-ran both at drillAMR=0 -- full
# level 13 everywhere, no regional ceiling -- and both died again in the same
# window, which disproves the regional-under-resolution diagnosis. The neck is
# not thin at death; the blow-up is abrupt and local, and no projection guard
# or solver warning fires.
#
# On this repository's Stack 1 (`FILTERED` + `two-phase.h` +
# `navier-stokes/double-projection.h`) the same Oh at the same Delta runs
# through the window with dt pinned at its ceiling and pinches a drop at
# t ~ 0.61 (crossdriver case 6506).
#
# WHAT THIS JOB SEPARATES
#
# The two campaigns differ in two ways at once, so the comparison above does
# not by itself assign a cause:
#
#   1. the solver stack, and
#   2. the gas-to-liquid viscosity ratio -- 2e-2 on the drill ladder against
#      the 1e-2 the manuscript states and this solver defaults to.
#
# Four cases, one node, horizon t=0.75. That clears the death window at
# t ~ 0.49-0.51 and the first pinch at t ~ 0.61, so the job answers both
# "does it survive" and "does it pinch at the same time". Measured on the
# running crossdriver group B, 48 threads per case advance ~0.043 t/hour at
# this Oh and MAXlevel, so t=0.75 is about 17 h; the 36 h walltime is the
# margin, not the estimate.
#
#   6551  Oh=0.0325  Stack 1   MuRin=2e-2   the decisive case: identical
#                                           physics to the dead 6330, only
#                                           the stack differs
#   6552  Oh=0.0280  Stack 1   MuRin=2e-2   the same for 6326
#   6553  Oh=0.0280  Stack 1S  MuRin=1e-2   is double projection load-bearing?
#   6554  Oh=0.0280  Stack 2   MuRin=1e-2   is `conserving.h` itself the
#                                           problem at MAXlevel 13, or only
#                                           its illegal pairing with FILTERED?
#
# 6552 also pairs with the running crossdriver 6506 (Oh=0.0280, Stack 1,
# MuRin=1e-2) to give the gas-viscosity sensitivity at fixed stack, which is
# needed before any drill-era number is compared with a Stack 1 number.
#
# Expected outcomes and what each would mean:
#   - 6551 and 6552 survive  -> the stack is the cause; the ladder re-run on
#     Stack 1 is the correct remediation and no drill-era point is salvageable.
#   - 6551 or 6552 dies      -> gas viscosity, not the stack, carries the
#     failure; Stack 1 at MuRin=1e-2 is merely a milder problem and the
#     re-measured ladder needs a resolution study before it is trusted.
#   - 6554 dies and 6553 survives -> `conserving.h` at MAXlevel 13, exactly as
#     OPERATIONAL-NOTES.md section 2 already recorded.

set -euo pipefail

MODE="${MODE:-production}"
MAXLEVEL="${MAXLEVEL:-13}"
ZWALL="${ZWALL:-4.0}"          # Ldomain = min(zWall+6, 16) = 10 -> matched Delta
WALLCLEARANCE="${WALLCLEARANCE:--1}"
THREADS="${SLURM_CPUS_PER_TASK:-48}"

CASE_IDS=(6551      6552      6553      6554)
OH_VALUES=(0.0325   0.0280    0.0280    0.0280)
MURIN_VALUES=(2e-2  2e-2      1e-2      1e-2)
STACK_FLAGS=(""     ""        "-DSINGLE_PROJECTION=1" "-DUSE_CONSERVING=1")
STACK_NAMES=("filtered+double-projection" "filtered+double-projection" \
             "filtered+single-projection" "conserving")

case "$MODE" in
  preflight)
    TMAX=0.15; SNAPINT=0.05
    for i in 0 1 2 3; do CASE_IDS[$i]=$(( ${CASE_IDS[$i]} + 80 )); done
    ;;
  production)
    TMAX="${TMAX:-0.75}"; SNAPINT="${SNAPINT:-0.05}"
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
test -x "${BASILISK}/qcc"
test -f "${BASILISK_ROOT}/.comphy-lock"
grep -qx "ref=${EXPECTED_BASILISK_REF}" "${BASILISK_ROOT}/.comphy-lock"

echo "mode=$MODE tmax=$TMAX MAXlevel=$MAXLEVEL zWall=$ZWALL threads_per_case=$THREADS"
echo "cases=${CASE_IDS[*]}"

cd simulationCases

pids=()
for index in 0 1 2 3; do
  case_id="${CASE_IDS[$index]}"
  oh="${OH_VALUES[$index]}"
  murin="${MURIN_VALUES[$index]}"
  flags="${STACK_FLAGS[$index]}"
  stack="${STACK_NAMES[$index]}"
  # Per case *and* per job: the case ids are unique across the design,
  # but two probe jobs sharing a staging root would still collide.
  binary="burstingBubbleInfiniteRr-$case_id-${SLURM_JOB_ID:-$$}"
  case_dir="$ROOT/simulationCases/$case_id"

  # One binary per case: the stacks are mutually exclusive at compile time,
  # so they cannot share an executable.
  # shellcheck disable=SC2086
  qcc -I. -I../src-local -O2 -Wall -disable-dimensions -fopenmp $flags \
    burstingBubbleInfiniteRr.c -o "$binary" -lm

  mkdir -p "$case_dir/intermediate"
  cp "$binary" "$case_dir/"
  ln -sfn "$ROOT/simulationCases/DataFiles" "$case_dir/DataFiles"

  coalescence_solver_args solver_default_lookup \
    "OhOut=$oh" "RhoIn=1e-3" "Rr=1000" "MAXlevel=$MAXLEVEL" \
    "tmax=$TMAX" "zWall=$ZWALL" \
    "dropRadiusMin=0.021005127" "dropPersistence=3" \
    "snapshotInterval=$SNAPINT" "drillAMR=0" \
    "geometryMode=halfspace" "wallClearance=$WALLCLEARANCE" \
    "interfaceFloor=1" "MuRin=$murin"
  solver_args=("${COALESCENCE_SOLVER_ARGS[@]}")

  {
    printf 'CaseNo=%s\nOh=%s\nMuRin=%s\nMAXlevel=%s\ntmax=%s\nzWall=%s\n' \
      "$case_id" "$oh" "$murin" "$MAXLEVEL" "$TMAX" "$ZWALL"
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

# A case that dies the way 6326 and 6330 died exits on its own through the ke
# gate, so this job needs no watchdog: the question it asks is precisely
# "does it die", and both answers are recorded.
rc=0
for pid in "${pids[@]}"; do
  wait "$pid" || rc=1
done

echo
echo "=== discriminator summary ==="
printf '%-6s %-8s %-7s %-28s %-9s %s\n' case Oh MuRin stack reached verdict
for index in 0 1 2 3; do
  case_id="${CASE_IDS[$index]}"
  d="$ROOT/simulationCases/$case_id"
  last_t="nodata"
  [[ -s "$d/jet.log" ]] && last_t=$(tail -1 "$d/jet.log" | cut -d, -f1)
  verdict="reached horizon"
  if grep -q 'kinetic energy blew up' "$d/runner.err" 2>/dev/null; then
    verdict="DIVERGED (ke gate)"
  elif [[ "$last_t" == "nodata" ]] || \
       awk -v t="$last_t" -v h="$TMAX" 'BEGIN {exit !(t < h - 0.02)}'; then
    verdict="stopped short"
  fi
  printf '%-6s %-8s %-7s %-28s %-9s %s\n' \
    "$case_id" "${OH_VALUES[$index]}" "${MURIN_VALUES[$index]}" \
    "${STACK_NAMES[$index]}" "$last_t" "$verdict"
done

# The job's own exit code is deliberately not tied to a case diverging: a
# divergence here is a result, not a failure.
exit 0
