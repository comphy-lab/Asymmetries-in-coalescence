#!/bin/bash
#SBATCH --job-name=sbbInfRr
#SBATCH --account=nctth181
#SBATCH --partition=genoa
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=48
#SBATCH --time=72:00:00
#SBATCH --output=slurm-infrr-%j.out
#SBATCH --error=slurm-infrr-%j.err

# Matched-Delta r_d(Oh) / v_jet(Oh) ladder on the unified driver
# simulationCases/burstingBubbleInfiniteRr.c (R_r -> infinity, arbitrary chi).
#
# The driver is SELF-CONTAINED in this repository (it wraps
# coalescenceBubble.c, Stack 1: FILTERED + double projection), so unlike the
# retired drop-map launcher there is no cross-repository bundle to stage and
# no foreign-file hash guard. Only the Basilisk ref is pinned.
#
# Cases run at a fixed adaptive ceiling (drillAMR=0 in this solver's contour
# machinery), matching the group-D remediation of the drill-solver ladder:
# the drill ceiling silently released to level 8 in the Oh=0.0100 case and
# regionally capped the neck in the two diverged cases.
#
# WHY THIS LAUNCHER SUPERSEDES THE DRILL-SOLVER LADDER (2026-08-30)
#
# The 2026-08-28/29 campaign (jobs 26162496-98, 26165434/35) ran the drill
# solver `burstingBubble-drillResolution.c`, whose stack is
# `FILTERED` + `two-phase.h` + `navier-stokes/conserving.h` with a single
# projection. `coalescenceBubble.c` rejects that pairing at compile time as
# PROJECT POLICY (Basilisk itself compiles it): conserving.h reconstructs
# u = (q1+q2)/rho(f) with the sharp volume fraction while FILTERED gives the
# projection and viscous operators smoothed-sf properties -- one momentum
# balance, two densities -- and no upstream example or test pairs them.
# Separately, conserving.h at MAXlevel 13 runs with dt pinned ~13x below
# Stack 1's (measured, case 6554). The ladder ran at MAXlevel 13 and lost
# four of twelve points; the group-D
# remediation (drillAMR=0, i.e. full level 13 everywhere) rescued the
# under-resolved point but reproduced both divergences, which disproves the
# regional-under-resolution diagnosis and leaves the stack as the cause.
# Group B of this launcher crossed the same Oh at the same Delta on Stack 1
# without incident, so the whole ladder is re-measured here.
#
# The drill ladder also ran a gas-to-liquid viscosity ratio of 2e-2 through
# its own parameter file, against the 1e-2 the manuscript states. MuRin is now
# argv 26 so the value is explicit in argv and in `case.params`.
#
# Four cases x 48 threads = 192 CPUs = one full genoa node per group.
# MODE=preflight runs the same quartet to a short horizon first.

set -euo pipefail

MODE="${MODE:-production}"
GROUP="${GROUP:-A}"
# bo0-bulk is the existing Stack-1 Bo=0 ladder (6501-6520).
# bo0001-bulk / bo001-bulk / bo010-bulk repeat that Oh list at Bo=0.001,
# 0.01 and 0.1 (Young-Laplace Bo%5.4f.dat, reduced.h G.x=-Bo).
# bo001-bulk is the production Bond (experiments rounded up).
# bo001-wall is the same Oh list at Bo=0.01 with wallClearance=0.027, so the
# small-bubble south pole sits next to the wall (solver-derived zWall~0.02).
# bo001-*-ml14 is the Bond=0.01 mesh-sensitivity rung: same geometries and
# pinned detector, MAXlevel 14. Compare 6641-44 to 6601/6606/6614/6620 and
# 6645-48 to 6621/6626/6634/6640. Same Ldomain family, Delta halved; do not
# mix bulk and wall points on one Delta curve. Snellius drops prefix env, so
# SERIES and GROUP must be passed as `sbatch --export=ALL,SERIES=...,GROUP=...`.
SERIES="${SERIES:-bo0-bulk}"
MAXLEVEL="${MAXLEVEL:-13}"
if [[ "$SERIES" == "bo0-bulk" && "${SLURM_JOB_NAME:-}" == *Bo001* ]]; then
  echo "ERROR: job name ${SLURM_JOB_NAME} looks like a Bo=0.01 run but SERIES defaulted to bo0-bulk. Pass sbatch --export=ALL,SERIES=...,GROUP=..." >&2
  exit 2
fi
MURIN="${MURIN:-1e-2}"         # Oh_gas/Oh_liquid; 1e-2 is the manuscript value
THREADS="${SLURM_CPUS_PER_TASK:-48}"

case "$SERIES" in
  bo0-bulk)
    BOND="${BOND:-0}"
    ZWALL="${ZWALL:-4.0}"
    WALLCLEARANCE="${WALLCLEARANCE:--1}"
    case "$GROUP" in
      A) CASE_IDS=(6501 6502 6503 6504); OH_VALUES=(0.0100 0.0150 0.0200 0.0250) ;;
      B) CASE_IDS=(6505 6506 6507 6508); OH_VALUES=(0.0270 0.0280 0.0290 0.0300) ;;
      C) CASE_IDS=(6509 6510 6511 6512); OH_VALUES=(0.0310 0.0325 0.0333 0.0365) ;;
      # F closes the gap between C's 0.0333 and 0.0365, where the forward-injection
      # boundary Oh_c^(U) sits. The drill-solver bracket (0.0345, 0.0355] was
      # measured on the superseded stack and at MuRin=2e-2, so it is re-measured
      # rather than extended.
      F) CASE_IDS=(6513 6514 6515 6516); OH_VALUES=(0.0340 0.0345 0.0350 0.0355) ;;
      # G extends the boundary search upward: the completed 16-point Stack 1
      # ladder is forward-ejecting EVERYWHERE, with U_d = +17.9 still at
      # Oh = 0.0365, so Oh_c^(U) > 0.0365 on this stack at MuRin = 1e-2 and the
      # F-group window (designed around the retired drill-era bracket at 2e-2)
      # sits entirely below the boundary. The spread covers the plausible window
      # up to and past the Walls et al. (2015) low-Bo reading (~0.037-0.045 for
      # mu_g/mu_l = 1e-2..1e-3): if all four still eject forward, the boundary
      # exceeds 0.046 and that disagreement is itself the result.
      G) CASE_IDS=(6517 6518 6519 6520); OH_VALUES=(0.0375 0.0400 0.0430 0.0460) ;;
      # H/I extend Bo=0 past 0.0460, where every completed rung still ejects.
      # E remains the turning-point refinement quartet (6521-6524).
      H) CASE_IDS=(6531 6532 6533 6534); OH_VALUES=(0.0500 0.0550 0.0600 0.0650) ;;
      I) CASE_IDS=(6535 6536 6537 6538); OH_VALUES=(0.0700 0.0800 0.0900 0.1000) ;;
      # E is reserved for the refinement quartet around the measured turning
      # point; fill the Oh values in once the ladder above has located it.
      E) CASE_IDS=(6521 6522 6523 6524); OH_VALUES=("${OH_E1:?set OH_E1..OH_E4}" "${OH_E2:?}" "${OH_E3:?}" "${OH_E4:?}") ;;
      *) echo "ERROR: GROUP must be A, B, C, E, F, G, H or I" >&2; exit 2 ;;
    esac
    ;;
  bo0001-bulk)
    BOND="${BOND:-0.001}"
    ZWALL="${ZWALL:-4.0}"
    WALLCLEARANCE="${WALLCLEARANCE:--1}"
    case "$GROUP" in
      A) CASE_IDS=(6701 6702 6703 6704); OH_VALUES=(0.0100 0.0150 0.0200 0.0250) ;;
      B) CASE_IDS=(6705 6706 6707 6708); OH_VALUES=(0.0270 0.0280 0.0290 0.0300) ;;
      C) CASE_IDS=(6709 6710 6711 6712); OH_VALUES=(0.0310 0.0325 0.0333 0.0365) ;;
      F) CASE_IDS=(6713 6714 6715 6716); OH_VALUES=(0.0340 0.0345 0.0350 0.0355) ;;
      G) CASE_IDS=(6717 6718 6719 6720); OH_VALUES=(0.0375 0.0400 0.0430 0.0460) ;;
      *) echo "ERROR: SERIES=$SERIES GROUP must be A, B, C, F or G" >&2; exit 2 ;;
    esac
    ;;
  bo001-bulk)
    BOND="${BOND:-0.01}"
    ZWALL="${ZWALL:-4.0}"
    WALLCLEARANCE="${WALLCLEARANCE:--1}"
    case "$GROUP" in
      A) CASE_IDS=(6601 6602 6603 6604); OH_VALUES=(0.0100 0.0150 0.0200 0.0250) ;;
      B) CASE_IDS=(6605 6606 6607 6608); OH_VALUES=(0.0270 0.0280 0.0290 0.0300) ;;
      C) CASE_IDS=(6609 6610 6611 6612); OH_VALUES=(0.0310 0.0325 0.0333 0.0365) ;;
      F) CASE_IDS=(6613 6614 6615 6616); OH_VALUES=(0.0340 0.0345 0.0350 0.0355) ;;
      G) CASE_IDS=(6617 6618 6619 6620); OH_VALUES=(0.0375 0.0400 0.0430 0.0460) ;;
      *) echo "ERROR: SERIES=$SERIES GROUP must be A, B, C, F or G" >&2; exit 2 ;;
    esac
    ;;
  bo010-bulk)
    BOND="${BOND:-0.1}"
    ZWALL="${ZWALL:-4.0}"
    WALLCLEARANCE="${WALLCLEARANCE:--1}"
    case "$GROUP" in
      A) CASE_IDS=(6801 6802 6803 6804); OH_VALUES=(0.0100 0.0150 0.0200 0.0250) ;;
      B) CASE_IDS=(6805 6806 6807 6808); OH_VALUES=(0.0270 0.0280 0.0290 0.0300) ;;
      C) CASE_IDS=(6809 6810 6811 6812); OH_VALUES=(0.0310 0.0325 0.0333 0.0365) ;;
      F) CASE_IDS=(6813 6814 6815 6816); OH_VALUES=(0.0340 0.0345 0.0350 0.0355) ;;
      G) CASE_IDS=(6817 6818 6819 6820); OH_VALUES=(0.0375 0.0400 0.0430 0.0460) ;;
      *) echo "ERROR: SERIES=$SERIES GROUP must be A, B, C, F or G" >&2; exit 2 ;;
    esac
    ;;
  bo001-bulk-densify)
    BOND="${BOND:-0.01}"
    ZWALL="${ZWALL:-4.0}"
    WALLCLEARANCE="${WALLCLEARANCE:--1}"
    case "$GROUP" in
      D) CASE_IDS=(6651 6652 6653 6654); OH_VALUES=(0.0440 0.0445 0.0450 0.0455) ;;
      *) echo "ERROR: SERIES=$SERIES GROUP must be D" >&2; exit 2 ;;
    esac
    ;;
  bo001-wall)
    BOND="${BOND:-0.01}"
    # argv zWall is a placeholder: wallClearance>0 resets origin from the
    # Young-Laplace south pole so the film gap is the physical 0.027 R0.
    ZWALL="${ZWALL:-4.0}"
    WALLCLEARANCE="${WALLCLEARANCE:-0.027}"
    case "$GROUP" in
      A) CASE_IDS=(6621 6622 6623 6624); OH_VALUES=(0.0100 0.0150 0.0200 0.0250) ;;
      B) CASE_IDS=(6625 6626 6627 6628); OH_VALUES=(0.0270 0.0280 0.0290 0.0300) ;;
      C) CASE_IDS=(6629 6630 6631 6632); OH_VALUES=(0.0310 0.0325 0.0333 0.0365) ;;
      F) CASE_IDS=(6633 6634 6635 6636); OH_VALUES=(0.0340 0.0345 0.0350 0.0355) ;;
      G) CASE_IDS=(6637 6638 6639 6640); OH_VALUES=(0.0375 0.0400 0.0430 0.0460) ;;
      *) echo "ERROR: SERIES=$SERIES GROUP must be A, B, C, F or G" >&2; exit 2 ;;
    esac
    ;;
  bo001-bulk-ml14)
    # Same bulk half-space as bo001-bulk (Ldomain=10). Forcing MAXlevel=14 so
    # a leftover MAXLEVEL=13 in the batch environment cannot silently rerun
    # the production rung. Oh span: vigorous jet, Bo=0 Delta-anchor, near the
    # Bo=0 Oh_c^(U) window, and the Bo=0.01 probe that returned the jet.
    BOND="${BOND:-0.01}"
    MAXLEVEL=14
    ZWALL="${ZWALL:-4.0}"
    WALLCLEARANCE="${WALLCLEARANCE:--1}"
    case "$GROUP" in
      M) CASE_IDS=(6641 6642 6643 6644); OH_VALUES=(0.0100 0.0280 0.0345 0.0460) ;;
      *) echo "ERROR: SERIES=$SERIES GROUP must be M" >&2; exit 2 ;;
    esac
    ;;
  bo001-wall-ml14)
    BOND="${BOND:-0.01}"
    MAXLEVEL=14
    ZWALL="${ZWALL:-4.0}"
    WALLCLEARANCE="${WALLCLEARANCE:-0.027}"
    case "$GROUP" in
      M) CASE_IDS=(6645 6646 6647 6648); OH_VALUES=(0.0100 0.0280 0.0345 0.0460) ;;
      *) echo "ERROR: SERIES=$SERIES GROUP must be M" >&2; exit 2 ;;
    esac
    ;;
  *) echo "ERROR: SERIES must be bo0-bulk, bo0001-bulk, bo001-bulk, bo001-bulk-densify, bo010-bulk, bo001-wall, bo001-bulk-ml14 or bo001-wall-ml14" >&2; exit 2 ;;
esac

case "$MODE" in
  preflight)
    TMAX=0.15; SNAPINT=0.05
    for i in 0 1 2 3; do CASE_IDS[$i]=$(( ${CASE_IDS[$i]} + 80 )); done
    ;;
  production)
    TMAX="${TMAX:-1.5}"; SNAPINT="${SNAPINT:-0.05}"
    ;;
  *) echo "ERROR: MODE must be preflight or production" >&2; exit 2 ;;
esac

ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
cd "$ROOT"

# Canonical 25-argument solver contract; compiled defaults fill every key a
# ladder case does not override.
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

echo "mode=$MODE series=$SERIES group=$GROUP tmax=$TMAX snapshotInterval=$SNAPINT MAXlevel=$MAXLEVEL zWall=$ZWALL wallClearance=$WALLCLEARANCE Bond=$BOND MuRin=$MURIN threads_per_case=$THREADS"
echo "cases=${CASE_IDS[*]} oh=${OH_VALUES[*]}"

for header in adapt_wavelet_limited.h jetFoot.h; do
  if [[ ! -f "$ROOT/src-local/$header" ]]; then
    echo "ERROR: staged tree is missing src-local/$header" >&2
    exit 1
  fi
done

cd simulationCases

# The build output is named per job. Groups share one staging root, so two
# groups submitted together compile the same source into the same directory:
# with a fixed name, one job's `qcc` truncates the file the other job's `cp`
# is reading, and the second dies with "cannot open ... for reading". That is
# exactly how group C failed 32 s in on 2026-08-30 while A and F survived.
BINARY="burstingBubbleInfiniteRr-${SERIES}-${GROUP}-${SLURM_JOB_ID:-$$}"
qcc -I. -I../src-local -O2 -Wall -disable-dimensions -fopenmp \
  burstingBubbleInfiniteRr.c -o "$BINARY" -lm

pids=()
for index in 0 1 2 3; do
  case_id="${CASE_IDS[$index]}"
  oh="${OH_VALUES[$index]}"
  case_dir="$ROOT/simulationCases/$case_id"
  mkdir -p "$case_dir/intermediate"
  cp "$BINARY" "$case_dir/"
  ln -sfn "$ROOT/simulationCases/DataFiles" "$case_dir/DataFiles"

  # Build argv 1..25 through the canonical contract in
  # src-local/solver_args.sh: compiled defaults for everything the ladder
  # does not set, explicit overrides for what it does. The pinned detector
  # radius is NEVER 0 (0 makes the solver derive a mesh-tied value); Rr is a
  # placeholder in halfspace geometry.
  coalescence_solver_args solver_default_lookup \
    "OhOut=$oh" "RhoIn=1e-3" "Rr=1000" "MAXlevel=$MAXLEVEL" \
    "tmax=$TMAX" "zWall=$ZWALL" \
    "dropRadiusMin=0.021005127" "dropPersistence=3" \
    "snapshotInterval=$SNAPINT" "drillAMR=0" \
    "geometryMode=halfspace" "wallClearance=$WALLCLEARANCE" \
    "interfaceFloor=1" "MuRin=$MURIN" "Bond=$BOND"
  solver_args=("${COALESCENCE_SOLVER_ARGS[@]}")

  # Record the resolved argv beside the outputs. The drill ladder's gas
  # viscosity differed from the manuscript's and nothing in the run directory
  # said so; this file is the fix.
  {
    printf 'CaseNo=%s\nOh=%s\nBond=%s\nMuRin=%s\nMAXlevel=%s\ntmax=%s\nzWall=%s\n' \
      "$case_id" "$oh" "$BOND" "$MURIN" "$MAXLEVEL" "$TMAX" "$ZWALL"
    printf 'drillAMR=0\ngeometryMode=halfspace\nwallClearance=%s\n' "$WALLCLEARANCE"
    printf 'series=%s\n' "$SERIES"
    printf 'initialShape=Bo%5.4f.dat\n' "$BOND"
    printf 'dropRadiusMin=0.021005127\ndropPersistence=3\nsnapshotInterval=%s\n' "$SNAPINT"
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

# Sleep for $2 seconds, but wake as soon as the case exits. A plain
# `sleep 600` inside a backgrounded shell function cannot be interrupted --
# bash defers the signal until the foreground `sleep` returns -- so the
# watchdogs would hold the job (and its stdout) open for up to ten minutes
# after the last case finished.
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

# dt/progress watchdog: terminate a case only when t is unchanged AND
# dt < 1e-9 on three consecutive samples ten minutes apart. A converging run
# advances t every step and cannot trip this. The log's columns are
# `i dt t ...`, same positions as the drill solver's.
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
    # Unlike the drill-solver ladder there is no keStop gate and
    # DETACH_STOP=0 removes the classifier stop, so a healthy case has no
    # legitimate way to end early: a short stop is a projection-guard
    # failure, a crash or a watchdog kill, and fails the job.
    awk -v t="$last_t" -v h="$TMAX" 'BEGIN {exit !(t < h - 0.02)}' \
      && { echo "  ERROR: $case_id stopped at t=$last_t before its $TMAX horizon" >&2
           [[ -s "$d/projection_guard.status" ]] && head -2 "$d/projection_guard.status" >&2
           rc=1; }
  fi
  [[ -s "$d/WATCHDOG" ]] && { cat "$d/WATCHDOG"; rc=1; }
  [[ -s "$d/classification.status" ]] && head -2 "$d/classification.status"
done

exit "$rc"
