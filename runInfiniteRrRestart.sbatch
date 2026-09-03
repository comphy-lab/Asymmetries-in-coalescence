#!/bin/bash
#SBATCH --job-name=sbbRestart
#SBATCH --account=nctth181
#SBATCH --partition=genoa
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --time=72:00:00
#SBATCH --output=slurm-restart-%j.out
#SBATCH --error=slurm-restart-%j.err

# Pinned single-case continuation of an existing run directory.
#
# WHY THIS EXISTS
#
# The group launchers start fresh quartets only. When one case needs to run
# further -- the first user is 6564, the ML14 Delta-anchor rung, which hit
# walltime near t~0.5 with its burst window unfinished -- an improvised
# `sbatch --wrap` would bypass the recorded launcher/configuration contract
# and cannot satisfy the dispatch gate's "pinned code/configuration plus
# explicitly resumable run directory" requirement. This launcher IS that
# contract for continuations.
#
# HOW IT PINS THE CONFIGURATION
#
# It does not rebuild the argument list. It replays the EXACT argv recorded
# in the case's own `case.params` at launch time, with exactly one field
# replaced: argv 5, tmax. The compile flags are likewise taken from the
# recorded `qccFlags`. Nothing else is recomputable, so nothing else can
# drift. The solver's own init event restores the full state from the
# case's `restart` dump (checkpoints are self-contained: nodump cleared,
# `ufc` face-velocity mirror carried), and the TLOG appenders extend rather
# than truncate the histories.
#
# SUBMISSION CONTRACT
#
#   CASE=6564 TMAX=0.80 sbatch --dependency=afterany:<owning-job-id> \
#     runInfiniteRrRestart.sbatch
#
# The dependency is mandatory whenever the owning job could still be
# running: two writers on one `restart` file corrupt the continuation. The
# in-script mtime guard below is a backstop, not a substitute.

set -euo pipefail

CASE="${CASE:?set CASE, e.g. CASE=6564}"
TMAX="${TMAX:?set TMAX, the new horizon, e.g. TMAX=0.80}"
# The continuation is a success once the physical burst window is complete,
# even if the last stretch to TMAX hits walltime. At ML13 the first
# confirmation lands near t~0.57 and the 0.15 burst window closes by ~0.72.
BURST_OK="${BURST_OK:-0.72}"
THREADS="${SLURM_CPUS_PER_TASK:-48}"

ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
cd "$ROOT"

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

case_dir="$ROOT/simulationCases/$CASE"
params="$case_dir/case.params"

# ---- resumability gates ---------------------------------------------------
if [[ ! -f "$params" ]]; then
  echo "ERROR: $params missing; this run directory is not resumable under the contract" >&2
  exit 1
fi
if [[ ! -s "$case_dir/restart" ]]; then
  echo "ERROR: $case_dir/restart missing or empty; nothing to restore from" >&2
  exit 1
fi
# Backstop against a still-running owner (the dependency is the real guard):
# a restart file younger than ten minutes is still being written to.
if [[ -n "$(find "$case_dir/restart" -mmin -10 2>/dev/null)" ]]; then
  echo "ERROR: $case_dir/restart modified within the last 10 minutes;" >&2
  echo "       the owning job may still be running. Submit with" >&2
  echo "       --dependency=afterany:<owning-job-id>." >&2
  exit 1
fi

# `|| true`: under pipefail a missing argv= line would kill the script at
# this pipeline with no diagnostic, and the intended error below would
# never print.
recorded_argv_line=$(grep -m1 '^argv=' "$params" | cut -d= -f2- || true)
if [[ -z "$recorded_argv_line" ]]; then
  echo "ERROR: no argv= record in $params; cannot pin the configuration" >&2
  exit 1
fi
read -r -a solver_args <<< "$recorded_argv_line"
if [[ "${#solver_args[@]}" -ne 26 && "${#solver_args[@]}" -ne 27 ]]; then
  echo "ERROR: recorded argv has ${#solver_args[@]} fields, expected 26 or 27" >&2
  exit 1
fi

old_tmax="${solver_args[4]}"
if ! awk -v new="$TMAX" -v old="$old_tmax" 'BEGIN {exit !(new > old)}'; then
  echo "ERROR: TMAX=$TMAX does not extend the recorded tmax=$old_tmax" >&2
  exit 1
fi
solver_args[4]="$TMAX"

flags=$(grep -m1 '^qccFlags=' "$params" | cut -d= -f2- || true)
stack=$(grep -m1 '^solverStack=' "$params" | cut -d= -f2- || true)

echo "case=$CASE stack='${stack:-filtered+double-projection}' qccFlags='$flags'"
echo "resuming from restart (last dump), tmax $old_tmax -> $TMAX, burst-ok gate $BURST_OK"
echo "argv: ${solver_args[*]}"

# ---- build (recorded flags, per-job binary) --------------------------------
cd simulationCases
binary="burstingBubbleInfiniteRr-$CASE-restart-${SLURM_JOB_ID:-$$}"
# shellcheck disable=SC2086
qcc -I. -I../src-local -O2 -Wall -disable-dimensions -fopenmp $flags \
  burstingBubbleInfiniteRr.c -o "$binary" -lm
cp "$binary" "$case_dir/"

# Continuation provenance under DISTINCT keys: duplicating tmax= in
# case.params would make its value parser-dependent (first-match awk vs
# last-wins dict have already coexisted in this project).
{
  printf 'restartJob=%s\n' "${SLURM_JOB_ID:-local}"
  printf 'restartAt=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf 'restartTmax=%s\n' "$TMAX"
  printf 'restartArgv=%s\n' "${solver_args[*]}"
} >> "$params"

# ---- run -------------------------------------------------------------------
cd "$case_dir"
export OMP_NUM_THREADS="$THREADS"
export OMP_PROC_BIND=close
export OMP_PLACES=cores
set +e
srun --exact -n1 -c "$THREADS" --cpu-bind=cores \
  "./$binary" "${solver_args[@]}" >> runner.out 2>> runner.err
worker=$?
set -e

# ---- verdict ---------------------------------------------------------------
last_t="nodata"
[[ -s jet.log ]] && last_t=$(tail -1 jet.log | cut -d, -f1)
echo
echo "=== restart summary: case $CASE reached t=$last_t (target $TMAX) ==="

if grep -q 'kinetic energy blew up' runner.err 2>/dev/null; then
  echo "verdict: DIVERGED (ke gate) — a result, not a launcher failure"
  exit 0
fi
if [[ "$last_t" != "nodata" ]] && \
   awk -v t="$last_t" -v h="$TMAX" 'BEGIN {exit !(t >= h - 0.02)}'; then
  echo "verdict: reached horizon"
  exit 0
fi
if [[ "$last_t" != "nodata" ]] && \
   awk -v t="$last_t" -v b="$BURST_OK" 'BEGIN {exit !(t >= b)}'; then
  echo "verdict: stopped short of $TMAX but past the burst-ok gate $BURST_OK — usable"
  exit 0
fi
echo "verdict: FAILED — t=$last_t is short of the burst-ok gate $BURST_OK (worker exit $worker)"
exit 1
