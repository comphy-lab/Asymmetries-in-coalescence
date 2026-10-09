#!/bin/bash
#SBATCH --job-name=dijCamp
#SBATCH --account=nctth181
#SBATCH --partition=genoa
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=5-00:00:00
#SBATCH --output=slurm-camp-%j.out
#SBATCH --error=slurm-camp-%j.err

# Packed, self-resuming campaign launcher for burstingBubbleInfiniteRr.c.
#
# WHY THIS EXISTS (2026-10-09)
#
# The group launcher runBurstingBubbleInfiniteRr.sbatch hard-codes four cases
# x 48 OpenMP threads per node and fresh starts only. ML14 cases need more
# than the 5-day genoa limit, so every continuation went through the
# single-case restart launcher by hand, and the September campaigns finished
# in scratch-shared, which SURF purged before anyone harvested them. This
# launcher:
#
# - reads a manifest (one case per line, its own thread count), so a node can
#   be packed with any number of cases that fit in 192 cores;
# - resumes every unfinished case from its own `restart` dump when the same
#   manifest is resubmitted (chain jobs with --dependency=afterany), and skips
#   finished ones, so a 5-day wall is a checkpoint, not a failure;
# - stops each solver 20 minutes before the wall so the job ends cleanly;
# - refuses to restore a `restart` that a kill may have truncated: `restart`
#   and `intermediate/snapshot-<t>` are written back to back at each snapshot,
#   so a size mismatch means the pair is incomplete and the previous snapshot
#   is restored instead. Every resume is appended to RESUMES (job, restore t),
#   because the TLOG appenders extend the histories and rows after the restore
#   time are repeated once; post-processing keeps the last segment.
#
# MANIFEST columns (whitespace separated, '#' comments allowed):
#   case_id series bond oh maxlevel zwall wallclearance ldomain threads
# ldomain=10 is the solver default; any other value is exported as LDOMAIN
# (far-floor half-space). Sum of threads must not exceed 192.
#
# Physics, detector and argv contract are identical to the group launcher:
# Stack 1 (FILTERED + double projection), pinned dropRadiusMin=0.021005127,
# dropPersistence=3, MuRin=1e-2, drillAMR=0, halfspace geometry.
#
# SUBMISSION
#   sbatch --export=ALL,MANIFEST=campaign.tsv runInfiniteRrCampaign.sbatch
#   sbatch --dependency=afterany:<jobid> --export=ALL,MANIFEST=campaign.tsv \
#     runInfiniteRrCampaign.sbatch
# PILOT_SECONDS=<s> caps every case at that wall time (thread-scaling pilot)
# and samples progress every 20 s into pilot_samples.tsv.

set -euo pipefail

MANIFEST="${MANIFEST:?set MANIFEST, e.g. MANIFEST=campaign.tsv}"
TMAX="${TMAX:-1.5}"
SNAPINT="${SNAPINT:-0.05}"
MURIN="${MURIN:-1e-2}"
STOP_MARGIN="${STOP_MARGIN:-1200}"
PILOT_SECONDS="${PILOT_SECONDS:-0}"
JOB="${SLURM_JOB_ID:-$$}"

ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
cd "$ROOT"
[[ -f "$MANIFEST" ]] || { echo "ERROR: manifest $MANIFEST not found in $ROOT" >&2; exit 2; }

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
if ! grep -qx "ref=${EXPECTED_BASILISK_REF}" "${BASILISK_ROOT}/.comphy-lock" 2>/dev/null; then
  echo "ERROR: Basilisk install is not ref-locked to ${EXPECTED_BASILISK_REF}" >&2
  exit 1
fi

# ---- read the manifest -------------------------------------------------------
CASES=(); SERIESES=(); BONDS=(); OHS=(); LEVELS=(); ZWALLS=(); CLEARS=(); LDOMS=(); THREADS=()
total_threads=0
while read -r c s b o l z w d n _; do
  [[ -z "${c:-}" || "$c" == \#* ]] && continue
  [[ -n "${n:-}" ]] || { echo "ERROR: manifest line for case $c has fewer than 9 columns" >&2; exit 2; }
  CASES+=("$c"); SERIESES+=("$s"); BONDS+=("$b"); OHS+=("$o"); LEVELS+=("$l")
  ZWALLS+=("$z"); CLEARS+=("$w"); LDOMS+=("$d"); THREADS+=("$n")
  total_threads=$(( total_threads + n ))
done < "$MANIFEST"
(( ${#CASES[@]} > 0 )) || { echo "ERROR: manifest $MANIFEST lists no cases" >&2; exit 2; }
if (( total_threads > ${SLURM_CPUS_ON_NODE:-192} )); then
  echo "ERROR: manifest asks for $total_threads threads on a ${SLURM_CPUS_ON_NODE:-192}-core node" >&2
  exit 2
fi

# ---- wall-clock budget per solver --------------------------------------------
now=$(date +%s)
if [[ -n "${SLURM_JOB_END_TIME:-}" ]]; then
  end_epoch="$SLURM_JOB_END_TIME"
else
  left=$(squeue -h -j "$JOB" -o %L 2>/dev/null || echo "5-00:00:00")
  end_epoch=$(( now + $(awk -F'[-:]' '{ if (NF==4) print $1*86400+$2*3600+$3*60+$4;
    else if (NF==3) print $1*3600+$2*60+$3; else print $1*60+$2 }' <<< "$left") ))
fi
RUN_SECONDS=$(( end_epoch - now - STOP_MARGIN ))
if (( PILOT_SECONDS > 0 && PILOT_SECONDS < RUN_SECONDS )); then RUN_SECONDS=$PILOT_SECONDS; fi
(( RUN_SECONDS > 600 )) || { echo "ERROR: only ${RUN_SECONDS}s of wall left" >&2; exit 1; }

SOURCE_COMMIT=$(cat "$ROOT/SOURCE_COMMIT.txt" 2>/dev/null || echo unknown)
echo "job=$JOB manifest=$MANIFEST cases=${#CASES[@]} threads=$total_threads tmax=$TMAX run_seconds=$RUN_SECONDS source=$SOURCE_COMMIT pilot=$PILOT_SECONDS"

cd simulationCases
BINARY="burstingBubbleInfiniteRr-camp-${JOB}"
qcc -I. -I../src-local -O2 -Wall -disable-dimensions -fopenmp \
  burstingBubbleInfiniteRr.c -o "$BINARY" -lm

last_t_of() {  # last physical time in jet.log, or empty
  [[ -s "$1/jet.log" ]] || return 0
  tail -1 "$1/jet.log" | cut -d, -f1
}

param_of() {  # value of KEY= in case.params
  grep -m1 "^$2=" "$1/case.params" 2>/dev/null | cut -d= -f2- || true
}

# Choose the dump to restore from; echoes the restore time ("fresh" if none).
prepare_restart() {
  local dir="$1" snaps latest prev rsize lsize
  mapfile -t snaps < <(ls "$dir/intermediate" 2>/dev/null | grep '^snapshot-' \
    | sort -t- -k2 -g | sed "s|^|$dir/intermediate/|")
  if (( ${#snaps[@]} == 0 )); then
    rm -f "$dir/restart"
    echo fresh; return 0
  fi
  latest="${snaps[-1]}"
  rsize=$(stat -c %s "$dir/restart" 2>/dev/null || echo 0)
  lsize=$(stat -c %s "$latest")
  if (( rsize > 0 && rsize == lsize )); then
    echo "${latest##*snapshot-}"; return 0
  fi
  if (( ${#snaps[@]} >= 2 )); then
    prev="${snaps[-2]}"
    cp "$prev" "$dir/restart"
    echo "${prev##*snapshot-}"; return 0
  fi
  rm -f "$dir/restart"
  echo fresh
}

pids=(); run_dirs=(); next_core=0; drift_skipped=0
for i in "${!CASES[@]}"; do
  case_id="${CASES[$i]}"; series="${SERIESES[$i]}"; bond="${BONDS[$i]}"
  oh="${OHS[$i]}"; level="${LEVELS[$i]}"; zwall="${ZWALLS[$i]}"
  clear="${CLEARS[$i]}"; ldom="${LDOMS[$i]}"; nthreads="${THREADS[$i]}"
  case_dir="$ROOT/simulationCases/$case_id"

  if [[ -f "$case_dir/DONE" ]]; then
    echo "case=$case_id already DONE; skipped"; continue
  fi
  last_t=$(last_t_of "$case_dir")
  if [[ -n "$last_t" ]] && awk -v t="$last_t" -v h="$TMAX" 'BEGIN {exit !(t >= h - 0.02)}'; then
    echo "t=$last_t" > "$case_dir/DONE"
    echo "case=$case_id reached t=$last_t; marked DONE"; continue
  fi

  mkdir -p "$case_dir/intermediate"
  ln -sfn "$ROOT/simulationCases/DataFiles" "$case_dir/DataFiles"
  cp "$BINARY" "$case_dir/"

  if [[ -f "$case_dir/case.params" ]]; then
    # Continuation: the recorded argv and source are the contract; refuse any
    # drift, including a run root restaged at another commit.
    for pair in "Oh=$oh" "Bond=$bond" "MAXlevel=$level" "Ldomain=$ldom" \
                "sourceCommit=$SOURCE_COMMIT"; do
      key="${pair%%=*}"; want="${pair#*=}"
      have=$(param_of "$case_dir" "$key")
      if [[ "$have" != "$want" || "$want" == unknown ]]; then
        echo "ERROR: case $case_id records $key=${have:-<none>} but this job has $want; not resumed" >&2
        drift_skipped=1
        continue 2
      fi
    done
    read -r -a solver_args <<< "$(param_of "$case_dir" argv)"
    restore_t=$(prepare_restart "$case_dir")
    echo "job=$JOB restore=$restore_t" >> "$case_dir/RESUMES"
    echo "case=$case_id resume from t=$restore_t threads=$nthreads"
  else
    coalescence_solver_args solver_default_lookup \
      "OhOut=$oh" "RhoIn=1e-3" "Rr=1000" "MAXlevel=$level" \
      "tmax=$TMAX" "zWall=$zwall" \
      "dropRadiusMin=0.021005127" "dropPersistence=3" \
      "snapshotInterval=$SNAPINT" "drillAMR=0" \
      "geometryMode=halfspace" "wallClearance=$clear" \
      "interfaceFloor=1" "MuRin=$MURIN" "Bond=$bond"
    solver_args=("${COALESCENCE_SOLVER_ARGS[@]}")
    {
      printf 'CaseNo=%s\nOh=%s\nBond=%s\nMuRin=%s\nMAXlevel=%s\ntmax=%s\nzWall=%s\n' \
        "$case_id" "$oh" "$bond" "$MURIN" "$level" "$TMAX" "$zwall"
      printf 'Ldomain=%s\ndrillAMR=0\ngeometryMode=halfspace\nwallClearance=%s\n' "$ldom" "$clear"
      printf 'series=%s\ninitialShape=Bo%5.4f.dat\n' "$series" "$bond"
      printf 'dropRadiusMin=0.021005127\ndropPersistence=3\nsnapshotInterval=%s\n' "$SNAPINT"
      printf 'solverStack=filtered+double-projection\nqccFlags=\nsourceCommit=%s\n' "$SOURCE_COMMIT"
      printf 'launcher=runInfiniteRrCampaign.sbatch\nargv=%s\n' "${solver_args[*]}"
    } > "$case_dir/case.params"
    echo "case=$case_id fresh Oh=$oh Bond=$bond ML=$level Ldomain=$ldom threads=$nthreads"
  fi

  # Each case owns a contiguous core block (genoa has no SMT; logical CPU k is
  # core k), pinned with taskset and mirrored in OMP_PLACES so the threads
  # cannot migrate onto a neighbour's cores.
  lo=$next_core; hi=$(( next_core + nthreads - 1 )); next_core=$(( hi + 1 ))
  echo "  cores $lo-$hi" >> "$case_dir/RESUMES"
  (
    cd "$case_dir"
    export OMP_NUM_THREADS="$nthreads" OMP_PROC_BIND=close OMP_PLACES="{$lo}:$nthreads"
    if [[ "$ldom" != "10" ]]; then export LDOMAIN="$ldom"; else unset LDOMAIN; fi
    taskset -c "$lo-$hi" timeout --signal=KILL "$RUN_SECONDS" \
      "./$BINARY" "${solver_args[@]}" > "runner-$JOB.out" 2> "runner-$JOB.err" &
    echo $! > .srun.pid
    wait $! || true
    rm -f .srun.pid
  ) &
  pids+=("$!"); run_dirs+=("$case_dir")
done

if (( ${#pids[@]} == 0 )); then
  if (( drift_skipped )); then
    echo "ERROR: nothing launched and at least one case was refused for drift" >&2; exit 3
  fi
  echo "nothing to run: every manifest case is DONE"; exit 0
fi

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

# Same rule as the group launcher: kill only when t is frozen and dt < 1e-9 on
# three consecutive samples ten minutes apart.
watchdog() {
  local dir="$1" last_t="" stalls=0 row t dt pid
  while :; do
    watchdog_nap "$dir" 600 || return 0
    pid=$(cat "$dir/.srun.pid" 2>/dev/null) || true
    [[ -n "${pid:-}" ]] && kill -0 "$pid" 2>/dev/null || return 0
    row=$(grep -E '^[0-9]+ ' "$dir/log" 2>/dev/null | tail -1) || continue
    [[ -n "$row" ]] || continue
    t=$(awk '{print $3}' <<< "$row"); dt=$(awk '{print $2}' <<< "$row")
    if [[ "$t" == "$last_t" ]] && awk -v d="$dt" 'BEGIN{exit !(d+0 < 1e-9)}'; then
      stalls=$((stalls + 1))
      echo "WATCHDOG: $(basename "$dir") no progress at t=$t dt=$dt (${stalls}/3)" >&2
      if (( stalls >= 3 )); then
        echo "watchdog=stalled t=$t dt=$dt job=$JOB" > "$dir/WATCHDOG"
        kill "$pid" 2>/dev/null || true
        return 0
      fi
    else
      stalls=0
    fi
    last_t="$t"
  done
}

sampler() {  # pilot only: wall-clock progress of every case every 20 s
  local out="$ROOT/pilot_samples-$JOB.tsv" d row
  printf 'epoch\tcase\ti\tt\n' > "$out"
  while :; do
    sleep 20
    for d in "${run_dirs[@]}"; do
      row=$(grep -E '^[0-9]+ ' "$d/log" 2>/dev/null | tail -1) || true
      [[ -n "$row" ]] && printf '%s\t%s\t%s\t%s\n' "$(date +%s)" "$(basename "$d")" \
        "$(awk '{print $1}' <<< "$row")" "$(awk '{print $3}' <<< "$row")" >> "$out"
    done
  done
}

wpids=()
for d in "${run_dirs[@]}"; do watchdog "$d" & wpids+=("$!"); done
if (( PILOT_SECONDS > 0 )); then sampler & wpids+=("$!"); fi

for pid in "${pids[@]}"; do wait "$pid" || true; done
for wpid in "${wpids[@]}"; do kill "$wpid" 2>/dev/null || true; done

rc=0
for d in "${run_dirs[@]}"; do
  case_id=$(basename "$d")
  last_t=$(last_t_of "$d")
  echo "=== case=$case_id reached t=${last_t:-none} ==="
  if [[ -n "$last_t" ]] && awk -v t="$last_t" -v h="$TMAX" 'BEGIN {exit !(t >= h - 0.02)}'; then
    echo "t=$last_t job=$JOB" > "$d/DONE"
  elif (( PILOT_SECONDS == 0 )); then
    echo "  unfinished: resubmit the same manifest to continue"
  fi
  [[ -s "$d/WATCHDOG" ]] && { cat "$d/WATCHDOG"; rc=1; }
  [[ -s "$d/projection_guard.status" ]] && { head -2 "$d/projection_guard.status"; rc=1; }
  [[ -s "$d/classification.status" ]] && sed -n 2p "$d/classification.status"
  tail -1 "$d/runner-$JOB.out" 2>/dev/null | grep '^#' || true
done
(( drift_skipped )) && { echo "ERROR: at least one case was refused for drift" >&2; rc=1; }
exit "$rc"
