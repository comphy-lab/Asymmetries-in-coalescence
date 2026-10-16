#!/bin/bash
#SBATCH --job-name=dijCoal
#SBATCH --account=nctth181
#SBATCH --partition=genoa
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --time=5-00:00:00
#SBATCH --output=slurm-camp-%j.out
#SBATCH --error=slurm-camp-%j.err

# Packed, self-resuming campaign launcher for both coalescence geometries:
#   halfspace  burstingBubbleInfiniteRr.c   (R_r = infinity, Bo0-style IC)
#   finite     coalescenceBubbleFiniteRr.c  (two spheres, R_r finite)
#
# Generalises runInfiniteRrCampaign.sbatch (2026-10-09) for the two-group
# campaign of 2026-10-16 (scratch: experiments-bond-review/CAMPAIGN-DESIGN):
# the restart, watchdog, wall-margin and DONE logic is unchanged; what is new
#
# - a geometry column, and the real R_r for finite cases;
# - one initial condition per case, named in DataFiles/ and passed with
#   INITIAL_SHAPE, its sha256 recorded in case.params and checked on resume;
# - LDOMAIN always exported (no "10" sentinel), so every case states its
#   domain; keep Delta = 10/2^14 with L0 = 10*2^k at MAXlevel 14+k;
# - the far-field refinement floor (FAR_FLOOR_RADIUS / FAR_FLOOR_DROP, 0 0 =
#   off) per case;
# - both wrappers write components.log / jet.log (src-local/dropCensus.h), so
#   DONE is decided from jet.log for either geometry.
#
# MANIFEST columns (whitespace separated, '#' comments allowed):
#   case_id series geometry rr bond oh maxlevel zwall wallclearance ldomain
#   ic_file farfloor_radius farfloor_drop min_level threads
# min_level 0 = solver default (MAXlevel-6); farfloor 0 0 = off.
# rr is ignored for halfspace (passed as the 1000 placeholder). Sum of threads
# must not exceed 192.
#
# Physics and detector: Stack 1 (FILTERED + double projection), pinned
# dropRadiusMin=0.021005127, dropPersistence=3, MuRin=1e-2, drillAMR=0,
# interfaceFloor=1, IC band refinement on (solver default).
#
# SUBMISSION
#   sbatch --export=ALL,MANIFEST=campaign.tsv runCoalescenceCampaign.sbatch
#   sbatch --dependency=afterany:<jobid> --export=ALL,MANIFEST=campaign.tsv \
#     runCoalescenceCampaign.sbatch
# PILOT_SECONDS=<s> caps every case at that wall time (scaling pilot).

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
CASES=(); SERIESES=(); GEOMS=(); RRS=(); BONDS=(); OHS=(); LEVELS=(); ZWALLS=(); CLEARS=()
LDOMS=(); ICS=(); FFRS=(); FFDS=(); MINLS=(); THREADS=()
total_threads=0; need_half=0; need_finite=0
while read -r c s g r b o l z w d ic ffr ffd ml n _; do
  [[ -z "${c:-}" || "$c" == \#* ]] && continue
  [[ -n "${n:-}" ]] || { echo "ERROR: manifest line for case $c has fewer than 15 columns" >&2; exit 2; }
  case "$g" in
    halfspace) need_half=1 ;;
    finite)    need_finite=1 ;;
    *) echo "ERROR: case $c geometry must be halfspace or finite, got $g" >&2; exit 2 ;;
  esac
  [[ -s "$ROOT/simulationCases/DataFiles/$ic" ]] || { echo "ERROR: case $c IC DataFiles/$ic missing" >&2; exit 2; }
  CASES+=("$c"); SERIESES+=("$s"); GEOMS+=("$g"); RRS+=("$r"); BONDS+=("$b"); OHS+=("$o")
  LEVELS+=("$l"); ZWALLS+=("$z"); CLEARS+=("$w"); LDOMS+=("$d"); ICS+=("$ic")
  FFRS+=("$ffr"); FFDS+=("$ffd"); MINLS+=("$ml"); THREADS+=("$n")
  total_threads=$(( total_threads + n ))
done < "$MANIFEST"
(( ${#CASES[@]} > 0 )) || { echo "ERROR: manifest $MANIFEST lists no cases" >&2; exit 2; }
if (( total_threads > ${SLURM_CPUS_ON_NODE:-192} )); then
  echo "ERROR: manifest asks for $total_threads threads on a ${SLURM_CPUS_ON_NODE:-192}-core node" >&2
  exit 2
fi

# ---- CPU set ----------------------------------------------------------------
# Cases are pinned to contiguous blocks counted from CPU 0, which is correct
# only when the job owns CPUs 0..N-1 (an exclusive genoa node). Refuse anything
# else rather than pinning cases onto CPUs the job does not own.
allowed=$(awk '/^Cpus_allowed_list:/ {print $2}' /proc/self/status)
if [[ ! "$allowed" =~ ^0-([0-9]+)$ ]] || (( BASH_REMATCH[1] + 1 < total_threads )); then
  echo "ERROR: job CPU set '$allowed' is not 0-N covering $total_threads threads" >&2
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
BIN_HALF="burstingBubbleInfiniteRr-camp-${JOB}"
BIN_FINITE="coalescenceBubbleFiniteRr-camp-${JOB}"
build() {  # build <source.c> <binary>; a qcc warning is treated as a failure
  local log="build-$2.log"
  qcc -I. -I../src-local -O2 -Wall -disable-dimensions -fopenmp "$1" -o "$2" -lm > "$log" 2>&1 \
    || { cat "$log" >&2; echo "ERROR: build of $1 failed" >&2; exit 1; }
  if grep -qiE "warning|error" "$log"; then
    cat "$log" >&2; echo "ERROR: build of $1 reported warnings" >&2; exit 1
  fi
  echo "built $2 sha256=$(sha256sum "$2" | cut -d' ' -f1)"
}
(( need_half ))   && build burstingBubbleInfiniteRr.c "$BIN_HALF"
(( need_finite )) && build coalescenceBubbleFiniteRr.c "$BIN_FINITE"

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
  case_id="${CASES[$i]}"; series="${SERIESES[$i]}"; geom="${GEOMS[$i]}"; rr="${RRS[$i]}"
  bond="${BONDS[$i]}"; oh="${OHS[$i]}"; level="${LEVELS[$i]}"; zwall="${ZWALLS[$i]}"
  clear="${CLEARS[$i]}"; ldom="${LDOMS[$i]}"; ic="${ICS[$i]}"
  ffr="${FFRS[$i]}"; ffd="${FFDS[$i]}"; minl="${MINLS[$i]}"; nthreads="${THREADS[$i]}"
  if [[ "$geom" == halfspace ]]; then BINARY="$BIN_HALF"; rr_arg=1000; else BINARY="$BIN_FINITE"; rr_arg="$rr"; fi
  ic_sha=$(sha256sum "$ROOT/simulationCases/DataFiles/$ic" | cut -d' ' -f1)
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
                "geometryMode=$geom" "Rr=$rr_arg" "initialShape=$ic" "initialShapeSha256=$ic_sha" \
                "farFloorRadius=$ffr" "farFloorDrop=$ffd" "minLevel=$minl" "sourceCommit=$SOURCE_COMMIT"; do
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
      "OhOut=$oh" "RhoIn=1e-3" "Rr=$rr_arg" "MAXlevel=$level" \
      "tmax=$TMAX" "zWall=$zwall" \
      "dropRadiusMin=0.021005127" "dropPersistence=3" \
      "snapshotInterval=$SNAPINT" "drillAMR=0" \
      "geometryMode=$geom" "wallClearance=$clear" \
      "interfaceFloor=1" "MuRin=$MURIN" "Bond=$bond"
    solver_args=("${COALESCENCE_SOLVER_ARGS[@]}")
    {
      printf 'CaseNo=%s\nOh=%s\nBond=%s\nMuRin=%s\nMAXlevel=%s\ntmax=%s\nzWall=%s\n' \
        "$case_id" "$oh" "$bond" "$MURIN" "$level" "$TMAX" "$zwall"
      printf 'Ldomain=%s\ndrillAMR=0\ngeometryMode=%s\nRr=%s\nwallClearance=%s\n' "$ldom" "$geom" "$rr_arg" "$clear"
      printf 'series=%s\ninitialShape=%s\ninitialShapeSha256=%s\n' "$series" "$ic" "$ic_sha"
      printf 'farFloorRadius=%s\nfarFloorDrop=%s\nminLevel=%s\nicBandRefine=1\n' "$ffr" "$ffd" "$minl"
      printf 'dropRadiusMin=0.021005127\ndropPersistence=3\nsnapshotInterval=%s\n' "$SNAPINT"
      printf 'solverStack=filtered+double-projection\nqccFlags=\nsourceCommit=%s\n' "$SOURCE_COMMIT"
      printf 'launcher=runCoalescenceCampaign.sbatch\nargv=%s\n' "${solver_args[*]}"
    } > "$case_dir/case.params"
    echo "case=$case_id fresh geometry=$geom Rr=$rr_arg Oh=$oh Bond=$bond ML=$level Ldomain=$ldom ic=$ic farFloor=$ffr/$ffd minLevel=$minl threads=$nthreads"
  fi

  # Each case owns a contiguous core block (genoa has no SMT; logical CPU k is
  # core k), pinned with taskset and mirrored in OMP_PLACES so the threads
  # cannot migrate onto a neighbour's cores.
  lo=$next_core; hi=$(( next_core + nthreads - 1 )); next_core=$(( hi + 1 ))
  echo "  cores $lo-$hi" >> "$case_dir/RESUMES"
  (
    cd "$case_dir"
    export OMP_NUM_THREADS="$nthreads" OMP_PROC_BIND=close OMP_PLACES="{$lo}:$nthreads"
    export LDOMAIN="$ldom" INITIAL_SHAPE="$ic"
    if awk -v r="$ffr" 'BEGIN{exit !(r+0 > 0)}'; then
      export FAR_FLOOR_RADIUS="$ffr" FAR_FLOOR_DROP="$ffd"
    else
      unset FAR_FLOOR_RADIUS FAR_FLOOR_DROP
    fi
    if (( minl > 0 )); then export MIN_LEVEL="$minl"; else unset MIN_LEVEL; fi
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
