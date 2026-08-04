#!/usr/bin/env bash

# Drive a fully local Bayesian contour campaign without an interactive agent.
# The controller owns proposals and collection; simulations run synchronously
# as bounded OpenMP batches inside this detached supervisor. Terminal numerical
# failures move to adjacent Oh values for a bounded number of attempts.

set -euo pipefail

campaign_root=${1:?usage: runContourCampaignLoop.sh CAMPAIGN_ROOT PREDICTOR_ROOT [POLL_SECONDS]}
predictor_root=${2:?usage: runContourCampaignLoop.sh CAMPAIGN_ROOT PREDICTOR_ROOT [POLL_SECONDS]}
poll_seconds=${3:-120}
project_root=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
controller="${project_root}/contourWorkflow/contour_campaign.py"
state_file="${campaign_root}/campaign-state.json"
heartbeat_file="${campaign_root}/loop-heartbeat.json"
max_attempts=${CONTOUR_MAX_ATTEMPTS:-6}

if [[ ! "$poll_seconds" =~ ^[1-9][0-9]*$ ]]; then
  echo "POLL_SECONDS must be a positive integer" >&2
  exit 2
fi
if [[ ! "$max_attempts" =~ ^[1-9][0-9]*$ ]]; then
  echo "CONTOUR_MAX_ATTEMPTS must be a positive integer" >&2
  exit 2
fi

command=(
  python3 "$controller"
  --campaign-root "$campaign_root"
  --project-root "$project_root"
  --predictor-root "$predictor_root"
  --backend inline
)

write_heartbeat() {
  local loop_state=$1
  local detail=${2:-}
  python3 - "$heartbeat_file" "$loop_state" "$detail" <<'PY'
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

path = Path(sys.argv[1])
payload = {
    "updated_at": datetime.now(timezone.utc).isoformat(),
    "loop_state": sys.argv[2],
    "detail": sys.argv[3],
    "pid": os.getppid(),
}
temporary = path.with_suffix(".json.tmp")
temporary.write_text(json.dumps(payload, indent=2) + "\n")
temporary.replace(path)
PY
}

on_exit() {
  rc=$?
  if (( rc == 0 )); then
    write_heartbeat stopped "campaign loop exited normally"
  else
    write_heartbeat failed "campaign loop exited with code ${rc}"
  fi
}
trap on_exit EXIT

while true; do
  write_heartbeat advancing
  "${command[@]}" advance --submit

  readarray -t state_fields < <(
    python3 - "$state_file" <<'PY'
import json
import sys

state = json.load(open(sys.argv[1]))
print(state.get("state", "unknown"))
print(state.get("iteration", state.get("last_completed_iteration", 0)))
print(state.get("attempt", 0))
print(state.get("reason", ""))
PY
  )
  state=${state_fields[0]}
  iteration=${state_fields[1]}
  attempt=${state_fields[2]}
  reason=${state_fields[3]}
  write_heartbeat "$state" "iteration=${iteration} attempt=${attempt} ${reason}"

  case "$state" in
    complete)
      exit 0
      ;;
    running)
      sleep "$poll_seconds"
      ;;
    needs_attention)
      if (( attempt < max_attempts )); then
        write_heartbeat retrying "iteration=${iteration} completed attempt=${attempt}"
        "${command[@]}" retry --submit
        sleep "$poll_seconds"
      else
        echo "Campaign needs attention after ${attempt} attempts: ${reason}" >&2
        exit 3
      fi
      ;;
    ready)
      sleep 2
      ;;
    *)
      echo "Unexpected campaign state: ${state}" >&2
      exit 4
      ;;
  esac
done
