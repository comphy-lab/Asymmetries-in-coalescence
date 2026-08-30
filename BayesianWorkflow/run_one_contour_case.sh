#!/bin/bash

set -euo pipefail

case_dir=${1:?usage: run_one_contour_case.sh CASE_DIR BINARY}
binary=${2:?usage: run_one_contour_case.sh CASE_DIR BINARY}
params="${case_dir}/case.params"
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)

if [[ ! -f "$params" ]]; then
  echo "Missing parameter file: $params" >&2
  exit 2
fi

# shellcheck source=../src-local/solver_args.sh
source "${script_dir}/../src-local/solver_args.sh"

# A materialised case table carries every control explicitly, so anything the
# contour campaign is meant to pin must be present. The four keys that
# postdate the earliest case tables fall back to the solver's compiled default
# instead, which keeps older batches reproducible.
case_param() {
  local key=$1
  local value
  value=$(awk -F= -v key="$key" '$1 == key {print $2; exit}' "$params")
  if [[ -z "$value" ]]; then
    case "$key" in
      geometryMode | wallClearance | interfaceFloor | MuRin) return 0 ;;
      *)
        echo "Missing $key in $params" >&2
        return 1
        ;;
    esac
  fi
  printf '%s\n' "$value"
}

coalescence_solver_args case_param

cd "$case_dir"
rm -f classification.status classification.status.tmp runner.status runner.status.tmp
cp "$binary" ./coalescenceBubbleContour
chmod u+x ./coalescenceBubbleContour

tmp_status=runner.status.tmp
final_status=runner.status
printf 'state=running\nstarted_at=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" > "$tmp_status"
mv "$tmp_status" "$final_status"

set +e
./coalescenceBubbleContour "${COALESCENCE_SOLVER_ARGS[@]}"
rc=$?
set -e

classification_id=""
if [[ -f classification.status ]]; then
  classification_id=$(awk -F, 'NR == 2 {print $1}' classification.status)
fi

if [[ $rc -eq 0 && "$classification_id" =~ ^[01]$ ]]; then
  state=complete
else
  state=failed
fi

{
  printf 'state=%s\n' "$state"
  printf 'exit_code=%s\n' "$rc"
  printf 'classification_id=%s\n' "${classification_id:-unresolved}"
  printf 'finished_at=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
} > "$tmp_status"
mv "$tmp_status" "$final_status"

[[ "$state" == complete ]]
