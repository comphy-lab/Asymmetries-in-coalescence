#!/bin/bash

set -euo pipefail

case_dir=${1:?usage: run_one_contour_case.sh CASE_DIR BINARY}
binary=${2:?usage: run_one_contour_case.sh CASE_DIR BINARY}
params="${case_dir}/case.params"

if [[ ! -f "$params" ]]; then
  echo "Missing parameter file: $params" >&2
  exit 2
fi

get_param() {
  local key=$1
  local value
  value=$(awk -F= -v key="$key" '$1 == key {print $2; exit}' "$params")
  if [[ -z "$value" ]]; then
    echo "Missing $key in $params" >&2
    return 1
  fi
  printf '%s\n' "$value"
}

OhOut=$(get_param OhOut)
RhoIn=$(get_param RhoIn)
Rr=$(get_param Rr)
MAXlevel=$(get_param MAXlevel)
tmax=$(get_param tmax)
zWall=$(get_param zWall)
dropRadiusMin=$(get_param dropRadiusMin)
dropPersistence=$(get_param dropPersistence)
snapshotInterval=$(get_param snapshotInterval)
drillAMR=$(get_param drillAMR)
drillMaxlevelStart=$(get_param drillMaxlevelStart)
drillMaxlevelFocus=$(get_param drillMaxlevelFocus)
drillNcells=$(get_param drillNcells)
drillRegionMinX=$(get_param drillRegionMinX)
drillArmSteps=$(get_param drillArmSteps)
drillArmTime=$(get_param drillArmTime)
drillCoarsenTime=$(get_param drillCoarsenTime)
drillRegionMaxX=$(get_param drillRegionMaxX)
drillRegionRadius=$(get_param drillRegionRadius)
drillFireX=$(get_param drillFireX)
drillTipRadius=$(get_param drillTipRadius)
drillRegionalOnly=$(get_param drillRegionalOnly)

cd "$case_dir"
rm -f classification.status classification.status.tmp runner.status runner.status.tmp
cp "$binary" ./coalescenceBubbleContour
chmod u+x ./coalescenceBubbleContour

tmp_status=runner.status.tmp
final_status=runner.status
printf 'state=running\nstarted_at=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" > "$tmp_status"
mv "$tmp_status" "$final_status"

set +e
./coalescenceBubbleContour \
  "$OhOut" "$RhoIn" "$Rr" "$MAXlevel" "$tmax" "$zWall" \
  "$dropRadiusMin" "$dropPersistence" "$snapshotInterval" "$drillAMR" \
  "$drillMaxlevelStart" "$drillMaxlevelFocus" "$drillNcells" \
  "$drillRegionMinX" "$drillArmSteps" "$drillArmTime" "$drillCoarsenTime" \
  "$drillRegionMaxX" "$drillRegionRadius" "$drillFireX" "$drillTipRadius" \
  "$drillRegionalOnly"
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
