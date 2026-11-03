#!/usr/bin/env bash
# Strict reader shared by launchers which use the repository Basilisk pin.

basilisk_read_pin() {
    local file=$1 line ref="" count=0
    if [[ ! -f $file || ! -r $file ]]; then
        printf 'ERROR: missing or unreadable Basilisk pin: %s\n' "$file" >&2
        return 1
    fi
    while IFS= read -r line || [[ -n $line ]]; do
        line=${line%%#*}
        line="${line#"${line%%[![:space:]]*}"}"
        line="${line%"${line##*[![:space:]]}"}"
        [[ -n $line ]] || continue
        ref=$line
        count=$((count + 1))
    done < "$file"
    if [[ $count -ne 1 || ! $ref =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]]; then
        printf 'ERROR: Basilisk pin must contain exactly one valid ref: %s\n' "$file" >&2
        return 1
    fi
    printf '%s\n' "$ref"
}

basilisk_require_campaign_pin() {
    local root=$1 pin
    pin=$(basilisk_read_pin "$root/basilisk.ref") || return 1
    if [[ ${EXPECTED_BASILISK_REF+x} && ${EXPECTED_BASILISK_REF} != "$pin" ]]; then
        printf 'ERROR: EXPECTED_BASILISK_REF conflicts with the component pin %s\n' "$pin" >&2
        return 1
    fi
    EXPECTED_BASILISK_REF=$pin
}

basilisk_require_lock() {
    local lock=$1 expected=$2 line ref="" count=0
    if [[ ! -f $lock || ! -r $lock ]]; then
        printf 'ERROR: missing or unreadable Basilisk lock: %s\n' "$lock" >&2
        return 1
    fi
    while IFS= read -r line || [[ -n $line ]]; do
        case $line in
            ref=*) ref=${line#ref=}; count=$((count + 1)) ;;
        esac
    done < "$lock"
    if [[ $count -ne 1 || $ref != "$expected" ]]; then
        printf 'ERROR: expected Basilisk ref=%s in %s, with exactly one ref line\n' "$expected" "$lock" >&2
        return 1
    fi
}
