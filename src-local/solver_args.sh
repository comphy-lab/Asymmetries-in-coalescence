#!/bin/bash
# ==============================================================================
# solver_args.sh - Canonical command-line contract for coalescenceBubble
#
# `simulationCases/coalescenceBubble.c` reads 25 positional arguments. Every
# launcher used to spell that list out for itself, which is how four of them
# drifted back to the retired six-argument call while the solver grew nineteen
# more controls. The order and the per-argument defaults live here instead, so
# a new solver argument is added in exactly two places: the solver and this
# file.
#
# Usage:
#   source src-local/solver_args.sh
#
# After sourcing:
#   COALESCENCE_SOLVER_ARG_KEYS       Parameter names in argv order (1..25)
#   COALESCENCE_SOLVER_ARG_DEFAULTS   Solver's compiled default per optional key
#   COALESCENCE_SOLVER_ARGS           Output array, rebuilt by the builders below
#
#   coalescence_solver_args <lookup_fn> [KEY=VALUE ...]
#       Generic builder. `lookup_fn KEY` prints a value or nothing; an empty
#       result falls back to the compiled default, and a key with no default is
#       reported as missing. Trailing KEY=VALUE pairs override any lookup.
#
#   coalescence_solver_args_from_params [KEY=VALUE ...]
#       Builder for the `.params` workflow. Requires `parse_params.sh` to be
#       sourced and `parse_param_file` to have run.
#
# The defaults below mirror the initialisers in `coalescenceBubble.c`. Passing
# them explicitly is equivalent to omitting the argument, so a launcher that
# adopts the full list keeps its previous behaviour until a parameter file
# actually sets one of the newer controls.
#
# Author: Vatsal Sanjay
# Organization: CoMPhy Lab, Durham University
# ==============================================================================

# argv 1..25, in the order `main()` parses them.
COALESCENCE_SOLVER_ARG_KEYS=(
    OhOut               # argv 1  - required
    RhoIn               # argv 2  - required
    Rr                  # argv 3  - required
    MAXlevel            # argv 4  - required
    tmax                # argv 5  - required
    zWall               # argv 6  - required
    dropRadiusMin       # argv 7
    dropPersistence     # argv 8
    snapshotInterval    # argv 9
    drillAMR            # argv 10
    drillMaxlevelStart  # argv 11
    drillMaxlevelFocus  # argv 12
    drillNcells         # argv 13
    drillRegionMinX     # argv 14
    drillArmSteps       # argv 15
    drillArmTime        # argv 16
    drillCoarsenTime    # argv 17
    drillRegionMaxX     # argv 18
    drillRegionRadius   # argv 19
    drillFireX          # argv 20
    drillTipRadius      # argv 21
    drillRegionalOnly   # argv 22
    geometryMode        # argv 23
    wallClearance       # argv 24
    interfaceFloor      # argv 25
)

# Compiled defaults for argv 7..25. A key absent from this map is required:
# argv 1..6 have no meaningful default and the solver refuses to start without
# them.
declare -gA COALESCENCE_SOLVER_ARG_DEFAULTS=(
    [dropRadiusMin]="-1"
    [dropPersistence]="3"
    [snapshotInterval]="1e-2"
    [drillAMR]="0"
    [drillMaxlevelStart]="-1"
    [drillMaxlevelFocus]="-1"
    [drillNcells]="5"
    [drillRegionMinX]="-2.1"
    [drillArmSteps]="5"
    [drillArmTime]="0"
    [drillCoarsenTime]="0"
    [drillRegionMaxX]="3"
    [drillRegionRadius]="1.5"
    [drillFireX]="0.25"
    [drillTipRadius]="0.25"
    [drillRegionalOnly]="0"
    [geometryMode]="finite"
    [wallClearance]="-1"
    [interfaceFloor]="1"
)

# Fallbacks the sweep launchers have always applied to the six required
# arguments when a parameter file omits them.
declare -gA COALESCENCE_SOLVER_BASE_DEFAULTS=(
    [OhOut]="1e-2"
    [RhoIn]="1e-3"
    [Rr]="1.0"
    [MAXlevel]="10"
    [tmax]="40.0"
    [zWall]="0.01"
)

COALESCENCE_SOLVER_ARGS=()

# ==============================================================================
# Function: coalescence_solver_args
#
# Description:
#   Build the complete positional argument list for coalescenceBubble.
#
# Parameters:
#   $1   - Name of a lookup function taking one key and printing its value
#   $2.. - Optional KEY=VALUE overrides applied after the lookup
#
# Returns:
#   0 on success; 1 if a required argument resolves to nothing
#
# Side Effects:
#   Replaces the contents of COALESCENCE_SOLVER_ARGS
# ==============================================================================
coalescence_solver_args() {
    local lookup=$1
    shift

    local -A overrides=()
    local pair
    for pair in "$@"; do
        if [[ "$pair" != *=* ]]; then
            echo "ERROR: solver argument override must be KEY=VALUE, got '$pair'" >&2
            return 1
        fi
        overrides["${pair%%=*}"]="${pair#*=}"
    done

    COALESCENCE_SOLVER_ARGS=()
    local key value
    for key in "${COALESCENCE_SOLVER_ARG_KEYS[@]}"; do
        if [[ -n "${overrides[$key]+set}" ]]; then
            value=${overrides[$key]}
        else
            value=$("$lookup" "$key") || return 1
        fi
        if [[ -z "$value" ]]; then
            if [[ -z "${COALESCENCE_SOLVER_ARG_DEFAULTS[$key]+set}" ]]; then
                echo "ERROR: missing required solver argument $key" >&2
                return 1
            fi
            value=${COALESCENCE_SOLVER_ARG_DEFAULTS[$key]}
        fi
        COALESCENCE_SOLVER_ARGS+=("$value")
    done
}

# ==============================================================================
# Function: coalescence_solver_args_from_params
#
# Description:
#   Build the argument list from the parameters most recently loaded by
#   `parse_param_file`, applying the sweep launchers' historical fallbacks for
#   the six required arguments.
#
# Parameters:
#   $@ - Optional KEY=VALUE overrides, e.g. a staged tmax
#
# Returns:
#   0 on success, 1 on a missing required argument
# ==============================================================================
coalescence_solver_args_from_params() {
    if ! declare -F get_param >/dev/null; then
        echo "ERROR: source src-local/parse_params.sh before solver_args.sh helpers" >&2
        return 1
    fi
    coalescence_solver_args _coalescence_param_lookup "$@"
}

# Internal lookup used by coalescence_solver_args_from_params.
_coalescence_param_lookup() {
    local key=$1
    get_param "$key" "${COALESCENCE_SOLVER_BASE_DEFAULTS[$key]:-}"
}
