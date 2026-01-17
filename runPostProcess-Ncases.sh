#!/bin/bash
# runPostProcess-Ncases.sh - Run post-processing pipeline on multiple simulation cases
# Author: Vatsal Sanjay
# vatsal.sanjay@comphy-lab.org
# CoMPhy Lab - Durham University
# Last updated: Jan 2026

set -euo pipefail  # Exit on error, unset variables, pipeline failures

# ============================================================
# Configuration
# ============================================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Source parameter parsing library
if [ -f "${SCRIPT_DIR}/src-local/parse_params.sh" ]; then
    source "${SCRIPT_DIR}/src-local/parse_params.sh"
else
    echo "ERROR: src-local/parse_params.sh not found" >&2
    exit 1
fi

# Post-processing script paths
VIDEO_SCRIPT="${SCRIPT_DIR}/postProcess/Video-generic.py"

# C helper executables
HELPER_GETFACET="${SCRIPT_DIR}/postProcess/getFacet"
HELPER_GETDATA="${SCRIPT_DIR}/postProcess/getData-generic"
HELPER_GETCOM="${SCRIPT_DIR}/postProcess/getCOM"

# Case directory root
CASES_DIR="${SCRIPT_DIR}/simulationCases"

# ============================================================
# Usage Information
# ============================================================
usage() {
    cat <<EOF
Usage: $0 [OPTIONS] CASE_NO [CASE_NO ...]

Run post-processing pipeline on multiple bubble coalescence simulation cases.
For each case, generates video frames with strain-rate/velocity fields,
interface facets, and center of mass tracking.

Options:
    --CPUs N            Number of parallel workers (default: 4)
    --nGFS N            Number of snapshots to process (default: 4000)
    --tsnap F           Time interval between snapshots (default: 0.01)
    --GridsPerR N       Radial grid resolution for video (default: 256)
    --ZMIN F            Override minimum Z (auto: max(-2-zWall, -3.0))
    --ZMAX F            Override maximum Z (auto: 5*Rr-2.0)
    --RMAX F            Override maximum R (auto: 2.5*Rr)

    Note: Rr and zWall are auto-extracted from each case's case.params file.
    Domain bounds are calculated per-case unless overridden above.

    --skip-video-encode Skip ffmpeg video encoding after frame generation

    -n, --dry-run       Show what would run without executing
    -v, --verbose       Verbose output
    -h, --help          Show this help message

Arguments:
    CASE_NO             4-digit case numbers (1000-9999) or ranges (e.g., 3000-3010)

Examples:
    # Process multiple cases with default settings
    $0 3000 3001 3002

    # Process a range of cases
    $0 3000-3010

    # Mix individual cases and ranges
    $0 3000 3005-3007 3010

    # Process with 8 CPUs
    $0 --CPUs 8 3000 3001

    # Process first 100 snapshots only (for testing)
    $0 --nGFS 100 3000

    # Skip video encoding (only generate frames)
    $0 --skip-video-encode 3000 3001

    # Dry run to preview commands
    $0 --dry-run 3000

Output locations:
    simulationCases/<CaseNo>/Video/              # PNG frames
    simulationCases/<CaseNo>/<CaseNo>_COMData.csv # COM time series
    simulationCases/<CaseNo>/<CaseNo>.mp4        # Encoded video

For more information, see CLAUDE.md
EOF
}

# ============================================================
# Helper: Expand Case Argument (single number or range)
# ============================================================
expand_case_arg() {
    local arg="$1"

    # Check if it's a range (contains hyphen between two numbers)
    if [[ "$arg" =~ ^([0-9]{4})-([0-9]{4})$ ]]; then
        local start="${BASH_REMATCH[1]}"
        local end="${BASH_REMATCH[2]}"

        # Validate range bounds
        if [ "$start" -lt 1000 ] || [ "$start" -gt 9999 ]; then
            echo "ERROR: Range start out of bounds (1000-9999): $start" >&2
            return 1
        fi
        if [ "$end" -lt 1000 ] || [ "$end" -gt 9999 ]; then
            echo "ERROR: Range end out of bounds (1000-9999): $end" >&2
            return 1
        fi

        # Validate range order
        if [ "$start" -gt "$end" ]; then
            echo "ERROR: Invalid range (start > end): $arg" >&2
            return 1
        fi

        # Expand range using seq
        seq "$start" "$end"
        return 0
    fi

    # Check if it's a single 4-digit number
    if [[ "$arg" =~ ^[0-9]{4}$ ]]; then
        if [ "$arg" -lt 1000 ] || [ "$arg" -gt 9999 ]; then
            echo "ERROR: Case number out of range (1000-9999): $arg" >&2
            return 1
        fi
        echo "$arg"
        return 0
    fi

    # Invalid format
    echo "ERROR: Invalid case argument (expected 4-digit number or range): $arg" >&2
    return 1
}

# ============================================================
# Parse Command Line Options
# ============================================================
CPUS=4
NGFS=4000
TSNAP=0.01
GRIDS_PER_R=256
ZMIN=""
ZMAX=""
RMAX=""

SKIP_VIDEO_ENCODE=0
DRY_RUN=0
VERBOSE=0

RAW_CASE_ARGS=()
CASE_NUMBERS=()

while [[ $# -gt 0 ]]; do
    case $1 in
        --CPUs)
            CPUS="$2"
            if ! [[ "$CPUS" =~ ^[0-9]+$ ]] || [ "$CPUS" -lt 1 ]; then
                echo "ERROR: --CPUs requires a positive integer, got: $CPUS" >&2
                exit 1
            fi
            shift 2
            ;;
        --nGFS)
            NGFS="$2"
            if ! [[ "$NGFS" =~ ^[0-9]+$ ]] || [ "$NGFS" -lt 1 ]; then
                echo "ERROR: --nGFS requires a positive integer, got: $NGFS" >&2
                exit 1
            fi
            shift 2
            ;;
        --tsnap)
            TSNAP="$2"
            shift 2
            ;;
        --GridsPerR)
            GRIDS_PER_R="$2"
            if ! [[ "$GRIDS_PER_R" =~ ^[0-9]+$ ]] || [ "$GRIDS_PER_R" -lt 1 ]; then
                echo "ERROR: --GridsPerR requires a positive integer, got: $GRIDS_PER_R" >&2
                exit 1
            fi
            shift 2
            ;;
        --ZMIN)
            ZMIN="$2"
            shift 2
            ;;
        --ZMAX)
            ZMAX="$2"
            shift 2
            ;;
        --RMAX)
            RMAX="$2"
            shift 2
            ;;
        --skip-video-encode)
            SKIP_VIDEO_ENCODE=1
            shift
            ;;
        -n|--dry-run)
            DRY_RUN=1
            shift
            ;;
        -v|--verbose)
            VERBOSE=1
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        -*)
            echo "ERROR: Unknown option: $1" >&2
            usage
            exit 1
            ;;
        *)
            # Collect raw case arguments (numbers or ranges)
            RAW_CASE_ARGS+=("$1")
            shift
            ;;
    esac
done

# Expand raw case arguments into individual case numbers
for raw_arg in "${RAW_CASE_ARGS[@]}"; do
    expanded=$(expand_case_arg "$raw_arg") || exit 1
    while IFS= read -r case_no; do
        CASE_NUMBERS+=("$case_no")
    done <<< "$expanded"
done

# ============================================================
# Validation
# ============================================================

# Check at least one case argument provided
if [ ${#RAW_CASE_ARGS[@]} -eq 0 ]; then
    echo "ERROR: No case numbers provided" >&2
    usage
    exit 1
fi

# Note: Individual case number validation is handled by expand_case_arg()

# Check Python availability
if ! command -v python &> /dev/null; then
    echo "ERROR: python not found in PATH" >&2
    exit 1
fi

# Check Python script exists
if [ ! -f "$VIDEO_SCRIPT" ]; then
    echo "ERROR: Python script not found: $VIDEO_SCRIPT" >&2
    exit 1
fi

# Check C helpers exist
for helper in "$HELPER_GETFACET" "$HELPER_GETDATA" "$HELPER_GETCOM"; do
    if [ ! -x "$helper" ]; then
        helper_name=$(basename "$helper")
        echo "ERROR: Compiled helper not found or not executable: $helper" >&2
        echo "       Compile with: qcc -O2 -Wall postProcess/${helper_name}.c -o postProcess/${helper_name} -lm" >&2
        exit 1
    fi
done

# ============================================================
# Display Configuration
# ============================================================
echo "========================================="
echo "Bubble Coalescence - Post-Processing Pipeline"
echo "========================================="
echo "Cases to process: ${CASE_NUMBERS[*]}"
echo "Total cases: ${#CASE_NUMBERS[@]}"
echo ""
echo "Settings:"
echo "  CPUs:       $CPUS"
echo "  nGFS:       $NGFS"
echo "  tsnap:      $TSNAP"
echo "  GridsPerR:  $GRIDS_PER_R"
echo "  Domain bounds: auto-calculated from case.params (Rr, zWall)"
if [ -n "$ZMIN" ]; then echo "  ZMIN:       $ZMIN (override)"; fi
if [ -n "$ZMAX" ]; then echo "  ZMAX:       $ZMAX (override)"; fi
if [ -n "$RMAX" ]; then echo "  RMAX:       $RMAX (override)"; fi
echo ""
echo "Pipeline:"
[ $SKIP_VIDEO_ENCODE -eq 0 ] && echo "  [1] Video-generic.py (frames + video)" || echo "  [1] Video-generic.py (frames only, video SKIPPED)"
echo ""
[ $DRY_RUN -eq 1 ] && echo "Mode: DRY RUN (no execution)"
echo ""

# ============================================================
# Processing Functions
# ============================================================

run_video() {
    local case_no="$1"
    local case_dir="${CASES_DIR}/${case_no}"
    local video_dir="${case_dir}/Video"

    # Read Rr and zWall from case.params
    local case_params="${case_dir}/case.params"
    local case_rr="1.0"
    local case_zwall="4.0"  # default from default.params

    if [ -f "$case_params" ]; then
        parse_param_file "$case_params"
        case_rr=$(get_param "Rr" "1.0")
        case_zwall=$(get_param "zWall" "4.0")
    else
        echo "  WARNING: case.params not found, using defaults (Rr=1.0, zWall=4.0)"
    fi

    # Calculate domain bounds
    local calc_zmin=$(echo "-2 - $case_zwall" | bc -l)
    local calc_zmax=$(echo "5 * $case_rr - 2.0" | bc -l)
    local calc_rmax=$(echo "2.5 * $case_rr" | bc -l)

    # Apply ZMIN clamp: max(-2-zWall, -3.0)
    if (( $(echo "$calc_zmin < -3.0" | bc -l) )); then
        calc_zmin="-3.0"
    fi

    # Use command-line overrides if provided, otherwise use calculated values
    local use_zmin="${ZMIN:-$calc_zmin}"
    local use_zmax="${ZMAX:-$calc_zmax}"
    local use_rmax="${RMAX:-$calc_rmax}"

    if [ $VERBOSE -eq 1 ]; then
        echo "  Extracted: Rr=$case_rr, zWall=$case_zwall"
        echo "  Calculated: ZMIN=$use_zmin, ZMAX=$use_zmax, RMAX=$use_rmax"
    fi

    # Build command with calculated/override values
    local cmd_args=(
        "--caseToProcess" "${case_dir}"
        "--folderToSave" "${video_dir}"
        "--CPUs" "${CPUS}"
        "--nGFS" "${NGFS}"
        "--tsnap" "${TSNAP}"
        "--GridsPerR" "${GRIDS_PER_R}"
        "--Rr" "${case_rr}"
        "--ZMIN" "${use_zmin}"
        "--ZMAX" "${use_zmax}"
        "--RMAX" "${use_rmax}"
    )

    # Add skip flag if needed
    [ $SKIP_VIDEO_ENCODE -eq 1 ] && cmd_args+=("--skip-video-encode")

    if [ $VERBOSE -eq 1 ] || [ $DRY_RUN -eq 1 ]; then
        echo "  CMD: python ${VIDEO_SCRIPT} ${cmd_args[*]}"
    fi

    if [ $DRY_RUN -eq 0 ]; then
        python "${VIDEO_SCRIPT}" "${cmd_args[@]}"
    fi
}

# ============================================================
# Main Processing Loop
# ============================================================
echo "========================================="
echo "Processing Cases"
echo "========================================="

SUCCESSFUL_CASES=()
FAILED_CASES=()
FAILURE_REASONS=()

for case_no in "${CASE_NUMBERS[@]}"; do
    echo ""
    echo "-----------------------------------------"
    echo "Case $case_no"
    echo "-----------------------------------------"

    case_dir="${CASES_DIR}/${case_no}"
    intermediate_dir="${case_dir}/intermediate"

    # Validate case directory exists
    if [ ! -d "$case_dir" ]; then
        echo "  ERROR: Case directory not found: $case_dir"
        FAILED_CASES+=("$case_no")
        FAILURE_REASONS+=("Case directory not found")
        continue
    fi

    if [ ! -d "$intermediate_dir" ]; then
        echo "  ERROR: Intermediate directory not found: $intermediate_dir"
        FAILED_CASES+=("$case_no")
        FAILURE_REASONS+=("No intermediate/ snapshots")
        continue
    fi

    # Count snapshots
    snapshot_count=$(find "$intermediate_dir" -name "snapshot-*" 2>/dev/null | wc -l | tr -d ' ')
    echo "  Found $snapshot_count snapshots in intermediate/"

    if [ "$snapshot_count" -eq 0 ]; then
        echo "  ERROR: No snapshots found"
        FAILED_CASES+=("$case_no")
        FAILURE_REASONS+=("No snapshots in intermediate/")
        continue
    fi

    # Track step failures
    step_failed=0

    # Run video generation
    echo ""
    echo "  [1/1] Running Video-generic.py..."
    if ! run_video "$case_no"; then
        echo "  ERROR: Video generation failed"
        step_failed=1
    else
        [ $DRY_RUN -eq 0 ] && echo "  [1/1] Video generation complete"
    fi

    # Record result
    if [ $step_failed -eq 0 ]; then
        SUCCESSFUL_CASES+=("$case_no")
        echo ""
        echo "  Case $case_no: SUCCESS"
    else
        FAILED_CASES+=("$case_no")
        FAILURE_REASONS+=("Processing step failed")
        echo ""
        echo "  Case $case_no: FAILED"
    fi
done

# ============================================================
# Summary
# ============================================================
echo ""
echo "========================================="
echo "Post-Processing Complete"
echo "========================================="
echo "Total cases: ${#CASE_NUMBERS[@]}"
echo "Successful:  ${#SUCCESSFUL_CASES[@]}"
echo "Failed:      ${#FAILED_CASES[@]}"

if [ ${#FAILED_CASES[@]} -gt 0 ]; then
    echo ""
    echo "Failed cases:"
    for i in "${!FAILED_CASES[@]}"; do
        echo "  - Case ${FAILED_CASES[$i]}: ${FAILURE_REASONS[$i]}"
    done
fi

if [ ${#SUCCESSFUL_CASES[@]} -gt 0 ]; then
    echo ""
    echo "Output locations:"
    for case_no in "${SUCCESSFUL_CASES[@]}"; do
        echo "  ${CASES_DIR}/${case_no}/"
    done
fi

echo ""

# Exit with error if any cases failed
if [ ${#FAILED_CASES[@]} -gt 0 ]; then
    exit 1
fi

exit 0
