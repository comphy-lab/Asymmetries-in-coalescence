#!/bin/bash
#SBATCH --job-name=bubbleCoal-S1
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --partition=genoa
#SBATCH --mail-type=ALL
#SBATCH --mail-user=vatsal.sanjay@comphy-lab.org
#SBATCH --output=slurm-stage1-%j.out
#SBATCH --error=slurm-stage1-%j.err

# ============================================================
# Bubble Coalescence - Stage 1 Serial Runner (Snellius)
# ============================================================
# This script runs Stage 1 (initialization) for all cases in parallel.
# Each case runs on 1 core (serial, no OpenMP/MPI).
# Up to 48 cases run simultaneously.
#
# Stage 1 generates restart files needed for Stage 2 (MPI simulation).
#
# AFTER RUNNING:
# Submit Stage 2 with: sbatch runSweepSnellius.sbatch
#
# SBATCH Parameters to Customize:
#   --ntasks: Number of concurrent cases (currently 48)
#   --time: Wall time for all Stage 1 runs (currently 12 hours)
#   --partition: Compute partition (currently genoa)
#   --mail-user: Email for job notifications
# ============================================================

set -e  # Exit on error (but we'll handle case failures gracefully)

# ============================================================
# Configuration
# ============================================================
SCRIPT_DIR="${SLURM_SUBMIT_DIR}"
SWEEP_FILE="${SCRIPT_DIR}/sweep.params"

# Maximum concurrent jobs (matches SLURM allocation)
MAX_CONCURRENT=${SLURM_NTASKS:-48}

# Stage 1 tmax - short run to generate restart file
STAGE1_TMAX="1e-2"

# ============================================================
# Print Job Information
# ============================================================
echo "============================================="
echo "Bubble Coalescence - Stage 1 (Serial)"
echo "============================================="
echo "Job started at: $(date)"
echo "Running on node: $(hostname)"
echo "Working directory: $(pwd)"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Max concurrent cases: ${MAX_CONCURRENT}"
echo "Partition: ${SLURM_JOB_PARTITION}"
echo ""
echo "Mode: Stage 1 only (serial initialization)"
echo "      Generates restart files for Stage 2"
echo "      Each case runs on 1 core (no OpenMP/MPI)"
echo ""

# ============================================================
# Load Required Modules (minimal for serial)
# ============================================================
echo "Loading modules..."
module purge
module load 2024
module load GCC/13.3.0
echo "Modules loaded (serial mode - GCC only, no MPI)"
echo ""

# ============================================================
# Setup Basilisk Environment
# ============================================================
echo "Setting up Basilisk environment..."
if [ -f "${SCRIPT_DIR}/.project_config" ]; then
    source "${SCRIPT_DIR}/.project_config"
    echo "Basilisk environment loaded from .project_config"
    echo "BASILISK: $BASILISK"
else
    echo "ERROR: .project_config not found" >&2
    exit 1
fi
echo ""

# ============================================================
# Validate Environment
# ============================================================
if [ -f "${SCRIPT_DIR}/src-local/parse_params.sh" ]; then
    source "${SCRIPT_DIR}/src-local/parse_params.sh"
else
    echo "ERROR: src-local/parse_params.sh not found" >&2
    exit 1
fi

if [ ! -f "$SWEEP_FILE" ]; then
    echo "ERROR: Sweep file not found: $SWEEP_FILE" >&2
    exit 1
fi

echo "Sweep file: $SWEEP_FILE"
echo ""

# ============================================================
# Parse Sweep Configuration
# ============================================================
echo "Parsing sweep configuration..."

source "$SWEEP_FILE"

if [ -z "$BASE_CONFIG" ]; then
    echo "ERROR: BASE_CONFIG not defined in sweep file" >&2
    exit 1
fi

if [ -z "$CASE_START" ] || [ -z "$CASE_END" ]; then
    echo "ERROR: CASE_START and CASE_END must be defined in sweep file" >&2
    exit 1
fi

if [ "$CASE_START" -lt 1000 ] || [ "$CASE_START" -gt 9999 ]; then
    echo "ERROR: CASE_START must be 4-digit (1000-9999), got: $CASE_START" >&2
    exit 1
fi

if [ "$CASE_END" -lt "$CASE_START" ] || [ "$CASE_END" -gt 9999 ]; then
    echo "ERROR: CASE_END must be >= CASE_START and <= 9999, got: $CASE_END" >&2
    exit 1
fi

if [ ! -f "$BASE_CONFIG" ]; then
    echo "ERROR: Base configuration file not found: $BASE_CONFIG" >&2
    exit 1
fi

echo "Base configuration: $BASE_CONFIG"
echo "Case number range: $CASE_START to $CASE_END"
echo ""

# ============================================================
# Extract Sweep Variables
# ============================================================
SWEEP_VARS=()
SWEEP_VALUES=()

while IFS='=' read -r key value; do
    [[ "$key" =~ ^[[:space:]]*# ]] && continue
    [[ -z "$key" ]] && continue

    if [[ "$key" =~ ^[[:space:]]*SWEEP_([^=]+) ]]; then
        var_name="${BASH_REMATCH[1]}"
        value=$(echo "$value" | sed 's/#.*//' | xargs)
        SWEEP_VARS+=("$var_name")
        SWEEP_VALUES+=("$value")
    fi
done < "$SWEEP_FILE"

if [ ${#SWEEP_VARS[@]} -eq 0 ]; then
    echo "ERROR: No SWEEP_* variables found in $SWEEP_FILE" >&2
    exit 1
fi

echo "Sweep variables:"
for i in "${!SWEEP_VARS[@]}"; do
    echo "  ${SWEEP_VARS[$i]} = ${SWEEP_VALUES[$i]}"
done
echo ""

# ============================================================
# Generate Parameter Combinations
# ============================================================
echo "Generating parameter combinations..."

TEMP_DIR="${SCRIPT_DIR}/.sweep_tmp_$$"
mkdir -p "$TEMP_DIR" || {
    echo "ERROR: Failed to create temp directory: $TEMP_DIR" >&2
    exit 1
}
trap "rm -rf $TEMP_DIR" EXIT

CASE_NUM=$CASE_START
COMBINATION_COUNT=0
CASE_FILES=()

generate_combinations() {
    local depth=$1
    shift
    local current_values=("$@")

    if [ $depth -eq ${#SWEEP_VARS[@]} ]; then
        local case_file="${TEMP_DIR}/case_$(printf "%04d" $CASE_NUM).params"
        cp "$BASE_CONFIG" "$case_file"

        if grep -q "^CaseNo=" "$case_file"; then
            sed -i'.bak' "s|^CaseNo=.*|CaseNo=${CASE_NUM}|" "$case_file"
        else
            echo "CaseNo=${CASE_NUM}" >> "$case_file"
        fi
        rm -f "${case_file}.bak"

        for i in "${!SWEEP_VARS[@]}"; do
            local var="${SWEEP_VARS[$i]}"
            local val="${current_values[$i]}"

            if grep -q "^${var}=" "$case_file"; then
                sed -i'.bak' "s|^${var}=.*|${var}=${val}|" "$case_file"
            else
                echo "${var}=${val}" >> "$case_file"
            fi
            rm -f "${case_file}.bak"
        done

        CASE_FILES+=("$case_file")
        ((COMBINATION_COUNT++)) || true

        echo "Case $CASE_NUM:"
        for i in "${!SWEEP_VARS[@]}"; do
            echo "  ${SWEEP_VARS[$i]} = ${current_values[$i]}"
        done
        echo ""

        ((CASE_NUM++)) || true
        return
    fi

    local values="${SWEEP_VALUES[$depth]}"
    IFS=',' read -ra value_array <<< "$values"

    for val in "${value_array[@]}"; do
        val=$(echo "$val" | xargs)
        generate_combinations $((depth + 1)) ${current_values[@]+"${current_values[@]}"} "$val"
    done
}

generate_combinations 0

echo "Generated $COMBINATION_COUNT parameter combinations"

EXPECTED_COUNT=$((CASE_END - CASE_START + 1))
if [ $COMBINATION_COUNT -ne $EXPECTED_COUNT ]; then
    echo "WARNING: Generated $COMBINATION_COUNT combinations, but CASE_END suggests $EXPECTED_COUNT" >&2
fi

if [ $COMBINATION_COUNT -gt $EXPECTED_COUNT ]; then
    echo "ERROR: Too many combinations ($COMBINATION_COUNT) for range $CASE_START-$CASE_END" >&2
    exit 1
fi

echo ""

# ============================================================
# Prepare Status Tracking
# ============================================================
STATUS_DIR="${TEMP_DIR}/status"
mkdir -p "$STATUS_DIR"

# ============================================================
# Function to Run Single Case (Stage 1)
# ============================================================
run_stage1_case() {
    local param_file=$1
    local status_dir=$2
    local script_dir=$3
    local stage1_tmax=$4

    # Source parameter parsing in subshell
    source "${script_dir}/src-local/parse_params.sh"

    parse_param_file "$param_file"
    local CASE_NO=$(get_param "CaseNo")
    local OhOut=$(get_param "OhOut" "1e-2")
    local RhoIn=$(get_param "RhoIn" "1e-3")
    local Rr=$(get_param "Rr" "1.0")
    local MAXlevel=$(get_param "MAXlevel" "10")
    local zWall=$(get_param "zWall" "0.01")

    if [ -z "$CASE_NO" ]; then
        echo "FAILED" > "${status_dir}/${CASE_NO:-unknown}"
        return 1
    fi

    local CASE_DIR="${script_dir}/simulationCases/${CASE_NO}"
    local LOG_FILE="${CASE_DIR}/stage1.log"

    {
        echo "========================================="
        echo "Case $CASE_NO - Stage 1 ($(date))"
        echo "========================================="
        echo "Parameters: OhOut=$OhOut, RhoIn=$RhoIn, Rr=$Rr"
        echo "            MAXlevel=$MAXlevel, tmax=$stage1_tmax, zWall=$zWall"

        # Create case directory
        mkdir -p "$CASE_DIR"
        cd "$CASE_DIR"

        # Copy source file
        local SRC_FILE_ORIG="${script_dir}/simulationCases/coalescenceBubble.c"
        if [ ! -f "$SRC_FILE_ORIG" ]; then
            echo "ERROR: Source file not found"
            echo "FAILED" > "${status_dir}/${CASE_NO}"
            return 1
        fi
        cp "$SRC_FILE_ORIG" "coalescenceBubble.c"

        # Create symlink to DataFiles
        if [ ! -e "DataFiles" ]; then
            ln -s ../DataFiles DataFiles
        fi

        # Compile SERIAL (no OpenMP, no MPI)
        echo "Compiling (serial)..."
        if ! qcc -I../../src-local -Wall -O2 -disable-dimensions \
            coalescenceBubble.c -o coalescenceBubble -lm 2>&1; then
            echo "ERROR: Compilation failed"
            echo "FAILED" > "${status_dir}/${CASE_NO}"
            return 1
        fi
        echo "Compilation successful"

        # Run Stage 1 with small tmax
        echo "Running Stage 1 (tmax=$stage1_tmax)..."
        if ./coalescenceBubble $OhOut $RhoIn $Rr $MAXlevel $stage1_tmax $zWall; then
            echo "Stage 1 completed"

            # Verify restart file was created
            if [ -f "restart" ]; then
                echo "Restart file created successfully"
                echo "SUCCESS" > "${status_dir}/${CASE_NO}"
            else
                echo "WARNING: restart file not found after Stage 1"
                echo "FAILED" > "${status_dir}/${CASE_NO}"
                return 1
            fi
        else
            echo "ERROR: Stage 1 failed"
            echo "FAILED" > "${status_dir}/${CASE_NO}"
            return 1
        fi

        echo "Case $CASE_NO Stage 1 complete"
    } > "$LOG_FILE" 2>&1

    return 0
}

export -f run_stage1_case

# ============================================================
# Run All Cases in Parallel (Stage 1)
# ============================================================
echo "============================================="
echo "Running $COMBINATION_COUNT Cases (Stage 1)"
echo "============================================="
echo "Max concurrent: ${MAX_CONCURRENT}"
echo "Each case: serial execution, tmax=${STAGE1_TMAX}"
echo ""

RUNNING_JOBS=0
PIDS=()

for param_file in "${CASE_FILES[@]}"; do
    # Wait if we've hit the concurrency limit
    while [ $RUNNING_JOBS -ge $MAX_CONCURRENT ]; do
        # Wait for any background job to finish
        wait -n 2>/dev/null || true
        ((RUNNING_JOBS--)) || true
    done

    # Launch case in background
    run_stage1_case "$param_file" "$STATUS_DIR" "$SCRIPT_DIR" "$STAGE1_TMAX" &
    PIDS+=($!)
    ((RUNNING_JOBS++)) || true

    # Extract case number for logging
    parse_param_file "$param_file"
    CASE_NO=$(get_param "CaseNo")
    echo "Launched case $CASE_NO (PID: ${PIDS[-1]}, running: $RUNNING_JOBS)"
done

echo ""
echo "All cases launched, waiting for completion..."
echo ""

# Wait for all remaining jobs
wait

# ============================================================
# Collect Results
# ============================================================
SUCCESSFUL_CASES=0
FAILED_CASES=0

for status_file in "${STATUS_DIR}"/*; do
    if [ -f "$status_file" ]; then
        status=$(cat "$status_file")
        case_no=$(basename "$status_file")
        if [ "$status" = "SUCCESS" ]; then
            ((SUCCESSFUL_CASES++)) || true
        else
            ((FAILED_CASES++)) || true
            echo "FAILED: Case $case_no"
        fi
    fi
done

# ============================================================
# Final Summary
# ============================================================
echo ""
echo "============================================="
echo "Stage 1 Complete"
echo "============================================="
echo "Job completed at: $(date)"
echo "Total cases: $COMBINATION_COUNT"
echo "Successful: $SUCCESSFUL_CASES"
echo "Failed: $FAILED_CASES"
echo ""
echo "Restart files location: simulationCases/*/restart"
echo ""
if [ $SUCCESSFUL_CASES -gt 0 ]; then
    echo "Next step: Submit Stage 2 with:"
    echo "  sbatch runSweepSnellius.sbatch"
fi
echo "============================================="

if [ $FAILED_CASES -gt 0 ]; then
    echo "WARNING: $FAILED_CASES case(s) failed. Check logs in simulationCases/*/stage1.log" >&2
    exit 1
fi

exit 0
