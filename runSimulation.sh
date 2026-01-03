#!/bin/bash
# runSimulation.sh - Run single bubble coalescence simulation from root directory
# Creates case folder in simulationCases/<CaseNo>/ and runs simulation there
#
# IMPORTANT: This script uses TWO-STAGE EXECUTION because distance.h is
# incompatible with MPI. Stage 1 generates the initial condition dump file
# using OpenMP, then Stage 2 runs the full simulation with MPI.

set -euo pipefail  # Exit on error, unset variables, pipeline failures

# ============================================================
# Configuration
# ============================================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Source project configuration
if [ -f "${SCRIPT_DIR}/.project_config" ]; then
    source "${SCRIPT_DIR}/.project_config"
else
    echo "ERROR: .project_config not found" >&2
    exit 1
fi

# Source parameter parsing library
if [ -f "${SCRIPT_DIR}/src-local/parse_params.sh" ]; then
    source "${SCRIPT_DIR}/src-local/parse_params.sh"
else
    echo "ERROR: src-local/parse_params.sh not found" >&2
    exit 1
fi

# ============================================================
# Usage Information
# ============================================================
usage() {
    cat <<EOF
Usage: $0 [OPTIONS] [params_file]

Run single bubble coalescence simulation from root directory.
Creates case folder in simulationCases/<CaseNo>/ based on parameter file.

This script uses TWO-STAGE EXECUTION:
  Stage 1: OpenMP compilation and brief run to generate dump file
           (Required because distance.h is incompatible with MPI)
  Stage 2: MPI compilation and full simulation from dump file

Options:
    -c, --compile-only    Compile but don't run simulation
    -d, --debug           Compile with debug flags (-g -DTRASH=1)
    -m, --mpi             Enable MPI parallel execution (default: serial)
    --cores N             Number of MPI cores (default: 4, requires --mpi)
    --omp-threads N       Number of OpenMP threads for Stage 1 (default: 4)
    --skip-stage1         Skip Stage 1 (assume dump file exists)
    -v, --verbose         Verbose output
    -h, --help           Show this help message

Parameter file mode (default):
    $0 default.params

If no parameter file specified, uses default.params from current directory.

Environment variables:
    QCC_FLAGS     Additional qcc compiler flags

Examples:
    # Run with default parameters (serial, both stages)
    $0

    # Run with MPI parallel execution (4 cores)
    $0 --mpi

    # Run with MPI using 8 cores
    $0 --mpi --cores 8 default.params

    # Skip Stage 1 if dump already exists
    $0 --mpi --skip-stage1

    # Compile only (check for errors)
    $0 --compile-only

For more information, see README.md
EOF
}

# ============================================================
# Parse Command Line Options
# ============================================================
COMPILE_ONLY=0
DEBUG_FLAGS=""
VERBOSE=0
MPI_ENABLED=0
MPI_CORES=4
OMP_THREADS=4
SKIP_STAGE1=0
QCC_FLAGS="${QCC_FLAGS:-}"

while [[ $# -gt 0 ]]; do
    case $1 in
        -c|--compile-only)
            COMPILE_ONLY=1
            shift
            ;;
        -d|--debug)
            DEBUG_FLAGS="-g -DTRASH=1"
            shift
            ;;
        -m|--mpi)
            MPI_ENABLED=1
            shift
            ;;
        --cores)
            MPI_CORES="$2"
            if ! [[ "$MPI_CORES" =~ ^[0-9]+$ ]] || [ "$MPI_CORES" -lt 1 ]; then
                echo "ERROR: --cores requires a positive integer, got: $MPI_CORES" >&2
                exit 1
            fi
            shift 2
            ;;
        --omp-threads)
            OMP_THREADS="$2"
            if ! [[ "$OMP_THREADS" =~ ^[0-9]+$ ]] || [ "$OMP_THREADS" -lt 1 ]; then
                echo "ERROR: --omp-threads requires a positive integer, got: $OMP_THREADS" >&2
                exit 1
            fi
            shift 2
            ;;
        --skip-stage1)
            SKIP_STAGE1=1
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
            break
            ;;
    esac
done

# ============================================================
# Detect OS and Verify MPI
# ============================================================
OS_TYPE=$(uname -s)

# Verify MPI tools if MPI is enabled
if [ $MPI_ENABLED -eq 1 ]; then
    if ! command -v mpicc &> /dev/null; then
        echo "ERROR: mpicc not found. MPI compilation requires mpicc (OpenMPI or MPICH)." >&2
        echo "       Install MPI tools or run without --mpi flag for serial execution." >&2
        exit 1
    fi
    if ! command -v mpirun &> /dev/null; then
        echo "ERROR: mpirun not found. MPI execution requires mpirun (OpenMPI or MPICH)." >&2
        echo "       Install MPI tools or run without --mpi flag for serial execution." >&2
        exit 1
    fi
fi

# ============================================================
# Determine Parameter File
# ============================================================
PARAM_FILE="${1:-default.params}"

if [ ! -f "$PARAM_FILE" ]; then
    echo "ERROR: Parameter file not found: $PARAM_FILE" >&2
    exit 1
fi

[ $VERBOSE -eq 1 ] && echo "Parameter file: $PARAM_FILE"

# ============================================================
# Parse Parameters
# ============================================================
parse_param_file "$PARAM_FILE"

CASE_NO=$(get_param "CaseNo")
OhOut=$(get_param "OhOut" "1e-2")
RhoIn=$(get_param "RhoIn" "1e-3")
Rr=$(get_param "Rr" "1.0")
MAXlevel=$(get_param "MAXlevel" "10")
tmax=$(get_param "tmax" "40.0")
zWall=$(get_param "zWall" "0.01")

if [ -z "$CASE_NO" ]; then
    echo "ERROR: CaseNo not found in parameter file" >&2
    exit 1
fi

# Validate CaseNo is 4 digits
if ! [[ "$CASE_NO" =~ ^[0-9]{4}$ ]] || [ "$CASE_NO" -lt 1000 ] || [ "$CASE_NO" -gt 9999 ]; then
    echo "ERROR: CaseNo must be 4-digit (1000-9999), got: $CASE_NO" >&2
    exit 1
fi

CASE_DIR="simulationCases/${CASE_NO}"

echo "========================================="
echo "Bubble Coalescence Simulation"
echo "========================================="
echo "Case Number: $CASE_NO"
echo "Case Directory: $CASE_DIR"
echo "Parameter File: $PARAM_FILE"
echo ""
echo "Physical Parameters:"
echo "  OhOut=$OhOut, RhoIn=$RhoIn, Rr=$Rr"
echo "  MAXlevel=$MAXlevel, tmax=$tmax, zWall=$zWall"
echo ""
if [ $MPI_ENABLED -eq 1 ]; then
    echo "Execution Mode: MPI Parallel ($MPI_CORES cores)"
else
    echo "Execution Mode: Serial"
fi
echo ""

# ============================================================
# Create Case Directory
# ============================================================
if [ ! -d "$CASE_DIR" ]; then
    echo "Creating case directory: $CASE_DIR"
    mkdir -p "$CASE_DIR"
else
    echo "Case directory exists"
fi

# Copy parameter file to case directory for record keeping
cp "$PARAM_FILE" "$CASE_DIR/case.params"

# Change to case directory
cd "$CASE_DIR"
[ $VERBOSE -eq 1 ] && echo "Working directory: $(pwd)"

# ============================================================
# Source File Setup
# ============================================================
SRC_FILE_ORIG="../coalescenceBubble.c"
SRC_FILE_LOCAL="coalescenceBubble.c"
EXECUTABLE="coalescenceBubble"

# Check if source file exists
if [ ! -f "$SRC_FILE_ORIG" ]; then
    echo "ERROR: Source file $SRC_FILE_ORIG not found" >&2
    exit 1
fi

# Copy source file to case directory
cp "$SRC_FILE_ORIG" "$SRC_FILE_LOCAL"
echo "Copied source file to case directory"

# ============================================================
# Stage 1: Generate Dump File with OpenMP
# ============================================================
# This stage is required because distance.h (used for initial conditions)
# is incompatible with MPI but works with OpenMP.

if [ ! -f "dump" ] && [ $SKIP_STAGE1 -eq 0 ]; then
    echo ""
    echo "========================================="
    echo "Stage 1: Generate Initial Condition (OpenMP)"
    echo "========================================="
    echo "Compiling with OpenMP..."

    [ $VERBOSE -eq 1 ] && echo "Compiler: qcc"
    [ $VERBOSE -eq 1 ] && echo "Include paths: -I../../src-local"
    [ $VERBOSE -eq 1 ] && echo "Flags: -O2 -Wall -disable-dimensions -fopenmp $DEBUG_FLAGS $QCC_FLAGS"

    qcc -I../../src-local \
        -O2 -Wall -disable-dimensions -fopenmp \
        $DEBUG_FLAGS $QCC_FLAGS \
        "$SRC_FILE_LOCAL" -o "${EXECUTABLE}_omp" -lm

    if [ $? -ne 0 ]; then
        echo "ERROR: OpenMP compilation failed" >&2
        exit 1
    fi

    echo "OpenMP compilation successful"
    echo ""
    echo "Running briefly to generate dump file..."
    echo "  OMP_NUM_THREADS=$OMP_THREADS"
    echo "  Command: ./${EXECUTABLE}_omp $OhOut $RhoIn $Rr $MAXlevel 0.01 $zWall"

    export OMP_NUM_THREADS=$OMP_THREADS
    ./${EXECUTABLE}_omp $OhOut $RhoIn $Rr $MAXlevel 0.01 $zWall

    if [ ! -f "dump" ]; then
        echo "ERROR: Stage 1 failed - dump file was not created" >&2
        exit 1
    fi

    echo "Stage 1 complete: dump file created"
elif [ -f "dump" ]; then
    echo "Dump file already exists - skipping Stage 1"
elif [ $SKIP_STAGE1 -eq 1 ]; then
    echo "Stage 1 skipped by user (--skip-stage1)"
    if [ ! -f "dump" ]; then
        echo "WARNING: dump file does not exist - Stage 2 may fail" >&2
    fi
fi

# Exit if compile-only mode
if [ $COMPILE_ONLY -eq 1 ]; then
    echo ""
    echo "Compile-only mode: Stopping here"
    cd ../..
    exit 0
fi

# ============================================================
# Stage 2: Full Simulation
# ============================================================
echo ""
echo "========================================="
echo "Stage 2: Full Simulation"
echo "========================================="

if [ $MPI_ENABLED -eq 1 ]; then
    # MPI parallel compilation
    echo "Compiling with MPI..."

    if [ "$OS_TYPE" = "Darwin" ]; then
        # macOS
        [ $VERBOSE -eq 1 ] && echo "Compiler: CC99='mpicc -std=c99' qcc"
        [ $VERBOSE -eq 1 ] && echo "Include paths: -I../../src-local"
        [ $VERBOSE -eq 1 ] && echo "Flags: -Wall -O2 -D_MPI=1 -disable-dimensions $DEBUG_FLAGS $QCC_FLAGS"

        CC99='mpicc -std=c99' qcc -I../../src-local \
            -Wall -O2 -D_MPI=1 -disable-dimensions \
            $DEBUG_FLAGS $QCC_FLAGS \
            "$SRC_FILE_LOCAL" -o "$EXECUTABLE" -lm
    else
        # Linux
        [ $VERBOSE -eq 1 ] && echo "Compiler: CC99='mpicc -std=c99 -D_GNU_SOURCE=1' qcc"
        [ $VERBOSE -eq 1 ] && echo "Include paths: -I../../src-local"
        [ $VERBOSE -eq 1 ] && echo "Flags: -Wall -O2 -D_MPI=1 -disable-dimensions $DEBUG_FLAGS $QCC_FLAGS"

        CC99='mpicc -std=c99 -D_GNU_SOURCE=1' qcc -I../../src-local \
            -Wall -O2 -D_MPI=1 -disable-dimensions \
            $DEBUG_FLAGS $QCC_FLAGS \
            "$SRC_FILE_LOCAL" -o "$EXECUTABLE" -lm
    fi
else
    # Serial compilation with OpenMP
    echo "Compiling for serial execution (with OpenMP)..."

    [ $VERBOSE -eq 1 ] && echo "Compiler: qcc"
    [ $VERBOSE -eq 1 ] && echo "Include paths: -I../../src-local"
    [ $VERBOSE -eq 1 ] && echo "Flags: -O2 -Wall -disable-dimensions -fopenmp $DEBUG_FLAGS $QCC_FLAGS"

    qcc -I../../src-local \
        -O2 -Wall -disable-dimensions -fopenmp \
        $DEBUG_FLAGS $QCC_FLAGS \
        "$SRC_FILE_LOCAL" -o "$EXECUTABLE" -lm
fi

if [ $? -ne 0 ]; then
    echo "ERROR: Compilation failed" >&2
    exit 1
fi

echo "Compilation successful: $EXECUTABLE"

# ============================================================
# Execution
# ============================================================
echo ""
echo "Starting full simulation..."
echo "  Command args: $OhOut $RhoIn $Rr $MAXlevel $tmax $zWall"
echo "========================================="

# Run simulation
if [ $MPI_ENABLED -eq 1 ]; then
    [ $VERBOSE -eq 1 ] && echo "Command: mpirun -np $MPI_CORES ./$EXECUTABLE $OhOut $RhoIn $Rr $MAXlevel $tmax $zWall"
    mpirun -np $MPI_CORES ./$EXECUTABLE $OhOut $RhoIn $Rr $MAXlevel $tmax $zWall
else
    export OMP_NUM_THREADS=$OMP_THREADS
    [ $VERBOSE -eq 1 ] && echo "Command: ./$EXECUTABLE $OhOut $RhoIn $Rr $MAXlevel $tmax $zWall"
    ./$EXECUTABLE $OhOut $RhoIn $Rr $MAXlevel $tmax $zWall
fi

EXIT_CODE=$?

echo "========================================="
if [ $EXIT_CODE -eq 0 ]; then
    echo "Simulation completed successfully"
    echo "Output location: $CASE_DIR"
else
    echo "Simulation failed with exit code $EXIT_CODE"
fi
echo "========================================="

# Return to root directory
cd ../..

exit $EXIT_CODE
