# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a computational fluid dynamics project studying asymmetric bubble coalescence using the Basilisk C framework. The simulation focuses on size asymmetry while maintaining axial symmetry using the Volume-of-Fluid (VOF) method.

## Build Commands

### Standard Build
```bash
qcc -O2 -Wall -disable-dimensions coalescenceBubble.c -o coalescenceBubble -lm
./coalescenceBubble <OhOut> <RhoIn> <Rr> <MAXlevel> <tmax> <zWall>
```

### Interactive Build with Browser Visualization
```bash
cd simulationCases
CFLAGS=-DDISPLAY=-1 make coalescenceBubble.tst
# Open http://basilisk.fr/three.js/editor/index.html?ws://localhost:7100
```

### OpenMP Parallel Build
```bash
qcc -O2 -Wall -disable-dimensions coalescenceBubble.c -o coalescenceBubble -lm -fopenmp
export OMP_NUM_THREADS=8
./coalescenceBubble <parameters>
```

### MPI Build (requires dump file)
```bash
# First run with OpenMP to generate dump file
qcc -O2 -Wall -disable-dimensions coalescenceBubble.c -o coalescenceBubble -lm -fopenmp
export OMP_NUM_THREADS=8
./coalescenceBubble <parameters> # Run briefly with tmax=1e-2

# Then compile and run with MPI
CC99='mpicc -std=c99' qcc -Wall -O2 -D_MPI=1 -disable-dimensions coalescenceBubble.c -o coalescenceBubble -lm
mpirun -np 8 coalescenceBubble <parameters>
```

## Command Line Parameters

1. `OhOut`: Ohnesorge number for outer fluid (e.g., 1e-2)
2. `RhoIn`: Density ratio (e.g., 1e-3)
3. `Rr`: Radius ratio - asymmetry parameter (e.g., 1.0)
4. `MAXlevel`: Maximum refinement level (e.g., 12)
5. `tmax`: Maximum simulation time (e.g., 40.0)
6. `zWall`: Wall position (e.g., 0.01)

## Architecture

The codebase uses Basilisk's adaptive mesh refinement with two-phase flow tracking:

- **Main simulation**: `simulationCases/coalescenceBubble.c` - implements bubble coalescence physics
- **Custom headers**: Located in `src-local/`, particularly `two-phase-tag.h` for interface tracking
- **Basilisk path**: Set via `.project_config` file
- **Data flow**: Initial conditions → Simulation → Snapshots → Post-processing

Key physics components:
- Uses both `f` and `ftag` fields for tracking different fluid interfaces
- Adaptive mesh refinement based on VOF and velocity error tolerances
- Logs kinetic energy, center of mass, interface positions, and velocities

## Key Directories

- `simulationCases/`: Main simulation code and Makefile
- `src-local/`: Custom headers extending Basilisk functionality
- `DataFiles/`: Initial condition files (InitialConditionRr-*.dat)
- `intermediate/`: Output directory for simulation snapshots
- `postProcess/`: Python and C tools for analyzing results

## Testing

Run tests using Basilisk's testing framework:
```bash
cd simulationCases
make coalescenceBubble.tst
```

## HPC Submission

Use `runCode.sh` for SLURM job submission. The script runs parameter sweeps with OpenMP parallelization.

## Post-Processing

- Python scripts in `postProcess/` for facet and center-of-mass analysis
- Video generation: `ffmpeg -r 30 -f image2 -s 1920x1080 -i %*.jpeg -c:v h264 -crf 1 -pix_fmt yuv420p video.mp4`