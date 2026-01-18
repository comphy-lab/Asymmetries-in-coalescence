# Asymmetries-in-coalescence
Asymmetries in coalescence: size asymmetry. Still axially symmetric.

## Basilisk (Required)

First-time install (or reinstall):
```bash
curl -sL https://raw.githubusercontent.com/comphy-lab/basilisk-C/main/reset_install_basilisk-ref-locked.sh | bash -s -- --ref=v2026-01-13 --hard
```

Subsequent runs (reuses existing `basilisk/` if same ref):
```bash
curl -sL https://raw.githubusercontent.com/comphy-lab/basilisk-C/main/reset_install_basilisk-ref-locked.sh | bash -s -- --ref=v2026-01-13
```

> **Note**: Replace `v2026-01-13` with the [latest release tag](https://github.com/comphy-lab/basilisk-C/releases).

## Repository Structure

```
├── simulationCases/                 Main simulation code
│   ├── coalescenceBubble.c         Primary simulation (production runs)
│   ├── coalescenceBubble-tag.c     Extended version with shape tracking
├── src-local/                       Custom Basilisk headers
│   ├── two-phase-tag.h             Two-phase flow with interface tagging
│   └── parse_params.sh             Parameter file parsing library
├── postProcess/                     Post-processing tools
│   ├── getData-generic.c           Field extraction on structured grids
│   ├── getFacet.c                  Interface geometry extraction
│   ├── getCOM.c                    Center of mass extraction
│   └── Video-generic.py            Frame-by-frame visualization pipeline
├── runSimulation.sh                 Single case runner (OpenMP/MPI)
├── runParameterSweep.sh             Parameter sweep runner
├── runPostProcess-Ncases.sh         Batch post-processing runner
├── default.params                   Single-case configuration
├── sweep.params                     Sweep configuration template
├── runSweepSnellius.sbatch          SLURM script for Snellius HPC
└── runSweepHamilton.sbatch          SLURM script for Hamilton HPC
```

## Simulation Files

This project contains two simulation files:

### coalescenceBubble.c (Primary)
The main simulation file used for all production runs. Outputs:
- `i dt t ke Xc Vcm` (6 columns)

### coalescenceBubble-tag.c (Optional)
An extended version with interface tagging for detailed shape tracking. Uses `tag.h` to identify and track only the largest connected bubble region (filters out satellite droplets). Outputs additional geometric measurements:
- `i dt t ke Xc Vcm Re ZNp ZSp` (9 columns)
- `Re`: Equatorial radius at center of mass
- `ZNp`: North pole position (positive x on axis)
- `ZSp`: South pole position (negative x on axis)

**Note:** All cases in this project are run with `coalescenceBubble.c`. The `-tag.c` variant is provided as an optional alternative for cases requiring detailed shape tracking.

## Why coalescenceBubble.c (not coalescenceBubble-tag.c)

The running scripts use `coalescenceBubble.c` because:

1. **`distance.h` is incompatible with MPI** - The `distance.h` header (used for computing initial conditions from shape files) cannot be compiled with `-D_MPI=1`
2. **OpenMP is compatible** - `distance.h` works fine with OpenMP (`-fopenmp`)
3. **Two-stage execution** - We first run briefly with OpenMP to generate the restart file with initial conditions, then run the full simulation with MPI (which restores from the restart file, bypassing `distance.h`)

The `coalescenceBubble-tag.c` file additionally uses `tag.h` for tracking shape metrics, which adds complexity not needed for basic coalescence studies.

## Running with Scripts (Recommended)

### Single Case
```bash
# Serial execution (both stages: OpenMP for restart, then serial)
./runSimulation.sh default.params
```

```bash
# MPI execution (Stage 1: OpenMP, Stage 2: MPI)
./runSimulation.sh --mpi --cores 8 default.params
```

```bash
# Compile only (check for errors)
./runSimulation.sh --compile-only default.params
```

### Parameter Sweep
```bash
# Dry run (see parameter combinations)
./runParameterSweep.sh --dry-run
```

```bash
# Run all combinations (serial)
./runParameterSweep.sh sweep.params
```

```bash
# Run with MPI (8 cores per case)
./runParameterSweep.sh --mpi --cores 8 sweep.params
```

### HPC (Snellius)
```bash
# Submit parameter sweep to SLURM
sbatch runSweepSnellius.sbatch
```

### Post-Processing
```bash
# Process multiple cases with default settings
./runPostProcess-Ncases.sh 3000 3001 3002
```

```bash
# Process a range of cases
./runPostProcess-Ncases.sh 3000-3010
```

```bash
# Process with 8 CPUs and custom snapshot count
./runPostProcess-Ncases.sh --CPUs 8 --nGFS 100 3000
```

```bash
# Skip video encoding (only generate frames)
./runPostProcess-Ncases.sh --skip-video-encode 3000
```

```bash
# Dry run to preview commands
./runPostProcess-Ncases.sh --dry-run 3000
```

**C Helper Tools** (must be compiled before running):
```bash
qcc -O2 -Wall postProcess/getFacet.c -o postProcess/getFacet -lm
qcc -O2 -Wall postProcess/getData-generic.c -o postProcess/getData-generic -lm
qcc -O2 -Wall postProcess/getCOM.c -o postProcess/getCOM -lm
```

- `getFacet`: Extracts interface facets using PLIC reconstruction
- `getData-generic`: Samples velocity/strain-rate fields on structured grids
- `getCOM`: Computes center of mass position and velocity

**Output locations:**
- `simulationCases/<CaseNo>/Video/` - PNG frames
- `simulationCases/<CaseNo>/<CaseNo>_COMData.csv` - COM time series
- `simulationCases/<CaseNo>/<CaseNo>.mp4` - Encoded video

### Command Line Parameters
The simulation takes 6 arguments: `OhOut RhoIn Rr MAXlevel tmax zWall`
- `OhOut`: Ohnesorge number for outer fluid (e.g., 1e-2)
- `RhoIn`: Density ratio inner/outer (e.g., 1e-3)
- `Rr`: Radius ratio - available values: 1.00, 1.50, 2.00, 4.00, 8.00 (matching DataFiles)
- `MAXlevel`: Maximum refinement level (e.g., 10)
- `tmax`: Maximum simulation time (e.g., 40.0)
- `zWall`: Wall position (e.g., 0.01)

## Running Manually (Legacy)

There are two ways to run the codes:

1. Using the vanilla basilisk method:

```bash
qcc -O2 -Wall -disable-dimensions coalescenceBubble.c -o coalescenceBubble -lm 
./coalescenceBubble
```

2. Using the makefile (can be interactively run using bview browser):

```bash
CFLAGS=-DDISPLAY=-1 make coalescenceBubble.tst
```
Check the localhost on coalescenceBubble/display.html. something like: [https://basilisk.fr/three.js/editor/index.html?ws://localhost:7100](https://basilisk.fr/three.js/editor/index.html?ws://localhost:7100) and run interactively.

### To run using openMP, please use the flag -fopenmp

```bash
qcc -O2 -Wall -disable-dimensions coalescenceBubble.c -o coalescenceBubble -lm -fopenmp
export OMP_NUM_THREADS=8
./coalescenceBubble
```

**Note:** The code will not directly work with openmpi (with -D_MPI flag and mpirun). To do that, please follow the procedure we use: 

1. Run the following for a few timesteps (stop using tmax=1e-2 or so)
```bash
qcc -O2 -Wall -disable-dimensions coalescenceBubble.c -o coalescenceBubble -lm -fopenmp
export OMP_NUM_THREADS=8
./coalescenceBubble
```

This will generate a "restart" file

2. Do not delete the restart file. Now you can use mpirun, like:

```bash
CC99='mpicc -std=c99' qcc -Wall -O2 -D_MPI=1 -disable-dimensions coalescenceBubble.c -o coalescenceBubble -lm
mpirun -np 8 coalescenceBubble
```
