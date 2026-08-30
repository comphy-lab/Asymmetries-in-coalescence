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
# Open https://basilisk.fr/three.js/editor/index.html?ws://localhost:7100
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
- `docs/`: Authored public reports (Markdown, LaTeX, PDF). Written as
  eventually public; not the project website.
- `.github/docs/`: Generated GitHub Pages site. Built by
  `.github/scripts/build.sh` and published by `.github/workflows/deploy.yml`.
  Never hand-edit it — change the generator instead.

## Phase convention (repository standard)

`f = 1` is the BUBBLE/GAS phase in every code in this repository — solvers,
headers and post-processing alike. `comphy-lab/Bursting-Bubble` (and the
singular-bursting-bubbles project built on it) uses the opposite convention
(`f = 1` is liquid), so anything ported from there must have its VOF sense
flipped as part of the port and the flip marked at the port site.
Ported-and-flipped so far: `src-local/jetFoot.h` (from `getBase.c` /
`getJetFoot.c`; its former `JETFOOT_F_IS_LIQUID` toggle is removed so the
convention cannot fork) and the drop-map instrumentation in
`simulationCases/burstingBubbleInfiniteRr.c`. Every liquid measurement in
this repository is a `1 - f` measurement.

## Solver stack (repository standard)

`simulationCases/coalescenceBubble.c` offers three mutually exclusive stacks
and enforces the exclusions with `#error`:

| flags | stack |
|---|---|
| (none) | 1: `FILTERED` property averaging + `navier-stokes/double-projection.h` |
| `-DSINGLE_PROJECTION=1` | 1S: `FILTERED`, single projection |
| `-DUSE_CONSERVING=1` | 2: `navier-stokes/conserving.h`, no `FILTERED` |

Stack 1 is the default and the production choice for asymmetric coalescence.

Two exclusions, with different standing (verified against Basilisk
`v2026-07-20` source, 2026-08-30):

- **`conserving.h` + `double-projection.h` is inconsistent by construction.**
  `centered.h` declares the `vof` event (where `conserving.h` advects
  momentum under `stokes = true`) before `advection_term`, where
  `double-projection.h` snapshots its baseline `Af`. The advective update
  therefore never enters the second projection, and the centred pressure
  gradient is built from a pressure missing that term.
- **`FILTERED` + `conserving.h` is a project policy, not an upstream rule.**
  Basilisk compiles it, and `coalescenceBubble-tag.c` (swimming-bubbles)
  uses it. The objection: `conserving.h` reconstructs
  `u = (q1+q2)/rho(f)` with the sharp volume fraction, while `FILTERED`
  hands the projection and viscous operators properties built from the
  smoothed `sf` — one momentum balance, two densities. No upstream example
  or test pairs them. Empirically here, the drill solver (exactly this
  pairing, single projection) diverged at the capillary-focusing event on
  2 of 16 ladder rungs at MAXlevel 13 and measured 1.55–1.77× lower peak
  jet speed than Stack 1 at matched Oh and Δ.
  `-DALLOW_FILTERED_CONSERVING` overrides the guard for controlled
  experiments; the banner then reports
  `filtered+conserving (ALLOW_FILTERED_CONSERVING)`.

For calibration: *no* composite stack here is upstream-endorsed for two-phase
work — upstream pairs `FILTERED` only with plain `two-phase.h`
(`rising.c`, `oscillation.c`), `conserving.h` only unfiltered
(`axi_rising_bubble.c`, `rt.c`, `tangaroa.c`), and `double-projection.h`
only single-phase (`starting.c`). Stack choices are settled by matched-Δ
cross-checks and Δ-convergence, not by appeal to upstream authority.

Any solver imported from another project must be checked against the table
above before its results are put on the same axes as this repository's — the
drill solver `burstingBubble-drillResolution.c`, used for the 2026-08 bulk
drop-map ladder, is `FILTERED` + `conserving.h` with a single projection.

Cost note (measured 2026-08-30, case 6554): pure `conserving.h` at
MAXlevel 13 runs with dt pinned ≈13× below Stack 1's at the same cell size
and settings — timestep-limited, not NITERMAX-limited; no convergence
warnings. A conserving ladder at this resolution costs ~13× Stack 1's.

## Gas-to-liquid viscosity ratio

`MuRin` (equivalently $Oh_g/Oh_l$) is **argv 26**, defaulting to `1e-2`, the
value the manuscript states. It was a compiled constant until 2026-08-30,
which is how the drop-map ladder came to run at `2e-2` — through a separate
parameter file — without the difference appearing anywhere in the argv, the
run directory or the logs. Pass it explicitly on any campaign whose points
will be plotted against another campaign's.

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

## Development Workflow

### Running Single Simulation
```bash
cd simulationCases
./runCode.sh coalescenceBubble
```

### Visualization
The simulation outputs:
- `dump` file: Current state for restart
- `intermediate/snapshot-*.` files: Time series data for post-processing
- `log` file: Time evolution of kinetic energy, center of mass, and velocity

### Error Tolerances
- VOF error: `fErr = 1e-3`
- Velocity error: `VelErr = 1e-2`
- Position tolerance: `TOL = 1e-2`

## Simulation Physics

The code simulates coalescence of two bubbles with:
- Axisymmetric geometry (`axi.h`)
- Two-phase flow with surface tension
- Adaptive mesh refinement (AMR) with levels from `MAXlevel-6` to `MAXlevel`
- Domain size automatically calculated as `Ldomain = zWall + 2 + 2*Rr + 4.0`
- Origin shifted to `(-2.0-zWall, 0.0)`

Initial conditions are loaded from pre-computed data files based on the radius ratio `Rr`.