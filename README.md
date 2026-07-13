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
│   ├── Video-generic.py            Frame-by-frame visualization pipeline
│   └── render_contour_pulse.py     Render lightweight live contour interfaces
├── contourWorkflow/                 Bayesian contour campaign state machine
│   ├── contour_campaign.py         Propose, submit, collect, and gate iterations
│   ├── materialize_cases.py        Validate proposals against initial shapes
│   ├── result_quality.py            Compute KE/facet corruption evidence
│   └── run_one_contour_case.sh     Run and classify one OpenMP case
├── runSimulation.sh                 Single case runner (OpenMP/MPI)
├── runParameterSweep.sh             Parameter sweep runner
├── runPostProcess-Ncases.sh         Batch post-processing runner
├── default.params                   Single-case configuration
├── sweep.params                     Sweep configuration template
├── runSweepSnellius.sbatch          SLURM script for Snellius HPC
├── runSweepHamilton.sbatch          Legacy sequential MPI runner
├── runContourHamilton.sbatch        Packed 16-case Hamilton runner
└── runContourLocal.sh               Bounded local-systemd/OpenMP runner
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

### Bayesian contour campaign

The rearmable contour workflow uses batches of 16 simulations. Hamilton runs
one 16-case Slurm allocation; the local backend launches a user-systemd unit
and executes three cases at a time by default. Each case receives eight OpenMP
threads, so the workstation default is 24 concurrent threads. The local runner
enforces a hard 48-thread ceiling through `CONTOUR_MAX_THREADS` and uses no MPI.
The simulation detects the first persistent leading-tip detachment during
runtime, writes `classification.status`, and stops after a component above the
configured radius cutoff persists for three checks. It is labelled as an
injected drop only when its volume-weighted axial velocity is positive at
pinch-off. A zero or negative velocity is no-drop; later Rayleigh--Plateau
breakup along the jet is deliberately excluded. The contour-case observation
horizon is `t=1.0` by default.

The campaign controller requires a Bayesian Contour Predictor checkout with
the `--x-candidates` interface (AnjaliML/Bayesian-Contour-Predictor PR #3 or a
later release). Initialise from a canonical seed and exclude any known
configuration-confounded column before proposing:

```bash
module load python/3.10.8
python3 contourWorkflow/contour_campaign.py \
  --campaign-root /nobackup/$USER/drop-injection-confined \
  --project-root "$PWD" \
  --predictor-root /nobackup/$USER/Bayesian-Contour-Predictor \
  init --seed NumConfinementSweep-0.csv --exclude-x 8

python3 contourWorkflow/contour_campaign.py \
  --campaign-root /nobackup/$USER/drop-injection-confined \
  --project-root "$PWD" \
  --predictor-root /nobackup/$USER/Bayesian-Contour-Predictor \
  advance --submit
```

On a controlled Linux workstation, select the local backend explicitly:

```bash
python3 contourWorkflow/contour_campaign.py \
  --campaign-root /path/to/drop-injection-confined \
  --project-root "$PWD" \
  --predictor-root /path/to/Bayesian-Contour-Predictor \
  --backend local \
  advance --submit
```

`runContourLocal.sh` compiles against the ref-locked Basilisk sources under the
project and uses a host-built `qcc`. Set `CONTOUR_QCC` only when `qcc` is not on
the non-interactive `PATH`.

`advance --submit` is idempotent. It collects only complete 16-row result
tables, submits the next batch, stops for manual selection after iteration 8,
and stops permanently at the configured final iteration. Radius-ratio
proposals are confined to `simulationCases/DataFiles/`; post-hoc rounding is
rejected because it changes acquisition scores.

At the iteration-8 checkpoint, generate (but do not approve) the review file
with `propose-manual-batch`, edit it only for an explicit physics or
information-gain reason, then install it with `approve-manual-batch`. A failed
or timed-out allocation remains unresolved; after inspection, `retry --submit`
creates a clean `attempt-NN` directory containing only unresolved cases at
adjacent, non-colliding Oh values (1% steps by default). Every replacement is
recorded in `replacements.csv`; superseded evidence stays auditable but can
never be promoted under the new coordinate.
Collection merges resolved labels across attempts in the original 16-case order
without overwriting earlier evidence.

Contour runs can enable the feature-driven drill with
`CONTOUR_DRILL_AMR=1`. The controller begins at a configured lower level and
uses persistent target-region curvature demand to arm. A configurable
fixed-Lmax bootstrap protects the freshly joined neck before any coarsening.
Arming does not release Lmax through the singular focus: regional Lmax fires
only after the leading near-axis tip advances persistently. Full resolution is
then confined to the physical end-pinchoff band while the parent-bubble
exterior remains capped one level lower. The public defaults are conservative;
the drill is off unless a campaign explicitly enables it. Calibrate labels and first-drop radii against
fixed-level references before production. Rayleigh--Plateau component counts
are not drill triggers.

`CONTOUR_DRILL_REGIONAL_ONLY=1` is the conservative mode: the complete
wave/focus/jet band remains at `MAXlevel` for the whole run and only the parent
bubble exterior is capped. Use this when dynamic pre-jet coarsening changes a
boundary label or first-drop size.

For a bounded unattended workstation campaign, initialise a fresh campaign
with explicit numerical and acquisition settings, then run the shell driver:

```bash
python3 contourWorkflow/contour_campaign.py \
  --campaign-root /path/to/confined-l11 \
  --project-root "$PWD" \
  --predictor-root /path/to/Bayesian-Contour-Predictor \
  --backend inline init \
  --seed /path/to/resolved-seed.csv \
  --final-iteration 20 --no-manual-checkpoint \
  --case-id-start 6000 --n-new 14 --n-repeats 2 \
  --max-level 11 --drop-radius-min 0.015625 \
  --workers 3 --threads-per-case 8 --max-threads 48 \
  --unit-prefix dropinj-l11 --allow-unbracketed-edges \
  --posterior-samples 64

./runContourCampaignLoop.sh \
  /path/to/confined-l11 \
  /path/to/Bayesian-Contour-Predictor
```

Run the driver in a detached `tmux` session. Its inline backend avoids depending
on login-session-scoped user-systemd services. The driver permits one selective
retry of unresolved cases, then stops with
`needs_attention`. It never reruns resolved cases. Full per-batch measurements,
including first-persistent-detachment radius and volume, are retained under
`measurements/`; the four-column `completed/` tables remain predictor input.

New Hamilton result rows include `max_ke`, `facet_lines`, `quality_state`, and
`quality_reason`. A row fails the conservative corruption gate when its maximum
kinetic energy exceeds 1000, its latest interface has more than 8000 nonblank
facet lines, or either evidence source is missing/non-finite. Failed rows remain
auditable in their attempt but are not promoted or counted as resolved. Rows
without quality columns are backfilled from retained attempt case evidence;
archived legacy tables with no retained `cases/` tree remain valid unless
manually quarantined.

Manual review can quarantine one exact piece of evidence by creating
`quality-quarantine.csv` in the campaign root:

```csv
iteration,attempt,caseId,reason
1,1,5000,corrupt interface confirmed by visual review
```

The exclusion applies only to that iteration, attempt, and case; a clean retry
of the same case remains eligible for collection.

For a full-node launch test, override the production header without editing
the file:

```bash
sbatch -p test --time=00:15:00 runContourHamilton.sbatch \
  /path/to/campaign/iterations/iteration-01
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
qcc -O2 -Wall -disable-dimensions postProcess/getFacet.c -o postProcess/getFacet -lm
qcc -O2 -Wall -disable-dimensions postProcess/getData-generic.c -o postProcess/getData-generic -lm
qcc -O2 -Wall -disable-dimensions postProcess/getCOM.c -o postProcess/getCOM -lm
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
- `OhOut`: Ohnesorge number for outer fluid, based on the small-bubble radius `R_s` (e.g., 1e-2)
- `RhoIn`: Density ratio inner/outer (e.g., 1e-3)
- `Rr`: Radius ratio `R_l/R_s` with `R_s = 1` and `R_l = Rr`; available values: 1.00, 1.50, 2.00, 4.00, 8.00 (matching DataFiles)
- `MAXlevel`: Maximum refinement level (e.g., 10)
- `tmax`: Maximum simulation time (e.g., 40.0)
- `zWall`: Wall position (e.g., 0.01)

## Running Manually (Legacy)

There are two ways to run the codes:

1. Using the vanilla basilisk method:

```bash
cd simulationCases
qcc -I../src-local -O2 -Wall -disable-dimensions coalescenceBubble.c -o coalescenceBubble -lm
./coalescenceBubble 1e-2 1e-3 1.00 12 1.5 0.05
```

2. Using the makefile (can be interactively run using bview browser):

```bash
cd simulationCases
CFLAGS=-DDISPLAY=-1 make coalescenceBubble.tst
```
Check the localhost on coalescenceBubble/display.html. something like: [https://basilisk.fr/three.js/editor/index.html?ws://localhost:7100](https://basilisk.fr/three.js/editor/index.html?ws://localhost:7100) and run interactively.

### To run using openMP, please use the flag -fopenmp

```bash
cd simulationCases
qcc -I../src-local -O2 -Wall -disable-dimensions coalescenceBubble.c -o coalescenceBubble -lm -fopenmp
export OMP_NUM_THREADS=8
./coalescenceBubble 1e-2 1e-3 1.00 12 1.5 0.05
```

**Note:** The code will not directly work with openmpi (with -D_MPI flag and mpirun). To do that, please follow the procedure we use: 

1. Run the following for a few timesteps (stop using tmax=1e-2 or so)
```bash
cd simulationCases
qcc -I../src-local -O2 -Wall -disable-dimensions coalescenceBubble.c -o coalescenceBubble -lm -fopenmp
export OMP_NUM_THREADS=8
./coalescenceBubble 1e-2 1e-3 1.00 12 1e-2 0.05
```

This will generate a "restart" file

2. Do not delete the restart file. Now you can use mpirun, like:

```bash
cd simulationCases
CC99='mpicc -std=c99' qcc -I../src-local -Wall -O2 -D_MPI=1 -disable-dimensions coalescenceBubble.c -o coalescenceBubble -lm
mpirun -np 8 ./coalescenceBubble 1e-2 1e-3 1.00 12 1.5 0.05
```
