/**
# Coalescence of Bubbles

Axisymmetric simulation of two-bubble coalescence with size asymmetry using
the Volume-of-Fluid (VOF) method. The simulation tracks the interface between
gas (inside bubble, $f=1$) and liquid (outside, $f=0$) phases.

## Scientific Context

This model isolates the hydrodynamic event responsible for electrolyte
entrainment in gas-evolving electrochemical systems: coalescence between a
large "parent" bubble and a smaller bubble. The key physical pathway is
capillary-wave focusing on the smaller bubble that generates a Worthington
jet and droplet pinch-off into the parent bubble.

## Hydrodynamic Pathway (Captured Here)

1. First contact creates a neck that expands rapidly.
2. The neck launches capillary waves over the small bubble.
3. Waves focus at the small bubble's south pole.
4. The leading jet tip end-pinches; it counts as injection only when the
   detached component is still moving into the larger bubble. Later
   Rayleigh--Plateau breakup of the jet is not the classified event.

## Physical Setup

Two initially spherical bubbles of different sizes are placed along the axis
of symmetry. The smaller bubble has radius 1 (reference scale), while the
larger bubble has radius `Rr` (the radius ratio parameter). Surface tension
drives the coalescence dynamics.

## Governing Equations

The incompressible Navier-Stokes equations with surface tension:

$$\nabla \cdot \mathbf{u} = 0$$

$$\rho\left(\frac{\partial \mathbf{u}}{\partial t} + \mathbf{u} \cdot \nabla \mathbf{u}\right) = -\nabla p + \nabla \cdot (2\mu \mathbf{D}) + \sigma \kappa \delta_s \mathbf{n}$$

where $\mathbf{D}$ is the strain-rate tensor, $\sigma$ is surface tension,
$\kappa$ is interface curvature, and $\delta_s$ is the interface delta function.

## Usage

```
./coalescenceBubble <OhOut> <RhoIn> <Rr> <MAXlevel> <tmax> <zWall> \
  [dropRadiusMin] [dropPersistence] [snapshotInterval] [drillAMR] \
  [drillMaxlevelStart] [drillMaxlevelFocus] [drillNcells] \
  [drillRegionMinX] [drillArmSteps] [drillArmTime] [drillCoarsenTime] \
  [drillRegionMaxX] [drillRegionRadius] [drillFireX] [drillTipRadius] \
  [drillRegionalOnly] [geometryMode] [wallClearance] [interfaceFloor]
```

## Solver Stacks

Three solver stacks are supported. A single compile-time flag selects between
them.

| Build | Stack |
|---|---|
| (no flag) | 1: `FILTERED` property averaging + `navier-stokes/double-projection.h` |
| `-DSINGLE_PROJECTION=1` | 1S: `FILTERED` property averaging, single projection |
| `-DUSE_CONSERVING=1` | 2: `navier-stokes/conserving.h`, no `FILTERED`, no double projection |

Stack 1 is the default. Double projection (Almgren et al. 2000) decouples the
face-velocity projection pressure from the centered-gradient pressure, so the
latter does not inherit the divergence history of adaptive refinement.

Stack 1S drops the second Poisson solve and keeps everything else. The
$R_r=30$ post-mortem showed that the *update* projection -- the second solve,
which carries the entire surface-tension impulse across the 1000:1 density
jump because `double-projection.h` zeroes `a` for the face solve -- is the one
that dies, while the face projection converges in two V-cycles with three times
the residual headroom. The AMR-divergence-history argument for the second solve
is also largely moot here: the unconditional interfacial refinement floor keeps
`contractDeficit` at zero for the whole run, so there is very little regrid
noise for the second pressure to be protected from. 1S is therefore the
cheaper and, empirically, the more robust member of the `FILTERED` family.

Stack 2 is the momentum-conserving VOF advection. It is mutually exclusive with
Stack 1 by construction: `conserving.h` sets `stokes = true` and moves the
momentum advection into the `vof` event, ahead of `advection_term`. The
double-projection update field $A_f$ is assembled in `advection_term` and would
therefore omit the advective update entirely. `FILTERED` is likewise dropped
because `conserving.h` supplies its own consistent face density.

The legacy `-DDOUBLE_PROJECTION=1` flag is retained only as a no-op alias for
Stack 1. It is rejected when combined with `-DUSE_CONSERVING` or
`-DSINGLE_PROJECTION`.

### Command-line Parameters

- `OhOut`: Solvent Ohnesorge number based on the small bubble radius,
  $Oh_s = \mu_{out}/\sqrt{\rho_{out} \sigma R_s}$
- `RhoIn`: Density ratio $\rho_{b}/\rho_{s}$ (typically $10^{-3}$ for air-water)
- `Rr`: Radius ratio $R_l/R_s$ (large to small bubble). The geometry files
  use the small bubble as the unit length, so the larger bubble radius is `Rr`.
- `MAXlevel`: Maximum adaptive mesh refinement level
- `tmax`: Maximum simulation time
- `zWall`: Wall position (distance from origin to left boundary)
- `dropRadiusMin`: Optional minimum equivalent detached-drop radius. Supplying
  zero selects the resolution-aware default $2\Delta_{min}$; omitting the
  argument keeps the legacy behaviour with detection disabled.
- `dropPersistence`: Number of consecutive checks above `dropRadiusMin`
  required before labelling a case as a drop (default: 3).
- `snapshotInterval`: Legacy full-snapshot interval (default: 0.01). Contour
  campaigns suppress the full time series, checkpoint every 0.5 time units,
  and publish lightweight facets every 0.05.
- `drillAMR`: Enable feature-driven regional refinement (`0` by default).
  The drill ramps from `drillMaxlevelStart`, arms on persistent target-region
  curvature demand, then gives full `MAXlevel` only to the end-pinchoff side.
- `geometryMode`: `finite` loads `InitialConditionRr-*.dat`; `halfspace`
  loads the Bursting-Bubble `Bo0.0000.dat` sphere-plane geometry and represents
  the true $R_r\to\infty$ limit.
- `wallClearance`: Optional physical distance from the bubble south pole to
  the left wall. A negative value preserves the legacy nominal `zWall`
  placement. Use `0.027` to match the finite-map `zWall=0.05` clearance.
- `interfaceFloor`: Enforce the unconditional `MAXlevel` refinement floor on
  every interfacial cell after each adaptation (`1` by default, `0` disables).

## Nondimensional Mapping Used in This Code

- The small-bubble radius is the reference length ($R_s=1$).
- The large-bubble radius is `R_l = Rr`.
- The bubble viscosity ratio is fixed as `MuRin = 1e-2`, so
  $Oh_b/Oh_s = 0.01$ with $Oh_b = \mu_b/\sqrt{\rho_s \sigma R_s}$.
- Confinement is controlled via `zWall`: smaller `zWall` places the small
  bubble closer to the wall, corresponding to smaller $\chi = d/R_s$.

## Parameter Sweeps

Typical sweeps vary `OhOut` and `Rr` (bulk) or `zWall` (confined) to build
droplet/no-droplet regime maps and to quantify the first injected droplet size.

## Author

Vatsal Sanjay
vatsalsanjay@gmail.com
Physics of Fluids

Last updated: Jan 2026
*/

/**
Reject every invalid stack combination before anything is included. These must
precede the `FILTERED` definition below, which is part of Stack 1 itself. */
#if defined(USE_CONSERVING) && defined(DOUBLE_PROJECTION)
#error "USE_CONSERVING and DOUBLE_PROJECTION are mutually exclusive: conserving.h sets stokes=true and advects momentum in the vof event, so double-projection's Af would omit the advective update. Pick one stack."
#endif
#if defined(USE_CONSERVING) && defined(FILTERED)
#error "USE_CONSERVING excludes FILTERED: conserving.h supplies its own consistent face density. Do not define FILTERED for the conserving stack."
#endif
#if defined(USE_CONSERVING) && defined(SINGLE_PROJECTION)
#error "USE_CONSERVING already runs a single projection: conserving.h never includes double-projection.h. Do not combine it with SINGLE_PROJECTION; pick one stack."
#endif
#if defined(SINGLE_PROJECTION) && defined(DOUBLE_PROJECTION)
#error "SINGLE_PROJECTION and DOUBLE_PROJECTION are mutually exclusive: the first removes the second Poisson solve that the second requests. Pick one stack."
#endif

#include "axi.h"
#include "navier-stokes/centered.h"
#ifndef USE_CONSERVING
/**
Stack 1: the second Poisson solve per timestep is unconditional here.
Stack 1S (`-DSINGLE_PROJECTION`): same `FILTERED` property averaging, but the
second solve is omitted, so `mgp` is the face-velocity projection and the
centred gradient `g` is built from the same pressure. */
#ifndef SINGLE_PROJECTION
#include "navier-stokes/double-projection.h"
#endif
#define FILTERED 1
#endif
#include "two-phase.h"
#ifdef USE_CONSERVING
/**
Stack 2: momentum-conserving VOF advection. Must follow `two-phase.h`. */
#include "navier-stokes/conserving.h"
#endif
#include "tension.h"
#ifdef USE_CONSERVING
#define SOLVER_STACK "conserving"
#define DUAL_PROJECTION 0
#elif defined(SINGLE_PROJECTION)
#define SOLVER_STACK "filtered+single-projection"
#define DUAL_PROJECTION 0
#else
#define SOLVER_STACK "filtered+double-projection"
#define DUAL_PROJECTION 1
#endif
#include "tag.h"
#include "adapt_wavelet_limited.h"
#include <float.h>
#include <string.h>

#if !_MPI
#include "distance.h"
#endif

/**
## Simulation Parameters

- `MAXlevel`: Maximum refinement level (set via command line)
- `tsnap`: Time interval for saving full snapshots (0.01)
- `tsnap2`: Time interval for logging diagnostics (0.0001)
*/

int MAXlevel, maxlevelLocal;

#define tsnap (1e-2)
#define tsnap2 (1e-4)

/**
## Error Tolerances

Adaptive mesh refinement is controlled by these error thresholds:

- `fErr`: VOF volume fraction tolerance (determines interface resolution)
- `VelErr`: Velocity field tolerance (captures flow gradients)
*/

#define fErr (1e-3)
#define VelErr (1e-2)

/**
## Boundary Conditions

Left boundary (solid substrate; the axisymmetric axis is `bottom`):
- `f[left]`: Liquid at the substrate (`f=0`, since `f` denotes gas)
- `u.n[left]`, `u.t[left]`: No penetration and no slip

Right boundary (far-field outlet, matching Bursting-Bubble):
- `u.n[right]`: Zero normal gradient
- `p[right]`: Reference pressure
*/

f[left] = dirichlet(0.0);
u.n[left] = dirichlet(0.0);
u.t[left] = dirichlet(0.0);
u.n[right] = neumann(0.0);
p[right] = dirichlet(0.0);

/**
## Global Variables

Physical parameters and output configuration:
*/

double tmax, MuRin, OhOut, RhoIn;
double Rr, zWall;
double Ldomain;
double wallClearance = -1.;
double shapeSouthPole = 0.;
char geometryMode[16] = "finite";
char initialConditionFile[80];
double dropRadiusMin = -1.;
double snapshotInterval = tsnap;
int dropPersistence = 3;
int dropConsecutive = 0;
bool dropDetected = false;
bool simulationInitialised = false;
double largestDetachedVolume = 0.;
double largestDetachedRadius = 0.;
double detachedAxialPosition = 0.;
double detachedAxialVelocity = 0.;
int drillAMR = 0;
int drillMaxlevelStart = -1;
int drillMaxlevelFocus = -1;
int drillArmSteps = 5;
int drillDemandSteps = 0;
int drillTipSteps = 0;
bool drillArmed = false;
bool drillFired = false;
int drillRegionalOnly = 0;
double drillNcells = 5.;
double drillRegionMinX = -2.1;
double drillRegionMaxX = 3.;
double drillRegionRadius = 1.5;
double drillArmTime = 0.;
double drillCoarsenTime = 0.;
double drillFireX = 0.25;
double drillTipRadius = 0.25;
char nameOut[80], dumpFile[80];

/**
## Numerics Guard State

- `interfaceFloor`: runtime switch for the unconditional `MAXlevel` refinement
  floor on interfacial cells (see the `adapt` event).
- `interfaceFloorDeficit`: AMR contract count. The number of interfacial cells
  left below `MAXlevel` after adaptation. It must be zero on every step and is
  reported per step in `projectionStats.dat`, not on stdout.
- `dtCap`/`dtCapMax`: persistent timestep ceiling maintained by the
  projection-failure backoff, and the ceiling it relaxes back to. The cap is
  published through the global `DT`, which `centered.h`'s `set_dtmax` event
  applies as `dtmax = DT` at the start of every step; `tension.h` and the
  advective CFL then reduce `dtmax` further as usual. It is inert while
  `projectionFailureLimit` is 1, because the run stops on the first failure
  before a second step can use the cap; the machinery is kept only for
  deliberate backoff experiments.
- `projectionFailures`: consecutive failed-projection steps. Reaching
  `projectionFailureLimit` stops the run cleanly instead of cascading to NaNs.
  The limit is **one**: the post-mortem showed that the step after a failed
  projection integrates a corrupted velocity field (`resb` jumped from
  $2.7\times10^{43}$ to $2.2\times10^{190}$ in a single step), so there is no
  useful information past the first failure and every later dump is garbage.
- `mgpFaceProjection`: statistics of the *face-velocity* projection. Under
  Stack 1, `double-projection.h` overwrites `mgp` with the second (update)
  projection in `end_timestep`, so the first solve's statistics are captured by
  the overloading `end_timestep` event below.
- `mgpUpdateProjection`: statistics of the *update* projection, i.e. `mgp` as
  it stands once every `end_timestep` overload has run. Under Stacks 1S and 2
  there is no second solve and this is the same solve as
  `mgpFaceProjection`; under Stack 1 it is the second (`Af`) solve, which is
  the one that actually diverges at $R_r=30$. Both are guarded.
- `cflTarget`: the advective CFL number the run is meant to use. Assigning
  `CFL` in `main()` is dead code -- `centered.h`'s `defaults` event (`i = 0`)
  overwrites it with 0.8 and `vof.h`'s `stability` event then clamps it to 0.5,
  both *after* `main()` returns. It is therefore (re)published every step by
  the `numericsControl` event below, which is a non-`last` event and so runs
  ahead of `set_dtmax`/`stability`.
- `ufc`: cell-centred mirror of the face velocity. `dump()` skips face vectors
  outright and skips `p`/`pf`/`dp` because `centered.h` marks the pressures
  `nodump`, which is why the plain `restart` file is not restartable. The
  `nodump` flags are cleared in `init` and this mirror is refreshed before
  every dump, so a checkpoint now carries the full state.
*/
int interfaceFloor = 1;
int interfaceFloorDeficit = 0;
int projectionFailures = 0;
int projectionFailureLimit = 1;
double dtCap = HUGE, dtCapMax = HUGE;
double cflTarget = 1e-1;
mgstats mgpFaceProjection, mgpUpdateProjection;
vector ufc[];

/**
Refresh the cell-centred face-velocity mirror. Called immediately before every
`dump()` so that checkpoints are self-contained. */
static void sync_uf_mirror (void)
{
  foreach()
    foreach_dimension()
      ufc.x[] = (uf.x[] + uf.x[1])/2.;
}

/**
Defined with the contour outputs below; forward-declared because the
projection guard publishes a final interface before stopping. */
static void write_contour_pulse (void);

/**
Return the smallest axial coordinate in an initial-shape polyline. Case
directories expose `DataFiles` as a symlink, while direct runs may already
have copied the shape into the working directory; support both layouts.
*/
static double initial_shape_south_pole (const char * filename)
{
  char path[160];
  FILE * fp = fopen (filename, "r");
  if (!fp) {
    snprintf (path, sizeof(path), "DataFiles/%s", filename);
    fp = fopen (path, "r");
  }
  if (!fp) {
    fprintf (ferr, "Cannot read initial-condition geometry '%s'\n", filename);
    exit (2);
  }

  double x, y, xmin = DBL_MAX;
  int points = 0;
  while (fscanf (fp, "%lf %lf", &x, &y) == 2) {
    xmin = min (xmin, x);
    points++;
  }
  fclose (fp);
  if (points < 2 || xmin == DBL_MAX) {
    fprintf (ferr, "Initial-condition geometry '%s' is empty or malformed\n",
             filename);
    exit (2);
  }
  return xmin;
}

/**
### write_classification_status()

Atomically publish the current contour-classification state. An `id` of `-1`
means the case is still running, `0` means no detached drop was detected before
normal termination, and `1` means a persistent detached drop crossed the
configured radius threshold.
*/
static void write_classification_status (int id, const char * reason)
{
  if (pid() != 0)
    return;

  FILE * fp = fopen ("classification.status.tmp", "w");
  if (!fp) {
    fprintf (ferr, "Could not write classification.status.tmp\n");
    return;
  }
  fprintf (fp, "id,reason,t,drop_volume,drop_radius,drop_axial_position,"
           "drop_axial_velocity,threshold,consecutive\n");
  fprintf (fp, "%d,%s,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%d\n", id,
           reason, t, largestDetachedVolume, largestDetachedRadius,
           detachedAxialPosition, detachedAxialVelocity, dropRadiusMin,
           dropConsecutive);
  fclose (fp);
  rename ("classification.status.tmp", "classification.status");
}

/**
### detached_tip_component()

Tag connected liquid regions, identify the exterior liquid as the largest
component, and return the leading detached component above the radius cutoff.
The leading component is the operational end-pinchoff candidate. Smaller
downstream fragments are ignored so generic Rayleigh--Plateau breakup does not
define the regime-map label.

The component accumulation is deliberately serial. `tag()` already performs
the connected-component labelling, while a serial accumulation avoids an
OpenMP race on the dynamically sized component-volume array. `dv()` supplies
the axisymmetric metric without the azimuthal $2\pi$ factor, which is restored
explicitly below.
*/
typedef struct {
  double volume;
  double radius;
  double axial_position;
  double axial_velocity;
} DetachedComponent;

static DetachedComponent detached_tip_component (void)
{
  DetachedComponent candidate = {0., 0., 0., 0.};
  scalar liquid[];
  foreach()
    liquid[] = (1. - f[]) > 1e-4;

  int n = tag (liquid);
  if (n < 2)
    return candidate;

  double volumes[n], axial_moments[n], velocity_moments[n];
  for (int j = 0; j < n; j++) {
    volumes[j] = 0.;
    axial_moments[j] = 0.;
    velocity_moments[j] = 0.;
  }

  foreach (serial) {
    if (liquid[] > 0.) {
      int label = (int) liquid[] - 1;
      double weight = 2.*pi*clamp (1. - f[], 0., 1.)*dv();
      volumes[label] += weight;
      axial_moments[label] += x*weight;
      velocity_moments[label] += u.x[]*weight;
    }
  }

  int exterior = 0;
  for (int j = 1; j < n; j++)
    if (volumes[j] > volumes[exterior])
      exterior = j;

  double leading_position = -1e30;
  for (int j = 0; j < n; j++) {
    if (j == exterior || volumes[j] <= 0.)
      continue;
    double radius = cbrt (3.*volumes[j]/(4.*pi));
    double axial_position = axial_moments[j]/volumes[j];
    if (radius >= dropRadiusMin && axial_position > leading_position) {
      leading_position = axial_position;
      candidate.volume = volumes[j];
      candidate.radius = radius;
      candidate.axial_position = axial_position;
      candidate.axial_velocity = velocity_moments[j]/volumes[j];
    }
  }
  return candidate;
}

/**
## Main Function

Initialize simulation parameters from command line and configure the domain.
*/

int main(int argc, char const *argv[]) {
  if (argc < 7){
    fprintf(ferr, "Lack of command line arguments. Check! Need %d more arguments\n", 7-argc);
    return 1;
  }

  /**
  Parse command-line arguments: */
  MuRin = 1e-2;

  OhOut = atof(argv[1]);
  RhoIn = atof(argv[2]);
  Rr = atof(argv[3]);
  MAXlevel = atoi(argv[4]);
  tmax = atof(argv[5]);
  zWall = atof(argv[6]);
  if (argc > 7)
    dropRadiusMin = atof(argv[7]);
  if (argc > 8)
    dropPersistence = atoi(argv[8]);
  if (argc > 9)
    snapshotInterval = atof(argv[9]);
  if (argc > 10)
    drillAMR = atoi(argv[10]);
  if (argc > 11)
    drillMaxlevelStart = atoi(argv[11]);
  if (argc > 12)
    drillMaxlevelFocus = atoi(argv[12]);
  if (argc > 13)
    drillNcells = atof(argv[13]);
  if (argc > 14)
    drillRegionMinX = atof(argv[14]);
  if (argc > 15)
    drillArmSteps = atoi(argv[15]);
  if (argc > 16)
    drillArmTime = atof(argv[16]);
  if (argc > 17)
    drillCoarsenTime = atof(argv[17]);
  if (argc > 18)
    drillRegionMaxX = atof(argv[18]);
  if (argc > 19)
    drillRegionRadius = atof(argv[19]);
  if (argc > 20)
    drillFireX = atof(argv[20]);
  if (argc > 21)
    drillTipRadius = atof(argv[21]);
  if (argc > 22)
    drillRegionalOnly = atoi(argv[22]);
  if (argc > 23)
    snprintf (geometryMode, sizeof(geometryMode), "%s", argv[23]);
  if (argc > 24)
    wallClearance = atof(argv[24]);
  if (argc > 25)
    interfaceFloor = atoi(argv[25]);

  bool halfspace = strcmp (geometryMode, "halfspace") == 0;
  bool finite = strcmp (geometryMode, "finite") == 0;
  if (!finite && !halfspace) {
    fprintf (ferr, "geometryMode must be 'finite' or 'halfspace', got '%s'\n",
             geometryMode);
    return 1;
  }
  if (halfspace)
    snprintf (initialConditionFile, sizeof(initialConditionFile),
              "Bo0.0000.dat");
  else
    snprintf (initialConditionFile, sizeof(initialConditionFile),
              "InitialConditionRr-%3.2f.dat", Rr);
  shapeSouthPole = initial_shape_south_pole (initialConditionFile);

  if (drillMaxlevelStart < 0)
    drillMaxlevelStart = max (MAXlevel - 2, 1);
  if (drillMaxlevelFocus < 0)
    drillMaxlevelFocus = max (MAXlevel - 1, drillMaxlevelStart);
  if (dropPersistence < 1 || snapshotInterval <= 0. ||
      drillMaxlevelStart < 1 || drillMaxlevelStart > drillMaxlevelFocus ||
      drillMaxlevelFocus > MAXlevel || drillNcells <= 0. ||
      drillArmSteps < 1 || drillArmTime < 0. || drillCoarsenTime < 0. ||
      drillCoarsenTime > drillArmTime || drillRegionMinX >= drillRegionMaxX ||
      drillRegionRadius <= 0. || drillFireX <= drillRegionMinX ||
      drillFireX >= drillRegionMaxX || drillTipRadius <= 0. ||
      drillTipRadius > drillRegionRadius || wallClearance == 0. ||
      (drillRegionalOnly != 0 && drillRegionalOnly != 1) ||
      (interfaceFloor != 0 && interfaceFloor != 1)) {
    fprintf (ferr, "Invalid contour controls: dropRadiusMin=%g, "
             "dropPersistence=%d, snapshotInterval=%g, drillAMR=%d, "
             "drillStart=%d, drillFocus=%d, drillNcells=%g, "
             "drillRegionMinX=%g, drillArmSteps=%d, drillArmTime=%g, "
             "drillCoarsenTime=%g, drillRegionMaxX=%g, "
             "drillRegionRadius=%g, drillFireX=%g, drillTipRadius=%g, "
             "drillRegionalOnly=%d, interfaceFloor=%d\n",
             dropRadiusMin, dropPersistence, snapshotInterval, drillAMR,
             drillMaxlevelStart, drillMaxlevelFocus, drillNcells,
             drillRegionMinX, drillArmSteps, drillArmTime, drillCoarsenTime,
             drillRegionMaxX, drillRegionRadius, drillFireX, drillTipRadius,
             drillRegionalOnly, interfaceFloor);
    return 1;
  }

  /**
  Place the wall first. For the half-space geometry, retain the canonical
  Bursting-Bubble right boundary at x=4 rather than forcing L0=16. The latter
  under-resolves the close-wall film by a factor 16/(zWall + 6). */
  double originX = -2.0 - zWall;
  if (wallClearance > 0.) {
    originX = shapeSouthPole - wallClearance;
    zWall = -2.0 - originX;
  }
  else
    wallClearance = shapeSouthPole - originX;

  /**
  The finite domain fits both bubbles plus a buffer. The half-space domain is
  the Bo=0 Bursting-Bubble domain, L0=min(zWall+6,16).

  The finite-geometry branch is deliberately uncapped: the previous
  `fmin(...,16.)` ceiling truncates the domain for any Rr large enough that
  `zWall+2+2*Rr+4 > 16` -- e.g. Rr=30 needs Ldomain=66.05 to contain a
  radius-30 bubble (`InitialConditionRr-30.00.dat` spans x up to ~60, y up to
  ~30), but the capped formula forces Ldomain=16, truncating the large bubble
  to roughly half its radius. The half-space branch keeps its cap, since it is
  unaffected by this and the 16-unit domain there is intentional (see the
  comment above). */
  Ldomain = halfspace ? fmin(zWall + 6.0, 16.) :
    zWall + 2. + 2.*Rr + 4.0;
  if (argc > 7 && dropRadiusMin == 0.)
    // exp2 avoids the int shift (1 << MAXlevel), which is UB/overflow for large
    // MAXlevel; the result is the same 2^MAXlevel for the levels used here.
    dropRadiusMin = 2.*Ldomain/exp2((double) MAXlevel);
  // Preserve the fully resolved initial neck. The drill may coarsen only after
  // drillCoarsenTime, mirroring the validated two-stage drill protocol.
  maxlevelLocal = MAXlevel;

  /**
  Configure domain and fluid properties: */
  L0=Ldomain;
  origin(originX, 0.0);
  init_grid (1 << (10));

  fprintf(ferr, "Level %d, Ldomain %g, tmax %3.2f, MuRin %3.2e, OhOut %3.2e, Rho21 %4.3f, Rr %f, geometry %s, initialShape %s, shapeSouthPole %g, wallClearance %g, zWall %g, dropRadiusMin %g, dropPersistence %d, snapshotInterval %g, drillAMR %d, drillStart %d, drillFocus %d, drillNcells %g, drillRegionMinX %g, drillArmSteps %d, drillArmTime %g, drillCoarsenTime %g, drillRegionMaxX %g, drillRegionRadius %g, drillFireX %g, drillTipRadius %g, drillRegionalOnly %d, interfaceFloor %d, solverStack %s\n", MAXlevel, Ldomain, tmax, MuRin, OhOut, RhoIn, Rr, geometryMode, initialConditionFile, shapeSouthPole, wallClearance, zWall, dropRadiusMin, dropPersistence, snapshotInterval, drillAMR, drillMaxlevelStart, drillMaxlevelFocus, drillNcells, drillRegionMinX, drillArmSteps, drillArmTime, drillCoarsenTime, drillRegionMaxX, drillRegionRadius, drillFireX, drillTipRadius, drillRegionalOnly, interfaceFloor, SOLVER_STACK);

  /**
  Set fluid properties:
  - Phase 1 (bubble interior): low density, low viscosity
  - Phase 2 (surrounding liquid): reference density and viscosity
  - Surface tension coefficient normalized to 1 */
  rho1 = RhoIn; mu1 = MuRin*OhOut;
  rho2 = 1e0; mu2 = OhOut;
  f.sigma = 1.0;

  char comm[80];
  sprintf (comm, "mkdir -p intermediate");
  system(comm);
  sprintf (dumpFile, "restart");

  TOLERANCE = 1e-4;
  /**
  Do *not* assign `CFL` here: `centered.h`'s `defaults` event resets it to 0.8
  and `vof.h`'s `stability` event clamps it to 0.5, both after `main()`
  returns, so the historical `CFL = 1e-1` in this spot never took effect (the
  measured effective CFL was exactly 0.5). Record the target instead; the
  `numericsControl` event publishes it every step. */
  cflTarget = 1e-1;

  /**
  Baseline for the projection-failure backoff. `DT` is Basilisk's global
  timestep ceiling (`HUGE` unless set), so the cap is inert until a projection
  actually fails and it relaxes all the way back here afterwards. */
  dtCapMax = DT;

  run();
}

/**
## Initialization Event

Load initial bubble shapes from pre-computed data files. The initial condition
files contain the interface coordinates for different radius ratios.

For MPI runs, the simulation must be restarted from a dump file created by
a prior serial/OpenMP run (due to `distance.h` incompatibility with MPI).
*/

event init(t = 0){
  /**
  This event is the head of the `init` chain, so it runs after every
  `defaults` event and before `centered.h`'s own `init`. That is the first
  point at which `centered.h`'s `CFL = 0.8` and `p.nodump = pf.nodump = true`
  can be undone.

  Clearing `nodump` makes `dump()` write the pressure. `double-projection.h`
  calls `scalar_clone (dp, p)` every step, which copies the attributes, so
  `dp` becomes dumpable too. `uf` is a face vector and can never be dumped;
  `centered.h`'s `init` rebuilds it from the restored `u` with
  `uf.x[] = fm.x[]*face_value (u.x, 0)`, and `ufc` carries the pre-restart
  face velocity for verification and post-processing. */
  CFL = cflTarget;
  p.nodump = pf.nodump = false;
#if _MPI
  if (!restore(file = dumpFile)) {
    fprintf(ferr, "Cannot restore from dump file '%s'. Run Stage 1 first.\n",
            dumpFile);
    return 1;
  }
#else
  if (!restore (file = dumpFile)){
    char comm[160];
    snprintf (comm, sizeof(comm), "cp DataFiles/%s .", initialConditionFile);
    if (system(comm) != 0) {
      fprintf(ferr, "Failed to copy initial-condition file '%s' from DataFiles/.\n",
              initialConditionFile);
      return 1;
    }

    FILE * fp = fopen(initialConditionFile,"rb");
    if (fp == NULL){
      fprintf(ferr, "There is no file named %s\n", initialConditionFile);
      return 1;
    }
    coord* InitialShape;
    InitialShape = input_xy(fp);
    fclose (fp);
    scalar d[];
    distance (d, InitialShape);
    /**
    Build the initial condition at `MAXlevel` everywhere, unconditionally.
    Seeding the distance-field adaptation at `drillMaxlevelStart` under the
    drill left the high-curvature neck coarse from step zero; the parasitic
    capillary currents it seeded are what later destroy the face-velocity
    projection. The startup cost of a uniform `MAXlevel` target is paid once,
    and it is the only IC that satisfies the AMR contract enforced in the
    `adapt` event below. */
    int initialLevel = MAXlevel;
    while (adapt_wavelet ((scalar *){f, d}, (double[]){1e-8, 1e-8},
                          initialLevel).nf);
    if (drillAMR && drillRegionalOnly) {
      /**
      Retained as a no-op accelerator. With `initialLevel = MAXlevel` the
      wavelet pass already targets the same ceiling; this only tops up any
      band cell whose distance-field wavelet error happened to fall below
      tolerance. It can raise resolution, never lower it. */
      refine (level < MAXlevel &&
              x >= drillRegionMinX && x <= drillRegionMaxX &&
              y <= drillRegionRadius && fabs(d[]) < 3.*Delta);
    }
    /**
    The distance function is defined at the center of each cell, we have
    to calculate the value of this function at each vertex. */
    vertex scalar phi[];
    foreach_vertex(){
      phi[] = (d[] + d[-1] + d[0,-1] + d[-1,-1])/4.;
    }
    /**
    We can now initialize the volume fraction of the domain. */
    fractions (phi, f);
    foreach(){
      u.x[] = 0.0;
      u.y[] = 0.0;
    }
    sync_uf_mirror();
    dump (file = dumpFile);
    static FILE * fp2;
    fp2 = fopen("InitialConditionStatus.dat","w");
    fprintf(fp2, "Initial condition is written to %s\n", dumpFile);
    fclose(fp2);
  }
#endif
  simulationInitialised = true;
}

/**
## Per-step Numerics Control

Republish the advective CFL target ahead of the timestep selection.

`CFL` is a global that three different places write to: `centered.h`'s
`defaults` event sets 0.8 at `i = 0`, `vof.h`'s `stability` event clamps
anything above 0.5 down to 0.5 on *every* step, and the user is expected to set
it in `main()` -- which is precisely the one place that cannot win, since both
of the others run afterwards. Setting it in `init` alone would survive, but
only by accident of `vof.h` clamping downwards rather than upwards.

This event is declared without the `last` attribute, so `qcc` registers it in
the leading block of the event array, ahead of `centered.h`'s
`set_dtmax`/`stability` and of `vof.h`'s `stability`. The assignment is
therefore in force for the `timestep (uf, dtmax)` call of the same step, on
every step including `i = 0`.
*/
event numericsControl (i++) {
  CFL = cflTarget;
}

/**
## Adaptive Mesh Refinement

Refine the mesh based on:
- Interface position (VOF field `f`)
- Velocity gradients (`u.x`, `u.y`)

Refinement ranges from `MAXlevel-6` (coarse, far from interface) to
`MAXlevel` (fine, near interface and in high-gradient regions).
*/

int drill_level_at (double x, double y, double z) {
  (void) z;
  bool in_wave_band = x >= drillRegionMinX && x <= drillRegionMaxX &&
    y <= drillRegionRadius;

  if (drillRegionalOnly)
    return in_wave_band ? MAXlevel : drillMaxlevelFocus;

  /*
  Before ARM, keep the whole domain at the validated pre-inception level. After
  ARM, release MAXlevel only in the capillary-wave, focusing and jet band while
  leaving the distant parent-bubble exterior at the cheaper focus level.
  */
  if (in_wave_band)
    return drillArmed ? MAXlevel : drillMaxlevelStart;
  return drillArmed ? drillMaxlevelFocus : drillMaxlevelStart;
}

event adapt(i++){
  if (drillAMR && t >= drillCoarsenTime)
    adapt_wavelet_limited ((scalar *){f, u.x, u.y},
      (double[]){fErr, VelErr, VelErr}, drill_level_at, MAXlevel-6);
  else
    adapt_wavelet ((scalar *){f, u.x, u.y},
      (double[]){fErr, VelErr, VelErr}, maxlevelLocal, MAXlevel-6);

  /**
  ### Interfacial refinement floor

  `adapt_wavelet*()` treats `MAXlevel` as a local *ceiling*, not a guarantee.
  A cut cell whose volume-fraction wavelet error falls below `fErr` is left
  coarse, and a nearly planar interface is exactly such a cell. That is how the
  $R_r=30$ neck and the thin wall film end up under-resolved between adapt
  steps, which seeds the parasitic capillary currents that break the face
  projection. Enforce the separate production contract instead: every
  interfacial cell sits at `MAXlevel`.

  `refine()` (`grid/tree-common.h`) already wraps its body in a
  `do { ... } while (refined)` loop and refines one level per pass, so a single
  invocation drives the whole band up to `MAXlevel`; no outer loop is needed. */
  const double interfaceEps = 1e-6;
  if (interfaceFloor)
    refine (level < MAXlevel &&
            f[] > interfaceEps && f[] < 1. - interfaceEps);

  /**
  ### AMR contract check

  Count the interfacial cells still below `MAXlevel`. With the floor engaged
  this must be zero on every step. It is reported per step in
  `projectionStats.dat` rather than on stdout, so a violation is auditable
  without flooding the run log. */
  int deficit = 0;
  foreach (reduction(+:deficit))
    if (f[] > interfaceEps && f[] < 1. - interfaceEps && level < MAXlevel)
      deficit++;
  interfaceFloorDeficit = deficit;
}

/**
Feature-driven arm/fire drill controller. Curvature demand arms the controller
after the bootstrap and releases the capillary-wave, focusing and jet band from
the pre-inception level to MAXlevel. The distant exterior remains capped at the
focus level. FIRE records persistent jet advance for diagnostics.
*/
event drillProbe(i++) {
  if (!drillAMR)
    return 0;
  if (drillRegionalOnly) {
    maxlevelLocal = MAXlevel;
    return 0;
  }
  scalar KAPPA[];
  curvature (f, KAPPA);
  double kmax = -1., tipx = -HUGE;
  foreach (reduction(max:kmax) reduction(max:tipx)) {
    if (x >= drillRegionMinX && x <= drillRegionMaxX &&
        y <= drillRegionRadius && f[] > 1e-6 && f[] < 1. - 1e-6 &&
        KAPPA[] != nodata && fabs(KAPPA[]) > kmax)
      kmax = fabs(KAPPA[]);
    if (x >= drillRegionMinX && x <= drillRegionMaxX &&
        y <= drillTipRadius && f[] > 1e-6 && f[] < 1. - 1e-6 && x > tipx)
      tipx = x;
  }

  int demanded = drillMaxlevelStart;
  if (kmax > 0.) {
    double need = drillNcells*Ldomain*kmax;
    while (demanded < MAXlevel && (double)(1 << demanded) < need)
      demanded++;
  }
  drillDemandSteps = t >= drillArmTime && demanded >= MAXlevel ?
    drillDemandSteps + 1 : 0;
  if (!drillArmed && drillDemandSteps >= drillArmSteps) {
    drillArmed = true;
    fprintf (ferr, "DRILL armed at t=%g after %d curvature-demand steps\n",
             t, drillDemandSteps);
  }
  if (t >= drillArmTime && tipx >= drillFireX)
    drillArmed = true; // self-heal a restart that already contains the jet
  drillTipSteps = drillArmed && tipx >= drillFireX ? drillTipSteps + 1 : 0;
  if (!drillFired && drillTipSteps >= drillArmSteps) {
    drillFired = true;
    fprintf (ferr, "DRILL fired at t=%g: tipx=%g persisted for %d steps; "
             "regional L%d enabled for x=[%g,%g], y<=%g\n", t, tipx,
             drillTipSteps, MAXlevel, drillRegionMinX, drillRegionMaxX,
             drillRegionRadius);
  }
  maxlevelLocal = MAXlevel;
  return 0;
}

/**
## Projection Diagnostics and Failure Guard

Under Stack 1 it is the *update* projection -- `double-projection.h`'s second
solve, `project (Af, p, ...)`, the one that produces the centred gradient `g`
-- that dies at $R_r=30$, not the face solve. Per-step dual instrumentation of
the crash showed the update solve running at a marginal fixed budget of about
six V-cycles at roughly $2.5\times$ reduction each, with a scaled residual
pinned at $9.7\times10^{-5}$ against `TOLERANCE` $=10^{-4}$ for hundreds of
consecutive steps, while the face solve converged in two cycles with three
times the headroom. It then failed to converge at all in one step, and the
resulting garbage `g` (of order $10^{15}$) destroyed the field. The reason the
update solve carries the load is structural: `double-projection.h` moves the
acceleration into $A_f$ and sets `a = zerof` for the face solve, so the entire
surface-tension impulse crosses the 1000:1 `alpha` jump in the second solve
alone.

Both solves are therefore recorded and both are guarded. Diagnostics go to
`projectionStats.dat`, one line per step. The main `log` column layout is
deliberately left untouched: `BayesianWorkflow/result_quality.py` parses `log`
positionally and reads `fields[3]` as the kinetic energy, so any numeric line
added there would be mistaken for a physical sample. New `projectionStats.dat`
columns are appended at the end for the same reason.
*/

/**
`double-projection.h` overwrites `mgp` with the second (update) projection in
its own `end_timestep`. Basilisk event inheritance runs the most recently
defined overload first, so this event -- defined last -- captures the
face-velocity projection statistics before they are lost. Under Stacks 1S and
2 there is no second projection and `mgp` is already the only solve, so the
same assignment is correct for all three stacks.
*/
event end_timestep (i++) {
  mgpFaceProjection = mgp;
}

event projectionStats (i++, last) {

  /**
  ### Update-projection capture

  This is a new event name, so `qcc` registers it after the whole
  `end_timestep` chain (including `double-projection.h`'s overload) and after
  `adapt`. `mgp` therefore holds the *last* projection performed this step:
  the update solve under Stack 1, and the single face solve under Stacks 1S
  and 2. */

  mgpUpdateProjection = mgp;

  /**
  ### Unweighted divergence monitor

  `project()` measures the residual of the *metric-weighted* discrete
  divergence. In axisymmetry `uf` already carries the face metric $f_m = y$, so
  the quantity the multigrid drives to zero is

  $$ D_w = \frac{1}{\Delta}\sum_d \left(u_{f,d}^{+} - u_{f,d}^{-}\right)
     = y\,\partial_x u_x + \partial_y (y\,u_y) $$

  which is $O(y)$ and therefore vanishes on the axis however bad the velocity
  field is there. Dividing by the cell metric $c_m = y$ recovers the honest
  divergence

  $$ \nabla\cdot\mathbf{u} = \partial_x u_x + \frac{1}{y}\partial_y (y\,u_y)
     = \frac{D_w}{c_m} $$

  and its maximum is what is reported. `cm[]` is at least $\Delta/2$ on the
  axis row so the division is safe; the floor is belt and braces. The field is
  sampled after adaptation, i.e. exactly as it is handed to the next step, so
  regrid-induced divergence is included by design. */

  double divMax = 0., fMin = 1., fMax = 0.;
  foreach (reduction(max:divMax) reduction(min:fMin) reduction(max:fMax)) {
    double divLocal = 0.;
    foreach_dimension()
      divLocal += uf.x[1] - uf.x[];
    divLocal = fabs(divLocal)/(Delta*max(cm[], 1e-30));
    if (divLocal > divMax)
      divMax = divLocal;
    if (f[] < fMin)
      fMin = f[];
    if (f[] > fMax)
      fMax = f[];
  }

  /**
  ### Measured advective CFL

  `timestep()` sets the advective ceiling to
  $\mathrm{CFL}\,\min_{\mathrm{faces}} \Delta/|u_f|$. Reporting that minimum
  and the achieved `dt` alongside the live value of `CFL` makes the effective
  CFL number an observable rather than an assumption -- which is exactly how
  the dead `CFL = 1e-1` assignment in `main()` went unnoticed. `dt/dtAdvLimit`
  equals `CFL` whenever the step is advection-limited and is smaller when the
  capillary or event-scheduling limits bind first. */

  double dtAdvLimit = HUGE;
  foreach_face (reduction(min:dtAdvLimit))
    if (uf.x[] != 0.) {
      double dtLocal = Delta*fm.x[]/fabs(uf.x[]);
      if (dtLocal < dtAdvLimit)
        dtAdvLimit = dtLocal;
    }

  /**
  ### Failure detection

  `project()` solves to the raw Poisson tolerance `TOLERANCE/dt^2`, so a fixed
  threshold on `resa` is dimensionally wrong and would reject healthy small
  timesteps. The scaled residual `resa*dt^2` is the dimensionless quantity the
  solver is actually asked to bring below `TOLERANCE`, so
  `resa*dt^2 > TOLERANCE` means precisely "the solve missed its bar", which can
  only happen by exhausting `NITERMAX`. Non-finite residuals count as failures
  too. */

  double faceScaledResidual = mgpFaceProjection.resa*sq(dt);
  double updScaledResidual = mgpUpdateProjection.resa*sq(dt);
  bool faceFailed = !isfinite (mgpFaceProjection.resa) ||
    !isfinite (faceScaledResidual) || faceScaledResidual > TOLERANCE;
  bool updFailed = !isfinite (mgpUpdateProjection.resa) ||
    !isfinite (updScaledResidual) || updScaledResidual > TOLERANCE;
  bool projectionFailed = faceFailed || updFailed;
  /**
  Report the worse of the two, so the guard status and the stdout line always
  name the solve that actually broke. */
  double scaledResidual = updScaledResidual > faceScaledResidual ||
    !isfinite (updScaledResidual) ? updScaledResidual : faceScaledResidual;
  const char * culprit = faceFailed && updFailed ? "both" :
    (updFailed ? "update" : (faceFailed ? "face" : "none"));

  /**
  ### First-failure stop

  The historical guard did three things wrong and all three are removed here.
  It watched only the face projection, so it never saw the update solve die.
  Its `DT` backoff took effect a step *late*: `centered.h`'s
  `set_dtmax`/`stability` are declared before this event, and on the failure
  step `timestep()` skips non-finite faces, so `dtnext()` returned the whole
  remaining logging interval and the measured `dt` doubled or tripled at the
  moment the cap was supposed to bite. And it continued integrating a corrupted
  field, whose `resb` reached $2.7\times10^{43}$ on the next step.

  The policy is now: the first failure of either projection stops the run. The
  `dtCap`/relax machinery is retained behind `projectionFailureLimit > 1` for
  deliberate backoff experiments, but it is no longer the production path. */

  if (projectionFailed) {
    projectionFailures++;
    if (projectionFailureLimit > 1) {
      dtCap = 0.5*(dtCap < dtCapMax ? dtCap : dt);
      DT = dtCap;
    }
    /**
    The failure dump is named with both the step number and a six-decimal
    time. The old `%5.4f` format collided whenever two failures fell inside
    the same $10^{-4}$ window, and the informative first dump was overwritten
    by the all-NaN second one. */
    char failName[96];
    snprintf (failName, sizeof(failName),
              "dump-projection-failure-i%06d-t%.6f", i, t);
    sync_uf_mirror();
    dump (file = failName);
    if (pid() == 0)
      fprintf (ferr, "PROJECTION-FAILURE %d/%d at i=%d t=%g dt=%g culprit=%s: "
               "face.i=%d/%d resb=%g resa=%g scaled=%g | upd.i=%d/%d resb=%g "
               "resa=%g scaled=%g | tolerance=%g, maxDivU=%g. Snapshot '%s'.\n",
               projectionFailures, projectionFailureLimit, i, t, dt, culprit,
               mgpFaceProjection.i, NITERMAX, mgpFaceProjection.resb,
               mgpFaceProjection.resa, faceScaledResidual,
               mgpUpdateProjection.i, NITERMAX, mgpUpdateProjection.resb,
               mgpUpdateProjection.resa, updScaledResidual, TOLERANCE, divMax,
               failName);
  }
  else {
    projectionFailures = 0;
    if (dtCap < dtCapMax) {
      dtCap = min (1.05*dtCap, dtCapMax);
      DT = dtCap;
    }
  }

  static FILE * fps;
  if (pid() == 0) {
    if (i == 0) {
      fps = fopen ("projectionStats.dat", "w");
      if (fps)
        fprintf (fps, "# solverStack %s, MAXlevel %d, TOLERANCE %g, "
                 "NITERMAX %d, interfaceFloor %d, projectionFailureLimit %d, "
                 "cflTarget %g, dualProjection %d\n"
                 "# i t dt mgp_i mgp_resb mgp_resa mgp_scaled_resa mgp_nrelax "
                 "maxDivU fMin fMax contractDeficit dtCap projectionFailed "
                 "mgu_i mgu_resa "
                 "upd_i upd_resb upd_resa upd_scaled_resa upd_nrelax "
                 "CFL dtAdvLimit\n",
                 SOLVER_STACK, MAXlevel, TOLERANCE, NITERMAX, interfaceFloor,
                 projectionFailureLimit, cflTarget, DUAL_PROJECTION);
    }
    else
      fps = fopen ("projectionStats.dat", "a");
    if (fps) {
      /**
      `mgu` (the viscous solve) is reported but not guarded. The conserving
      stack exhausts `NITERMAX` there long before the pressure projection is in
      trouble, so a projection-only diagnostic would report a healthy step
      while the run is in fact grinding. Two extra columns make that visible
      without changing any behaviour.

      The seven trailing columns are appended, never inserted: existing
      readers index `projectionStats.dat` positionally. `mgp_*` remains the
      face projection and `upd_*` is the update projection; under Stacks 1S
      and 2 the two blocks are identical by construction. */
      fprintf (fps, "%d %.8g %.8g %d %.8g %.8g %.8g %d %.8g %.8g %.8g %d "
               "%.8g %d %d %.8g %d %.8g %.8g %.8g %d %.8g %.8g\n",
               i, t, dt, mgpFaceProjection.i,
               mgpFaceProjection.resb, mgpFaceProjection.resa,
               faceScaledResidual,
               mgpFaceProjection.nrelax, divMax, fMin, fMax,
               interfaceFloorDeficit, dtCap, (int) projectionFailed,
               mgu.i, mgu.resa,
               mgpUpdateProjection.i, mgpUpdateProjection.resb,
               mgpUpdateProjection.resa, updScaledResidual,
               mgpUpdateProjection.nrelax, CFL, dtAdvLimit);
      fclose (fps);
    }
  }

  /**
  ### Clean stop

  Stop with a distinct status rather than letting the residual cascade into
  NaNs, and leave every dump in place. `classification.status` is finalised
  here rather than in `event end`, because returning a non-zero value from an
  event makes `events()` return immediately and the `end` event never runs --
  which is how the original crash left a stale `running` classification behind.
  No physical verdict is fabricated. */

  if (projectionFailures >= projectionFailureLimit) {
    if (pid() == 0) {
      fprintf (ferr, "PROJECTION-GUARD: %d consecutive failed projections "
               "(culprit=%s) at i=%d t=%g dt=%g (worst scaled resa=%g > "
               "tolerance=%g). Stopping cleanly; dumps retained; no physical "
               "verdict.\n",
               projectionFailures, culprit, i, t, dt, scaledResidual,
               TOLERANCE);
      FILE * gfp = fopen ("projection_guard.status", "w");
      if (gfp) {
        fprintf (gfp, "status=projection_failure_first_failure_stop\n"
                 "solverStack=%s\nculprit=%s\ni=%d\nt=%.8g\ndt=%.8g\n"
                 "consecutive=%d\nlimit=%d\n"
                 "face_i=%d\nface_resb=%.8g\nface_resa=%.8g\n"
                 "face_scaled_resa=%.8g\n"
                 "upd_i=%d\nupd_resb=%.8g\nupd_resa=%.8g\n"
                 "upd_scaled_resa=%.8g\n"
                 "tolerance=%.8g\nnitermax=%d\nmaxDivU=%.8g\n"
                 "contractDeficit=%d\ndtCap=%.8g\nCFL=%.8g\ndtAdvLimit=%.8g\n",
                 SOLVER_STACK, culprit, i, t, dt,
                 projectionFailures, projectionFailureLimit,
                 mgpFaceProjection.i, mgpFaceProjection.resb,
                 mgpFaceProjection.resa, faceScaledResidual,
                 mgpUpdateProjection.i, mgpUpdateProjection.resb,
                 mgpUpdateProjection.resa, updScaledResidual,
                 TOLERANCE, NITERMAX, divMax, interfaceFloorDeficit,
                 dtCap, CFL, dtAdvLimit);
        fclose (gfp);
      }
    }
    write_contour_pulse();
    write_classification_status (-1, "abnormal_termination_unclassified");
    return 1;
  }
  return 0;
}

/**
## Output Files

Save simulation snapshots at regular intervals:
- `restart`: Current state for checkpoint/restart
- `intermediate/snapshot-<time>`: Time series for post-processing
*/

event writingFiles (t = 0; t += snapshotInterval; t <= tmax + snapshotInterval) {
  /**
  The restart checkpoint is unconditional. Contour mode previously returned
  before any `dump()`, so a crash between the 0.5-interval `contourCheckpoint`
  dumps left nothing to restart from -- including at $t=0$, which made every
  early-death case unrecoverable and unreproducible.

  The snapshot *time series* is now also unconditional. Suppressing it in
  contour mode was a filesystem-load decision taken when `snapshotInterval`
  was the dense legacy $10^{-2}$ series; the campaigns actually pass a sparse
  interval, and the post-mortem had to be run on a bespoke instrumented build
  purely because the production contour runs kept no intermediate state at all
  between the 0.5-unit checkpoints. Control the load with `snapshotInterval`,
  not by discarding the series. */
  sync_uf_mirror();
  dump (file = dumpFile);
  sprintf (nameOut, "intermediate/snapshot-%5.4f", t);
  dump(file=nameOut);
}

event contourCheckpoint (t = 0.5; t += 0.5; t < tmax) {
  if (dropRadiusMin >= 0.) {
    sync_uf_mirror();
    dump (file = dumpFile);
  }
  return 0;
}

/**
## Contour-Campaign End-Pinchoff Detection

When `dropRadiusMin >= 0`, inspect connected liquid components every 0.02 time
units. The leading detached component must exceed the equivalent spherical
radius threshold on `dropPersistence` consecutive checks. It is a drop only
when its volume-weighted axial velocity is positive at confirmed pinch-off.
A zero or negative velocity is a resolved no-drop event; later
Rayleigh--Plateau fragments do not define injection.
*/
static void write_contour_pulse (void)
{
  if (dropRadiusMin < 0. || pid() != 0)
    return;

  FILE * fp = fopen ("interface-latest.dat.tmp", "w");
  if (fp) {
    output_facets (f, fp);
    fclose (fp);
    rename ("interface-latest.dat.tmp", "interface-latest.dat");
  }
  FILE * meta = fopen ("interface-latest.meta.tmp", "w");
  if (meta) {
    fprintf (meta, "t=%.8g\nRr=%.8g\nOh=%.8g\nzWall=%.8g\n"
             "geometry=%s\nwallClearance=%.8g\n",
             t, Rr, OhOut, zWall, geometryMode, wallClearance);
    fclose (meta);
    rename ("interface-latest.meta.tmp", "interface-latest.meta");
  }
}

event detectDetachedDrop (t = 0.05; t += 2.*tsnap; t <= tmax + tsnap) {
  if (dropRadiusMin < 0.)
    return 0;

  DetachedComponent candidate = detached_tip_component();
  largestDetachedVolume = candidate.volume;
  largestDetachedRadius = candidate.radius;
  detachedAxialPosition = candidate.axial_position;
  detachedAxialVelocity = candidate.axial_velocity;

  if (largestDetachedRadius >= dropRadiusMin)
    dropConsecutive++;
  else
    dropConsecutive = 0;

  write_classification_status (-1, "running");
  if (dropConsecutive >= dropPersistence) {
    dropDetected = true;
    sync_uf_mirror();
    dump (file = dumpFile);
    write_contour_pulse();
    int injected = detachedAxialVelocity > 0.;
    write_classification_status (injected,
      injected ? "ejected_end_pinchoff_drop" : "end_pinchoff_not_ejected");
    fprintf (ferr, "End pinch-off classified at t=%g: id=%d, volume=%g, "
             "radius=%g, x=%g, ux=%g (threshold=%g, consecutive=%d). "
             "Stopping.\n", t, injected, largestDetachedVolume,
             largestDetachedRadius, detachedAxialPosition,
             detachedAxialVelocity, dropRadiusMin, dropConsecutive);
    return 1;
  }
  return 0;
}

/**
## Lightweight Pulse Output

Contour runs atomically refresh `interface-latest.dat` every 0.05 time units.
The monitor renders this small PLIC-facet file locally, avoiding snapshot
post-processing or compilation on the Hamilton login node.
*/
event contourPulse (t = 0.; t += 5.*tsnap; t <= tmax + tsnap) {
  write_contour_pulse();
  return 0;
}

/**
Stop contour campaigns only at their validated observation horizon. The legacy
bubble-phase kinetic-energy criterion is not a safe no-drop classifier: a quiet
bubble does not prove that the surrounding liquid jet cannot later pinch off.
*/
event stopAtObservationHorizon (t = tmax) {
  if (dropRadiusMin >= 0.)
    return 1;
  return 0;
}

event end (t = end) {
  fprintf(ferr, "Level %d, Ldomain %g, tmax %3.2f, MuRin %3.2e, Oh2 %3.2e, Rho21 %4.3f, Rr %f\n", MAXlevel, Ldomain, tmax, MuRin, OhOut, RhoIn, Rr);
  if (dropRadiusMin >= 0. && !dropDetected) {
    write_contour_pulse();
    if (simulationInitialised && t >= tmax - 1e-8)
      write_classification_status (0, "observation_horizon_without_drop");
    else
      write_classification_status (-1, "abnormal_termination_unclassified");
  }
}

/**
## Logging Event

Track simulation diagnostics at fine time intervals. Computed quantities:

- `ke`: Total kinetic energy in the bubble phase
  $$KE = \int_V \frac{1}{2} f \rho |\mathbf{u}|^2 \, dV$$

- `xCOM`: Axial position of center of mass
  $$x_{COM} = \frac{\int_V f \cdot x \, dV}{\int_V f \, dV}$$

- `Vcm`: Axial velocity of center of mass
  $$V_{cm} = \frac{\int_V f \cdot u_x \, dV}{\int_V f \, dV}$$

Legacy simulations stop when bubble-phase kinetic energy drops below
$10^{-5}$. Contour campaigns ignore that early-exit criterion and run to the
validated observation horizon unless a persistent detached drop is detected.
*/

event logWriting (t = 0; t += tsnap2; t <= tmax+tsnap) {

  double ke = 0., wt = 0., xCOM = 0., Vcm = 0.;

  foreach (reduction(+:ke) reduction(+:wt) reduction(+:xCOM) reduction(+:Vcm)){
    ke += 2*pi*y*(0.5*clamp(f[], 0., 1.)*(sq(u.x[]) + sq(u.y[])))*sq(Delta);
    xCOM += 2*pi*y*x*clamp(f[], 0.0, 1.0)*sq(Delta);
    Vcm += 2*pi*y*u.x[]*clamp(f[], 0.0, 1.0)*sq(Delta);
    wt += 2*pi*y*clamp(f[], 0.0, 1.0)*sq(Delta);
  }
  xCOM /= wt;

  static FILE * fp;

  if (pid() == 0) {
    if (i == 0) {
      fprintf (ferr, "i dt t ke Xc Vcm maxlevel drillArmed drillFired\n");
      fp = fopen ("log", "w");
      fprintf(fp, "Level %d, Ldomain %g, tmax %3.2f, MuRin %3.2e, OhOut %3.2e, Rho21 %4.3f, Rr %f\n", MAXlevel, Ldomain, tmax, MuRin, OhOut, RhoIn, Rr);
      fprintf (fp, "i dt t ke Xc Vcm maxlevel drillArmed drillFired\n");
      fprintf (fp, "%d %g %g %g %g %g %d %d %d\n", i, dt, t, ke,
               xCOM, Vcm/wt, maxlevelLocal, drillArmed, drillFired);
      fclose(fp);
    } else {
      fp = fopen ("log", "a");
      fprintf (fp, "%d %g %g %g %g %g %d %d %d\n", i, dt, t, ke,
               xCOM, Vcm/wt, maxlevelLocal, drillArmed, drillFired);
      fclose(fp);
    }
    fprintf (ferr, "%d %g %g %g %g %g %d %d %d\n", i, dt, t, ke,
             xCOM, Vcm/wt, maxlevelLocal, drillArmed, drillFired);
  }

  /**
  `assert (ke > -1e-10)` was the process's actual cause of death in the
  $R_r=30$ crashes: once the projection produced a NaN velocity field, `ke`
  became NaN, the comparison was false, and `abort()` fired before anything
  could be written. That left no `projection_guard.status` and a stale
  `running` `classification.status` on disk, so the case looked hung rather
  than failed. With the first-failure projection guard above this path should
  now be unreachable, but if it is ever reached it must terminate the same way
  the guard does: an auditable status file and a clean return. */
  if (!isfinite (ke) || ke <= -1e-10) {
    if (pid() == 0) {
      fprintf (ferr, "NUMERICS-GUARD: non-physical kinetic energy ke=%g at "
               "i=%d t=%g dt=%g. Stopping cleanly.\n", ke, i, t, dt);
      FILE * gfp = fopen ("numerics_guard.status", "w");
      if (gfp) {
        fprintf (gfp, "status=nonphysical_kinetic_energy\nsolverStack=%s\n"
                 "i=%d\nt=%.8g\ndt=%.8g\nke=%.8g\n", SOLVER_STACK, i, t, dt,
                 ke);
        fclose (gfp);
      }
    }
    write_classification_status (-1, "abnormal_termination_unclassified");
    return 1;
  }
  if (ke < 1e-5 && i > 1000 && dropRadiusMin < 0.){
    if (pid() == 0){
      fprintf(ferr, "kinetic energy too small now! Stopping!\n");
      fp = fopen ("log", "a");
      fprintf(fp, "kinetic energy too small now! Stopping!\n");
      fclose(fp);
    }
    return 1;
  }
}
