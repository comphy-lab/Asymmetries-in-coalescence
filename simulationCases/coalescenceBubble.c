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
  [drillRegionalOnly] [geometryMode] [wallClearance]
```

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

#include "axi.h"
#include "navier-stokes/centered.h"
#define FILTERED 1
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
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
      (drillRegionalOnly != 0 && drillRegionalOnly != 1)) {
    fprintf (ferr, "Invalid contour controls: dropRadiusMin=%g, "
             "dropPersistence=%d, snapshotInterval=%g, drillAMR=%d, "
             "drillStart=%d, drillFocus=%d, drillNcells=%g, "
             "drillRegionMinX=%g, drillArmSteps=%d, drillArmTime=%g, "
             "drillCoarsenTime=%g, drillRegionMaxX=%g, "
             "drillRegionRadius=%g, drillFireX=%g, drillTipRadius=%g, "
             "drillRegionalOnly=%d\n",
             dropRadiusMin, dropPersistence, snapshotInterval, drillAMR,
             drillMaxlevelStart, drillMaxlevelFocus, drillNcells,
             drillRegionMinX, drillArmSteps, drillArmTime, drillCoarsenTime,
             drillRegionMaxX, drillRegionRadius, drillFireX, drillTipRadius,
             drillRegionalOnly);
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
  the Bo=0 Bursting-Bubble domain, L0=min(zWall+6,16). */
  Ldomain = halfspace ? fmin(zWall + 6.0, 16.) :
    fmin(zWall + 2. + 2.*Rr + 4.0, 16.);
  if (argc > 7 && dropRadiusMin == 0.)
    dropRadiusMin = 2.*Ldomain/(1 << MAXlevel);
  // Preserve the fully resolved initial neck. The drill may coarsen only after
  // drillCoarsenTime, mirroring the validated two-stage drill protocol.
  maxlevelLocal = MAXlevel;

  /**
  Configure domain and fluid properties: */
  L0=Ldomain;
  origin(originX, 0.0);
  init_grid (1 << (6));

  fprintf(ferr, "Level %d, Ldomain %g, tmax %3.2f, MuRin %3.2e, OhOut %3.2e, Rho21 %4.3f, Rr %f, geometry %s, initialShape %s, shapeSouthPole %g, wallClearance %g, zWall %g, dropRadiusMin %g, dropPersistence %d, snapshotInterval %g, drillAMR %d, drillStart %d, drillFocus %d, drillNcells %g, drillRegionMinX %g, drillArmSteps %d, drillArmTime %g, drillCoarsenTime %g, drillRegionMaxX %g, drillRegionRadius %g, drillFireX %g, drillTipRadius %g, drillRegionalOnly %d\n", MAXlevel, Ldomain, tmax, MuRin, OhOut, RhoIn, Rr, geometryMode, initialConditionFile, shapeSouthPole, wallClearance, zWall, dropRadiusMin, dropPersistence, snapshotInterval, drillAMR, drillMaxlevelStart, drillMaxlevelFocus, drillNcells, drillRegionMinX, drillArmSteps, drillArmTime, drillCoarsenTime, drillRegionMaxX, drillRegionRadius, drillFireX, drillTipRadius, drillRegionalOnly);

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
  CFL = 1e-1;


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
    int initialLevel = drillAMR ? drillMaxlevelStart : MAXlevel;
    while (adapt_wavelet ((scalar *){f, d}, (double[]){1e-8, 1e-8},
                          initialLevel).nf);
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
## Output Files

Save simulation snapshots at regular intervals:
- `restart`: Current state for checkpoint/restart
- `intermediate/snapshot-<time>`: Time series for post-processing
*/

event writingFiles (t = 0; t += snapshotInterval; t <= tmax + snapshotInterval) {
  // Contour campaigns use lightweight facets and coarse restart checkpoints.
  // Writing two full dumps every 0.05 time units for 16 simultaneous cases
  // overwhelms the shared filesystem without improving classification.
  if (dropRadiusMin >= 0.)
    return 0;
  dump (file = dumpFile);
  sprintf (nameOut, "intermediate/snapshot-%5.4f", t);
  dump(file=nameOut);
}

event contourCheckpoint (t = 0.5; t += 0.5; t < tmax) {
  if (dropRadiusMin >= 0.)
    dump (file = dumpFile);
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

  assert(ke > -1e-10);
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
