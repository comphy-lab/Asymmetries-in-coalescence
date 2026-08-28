/**
# Bursting-Bubble Drop Map

Instrumentation wrapper around the ref-locked bursting-bubble production
solver, used to measure the first-drop radius $r_d$ and the jet velocity
$v(t)$ as functions of the Ohnesorge number at a single, matched grid
spacing.

## Scientific Context

The first ejected drop from a bursting bubble does not shrink monotonically
with viscosity. Increasing $Oh$ from small values sharpens the capillary-wave
focusing and thins the jet, so the drop shrinks; this is the branch described
by the Blanco-Rodriguez & Gordillo drop-size law. At the Ohnesorge number
where the viscous cutoff catches up with the inviscid jet radius, the trend
turns: viscous regularisation thickens the jet again and the drop grows.
The turning point should therefore coincide with the minimum jet radius and
the maximum jet speed.

Resolving that turning point requires a drop-size curve taken at one grid
spacing with a resolution-independent detection radius, which is what this
driver exists to produce.

## Why a Separate Driver

Three properties are needed that a first-drop classifier deliberately does
not have.

1. **Every component, not the largest.** The first topological pinch is
   routinely only a handful of cells across and is frequently not the drop
   that should be classified. A single frame can hold the main body, a
   five-cell fragment and the physically meaningful drop at once. Selecting
   "the first detached component" silently picks the fragment, so this driver
   records every liquid component together with its cell count and defers
   selection to post-processing.
2. **A jet time series.** Tip position, tip-cap velocity, cap radius, cap
   volume and kinetic energy are written alongside the components, so the
   drop-size turning point can be tested against the jet-speed maximum
   independently.
3. **No early termination.** Each case runs to its declared horizon. Stopping
   at the first confirmed event hides the rest of the drop sequence, and later
   Rayleigh--Plateau fragments must be excluded by rule rather than by the run
   having already stopped.

The pinned detector is still evaluated every check and written to
`classification.status`, but purely as a diagnostic: it never ends a run.

## Dependency

This is a wrapper, not a standalone solver. It includes
`burstingBubble-drillResolution.c`, which is owned by the
`singular-bursting-bubbles` project and will not compile on its own.

That solver in turn includes `params.h` and `adapt_wavelet_limited.h`. Both
must come from the same bursting-bubble tree: this repository's `src-local`
has no `params.h` and carries a *different* `adapt_wavelet_limited.h`, so a
build that resolves `-I../src-local` to this checkout would either fail or
silently use different regional-AMR behaviour. The launcher hash-guards all
three files for that reason. Stage the bundle; do not run from a bare
checkout of this repository.

## Outputs

| file | contents |
|---|---|
| `components.log` | `t`, component index, `is_main`, cell count, volume, equivalent radius, axial centroid, mean axial velocity, axial extent |
| `jet.log` | `t`, tip position, tip-cap velocity, cap radius, cap volume, kinetic energy, component count, and the largest detached component |
| `classification.status` | pinned detector verdict, diagnostic only |

Both logs are appended at a fixed cadence (`TLOG`), so a restarted case
extends rather than truncates its history.

## Usage

~~~bash
# from the staged bundle, where src-local is the bursting-bubble one
qcc -O2 -Wall -disable-dimensions -fopenmp -I../src-local \
  burstingBubbleDropMap.c -o bursting-bubble-dropmap -lm
./bursting-bubble-dropmap case.params
~~~

No `-DFILTERED`. The comparator cases 6211/6212 were built without it, and
matching them is the point of the ladder.

Run with `drillRelaxLevel=-1` so detached components and the rounded tip stay
at production resolution; relaxing the mesh after pinch changes the measured
drop radius.
*/

#include "burstingBubble-drillResolution.c"

#define VISIBLE_RADIUS 0.021005127
#define VISIBLE_PERSISTENCE 3
#define TLOG 0.005          /* component + jet cadence */
#define RTIP 0.50           /* near-axis band defining the jet tip */
#define CAPDEPTH 0.05       /* axial depth of the tip cap used for v(t) */

static int visible_consecutive = 0;
static FILE * fp_components = NULL;
static FILE * fp_jet = NULL;

static void dropmap_open (void)
{
  if (fp_components)
    return;
  fp_components = fopen ("components.log", "a");
  if (fp_components && ftell (fp_components) == 0)
    fprintf (fp_components,
             "t,component,is_main,cells,volume,r_eq,x_mean,u_x_mean,x_min,x_max\n");
  fp_jet = fopen ("jet.log", "a");
  if (fp_jet && ftell (fp_jet) == 0)
    fprintf (fp_jet,
             "t,z_tip,u_tip_cap,r_cap_max,vol_cap,ke,n_components,"
             "r_largest_detached,u_largest_detached,x_largest_detached\n");
}

event dropmap (t = 0.; t += TLOG; t <= tmax + 1e-9)
{
  dropmap_open();

  /* ---- every liquid component ---------------------------------------- */
  scalar component[];
  foreach()
    component[] = f[] > 1e-4;
  int n = tag (component);

  double * volume = NULL, * momentum = NULL, * xmoment = NULL;
  double * cells = NULL, * xmin = NULL, * xmax = NULL;
  int main_component = 0;
  double r_big = 0., u_big = 0., x_big = -1000.;

  if (n > 0) {
    volume   = calloc (n, sizeof(double));
    momentum = calloc (n, sizeof(double));
    xmoment  = calloc (n, sizeof(double));
    cells    = calloc (n, sizeof(double));
    xmin     = calloc (n, sizeof(double));
    xmax     = calloc (n, sizeof(double));
    for (int j = 0; j < n; j++) {
      xmin[j] =  HUGE;
      xmax[j] = -HUGE;
    }
    foreach (serial)
      if (component[] > 0.) {
        int j = (int) component[] - 1;
        double dV = 2.*pi*y*sq(Delta)*clamp (f[], 0., 1.);
        volume[j]   += dV;
        momentum[j] += dV*u.x[];
        xmoment[j]  += dV*x;
        cells[j]    += 1.;
        if (x < xmin[j]) xmin[j] = x;
        if (x > xmax[j]) xmax[j] = x;
      }
    for (int j = 1; j < n; j++)
      if (volume[j] > volume[main_component])
        main_component = j;

    if (pid() == 0 && fp_components)
      for (int j = 0; j < n; j++)
        if (volume[j] > 0.)
          fprintf (fp_components,
                   "%.10g,%d,%d,%g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g\n",
                   t, j + 1, j == main_component ? 1 : 0, cells[j], volume[j],
                   cbrt (3.*volume[j]/(4.*pi)), xmoment[j]/volume[j],
                   momentum[j]/volume[j], xmin[j], xmax[j]);

    for (int j = 0; j < n; j++)
      if (j != main_component && volume[j] > 0.) {
        double radius = cbrt (3.*volume[j]/(4.*pi));
        if (radius > r_big) {
          r_big = radius;
          u_big = momentum[j]/volume[j];
          x_big = xmoment[j]/volume[j];
        }
      }
  }

  /* ---- jet tip and its cap velocity ----------------------------------- */
  double ztip = -HUGE;
  foreach (reduction(max:ztip))
    if (f[] > 0.5 && y < RTIP && x > ztip)
      ztip = x;

  double vol_cap = 0., mom_cap = 0., r_cap = 0.;
  if (ztip > -HUGE) {
    foreach (reduction(+:vol_cap) reduction(+:mom_cap) reduction(max:r_cap))
      if (f[] > 0.5 && y < RTIP && x > ztip - CAPDEPTH) {
        double dV = 2.*pi*y*sq(Delta)*clamp (f[], 0., 1.);
        vol_cap += dV;
        mom_cap += dV*u.x[];
        if (y > r_cap)
          r_cap = y;
      }
  }

  double kinetic = 0.;
  foreach (reduction(+:kinetic))
    kinetic += 2.*pi*y*sq(Delta)*clamp (f[], 0., 1.)*
               (sq(u.x[]) + sq(u.y[]))/2.;

  if (pid() == 0 && fp_jet) {
    fprintf (fp_jet, "%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%d,%.10g,%.10g,%.10g\n",
             t, ztip > -HUGE ? ztip : nodata,
             vol_cap > 0. ? mom_cap/vol_cap : nodata,
             r_cap, vol_cap, kinetic, n, r_big, u_big, x_big);
    fflush (fp_jet);
    if (fp_components)
      fflush (fp_components);
  }

  /* ---- pinned classifier, recorded as a diagnostic only --------------- */
  if (r_big >= VISIBLE_RADIUS)
    visible_consecutive++;
  else
    visible_consecutive = 0;

  if (pid() == 0) {
    FILE * fp = fopen ("classification.status.tmp", "w");
    if (fp) {
      int decided = visible_consecutive >= VISIBLE_PERSISTENCE;
      fprintf (fp,
               "classification=%d\nreason=%s\nt=%.12g\n"
               "drop_radius=%.12g\ndrop_axial_position=%.12g\n"
               "drop_axial_velocity=%.12g\nthreshold=%.12g\nconsecutive=%d\n",
               decided ? (u_big > 0.) : -1,
               decided ? (u_big > 0. ? "ejected_end_pinchoff_drop"
                                     : "end_pinchoff_not_ejected")
                       : "running",
               t, r_big, x_big, u_big, VISIBLE_RADIUS, visible_consecutive);
      fclose (fp);
      rename ("classification.status.tmp", "classification.status");
    }
  }

  if (n > 0) {
    free (volume); free (momentum); free (xmoment);
    free (cells); free (xmin); free (xmax);
  }
  /* deliberately no early return: every case runs to its horizon */
}

event dropmap_close (t = end)
{
  if (fp_components) { fclose (fp_components); fp_components = NULL; }
  if (fp_jet)        { fclose (fp_jet);        fp_jet = NULL; }
}
