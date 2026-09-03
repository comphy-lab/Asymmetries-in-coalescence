/**
# burstingBubbleInfiniteRr.c — the $R_r \to \infty$ driver, arbitrary $\chi$

One driver for every infinite-radius-ratio (bursting-bubble) run in this
repository, from the near-wall half-space limit ($\chi \ll 1$) to the bulk
far-from-wall limit ($\chi \to \infty$). It supersedes and replaces two
earlier files:

- `halfspaceJet.c` — the $\chi \ll 1$ jet-foot build (its `foot.dat`
  instrumentation is kept, unchanged);
- `burstingBubbleDropMap.c` — the drop-map wrapper around the
  `singular-bursting-bubbles` drill solver (its `components.log` / `jet.log`
  instrumentation is kept and extended).

The physics core is this repository's own `coalescenceBubble.c`: `axi.h` +
`navier-stokes/centered.h` + `#define FILTERED 1` + `two-phase.h` +
`double-projection.h` (Stack 1). That stack is the empirically most stable
configuration for the asymmetric-coalescence family (per Vatsal, 2026-08-29,
and the solver-stack discussion in `coalescenceBubble.c` itself), and using
the in-repo core removes the cross-repository staged dependency the drop-map
wrapper needed: this file compiles from a bare checkout.

Byte-identity with the `singular-bursting-bubbles` drill solver is
deliberately **not** claimed. Runs made with this driver are a different
numerical realisation of the same physics; cross-driver consistency is a
thing to test (matched-$\Delta$ spot checks against the drill-solver ladder),
not to assume.

## Scientific context — the two critical Ohnesorge numbers

The first ejected drop does not shrink monotonically with viscosity. On the
falling branch, increasing $Oh$ sharpens the capillary-wave focusing and
thins the jet, so the drop shrinks (the Gordillo & Blanco-Rodríguez regime).
Where the viscous cutoff catches up with the inviscid jet radius the trend
turns: past $Oh_{\mathrm{opt}}$ the jet fattens and the drop grows, until at
$Oh_c^{(U)}$ the drop velocity falls through zero and forward injection
stops. Locating both requires, at one grid spacing:

1. $r_d(Oh)$ under a resolution-independent detector, with EVERY liquid
   component logged so the drop is selected in post-processing;
2. the jet time series $v(t)$, plus the jet radius at the undisturbed
   free-surface station, so $r_{jet}(Oh)$, $v_{jet,0}(Oh)$ and
   $\max_t v(t)\,(Oh)$ can be measured independently of any drop detector;
3. no early termination, so the full drop sequence is observed.

## Outputs

Everything `coalescenceBubble.c` writes (`log`, snapshots, restart,
classification), plus:

| file | contents |
|---|---|
| `foot.dat` | `i dt t ke maxlevel r_b z_b r_base z_base q_jet q_l` at `tsnap2` cadence (jet-foot protocol; `ENABLE_JET_FOOT`) |
| `components.log` | `t`, component index, `is_main`, cell count, volume, equivalent radius, axial centroid, mean axial velocity, axial extent — EVERY component, every `TLOG` |
| `jet.log` | `t`, tip position, tip-cap velocity, cap radius, cap volume, `r_jet_z0`, kinetic energy, component count, largest detached component |

`r_jet_z0` is new relative to the drop-map wrapper: the minimum radial
position of interface cells in the thin axial slab $|x| \le Z0_SLAB$ within
the near-axis band. Once the jet tip has crossed the undisturbed
free-surface plane $x = 0$ this is the jet radius at that station — the
$r_{jet}$ of the Gordillo & Blanco-Rodríguez protocol, where $v_{jet,0}$ is
read at the same crossing instant (interpolated from `z_tip(t)` offline).
Before the crossing the slab holds no near-axis interface and the column is
`nodata`.

**Phase convention.** This repository standardises on `f = 1` as the
bubble/gas phase in every code (see `AGENTS.md`). The drop-map
instrumentation below was ported from the `comphy-lab/Bursting-Bubble`
drill-solver wrapper, where `f = 1` is liquid, and the VOF sense was
flipped as part of that port: every liquid measurement uses `1 - f`. The
flip was caught by the first smoke run of this driver — with the source
convention the "main liquid component" was the atmosphere, volume 1260.8.

`classification.status` is owned by the base solver's pinned detector
(`dropRadiusMin` argv 7, persistence argv 8). `DETACH_STOP=0` freezes the
first verdict and keeps the run going to its horizon.

## Compile

From `simulationCases/`, self-contained:

~~~bash
qcc -I. -I../src-local -O2 -Wall -disable-dimensions -fopenmp \
  burstingBubbleInfiniteRr.c -o burstingBubbleInfiniteRr -lm
~~~

## Usage

Same argv contract as `coalescenceBubble.c`. `geometryMode` is forced to
`halfspace` (the $R_r \to \infty$ initial condition, `Bo%5.4f.dat`;
default $Bo=0$ is `Bo0.0000.dat`).
$\chi$ is set at run time:

- **bulk / far from wall**: leave `wallClearance` (argv 24) unset (`-1`);
  `zWall` (argv 6) alone sets the domain, `Ldomain = min(zWall+6, 16)`
  (zWall=4 gives Ldomain=10, $\Delta = 10/2^{13}$ at MAXlevel 13 — the
  matched-$\Delta$ ladder spacing);
- **near wall** ($\chi \ll 1$): pass the physical clearance, e.g.
  `wallClearance = 0.027`, as in the $K_q$ measurement.

Always pass the pinned detector radius explicitly (argv 7
`= 0.021005127`, never 0: zero makes the solver derive a mesh-tied value).

~~~bash
./burstingBubbleInfiniteRr <OhOut> <RhoIn> <Rr-placeholder> <MAXlevel> \
  <tmax> <zWall> 0.021005127 3 <snapshotInterval> <drillAMR> ... \
  halfspace <wallClearance|-1> <interfaceFloor> <MuRin> <Bond>
~~~
*/

#define ENABLE_JET_FOOT 1
#define FORCE_GEOMETRY_HALFSPACE 1
#define DETACH_STOP 0
#include "coalescenceBubble.c"

/* The dropmap census below accumulates per-component sums in a
   foreach(serial) loop, which is rank-local: correct under OpenMP (shared
   memory), silently wrong under MPI. Fail the build rather than the
   science. */
#if _MPI
#error "burstingBubbleInfiniteRr.c instrumentation is OpenMP-only; compile without -D_MPI"
#endif

#define TLOG 0.005          /* component + jet cadence */
#define RTIP 0.50           /* near-axis band defining the jet tip */
#define CAPDEPTH 0.05       /* axial depth of the tip cap used for v(t) */
#define Z0_SLAB 0.01        /* half-width of the r_jet station at x = 0 */

static FILE * fp_components = NULL;
static FILE * fp_jet = NULL;

static void dropmap_open (void)
{
  if (fp_components)
    return;
  fp_components = fopen ("components.log", "a");
  fp_jet = fopen ("jet.log", "a");
  if (!fp_components || !fp_jet) {
    /* a case without its diagnostics is not a result: fail loudly rather
       than burn node-hours producing an unusable run */
    fprintf (ferr, "dropmap: cannot open %s in the case directory\n",
             !fp_components ? "components.log" : "jet.log");
    exit (1);
  }
  if (ftell (fp_components) == 0)
    fprintf (fp_components,
             "t,component,is_main,cells,volume,r_eq,x_mean,u_x_mean,x_min,x_max\n");
  if (ftell (fp_jet) == 0)
    fprintf (fp_jet,
             "t,z_tip,u_tip_cap,r_cap_max,vol_cap,r_jet_z0,ke,n_components,"
             "r_largest_detached,u_largest_detached,x_largest_detached\n");
}

event dropmap (t = 0.; t += TLOG; t <= tmax + 1e-9)
{
  dropmap_open();

  /* ---- every liquid component ---------------------------------------- */
  scalar component[];
  foreach()
    component[] = (1. - f[]) > 1e-4;   /* f = 1 is GAS in this solver */
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
        double dV = 2.*pi*y*sq(Delta)*clamp (1. - f[], 0., 1.);
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
    if ((1. - f[]) > 0.5 && y < RTIP && x > ztip)
      ztip = x;

  double vol_cap = 0., mom_cap = 0., r_cap = 0.;
  if (ztip > -HUGE) {
    foreach (reduction(+:vol_cap) reduction(+:mom_cap) reduction(max:r_cap))
      if ((1. - f[]) > 0.5 && y < RTIP && x > ztip - CAPDEPTH) {
        double dV = 2.*pi*y*sq(Delta)*clamp (1. - f[], 0., 1.);
        vol_cap += dV;
        mom_cap += dV*u.x[];
        if (y > r_cap)
          r_cap = y;
      }
  }

  /* ---- jet radius at the undisturbed-surface station x = 0 ------------ */
  /* Near the axis at x = 0 the only interface is the jet surface (before
     jet crossing the slab is inside the gas cavity and holds none), so the
     minimum interface radius in the slab IS the jet radius there. The
     f-in-(0.01,0.99) interface test is symmetric under the gas/liquid
     convention flip, so it needs no 1-f. */
  double rjet0 = HUGE;
  foreach (reduction(min:rjet0))
    if (f[] > 0.01 && f[] < 0.99 && fabs(x) <= Z0_SLAB && y < RTIP &&
        y < rjet0)
      rjet0 = y;

  double kinetic = 0.;
  foreach (reduction(+:kinetic))
    kinetic += 2.*pi*y*sq(Delta)*clamp (1. - f[], 0., 1.)*
               (sq(u.x[]) + sq(u.y[]))/2.;

  if (pid() == 0 && fp_jet) {
    fprintf (fp_jet,
             "%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%d,%.10g,%.10g,%.10g\n",
             t, ztip > -HUGE ? ztip : nodata,
             vol_cap > 0. ? mom_cap/vol_cap : nodata,
             r_cap, vol_cap, rjet0 < HUGE ? rjet0 : nodata,
             kinetic, n, r_big, u_big, x_big);
    fflush (fp_jet);
    if (fp_components)
      fflush (fp_components);
  }

  if (n > 0) {
    free (volume); free (momentum); free (xmoment);
    free (cells); free (xmin); free (xmax);
  }
  /* deliberately no early return: every case runs to its horizon
     (DETACH_STOP=0 also disables the base classifier's stop) */
}

event dropmap_close (t = end)
{
  if (fp_components) { fclose (fp_components); fp_components = NULL; }
  if (fp_jet)        { fclose (fp_jet);        fp_jet = NULL; }
}
