/**
# Drop-map census (shared by the half-space and finite-$R_r$ wrappers)

Writes `components.log` (every liquid component every `TLOG`) and `jet.log`
(jet tip, cap velocity, jet radius at the station $x = x_0$), so that any
drop-radius gate can be re-applied offline. The half-space wrapper
`burstingBubbleInfiniteRr.c` and the finite wrapper
`coalescenceBubbleFiniteRr.c` include this header after
`coalescenceBubble.c`, so both geometries log in the same format.

Coordinates: in both geometries the parent interface's undisturbed south pole
sits at $x = 0$ and the small bubble's south pole near $x = -2-\delta-R_f$.
For finite $R_r$ this is the convention of
`postProcess/generate_finiteRr_IC.py --align-parent-pole`.

`CENSUS_FINITE` (default 0) changes one thing: on the axis the outer liquid
lies beyond the parent's north pole, so the jet-tip search is restricted to
$x$ below the top of the near-axis gas column (the parent's north pole). In the
half-space geometry the region above the pool is gas and no restriction is
needed, so the half-space output is unchanged.

The census accumulates per-component sums in a `foreach(serial)` loop, which
is rank-local: correct under OpenMP, wrong under MPI. The build fails under
MPI rather than producing wrong numbers.
*/

#if _MPI
#error "dropCensus.h instrumentation is OpenMP-only; compile without -D_MPI"
#endif

#ifndef CENSUS_FINITE
#define CENSUS_FINITE 0
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
  /* x_far: upper bound of the tip search. Half-space: none (the region above
     the pool is gas up to the domain edge). Finite: the top of the near-axis
     gas column, i.e. the parent's north pole. The jet and any drops are inside
     that column, so the bound does not depend on how far the jet has gone;
     beyond it lies only the outer liquid. */
  double x_far = HUGE;
#if CENSUS_FINITE
  double x_gas_top = -HUGE;
  foreach (reduction(max:x_gas_top))
    if (f[] > 0.5 && y < RTIP && x > x_gas_top)
      x_gas_top = x;
  /* The parent's cap sags below its axial north pole by about
     y^2/(2 Rr) at radius y, so within the RTIP band outer liquid sits as low as
     x_gas_top - RTIP^2/(2 Rr). Keep the bound below that, with a few cells of
     margin; a jet tip that close to the far wall is outside the census. */
  if (x_gas_top > -HUGE)
    x_far = x_gas_top - sq(RTIP)/Rr - 4.*L0/exp2 ((double) MAXlevel);
#endif
  double ztip = -HUGE;
  foreach (reduction(max:ztip))
    if ((1. - f[]) > 0.5 && y < RTIP && x < x_far && x > ztip)
      ztip = x;

  double vol_cap = 0., mom_cap = 0., r_cap = 0.;
  if (ztip > -HUGE) {
    foreach (reduction(+:vol_cap) reduction(+:mom_cap) reduction(max:r_cap))
      if ((1. - f[]) > 0.5 && y < RTIP && x > ztip - CAPDEPTH && x < x_far) {
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
