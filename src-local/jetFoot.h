/**
# jetFoot.h — outer-surface jet base and base fluxes

Ported from `comphy-lab/Bursting-Bubble` (`getBase.c` protocol plus the
`getJetFoot.c` flux definition).

**VOF sense flipped as part of that port.** The source repository uses
$f=1$ as liquid; this repository standardises on $f=1$ as the bubble/gas
phase in every code (see `AGENTS.md`), so the liquid/gas tags below are
inverted relative to the original. There is no compile-time switch: a
liquid-is-$f=1$ caller does not exist in this repository, and the former
`JETFOOT_F_IS_LIQUID` toggle was removed so the convention cannot fork
silently.

Protocol (Vatsal / Bursting-Bubble `getBase.c`):

1. Tag pure-liquid cells → `MainLiq` = largest component (drops detached
   droplets).
2. Tag pure-gas cells → `MainGas` = largest component (atmosphere + cavity;
   entrained bubbles are separate).
3. An interfacial cell is on the outer free surface iff a face neighbour is
   `MainGas` and a face neighbour is `MainLiq`.
4. Base = lowest such cell (min axial $x$) with $y < $ `JETFOOT_RCAV`.
5. Fluxes on the single cell-row at $z_{\mathrm{base}}$, $0 < y < r_{\mathrm{base}}$,
   no $2\pi$, no $f$-weighting:
   $$q_{\mathrm{jet}} = \sum u_x\, y\,\Delta,\qquad
     q_l = \sum u_x\,\Delta.$$

Paper conversion (not applied here; see `postProcess/fitJetKq.py`):
$Q_j = 2\pi q_{\mathrm{jet}}$, $q_j = 2 q_{\mathrm{jet}}/r_j$,
$q_j = K_q r_j^{\nu}$.

Sentinels are $-1000$ when no outer-surface base exists.

MPI-safe: size tallies are `Allreduce`'d; the base uses Basilisk
`reduction(min:)`. Call from every rank; only rank 0 should write files.
*/

#ifndef JETFOOT_H
#define JETFOOT_H

#ifndef JETFOOT_RCAV
#define JETFOOT_RCAV 1.20
#endif
#ifndef JETFOOT_R_TIP
#define JETFOOT_R_TIP 0.25
#endif
typedef struct {
  double r_base, z_base, q_jet, q_l;
  double r_tip, z_tip;
} JetFoot;

static int jetfoot_largest_tag (scalar s, int n)
{
  int main_id = 0;
  if (n <= 0)
    return 0;
  double * sz = calloc (n, sizeof(double));
  foreach (serial)
    if (s[] > 0)
      sz[(int) s[] - 1] += 1.;
#if _MPI
  MPI_Allreduce (MPI_IN_PLACE, sz, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
  double sm = -1.;
  for (int j = 0; j < n; j++)
    if (sz[j] > sm) {
      sm = sz[j];
      main_id = j + 1;
    }
  free (sz);
  return main_id;
}

static JetFoot measure_jet_foot (void)
{
  JetFoot foot = {-1000., -1000., -1000., -1000., -1000., -1000.};

  scalar dl[], dg[];
  foreach() {
    /* repository convention: f = 1 is the bubble/gas phase (flipped from
       the Bursting-Bubble original as part of the port) */
    dl[] = (f[] < 1e-4);
    dg[] = (f[] > 1. - 1e-4);
  }
  int nliq = tag (dl), ngas = tag (dg);
  int MainLiq = jetfoot_largest_tag (dl, nliq);
  int MainGas = jetfoot_largest_tag (dg, ngas);
  boundary ((scalar *){dl, dg});

  double zbase = HUGE, ztip = -1000.;
  foreach (reduction(min:zbase) reduction(max:ztip)) {
    if (f[] <= 1e-6 || f[] >= 1. - 1e-6 || y > JETFOOT_RCAV)
      continue;
    if (y < JETFOOT_R_TIP && x > ztip)
      ztip = x;
    bool touchGas = ((int) dg[1,0] == MainGas) || ((int) dg[-1,0] == MainGas) ||
                    ((int) dg[0,1] == MainGas) || ((int) dg[0,-1] == MainGas);
    bool touchLiq = ((int) dl[1,0] == MainLiq) || ((int) dl[-1,0] == MainLiq) ||
                    ((int) dl[0,1] == MainLiq) || ((int) dl[0,-1] == MainLiq);
    if (touchGas && touchLiq && x < zbase)
      zbase = x;
  }

  double rbase = HUGE, rtip = HUGE;
  foreach (reduction(min:rbase) reduction(min:rtip)) {
    if (f[] <= 1e-6 || f[] >= 1. - 1e-6 || y > JETFOOT_RCAV)
      continue;
    if (ztip > -900. && y < JETFOOT_R_TIP && x == ztip && y < rtip)
      rtip = y;
    bool touchGas = ((int) dg[1,0] == MainGas) || ((int) dg[-1,0] == MainGas) ||
                    ((int) dg[0,1] == MainGas) || ((int) dg[0,-1] == MainGas);
    bool touchLiq = ((int) dl[1,0] == MainLiq) || ((int) dl[-1,0] == MainLiq) ||
                    ((int) dl[0,1] == MainLiq) || ((int) dl[0,-1] == MainLiq);
    if (touchGas && touchLiq && zbase != HUGE && x == zbase && y < rbase)
      rbase = y;
  }

  double qjet = 0., ql = 0.;
  int haveBase = (zbase != HUGE && rbase != HUGE && rbase > 0.);
  if (haveBase) {
    foreach (reduction(+:qjet) reduction(+:ql)) {
      if (fabs (x - zbase) < 0.5*Delta && y > 0. && y < rbase) {
        qjet += u.x[]*y*Delta;
        ql   += u.x[]*Delta;
      }
    }
    foot.z_base = zbase;
    foot.r_base = rbase;
    foot.q_jet = qjet;
    foot.q_l = ql;
  }
  if (ztip > -900. && rtip != HUGE) {
    foot.z_tip = ztip;
    foot.r_tip = rtip;
  }
  return foot;
}

#endif
