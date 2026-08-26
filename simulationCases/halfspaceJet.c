/**
# halfspaceJet.c — χ → 0 half-space jet-foot solver

Dedicated χ ≪ 1 (near-wall half-space) build of the drop-injection stack.
Numerics are identical to `coalescenceBubble.c` (FILTERED + double-projection
by default). Geometry is forced to `halfspace`. Every logging interval
(`tsnap2 = 10^{-4}`) the bursting-bubble outer-surface base and fluxes are
written to `foot.dat`.

This file exists because production `coalescenceBubble.c` does not log
$q_{\mathrm{jet}}$. Do not use this binary for finite-$R_r$ Oh$_c$ ladders.

## Observables (`foot.dat`)

Columns match singular-bursting-bubbles `log` /
`plotJetMetricsTheory.py` / `make_fig2_flux_scalings.py`:

```
i dt t ke maxlevel r_b z_b r_base z_base q_jet q_l
```

- `r_base`, `z_base`: outer-surface jet base (`getBase.c` protocol).
- `q_jet = \sum u_x y \Delta`  [L³/T], no $2\pi$, no $f$-weighting.
- `q_l   = \sum u_x \Delta`    [L²/T].
- `r_b`, `z_b`: sentinels ($-1000$). This solver has no AMR drill probe;
  the science columns are `r_base` and `q_jet`.

Paper conversion (post-process, `postProcess/fitJetKq.py`):

$$Q_j = 2\pi q_{\mathrm{jet}},\qquad
  q_j = 2 q_{\mathrm{jet}}/r_j,\qquad
  q_j = K_q r_j^{\nu}.$$

The existing drop-injection `log` (ke, Xc, Vcm, classification) is unchanged.

## Compile

From `simulationCases/`:

```
qcc -I. -I../src-local -Wall -O2 -fopenmp -disable-dimensions \
  halfspaceJet.c -o halfspaceJet -lm
```

Dump probe (per snapshot):

```
qcc -I../src-local -O2 -Wall -disable-dimensions \
  ../postProcess/getJetFoot-halfspace.c -o ../postProcess/getJetFoot-halfspace -lm
```

The sweep runners copy `coalescenceBubble.c` only. For a campaign, compile
this file once and copy the binary, or copy `halfspaceJet.c`,
`coalescenceBubble.c` and `src-local/jetFoot.h` into the case directory.

## Usage

Same argv as `coalescenceBubble.c`. `geometryMode` is forced to `halfspace`
even if argv 23 says otherwise. Still pass the production half-space
controls:

- `wallClearance = 0.027` (χ ≪ 1)
- pinned `dropRadiusMin = 0.021005127` as argv 7 (never 0)
- `interfaceFloor = 1`

Placeholder `Rr` is ignored for the Bo=0 IC (`Bo0.0000.dat`).

```
./halfspaceJet <OhOut> <RhoIn> <Rr> <MAXlevel> <tmax> <zWall> \
  0.021005127 3 0.01 0  ...  halfspace 0.027 1
```

## Author

Vatsal Sanjay / CoMPhy Lab. Wrapper for the χ → 0 $K_q$ measurement.
*/

#define ENABLE_JET_FOOT 1
#define FORCE_GEOMETRY_HALFSPACE 1
#include "coalescenceBubble.c"
