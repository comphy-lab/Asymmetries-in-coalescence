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
  matched-$\Delta$ ladder spacing). South-pole-to-floor clearance is
  \(\approx 4\) at zWall=4. Doubling that gap uses `wallClearance=8`,
  `LDOMAIN=20`, MAXlevel 14 so \(\Delta\) is unchanged;
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

/* The census (components.log, jet.log) lives in src-local/dropCensus.h so
   the finite-Rr wrapper logs in exactly the same format. It refuses to
   build under MPI. */
#include "dropCensus.h"
