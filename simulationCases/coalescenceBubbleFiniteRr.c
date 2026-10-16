/**
# Finite-$R_r$ coalescence with the drop-map census

Two bubbles of radius 1 and $R_r$ coalesce through a regularised neck. This
wrapper is the finite-$R_r$ sibling of `burstingBubbleInfiniteRr.c`: same
solver (`coalescenceBubble.c`, Stack 1 by default), same classifier, same
census (`src-local/dropCensus.h`), and every case runs to its horizon.

- `DETACH_STOP 0`: the first classified drop does not stop the run, so the
  whole history is available for offline re-cuts with any gate.
- `CENSUS_FINITE 1`: the jet-tip search excludes the outer liquid beyond the
  parent's north pole.
- The jet-foot diagnostic (`foot.dat`) is a bursting-bubble quantity defined
  on an outer free surface; it is not built here.

## Initial condition and coordinates

Generate the polyline with
`postProcess/generate_finiteRr_IC.py --Rr <Rr> --delta <delta>
--align-parent-pole --ldomain <L0> --maxlevel <MAXlevel>`. With
`--align-parent-pole` the parent's south pole sits at $x = 0$ and the small
bubble's south pole at $-(2 + \delta + R_f)$, exactly the half-space
convention of `Bo0.0000.dat`. The census station $x = 0$ therefore means the
same thing in both geometries, and the $R_r \to \infty$ limit of this family
is `Bo0.0000.dat` itself. Pass the file by name with `INITIAL_SHAPE`.

## Usage

Same positional arguments as `coalescenceBubble.c` (argv 3 is the real
$R_r$; argv 23 must be `finite`). Set `LDOMAIN` = $10\cdot 2^k$ at
MAXlevel $14 + k$ to keep $\Delta = 10/2^{14}$, and `FAR_FLOOR_RADIUS` /
`FAR_FLOOR_DROP` to keep the parent surface away from the axis off the
`MAXlevel` floor.
*/

#define ENABLE_JET_FOOT 0
#define DETACH_STOP 0
#include "coalescenceBubble.c"

#define CENSUS_FINITE 1
#include "dropCensus.h"
