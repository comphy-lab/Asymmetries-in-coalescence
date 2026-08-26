/**
# getJetFoot-halfspace.c — per-snapshot jet base + fluxes (χ → 0)

Dump-probe twin of Bursting-Bubble `getJetFoot.c` / `getBase.c` for
drop-injection snapshots, where $f=1$ is gas.

Emits one stderr line:

```
t  z_base r_base  q_jet q_l  z_tip r_tip
```

Sentinel $-1000$ if a candidate is missing. Compile from `postProcess/`:

```
qcc -I../src-local -O2 -Wall -disable-dimensions \
  getJetFoot-halfspace.c -o getJetFoot-halfspace -lm
./getJetFoot-halfspace intermediate/snapshot-0.5000
```

The on-the-fly time series from `halfspaceJet.c` (`foot.dat`) is the
primary $K_q$ source; this tool is for post-hoc snapshot checks.
*/

#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "fractions.h"
#include "tag.h"
#include "jetFoot.h"

char filename[80];

int main (int a, char const * arguments[])
{
  if (a < 2) {
    fprintf (ferr, "Usage: %s <snapshot>\n", arguments[0]);
    return 1;
  }
  sprintf (filename, "%s", arguments[1]);
  restore (file = filename);
#if TREE
  f.prolongation = fraction_refine;
#endif
  boundary ((scalar *){f, u.x, u.y});

  JetFoot foot = measure_jet_foot();
  fprintf (ferr, "%f %7.6e %7.6e %7.6e %7.6e %7.6e %7.6e\n",
           t, foot.z_base, foot.r_base, foot.q_jet, foot.q_l,
           foot.z_tip, foot.r_tip);
  fflush (ferr);
  fclose (ferr);
  return 0;
}
