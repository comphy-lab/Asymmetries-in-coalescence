/**
# Getting Center of Mass Data from Simulation Snapshot

Utility for extracting center of mass (COM) position and velocity from
Basilisk VOF snapshots. Uses axisymmetric volume integration weighted by
the VOF field `f[]`.

## Coordinate System

In axisymmetric geometry:
- `x`: axial coordinate (along symmetry axis)
- `y`: radial coordinate (perpendicular to axis)

The center of mass is computed for the axial direction only.

## Output Format

Outputs three space-separated values to stderr:

```
t zCOM uCOM
```

Where:
- `t`: simulation time
- `zCOM`: z-position (axial coordinate) of center of mass
- `uCOM`: z-velocity of center of mass

## Physics

The center of mass position and velocity are computed as volume-weighted
integrals over the VOF field:

$$x_{COM} = \frac{\int 2\pi y \cdot f \cdot x \, dA}{\int 2\pi y \cdot f \, dA}$$

$$u_{COM} = \frac{\int 2\pi y \cdot f \cdot u_x \, dA}{\int 2\pi y \cdot f \, dA}$$

where $f$ is the volume fraction ($f=1$ inside bubble, $f=0$ outside).

## Usage

```
./getCOM <snapshot-file>
```

## Author

Vatsal Sanjay
vatsal.sanjay@comphy-lab.org
CoMPhy Lab, Durham University

Last updated: Jan 2026
*/

#include "axi.h"
#include "navier-stokes/centered.h"
#include "fractions.h"

char filename[4096];
scalar f[];

int main(int a, char const *arguments[])
{
  if (a != 2) {
    fprintf(stderr, "Error: Expected 1 argument\n");
    fprintf(stderr, "Usage: %s <snapshot-file>\n", arguments[0]);
    return 1;
  }

  snprintf(filename, sizeof(filename), "%s", arguments[1]);

  /**
  Restore the simulation snapshot and set up proper boundary conditions: */
  restore (file = filename);
  f.prolongation = fraction_refine;
  boundary((scalar *){f, u.x, u.y});

  /**
  Calculate the center of mass position and velocity by integrating
  over all cells weighted by volume fraction and cell volume: */
  double xcom = 0.0, ucom = 0.0, wt = 0.0;

  foreach() {
    double fval = clamp(f[], 0.0, 1.0);
    double weight = 2*pi*y*fval*sq(Delta);
    xcom += weight * x;
    ucom += weight * u.x[];
    wt += weight;
  }

  if (wt > 0.0) {
    xcom /= wt;
    ucom /= wt;
  }

  fprintf(ferr, "%g %g %g", t, xcom, ucom);

  return 0;
}
