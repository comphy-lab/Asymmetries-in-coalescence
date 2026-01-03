/* Title: Getting Center of Mass Data from Simulation Snapshot
# Author: Vatsal Sanjay
# vatsal.sanjay@comphy-lab.org
# CoMPhy Lab
# Durham University
# Last updated: Jan 2026
*/

#include "axi.h"
#include "navier-stokes/centered.h"
#include "fractions.h"

/**
 * Utility for extracting center of mass (COM) position and velocity from
 * Basilisk VOF snapshots. Uses axisymmetric volume integration weighted by
 * the VOF field f[].
 *
 * Output format (space-separated to stderr):
 *   t zCOM uCOM
 *
 * Where:
 *   t     - simulation time
 *   zCOM  - z-position of center of mass (axial coordinate)
 *   uCOM  - z-velocity of center of mass
 *
 * Physics:
 *   xcom = integral(2*pi*y * f * x * dA) / integral(2*pi*y * f * dA)
 *   ucom = integral(2*pi*y * f * u.x * dA) / integral(2*pi*y * f * dA)
 *
 * Usage: ./getCOM <snapshot-file>
 */

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

  // Restore snapshot
  restore (file = filename);
  f.prolongation = fraction_refine;
  boundary((scalar *){f, u.x, u.y});

  // Calculate the center of mass of f[]
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
