/**
# Getting Interface Facets from Simulation Snapshot

Utility for extracting interface facets from Basilisk VOF snapshots using
piecewise linear interface reconstruction (PLIC/MYC approximation).

## Interface Reconstruction Method

The Volume-of-Fluid (VOF) method represents interfaces by storing volume
fractions in each cell. To visualize the interface, we use:

1. **PLIC** (Piecewise Linear Interface Calculation): Reconstructs a linear
   segment in each interfacial cell ($0 < f < 1$)
2. **MYC** (Mixed Youngs' Centered): Basilisk's default normal estimation
   algorithm for determining interface orientation

## Output Format

Outputs gnuplot-compatible line segments to stderr:

```
x1 y1
x2 y2
[blank line]
x3 y3
x4 y4
[blank line]
...
```

Each pair of coordinates defines one interface segment. Blank lines separate
disconnected segments, making the output directly usable with gnuplot:

```
gnuplot> plot 'facets.dat' with lines
```

## Usage

```
./getFacet <snapshot-file>
```

## Author

Vatsal Sanjay
vatsal.sanjay@comphy-lab.org
CoMPhy Lab, Durham University

Last updated: Jan 2026
*/

#include "utils.h"
#include "output.h"
#include "fractions.h"

scalar f[];
char filename[4096];

int main(int a, char const *arguments[])
{
  if (a != 2) {
    fprintf(stderr, "Error: Expected 1 argument\n");
    fprintf(stderr, "Usage: %s <snapshot-file>\n", arguments[0]);
    return 1;
  }

  snprintf(filename, sizeof(filename), "%s", arguments[1]);
  restore (file = filename);

  /**
  Set boundary condition: no fluid at left (axis of symmetry).
  Also configure proper VOF refinement for interface cells. */
  f[left] = dirichlet(0.);
  f.prolongation = fraction_refine;
  f.dirty = true;

  /**
  Output facets (interface segments where $0 < f < 1$): */
  FILE * fp = ferr;
  output_facets(f, fp);
  fflush (fp);
  fclose (fp);

  return 0;
}
