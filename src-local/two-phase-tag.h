/**
# Two-phase Interfacial Flows with Interface Tagging

This file extends the standard [two-phase.h](two-phase.h) setup for flows of
two fluids separated by an interface. It is typically used in combination with
a [Navier-Stokes solver](navier-stokes/centered.h).

## Overview

The interface between the fluids is tracked with a Volume-Of-Fluid method.
The volume fraction in fluid 1 is $f=1$ and $f=0$ in fluid 2. The densities
and dynamic viscosities for fluid 1 and 2 are `rho1`, `mu1`, `rho2`, `mu2`,
respectively.

## Interface Fields

This header declares **two** VOF fields:

- `f[]`: The primary VOF field used for all physical calculations (density,
  viscosity, surface tension). This is the field that appears in the
  Navier-Stokes equations.

- `ftag[]`: A secondary VOF field used **only for tracking purposes**. This
  field can be manipulated (e.g., filtered with `tag.h` to remove satellite
  droplets) without affecting the physical simulation.

## Typical Usage

In `coalescenceBubble-tag.c`, the `ftag` field is used with connected
component analysis (`tag.h`) to identify and track only the main bubble,
filtering out small satellite droplets that may form during coalescence.

```c
// In main():
f.sigma = 1.0;      // Physical surface tension
ftag.sigma = 0.0;   // No surface tension (tracking only)

// In logging event:
foreach() ftag[] = f[];           // Copy current interface
// ... apply tag() filtering to ftag ...
// ... compute diagnostics using ftag (filtered) ...
```

## When to Use This vs Standard two-phase.h

- Use **two-phase.h** for most simulations where you only need to track
  a single interface and don't need satellite droplet filtering.

- Use **two-phase-tag.h** when you need to:
  - Filter out satellite droplets for cleaner diagnostics
  - Track multiple disconnected regions separately
  - Compute geometric measurements on only the main fluid body

## See also

- [Two-phase interfacial flows](two-phase.h)
- [Two-phase interfacial flows with coupled VOF and levelset](two-phase-clsvof.h)
- [Two-phase interfacial flows with levelset](two-phase-levelset.h)
- [Connected component analysis](tag.h)
*/

#include "vof.h"

/**
## VOF Field Declarations

Both fields are added to the `interfaces` list for proper advection and
refinement handling by the VOF solver.
*/

scalar f[], ftag[], * interfaces = {f, ftag};

#include "two-phase-generic.h"
