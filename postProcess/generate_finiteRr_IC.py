#!/usr/bin/env python3
"""
generate_finiteRr_IC.py — resolution-aware finite radius-ratio two-bubble
coalescence initial condition, for the solver's `DataFiles/InitialConditionRr-<Rr>.dat`.

## Geometry (unchanged from the legacy generator)

Two tangent spherical bubbles sit on the axis of symmetry (the small one
nearest the wall). The small bubble has radius 1 (the reference length); the
large bubble has radius `Rr`. A point contact between two spheres cannot be
resolved, so the contact is regularised by a small circular fillet (a
capillary neck bridge), the same GetCircles convention used throughout this
project (see `generate_Bo0_IC.py`, `comphy-lab/Circle-Contacts-Line`).

  delta = neck-regularisation scale                 (--delta)
  phic1 = arcsin(2*delta)               (small-bubble tangent angle)
  Yfc   = (1+delta)*tan(phic1)          (fillet-circle centre, radial)
  Rf    = (1+delta)/cos(phic1) - 1      (fillet radius)
  cx3   = sqrt((Rr+Rf)**2 - Yfc**2)     (large-bubble centre, axial;
                                         fillet & large circle externally
                                         tangent: |C_fillet-C_large| = Rr+Rf)
  te    = atan2(-Yfc, cx3)              (polar angle of the tangent point on
                                         the fillet circle centred at (0,Yfc))

  1. small bubble  (centre (-(1+delta), 0), radius 1)
       phi1 in [pi, phic1];         X = -(1+delta)+cos(phi)   Y = sin(phi)
  2. fillet        (centre (0, Yfc), radius Rf)
       phi2 in [pi/2-phic1, -pi/2-te]; X = -Rf*sin(phi)       Y = Yfc-Rf*cos(phi)
  3. large bubble  (centre (cx3, 0), radius Rr)
       phi3 in [pi+te, 0];          X = cx3+Rr*cos(phi)       Y = Rr*sin(phi)

This is bit-compatible with the formulas in
`stokes:~/running-simulations/drop-injection-simulations/rr20-large-domain-validation/postProcess/generate_finiteRr_IC.py`
(verified: same `circle_params`, same per-arc parameterisation). Only the
SAMPLING changes — see below.

## Defects fixed relative to the legacy generator

1. **Fixed 1000-point-per-arc sampling.** At Rr=30, MAXlevel=15 (cell size
   Delta ~= 0.00201), the legacy generator's large-circle arc had chords of
   ~0.094 -- 47x the cell size -- which reads to Basilisk's `distance.h` /
   VOF fraction reconstruction as facet-kink curvature far above the true
   1/Rr curvature, injecting spurious capillary forcing. Simultaneously the
   fillet arc (radius ~delta) was ~500x over-sampled for no benefit.
   Fixed here by **arc-length-uniform sampling**: every arc gets exactly
   enough equal-angle points so its chord length is <= `h = Delta_target/2`,
   with `Delta_target = Ldomain/2**MAXlevel` computed from the consumer's own
   `Ldomain` formula (see `target_geometry()` below), and a floor of
   `--n-min` (default 64) points per arc regardless of how few the chord
   bound would otherwise demand.
2. **Exactly-duplicated junction points.** The legacy generator concatenates
   three independently-parameterised arcs whose shared tangency points are
   computed by different trigonometric expressions; they agree to the
   analytic value but differ at the ~1e-16 level, leaving a near-zero-length
   polyline edge at each of the two junctions. Basilisk's 2-D `input_xy` +
   `distance()` has no degenerate-segment filter, so these near-zero edges
   can poison the signed-distance computation. Fixed here by a general
   post-concatenation de-duplication pass: any consecutive point pair closer
   than 1e-12 collapses to a single point.

## Self-checks (run on every invocation; any failure aborts loudly)

- Every output point lies on its analytic arc to within 1e-12 (radius
  deviation from the arc's circle centre).
- Each arc's angular parametrization is strictly monotonic.
- The two endpoints hit the exact analytic tangency points (south pole of
  the small bubble, north pole of the large bubble).
- No consecutive-pair spacing in the final (de-duplicated) polyline is below
  1e-12, and none exceeds the target chord bound `h`.

## Output format

Identical to the legacy generator and to what `coalescenceBubble.c` reads via
`input_xy(fp)`: a plain two-column `x y` polyline, one point per line, no
header, ordered south pole -> small bubble -> fillet -> large bubble -> north
pole (`fscanf(fp, "%lf %lf", ...)` in a loop; consecutive-point polyline, NOT
explicit segment pairs). The consumer's fixed filename convention is
`InitialConditionRr-%3.2f.dat` (C `snprintf`), matched here by Python's
`f"{Rr:.2f}"` (verified identical for all Rr in the existing DataFiles/ tree,
including two-digit values like 30.00 and one-digit values like 5.00).

## Target-resolution formula (verified against coalescenceBubble.c)

For `geometryMode=finite` (the only mode this generator serves) with
`wallClearance > 0`:

    zWall   = delta + wallClearance          (coalescenceBubble.c derives this
                                               from shapeSouthPole = -(2+delta)
                                               and originX = shapeSouthPole -
                                               wallClearance, then
                                               zWall = -2 - originX)
    Ldomain = zWall + 2 + 2*Rr + 4           (finite branch; UNCAPPED)

Note the `halfspace` geometryMode branch (`Bo0.0000.dat`, not produced by
this script) uses a *different*, Rr-independent formula,
`Ldomain = min(zWall + 6, 16)` -- irrelevant here but flagged so nobody
reuses this script's Ldomain formula for that geometry by mistake.

## Usage

    python3 generate_finiteRr_IC.py --Rr 30 --delta 0.01 --maxlevel 15 \\
        --outdir DataFiles-M15-delta0010

    python3 generate_finiteRr_IC.py --Rr 12.5   # legacy-compatible defaults
"""
from __future__ import annotations

import argparse
import hashlib
import math
import os
import sys

import numpy as np

WALL_CLEARANCE_DEFAULT = 0.027  # matches the finite-map zWall=0.05 clearance
                                 # convention documented in coalescenceBubble.c


def circle_params(delta: float, Rr: float) -> dict:
    """Analytic circle parameters, bit-compatible with the legacy generator."""
    phic1 = np.arcsin(2 * delta)
    Yfc = (1 + delta) * np.tan(phic1)
    Rf = (1 + delta) / np.cos(phic1) - 1
    cx3 = np.sqrt((Rr + Rf) ** 2 - Yfc ** 2)
    te = np.arctan2(-Yfc, cx3)
    return dict(phic1=phic1, Yfc=Yfc, Rf=Rf, cx3=cx3, te=te)


def target_geometry(Rr: float, delta: float, maxlevel: int,
                     wall_clearance: float = WALL_CLEARANCE_DEFAULT) -> dict:
    """Reproduce coalescenceBubble.c's zWall / Ldomain / Delta_target for the
    finite geometry branch (wallClearance > 0)."""
    zWall = delta + wall_clearance
    Ldomain = zWall + 2.0 + 2.0 * Rr + 4.0
    Delta_target = Ldomain / (2.0 ** maxlevel)
    h = Delta_target / 2.0
    return dict(zWall=zWall, Ldomain=Ldomain, Delta_target=Delta_target, h=h)


def points_for_arc(R: float, phi_start: float, phi_end: float, h: float,
                    n_min: int = 64) -> int:
    """Minimum equal-angle point count so every chord on this arc is <= h,
    with a floor of n_min points."""
    span = abs(phi_end - phi_start)
    if span <= 0.0:
        return max(2, n_min)
    ratio = h / (2.0 * R)
    if ratio >= 1.0:
        n_seg = 1
    else:
        dphi_max = 2.0 * math.asin(ratio)
        n_seg = max(1, math.ceil(span / dphi_max))
    return max(n_seg + 1, n_min)


def build_arcs(Rr: float, delta: float, h: float, n_min: int) -> tuple[list, dict]:
    """Build the three arcs at arc-length-uniform (chord <= h) resolution."""
    p = circle_params(delta, Rr)
    phic1, Yfc, Rf, cx3, te = p["phic1"], p["Yfc"], p["Rf"], p["cx3"], p["te"]

    arcs = []

    # 1. small bubble: south pole -> neck tangent
    c1 = (-(1.0 + delta), 0.0)
    phi1_start, phi1_end = math.pi, float(phic1)
    n1 = points_for_arc(1.0, phi1_start, phi1_end, h, n_min)
    phi1 = np.linspace(phi1_start, phi1_end, n1)
    X1 = c1[0] + np.cos(phi1)
    Y1 = np.sin(phi1)
    arcs.append(dict(name="small", center=c1, R=1.0, phi=phi1, X=X1, Y=Y1))

    # 2. fillet / neck bridge
    c2 = (0.0, float(Yfc))
    phi2_start, phi2_end = math.pi / 2 - float(phic1), -math.pi / 2 - float(te)
    n2 = points_for_arc(float(Rf), phi2_start, phi2_end, h, n_min)
    phi2 = np.linspace(phi2_start, phi2_end, n2)
    X2 = -Rf * np.sin(phi2)
    Y2 = Yfc - Rf * np.cos(phi2)
    arcs.append(dict(name="fillet", center=c2, R=float(Rf), phi=phi2, X=X2, Y=Y2))

    # 3. large bubble: neck tangent -> north pole
    c3 = (float(cx3), 0.0)
    phi3_start, phi3_end = math.pi + float(te), 0.0
    n3 = points_for_arc(Rr, phi3_start, phi3_end, h, n_min)
    phi3 = np.linspace(phi3_start, phi3_end, n3)
    X3 = cx3 + Rr * np.cos(phi3)
    Y3 = Rr * np.sin(phi3)
    arcs.append(dict(name="large", center=c3, R=Rr, phi=phi3, X=X3, Y=Y3))

    return arcs, p


def dedup_polyline(X: np.ndarray, Y: np.ndarray, seg: np.ndarray,
                    tol: float = 1e-12) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Collapse any consecutive point pair closer than tol into one point.
    Chain-safe: compares each candidate against the last *kept* point, not
    the last *original* point, so runs of near-duplicates collapse fully."""
    keepX = [X[0]]
    keepY = [Y[0]]
    keepSeg = [seg[0]]
    for i in range(1, len(X)):
        dx = X[i] - keepX[-1]
        dy = Y[i] - keepY[-1]
        if math.hypot(dx, dy) < tol:
            continue
        keepX.append(X[i])
        keepY.append(Y[i])
        keepSeg.append(seg[i])
    return np.array(keepX), np.array(keepY), np.array(keepSeg)


def self_check(X: np.ndarray, Y: np.ndarray, seg: np.ndarray, arcs: list,
                delta: float, Rr: float, h: float) -> np.ndarray:
    """Run all built-in self-checks; raise AssertionError (fail loudly) on
    any violation. Returns the array of consecutive-point chord lengths."""
    # 1. max deviation of sampled points from their analytic arc
    max_dev = 0.0
    for idx, a in enumerate(arcs):
        m = seg == idx
        if not np.any(m):
            continue
        cx, cy = a["center"]
        r = np.hypot(X[m] - cx, Y[m] - cy)
        dev = np.abs(r - a["R"])
        max_dev = max(max_dev, float(dev.max()))
    assert max_dev < 1e-12, f"arc deviation too large: {max_dev:.3e} (limit 1e-12)"

    # 2. monotonic parametrization per arc
    for a in arcs:
        dphi = np.diff(a["phi"])
        assert np.all(dphi > 0) or np.all(dphi < 0), \
            f"non-monotonic phi parametrization in arc '{a['name']}'"

    # 3. endpoints hit the exact analytic tangency points
    assert abs(Y[0]) < 1e-9, f"south pole not on axis: Y[0]={Y[0]:.3e}"
    assert abs(Y[-1]) < 1e-9, f"north pole not on axis: Y[-1]={Y[-1]:.3e}"
    south_pole_x = -(2.0 + delta)
    north_pole_x = arcs[2]["center"][0] + Rr
    assert abs(X[0] - south_pole_x) < 1e-9, \
        f"south pole x mismatch: {X[0]:.10f} vs analytic {south_pole_x:.10f}"
    assert abs(X[-1] - north_pole_x) < 1e-9, \
        f"north pole x mismatch: {X[-1]:.10f} vs analytic {north_pole_x:.10f}"

    # 4/5. junction de-duplication + chord-bound checks on the final polyline
    d = np.hypot(np.diff(X), np.diff(Y))
    assert d.min() > 1e-12, \
        f"duplicate/near-duplicate points remain, min spacing {d.min():.3e}"
    assert d.max() <= h * (1.0 + 1e-6) + 1e-12, \
        f"chord exceeds target bound: max={d.max():.6e} > h={h:.6e}"

    return d


def generate(Rr: float, delta: float, maxlevel: int, outdir: str,
             wall_clearance: float = WALL_CLEARANCE_DEFAULT,
             n_min: int = 64, verbose: bool = True) -> dict:
    """Generate one InitialConditionRr-<Rr>.dat file. Returns a report dict."""
    # Reject degenerate geometry before any trigonometry or allocation: Rr <= 0
    # divides by zero in the large-arc parametrisation, delta <= 0 collapses the
    # fillet, and delta >= 0.5 puts arcsin(2*delta) outside its domain. Sampling
    # is arc-length-uniform, so an unbounded Rr would also size the point array
    # before anything else could complain.
    if not math.isfinite(Rr) or Rr <= 0.0:
        raise ValueError(f"Rr must be finite and greater than zero (got {Rr!r})")
    if not math.isfinite(delta) or not 0.0 < delta < 0.5:
        raise ValueError(
            f"delta must be finite and satisfy 0 < delta < 0.5 (got {delta!r})")

    tg = target_geometry(Rr, delta, maxlevel, wall_clearance)
    h = tg["h"]

    arcs, params = build_arcs(Rr, delta, h, n_min)
    X_raw = np.concatenate([a["X"] for a in arcs])
    Y_raw = np.concatenate([a["Y"] for a in arcs])
    seg_raw = np.concatenate([np.full(len(a["X"]), i) for i, a in enumerate(arcs)])

    X, Y, seg = dedup_polyline(X_raw, Y_raw, seg_raw)

    chords = self_check(X, Y, seg, arcs, delta, Rr, h)

    os.makedirs(outdir, exist_ok=True)
    filename = f"InitialConditionRr-{Rr:.2f}.dat"
    filepath = os.path.join(outdir, filename)
    with open(filepath, "w") as fh:
        for x, y in zip(X, Y):
            fh.write(f"{x:.17g} {y:.17g}\n")

    sha256 = hashlib.sha256()
    with open(filepath, "rb") as fh:
        sha256.update(fh.read())
    digest = sha256.hexdigest()
    size_bytes = os.path.getsize(filepath)

    report = dict(
        Rr=Rr, delta=delta, maxlevel=maxlevel, wall_clearance=wall_clearance,
        zWall=tg["zWall"], Ldomain=tg["Ldomain"], Delta_target=tg["Delta_target"],
        h=h, n_per_arc=[len(a["X"]) for a in arcs],
        n_points_raw=len(X_raw), n_points=len(X),
        n_deduped=len(X_raw) - len(X),
        chord_min=float(chords.min()), chord_max=float(chords.max()),
        filepath=filepath, sha256=digest, size_bytes=size_bytes,
        Rf=float(params["Rf"]),
    )

    if verbose:
        print(f"wrote {filepath}")
        print(f"  Rr={Rr:g}  delta={delta:g}  MAXlevel={maxlevel}  "
              f"wallClearance={wall_clearance:g}")
        print(f"  zWall={tg['zWall']:.6f}  Ldomain={tg['Ldomain']:.6f}  "
              f"Delta_target={tg['Delta_target']:.6e}  h={h:.6e}")
        print(f"  points per arc [small, fillet, large] = {report['n_per_arc']}  "
              f"-> {len(X_raw)} raw, {report['n_deduped']} de-duplicated, "
              f"{len(X)} written")
        print(f"  chord min={report['chord_min']:.6e}  max={report['chord_max']:.6e}  "
              f"(target h={h:.6e})")
        print(f"  south pole x={X[0]:.6f}  north pole x={X[-1]:.6f}")
        print(f"  size={size_bytes} bytes  sha256={digest}")

    return report


def main():
    ap = argparse.ArgumentParser(
        description="Resolution-aware finite radius-ratio coalescence IC generator.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--Rr", type=float, required=True,
                     help="radius ratio R_large / R_small (small bubble = unit)")
    ap.add_argument("--delta", type=float, default=0.023,
                     help="neck-regularisation scale (default 0.023, legacy default)")
    ap.add_argument("--maxlevel", type=int, default=15,
                     help="MAXlevel used to derive the target chord bound (default 15)")
    ap.add_argument("--outdir", type=str, default="simulationCases/DataFiles",
                     help="output directory; file is named InitialConditionRr-<Rr>.2f.dat "
                          "inside it (default simulationCases/DataFiles, legacy-compatible)")
    ap.add_argument("--wall-clearance", type=float, default=WALL_CLEARANCE_DEFAULT,
                     help=f"wallClearance used to derive zWall = delta + wallClearance "
                          f"(default {WALL_CLEARANCE_DEFAULT}, matches coalescenceBubble.c)")
    ap.add_argument("--n-min", type=int, default=64,
                     help="minimum points per arc regardless of chord bound (default 64)")
    a = ap.parse_args()

    try:
        generate(a.Rr, a.delta, a.maxlevel, a.outdir, a.wall_clearance, a.n_min)
    except (AssertionError, ValueError) as exc:
        print(f"SELF-CHECK FAILED: {exc}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
