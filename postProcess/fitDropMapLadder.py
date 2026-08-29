#!/usr/bin/env python3
"""Pre-registered estimation protocols for the two critical Ohnesorge numbers.

Operates on ``ladder_summary.csv`` from ``dropMapHarvest.py``. This module
deliberately contains NO reference values from the manuscript or the
literature: the protocols are frozen first, the measured numbers come out,
and only then are they compared against anything external. Keeping the
comparators out of the estimation code is what makes the comparison a test
rather than a calibration.

Protocols
---------
``Oh_opt`` (the turning point) is estimated three independent ways; their
agreement is itself a physics test (the mechanism says minimum drop size,
minimum jet radius and maximum jet speed coincide):

  P1  minimiser of r_d(Oh)          — parabola through the discrete minimum
                                      and its neighbours in (log Oh, log r);
  P2  maximiser of v_max(Oh)        — same parabola protocol on (log Oh, v);
  P3  intersection of the falling- and rising-branch power-law fits of
      r_d(Oh) (or r_jet0(Oh)).

``Oh_c^(r)`` (the virtual zero of the falling branch) is the fit parameter of

      r = C (Oh_c - Oh)^p            on falling-branch points only,

estimated by profiling Oh_c on a grid: for each candidate Oh_c the model is
linear in (log(Oh_c - Oh), log r), so C and p come from least squares and
Oh_c minimises the residual. Falling-branch membership is decided by the
discrete minimum, never by hand-picking.

``Oh_c^(U)`` (the forward-injection boundary) is the zero of U_d(Oh) on the
rising side: fit U_d = A (Oh_u - Oh)^q for points with U_d > 0 and
Oh > Oh_opt, same profiling scheme.

Every estimate is reported with its fit window, exponent and residual so a
reader can reject it.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np


# --------------------------------------------------------------------------
# building blocks
# --------------------------------------------------------------------------

def parabola_extremum(x: np.ndarray, y: np.ndarray) -> dict:
    """Vertex of the least-squares parabola through >= 3 points."""
    if len(x) < 3:
        return {"x_star": None, "y_star": None, "n": int(len(x))}
    a, b, c = np.polyfit(x, y, 2)
    if a == 0:
        return {"x_star": None, "y_star": None, "n": int(len(x))}
    xs = -b / (2 * a)
    return {"x_star": float(xs), "y_star": float(np.polyval([a, b, c], xs)),
            "curvature": float(a), "n": int(len(x))}


def profile_powerlaw_zero(oh: np.ndarray, y: np.ndarray, side: str,
                          n_grid: int = 2000) -> dict:
    """Fit y = C (oh_c - oh)^p (side='left', oh_c > max(oh)) or
    y = C (oh - oh_c)^p (side='right', oh_c < min(oh)) by profiling oh_c."""
    oh = np.asarray(oh, float)
    y = np.asarray(y, float)
    keep = np.isfinite(oh) & np.isfinite(y) & (y > 0)
    oh, y = oh[keep], y[keep]
    if len(oh) < 3:
        return {"oh_c": None, "p": None, "C": None, "rms": None,
                "n": int(len(oh))}
    if side == "left":
        lo, hi = oh.max() * (1 + 1e-6), oh.max() * 3

        def dist(ohc):
            return ohc - oh
    elif side == "right":
        lo, hi = oh.min() / 3, oh.min() * (1 - 1e-6)

        def dist(ohc):
            return oh - ohc
    else:
        raise ValueError(side)
    gaps = np.linspace(lo, hi, n_grid)
    best = None
    for ohc in gaps:
        d = dist(ohc)
        if np.any(d <= 0):
            continue
        X = np.log(d)
        Y = np.log(y)
        p, logC = np.polyfit(X, Y, 1)
        resid = Y - (p * X + logC)
        rms = float(np.sqrt(np.mean(resid ** 2)))
        if best is None or rms < best["rms"]:
            best = {"oh_c": float(ohc), "p": float(p),
                    "C": float(math.exp(logC)), "rms": rms, "n": int(len(oh))}
    if best is not None:
        # A minimiser at (or one grid step from) either search bound is a
        # saturated search, not a resolved virtual zero; flag it so the
        # reader can reject the estimate.
        step = (hi - lo) / max(n_grid - 1, 1)
        best["at_search_boundary"] = bool(
            best["oh_c"] <= lo + step or best["oh_c"] >= hi - step)
    return best or {"oh_c": None, "p": None, "C": None, "rms": None,
                    "n": int(len(oh))}


def split_branches(oh: np.ndarray, r: np.ndarray) -> tuple:
    """Falling/rising membership from the discrete minimum of r(Oh).

    Returns (falling_mask, rising_mask, i_min). The minimum point itself
    belongs to neither branch: near the turn both power laws break down.
    """
    i_min = int(np.nanargmin(r))
    falling = np.arange(len(oh)) < i_min
    rising = np.arange(len(oh)) > i_min
    return falling, rising, i_min


# --------------------------------------------------------------------------
# the estimation report
# --------------------------------------------------------------------------

def estimate(rows: list[dict]) -> dict:
    rows = sorted((r for r in rows if r.get("oh") not in (None, "")),
                  key=lambda r: float(r["oh"]))

    def col(name):
        out = []
        for r in rows:
            v = r.get(name)
            try:
                out.append(float(v))
            except (TypeError, ValueError):
                out.append(np.nan)
        return np.array(out)

    oh = col("oh")
    report: dict = {"n_cases": len(rows),
                    "oh_grid": [float(v) for v in oh]}

    for label, series in (("r_visible", col("r_drop_at_confirm")),
                          ("r_pinch_topo", col("r_pinch_topo")),
                          ("r_jet0", col("r_jet0"))):
        ok = np.isfinite(series) & (series > 0)
        if ok.sum() < 3:
            report[label] = {"skipped": f"only {int(ok.sum())} points"}
            continue
        o, s = oh[ok], series[ok]
        falling, rising, i_min = split_branches(o, s)
        blk = {"i_min": int(i_min), "oh_at_min": float(o[i_min]),
               "r_at_min": float(s[i_min])}
        # P1: parabola in log-log around the minimum
        w = slice(max(0, i_min - 2), min(len(o), i_min + 3))
        blk["P1_parabola"] = parabola_extremum(np.log(o[w]), np.log(s[w]))
        if blk["P1_parabola"].get("x_star") is not None:
            blk["P1_parabola"]["oh_opt"] = float(
                math.exp(blk["P1_parabola"]["x_star"]))
        # falling-branch virtual zero
        blk["falling_fit"] = profile_powerlaw_zero(o[falling], s[falling],
                                                   side="left")
        # rising branch (for P3)
        blk["rising_fit"] = profile_powerlaw_zero(o[rising], s[rising],
                                                  side="right")
        # P3: branch intersection
        f, g = blk["falling_fit"], blk["rising_fit"]
        if all(v is not None for v in
               (f["oh_c"], f["p"], f["C"], g["oh_c"], g["p"], g["C"])):
            # Both models must be defined at a candidate point: outside the
            # overlap inf - inf gives NaN and argmin would pick it.
            grid = np.linspace(o[0], o[-1], 20001)
            support = (grid < f["oh_c"]) & (grid > g["oh_c"])
            if support.any():
                gsub = grid[support]
                lhs = f["C"] * (f["oh_c"] - gsub) ** f["p"]
                rhs = g["C"] * (gsub - g["oh_c"]) ** g["p"]
                gap = np.log(lhs) - np.log(rhs)
                # a genuine intersection changes sign; the closest endpoint
                # of a non-crossing pair must not be reported as Oh_opt
                if np.any(gap > 0) and np.any(gap < 0):
                    k = int(np.argmin(np.abs(gap)))
                    blk["P3_branch_intersection_oh"] = float(gsub[k])
                    blk["P3_residual"] = float(abs(gap[k]))
                else:
                    blk["P3_branch_intersection_oh"] = None
                    blk["P3_note"] = ("branch fits do not cross inside the "
                                      "overlap")
            else:
                blk["P3_branch_intersection_oh"] = None
                blk["P3_note"] = "falling and rising supports do not overlap"
        report[label] = blk

    # P2: maximum of v_max(Oh)
    vmax = col("v_max")
    ok = np.isfinite(vmax)
    if ok.sum() >= 3:
        o, v = oh[ok], vmax[ok]
        i_max = int(np.nanargmax(v))
        w = slice(max(0, i_max - 2), min(len(o), i_max + 3))
        p2 = parabola_extremum(np.log(o[w]), v[w])
        if p2.get("x_star") is not None:
            p2["oh_opt"] = float(math.exp(p2["x_star"]))
        report["P2_vmax"] = {"oh_at_max": float(o[i_max]),
                             "v_at_max": float(v[i_max]), **p2}

    # Oh_c^(U): zero of U_d on the rising side
    ud = col("u_drop_at_confirm")
    ok = np.isfinite(ud) & (ud > 0)
    if ok.sum() >= 3:
        o, u = oh[ok], ud[ok]
        # rising side only: past the drop-size minimum if one is resolved
        rv = report.get("r_visible", {})
        oh_min = rv.get("oh_at_min")
        if oh_min is not None:
            sel = o > oh_min
            if sel.sum() >= 3:
                o, u = o[sel], u[sel]
        report["ohc_U"] = profile_powerlaw_zero(o, u, side="left")
        report["ohc_U"]["note"] = ("zero of U_d(Oh); the forward-injection "
                                   "boundary Oh_c^(U)")

    return report


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("summary", help="ladder_summary.csv from dropMapHarvest")
    ap.add_argument("--out", default="ladder_estimates.json")
    args = ap.parse_args(argv)
    with open(args.summary, newline="") as fh:
        rows = list(csv.DictReader(fh))
    report = estimate(rows)
    Path(args.out).write_text(json.dumps(report, indent=2))
    print(f"wrote {args.out}")
    for key in ("r_visible", "P2_vmax", "ohc_U"):
        if key in report:
            print(f"  {key}: {json.dumps(report[key], default=str)[:200]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
