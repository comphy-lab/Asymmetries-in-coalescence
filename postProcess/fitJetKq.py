#!/usr/bin/env python3
"""Fit K_q from a halfspaceJet foot.dat (bursting-bubble log columns).

Paper identities (singular-bursting-bubbles Fig. 2 / make_fig2_flux_scalings.py):

    Q_j = 2 π q_jet
    q_j = Q_j / (π r_j) = 2 q_jet / r_j
    q_j = K_q r_j^ν

Default inertial window is the bursting-bubble cone-fit band
r ∈ [0.005, 0.024]. Pass --nu to freeze the exponent (χ→0: ν* ≈ 0.447
from β* ≈ 43.09°) and fit only K_q; omit --nu to fit both.

Usage:
    python3 postProcess/fitJetKq.py --dat path/to/foot.dat --nu 0.447
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

COLUMNS = (
    "i",
    "dt",
    "t",
    "ke",
    "maxlevel",
    "r_b",
    "z_b",
    "r_base",
    "z_base",
    "q_jet",
    "q_l",
)


def read_foot(path: Path) -> list[dict[str, float]]:
    rows: dict[float, dict[str, float]] = {}
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("MAXlevel"):
                continue
            parts = line.split()
            if not parts or parts[0] == "i":
                continue
            if len(parts) < 11:
                continue
            try:
                values = [float(p) for p in parts[:11]]
            except ValueError:
                continue
            row = dict(zip(COLUMNS, values))
            if row["r_base"] <= -900 or row["z_base"] <= -900:
                continue
            rows[round(row["t"], 8)] = row
    return [rows[k] for k in sorted(rows)]


def q_j_of(row: dict[str, float]) -> float:
    r_j = row["r_base"]
    return 2.0 * row["q_jet"] / r_j


def fit_log(xs: list[float], ys: list[float]) -> tuple[float, float, float]:
    """Return (slope, intercept, R²) of log y = slope log x + intercept."""
    n = len(xs)
    lx = [math.log(x) for x in xs]
    ly = [math.log(y) for y in ys]
    mx = sum(lx) / n
    my = sum(ly) / n
    sxx = sum((x - mx) ** 2 for x in lx)
    sxy = sum((x - mx) * (y - my) for x, y in zip(lx, ly))
    syy = sum((y - my) ** 2 for y in ly)
    if sxx <= 0.0:
        raise SystemExit("fitJetKq: degenerate r_j window")
    slope = sxy / sxx
    intercept = my - slope * mx
    r2 = (sxy * sxy) / (sxx * syy) if syy > 0.0 else 1.0
    return slope, intercept, r2


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dat", required=True, type=Path)
    parser.add_argument("--rmin", type=float, default=0.005)
    parser.add_argument("--rmax", type=float, default=0.024)
    parser.add_argument(
        "--nu",
        type=float,
        default=None,
        help="freeze ν and fit only K_q (q_j = K_q r_j^ν)",
    )
    parser.add_argument(
        "--tmin",
        type=float,
        default=None,
        help="optional inception cut; default keeps all positive-q_jet rows in window",
    )
    args = parser.parse_args()
    if not args.dat.is_file():
        print(f"fitJetKq: missing {args.dat}", file=sys.stderr)
        return 2

    rows = read_foot(args.dat)
    selected: list[dict[str, float]] = []
    for row in rows:
        r_j = row["r_base"]
        qjet = row["q_jet"]
        if r_j <= 0.0 or qjet <= 0.0:
            continue
        if r_j < args.rmin or r_j > args.rmax:
            continue
        if args.tmin is not None and row["t"] <= args.tmin:
            continue
        if q_j_of(row) <= 0.0:
            continue
        selected.append(row)

    if len(selected) < 4:
        print(
            f"fitJetKq: only {len(selected)} points in r∈[{args.rmin},{args.rmax}] "
            f"(need ≥4). foot.dat has {len(rows)} valid base rows.",
            file=sys.stderr,
        )
        return 1

    r_j = [row["r_base"] for row in selected]
    q_j = [q_j_of(row) for row in selected]
    q_jet = [row["q_jet"] for row in selected]

    if args.nu is None:
        nu, intercept, r2 = fit_log(r_j, q_j)
        kq = math.exp(intercept)
        source = "fitted ν and K_q"
    else:
        nu = args.nu
        logs = [math.log(qj) - nu * math.log(r) for r, qj in zip(r_j, q_j)]
        kq = math.exp(sum(logs) / len(logs))
        pred = [kq * (r**nu) for r in r_j]
        mean = sum(q_j) / len(q_j)
        ss_res = sum((y - p) ** 2 for y, p in zip(q_j, pred))
        ss_tot = sum((y - mean) ** 2 for y in q_j)
        r2 = 1.0 - ss_res / ss_tot if ss_tot > 0.0 else 1.0
        source = f"frozen ν={nu:g}"

    qjet_exp, qjet_int, qjet_r2 = fit_log(r_j, q_jet)
    print(f"file          {args.dat}")
    print(f"n             {len(selected)}")
    print(f"r window      [{args.rmin:g}, {args.rmax:g}]")
    print(f"t range       {selected[0]['t']:.6f} → {selected[-1]['t']:.6f}")
    print(f"fit           {source}")
    print(f"nu            {nu:.6f}")
    print(f"K_q           {kq:.6e}     (q_j = K_q r_j^nu, q_j = 2 q_jet / r_j)")
    print(f"R^2(q_j)      {r2:.4f}")
    print(
        f"q_jet ~ r^{qjet_exp:.3f} (R^2={qjet_r2:.3f}); "
        f"prefactor {math.exp(qjet_int):.6e}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
