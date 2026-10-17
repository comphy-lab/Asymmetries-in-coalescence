#!/usr/bin/env python3
"""
generate_Bo0_IC.py — Bo=0 bursting-bubble initial condition (sphere-plane
contact), for the drill solver's DataFiles/Bo0.0000.dat.

At Bo=0 the bubble is an UNDEFORMED sphere of radius 1 intersecting the flat
free surface — the physically correct shape of a ~10 micron bubble, vs the
gravity-deformed Bo>0 shapes. A point contact cannot be resolved, so the
contact line is regularised by a fillet of radius ~delta (a small capillary
bridge), leaving an initial rupture "hole" of radius ~2*delta at the top.

Geometry from comphy-lab/Circle-Contacts-Line (GetCircles): a radius-1 cavity
arc from the bottom (on axis) up to the neck, a fillet, and the flat free
surface. That code emits (x_axial-into-liquid, y_radial) with the cavity in
+x and the surface at x=0. The drill solver's convention is the OPPOSITE
axial sign (cavity in -x, south pole near -2, free surface at x=0, then the
solver's init negates the distance field), so we flip x here. The polyline
orientation (cavity bottom -> surface far field) then matches the existing
DataFiles/Bo0.0010.dat exactly, so the same init() code fills liquid outside
the cavity and gas inside + above the surface.

Physical-neck mode (--waist W [--fillet RHO]). The delta family above rounds
the bridge with a fillet as wide as the neck (Rf ~ delta). Just after film
rupture the real bridge is much thinner: at waist W the liquid wedge between
a sphere touching the plane is only ~W^2/2 thick. This mode builds that
state directly. The sphere (radius 1) sits a gap g below the plane. A rim
circle of radius RHO is tangent to both, and its innermost point is the waist
W. RHO defaults to the tangent value (2 - W) - 2 sqrt(1 - W), where g = 0;
a larger RHO (e.g. a few cells, when the tangent rim is unresolved) sets
g > 0, i.e. a finite rupture gap. RHO below the tangent value would make the
sphere cross the plane and is refused.

Usage:
    python3 generate_Bo0_IC.py [--delta 0.01] [--rmax 32] [--out DataFiles/Bo0.0000.dat]
    python3 generate_Bo0_IC.py --waist 0.1 [--fillet RHO] [--rmax 32] --out ...
"""
import argparse
import numpy as np


def get_circles(delta, rmax, n=2000):
    X1c = -(1 + delta)
    phic1 = np.arcsin(2 * delta)
    phi1 = np.linspace(np.pi, phic1, n)
    X1 = X1c + np.cos(phi1)
    Y1 = np.sin(phi1)

    Yfc = (1 + delta) * np.tan(phic1)
    Rf = (1 + delta) / np.cos(phic1) - 1
    phif = np.linspace(np.pi / 2 - phic1, -np.pi / 2, n)
    Xf = -Rf * np.sin(phif)
    Yf = Yfc - Rf * np.cos(phif)

    # shift so the neck (fillet end) sits at x=0, and mirror
    X1 = -(X1 - Xf[-1])
    Xf = -(Xf - Xf[-1])

    X2 = Xf[-1] * np.ones(n)          # flat free surface at x=0
    Y2 = np.linspace(Yf[-1], rmax, n)  # extend radially past the domain

    X = np.concatenate([X1, Xf, X2])
    Y = np.concatenate([Y1, Yf, Y2])
    return X, Y, Rf


def tangent_fillet(waist):
    """Rim radius at which a rim of waist `waist` touches a sphere tangent to
    the plane (gap 0)."""
    # (2 - W) - 2 sqrt(1 - W) = (1 - sqrt(1 - W))^2, written without the
    # cancellation that loses small waists
    return (waist / (1. + np.sqrt(1. - waist)))**2


def get_physical_neck(waist, rho, rmax, n=2000):
    """Sphere-plane bridge with waist `waist` and rim radius `rho`.

    Built in the same frame as get_circles (cavity in +x, surface at x=0,
    ordered cavity bottom -> rim -> flat surface). Sphere centre (1+g, 0);
    rim centre (rho, waist+rho), tangent to the plane at (0, waist+rho) and
    externally tangent to the sphere. Returns X, Y and the gap g.
    """
    g = np.sqrt((1. + rho)**2 - (waist + rho)**2) - 1. + rho
    if g < -1e-12:
        raise ValueError("fillet %g is below the tangent value %g for waist %g"
                         % (rho, tangent_fillet(waist), waist))
    g = max(g, 0.)
    cx, fx, fy = 1. + g, rho, waist + rho
    if fy >= rmax:
        raise ValueError("rim extent %g reaches rmax %g" % (fy, rmax))
    # sphere: far pole (phi = 0, on the axis) to the rim tangent point
    phiT = np.arctan2(fy, fx - cx)
    phi = np.linspace(0., phiT, n)
    X1 = cx + np.cos(phi)
    Y1 = np.sin(phi)
    # rim: from the sphere tangent point, round the waist, to the plane
    th1 = np.arctan2(-fy, cx - fx)
    th = np.linspace(th1, -np.pi, n)
    Xf = fx + rho * np.cos(th)
    Yf = fy + rho * np.sin(th)
    X2 = np.zeros(n)
    Y2 = np.linspace(fy, rmax, n)
    X = np.concatenate([X1, Xf, X2])
    Y = np.concatenate([Y1, Yf, Y2])
    return X, Y, g


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--delta", type=float, default=0.01,
                    help="fillet/neck regularisation scale (initial hole radius ~2*delta)")
    ap.add_argument("--rmax", type=float, default=32.0,
                    help="radial extent of the flat free surface (>= domain size)")
    ap.add_argument("--waist", type=float, default=None,
                    help="physical-neck mode: waist radius of the opened bridge")
    ap.add_argument("--fillet", type=float, default=None,
                    help="physical-neck mode: rim radius (default: tangent value, gap 0)")
    ap.add_argument("--out", default="simulationCases/DataFiles/Bo0.0000.dat")
    a = ap.parse_args()

    if a.waist is None:
        if a.fillet is not None:
            ap.error("--fillet needs --waist")
        X, Y, Rf = get_circles(a.delta, a.rmax)
    else:
        if not (np.isfinite(a.waist) and 0. < a.waist < 1.):
            ap.error("--waist must lie in (0, 1)")
        rho = tangent_fillet(a.waist) if a.fillet is None else a.fillet
        if not (np.isfinite(rho) and rho > 0.):
            ap.error("--fillet must be a positive finite number")
        try:
            X, Y, gap = get_physical_neck(a.waist, rho, a.rmax)
        except ValueError as err:
            ap.error(str(err))
        Rf = rho
        print("physical neck: waist %g, rim radius %.6g (tangent %.6g), gap %.6g"
              % (a.waist, rho, tangent_fillet(a.waist), gap))
    # flip axial sign to the drill-solver convention (cavity in -x, south pole
    # near -2, surface at x=0); keeps the cavity-bottom -> surface ordering.
    Xbb = -X
    Ybb = np.maximum(Y, 1e-8)   # avoid exact 0 on the axis (matches Bo*.dat)

    with open(a.out, "w") as fh:
        for x, y in zip(Xbb, Ybb):
            fh.write("% .7e   % .7e\n" % (x, y))
    if a.waist is None:
        print("wrote %s  (delta=%g, Rf=%.5f, %d points)" % (a.out, a.delta, Rf, len(Xbb)))
    else:
        print("wrote %s  (waist=%g, rim=%.6g, %d points)" % (a.out, a.waist, Rf, len(Xbb)))
    print("  x_axial [%.4f, %.4f]  (south pole %.4f, surface 0)" % (Xbb.min(), Xbb.max(), Xbb.min()))
    print("  y_radial[%.4f, %.4f]" % (Ybb.min(), Ybb.max()))


if __name__ == "__main__":
    main()
