"""Recovery tests for fitDropMapLadder: synthetic ladders with known
parameters must round-trip through the pre-registered protocols."""

import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "postProcess"))

import fitDropMapLadder as fit  # noqa: E402


def synth_rows(oh_c_left=0.031, p_left=0.6, oh_c_right=0.024, p_right=1.2,
               oh_u=0.035, oh_cross=0.025):
    """Piecewise two-branch r(Oh): the falling power law up to oh_cross, the
    rising one beyond (the physical curve switches branch at the turn — the
    virtual zero of the falling branch is never attained). U_d falls through
    zero at oh_u."""
    ohs = [0.010, 0.015, 0.020, 0.024, 0.026, 0.028,
           0.030, 0.0315, 0.0325, 0.034]
    rows = []
    for oh in ohs:
        if oh <= oh_cross:
            r = 0.9 * (oh_c_left - oh) ** p_left
        else:
            r = 40.0 * (oh - oh_c_right) ** p_right
        ud = 5.0 * (oh_u - oh) if oh < oh_u else 0.0
        vmax = 30.0 - 4e4 * (oh - 0.028) ** 2   # peak at 0.028
        rows.append({"oh": oh, "r_drop_at_confirm": r,
                     "r_pinch_topo": r * 0.8, "r_jet0": r / 1.5,
                     "u_drop_at_confirm": ud, "v_max": vmax})
    return rows


def test_falling_branch_zero_recovered():
    rows = synth_rows()
    rep = fit.estimate(rows)
    f = rep["r_visible"]["falling_fit"]
    assert f["oh_c"] is not None
    assert abs(f["oh_c"] - 0.031) < 0.002
    assert abs(f["p"] - 0.6) < 0.1


def test_turning_point_protocols_agree():
    rows = synth_rows()
    rep = fit.estimate(rows)
    p1 = rep["r_visible"]["P1_parabola"].get("oh_opt")
    p2 = rep["P2_vmax"].get("oh_opt")
    p3 = rep["r_visible"].get("P3_branch_intersection_oh")
    for v in (p1, p2, p3):
        assert v is not None
        assert 0.024 < v < 0.032, v


def test_ohc_U_zero_recovered():
    rows = synth_rows(oh_u=0.035)
    rep = fit.estimate(rows)
    u = rep["ohc_U"]
    assert u["oh_c"] is not None
    assert abs(u["oh_c"] - 0.035) < 0.002
    assert abs(u["p"] - 1.0) < 0.15  # linear synthetic law


def test_missing_series_is_skipped_not_fatal():
    rows = [{"oh": o, "r_drop_at_confirm": r}
            for o, r in [(0.01, 0.1), (0.02, 0.05)]]  # only 2 points
    rep = fit.estimate(rows)
    assert "skipped" in rep["r_visible"]


def test_protocols_contain_no_external_comparators():
    """The estimation source must not embed manuscript or literature values
    (0.0325, 0.2845, 0.211 ...) — comparison happens after measurement."""
    src = (Path(__file__).resolve().parent.parent / "postProcess"
           / "fitDropMapLadder.py").read_text()
    for banned in ("0.0325", "0.2845", "0.211", "0.0326", "0.0398", "0.348"):
        assert banned not in src, f"external comparator {banned} in fitter"
