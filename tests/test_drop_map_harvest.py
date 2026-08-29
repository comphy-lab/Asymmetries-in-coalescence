"""Tests for postProcess/dropMapHarvest.py.

The decisive behaviours: tag() indices are NOT stable between frames, so
tracking must survive index permutation; the first topological pinch and the
first visible drop are different measurements; the jet-crossing interpolation
must use only the post-formation record.
"""

import csv
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "postProcess"))

import dropMapHarvest as h  # noqa: E402


def vol(r):
    return 4.0 / 3.0 * math.pi * r ** 3


def write_components(path, frames):
    """frames: list of (t, [(index, is_main, r_eq, x_mean, u_x)])."""
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["t", "component", "is_main", "cells", "volume", "r_eq",
                    "x_mean", "u_x_mean", "x_min", "x_max"])
        for t, comps in frames:
            for idx, main, r, x, u in comps:
                v = vol(r)
                w.writerow([t, idx, int(main), max(4, int(v / 1e-6)), v, r,
                            x, u, x - r, x + r])


def write_jet(path, rows, with_rjet=True):
    cols = ["t", "z_tip", "u_tip_cap", "r_cap_max", "vol_cap"]
    if with_rjet:
        cols.append("r_jet_z0")
    cols += ["ke", "n_components", "r_largest_detached",
             "u_largest_detached", "x_largest_detached"]
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(cols)
        for r in rows:
            w.writerow([r.get(c, 0.0) for c in cols])


def ladder_case(tmp_path, frames, jet_rows=None, with_rjet=True):
    d = tmp_path / "9001"
    d.mkdir()
    write_components(d / "components.log", frames)
    write_jet(d / "jet.log", jet_rows or [], with_rjet=with_rjet)
    return d


def test_tracking_survives_index_permutation(tmp_path):
    # main body always present; a drop of r=0.05 whose tag() index flips
    # between 2 and 1 across frames must still confirm as ONE track
    frames = []
    for k in range(5):
        t = 0.5 + 0.005 * k
        x = 1.0 + 0.05 * k
        if k % 2 == 0:
            comps = [(1, True, 0.8, -1.0, 0.1), (2, False, 0.05, x, 10.0)]
        else:
            comps = [(2, True, 0.8, -1.0, 0.1), (1, False, 0.05, x, 10.0)]
        frames.append((t, comps))
    d = ladder_case(tmp_path, frames)
    row = h.harvest_case(d, 0.03)
    assert row["n_confirmed_drops"] == 1
    assert abs(row["r_drop_at_confirm"] - 0.05) < 1e-12
    # the drop was visible from its first frame, so confirmation lands on
    # the PERSISTENCE-th frame
    assert abs(row["t_drop_confirmed"] - (0.5 + 0.005 * 2)) < 1e-9


def test_first_pinch_differs_from_first_visible(tmp_path):
    # a 3-cell fragment appears at t=0.50 (below threshold), the real drop
    # (r=0.06) appears at t=0.52: definitions must separate
    frames = [
        (0.500, [(1, True, 0.8, -1.0, 0.1), (2, False, 0.004, 0.9, 5.0)]),
        (0.505, [(1, True, 0.8, -1.0, 0.1)]),
        (0.510, [(1, True, 0.8, -1.0, 0.1)]),
        (0.520, [(1, True, 0.8, -1.0, 0.1), (2, False, 0.06, 1.0, 8.0)]),
        (0.525, [(1, True, 0.8, -1.0, 0.1), (2, False, 0.06, 1.05, 8.0)]),
        (0.530, [(1, True, 0.8, -1.0, 0.1), (2, False, 0.06, 1.10, 8.0)]),
    ]
    d = ladder_case(tmp_path, frames)
    row = h.harvest_case(d, 0.02)
    assert abs(row["t_pinch_topo"] - 0.500) < 1e-12          # fragment
    assert abs(row["r_pinch_topo"] - 0.004) < 1e-12
    assert abs(row["t_drop_first_seen"] - 0.520) < 1e-12     # real drop
    assert abs(row["r_drop_at_confirm"] - 0.06) < 1e-12


def test_backward_moving_blob_never_confirms(tmp_path):
    frames = [(0.5 + 0.005 * k,
               [(1, True, 0.8, -1.0, 0.1),
                (2, False, 0.05, 1.0 - 0.01 * k, -3.0)]) for k in range(5)]
    d = ladder_case(tmp_path, frames)
    row = h.harvest_case(d, 0.03)
    assert row["n_confirmed_drops"] == 0
    assert row["t_pinch_topo"] is not None  # it still pinched topologically


def test_jet_crossing_and_vmax(tmp_path):
    # tip rises through 0 between t=0.50 and 0.51; velocity peaks later
    # rim hovers at ~0 early (must NOT count as the crossing), collapse
    # takes z_tip to -1, then the jet rises through 0
    rows = [
        {"t": 0.01, "z_tip": +0.005, "u_tip_cap": -0.01, "vol_cap": 0.03,
         "ke": 1, "n_components": 1},
        {"t": 0.02, "z_tip": -0.005, "u_tip_cap": -0.02, "vol_cap": 0.03,
         "ke": 1, "n_components": 1},
        {"t": 0.03, "z_tip": +0.004, "u_tip_cap": -0.01, "vol_cap": 0.03,
         "ke": 1, "n_components": 1},
        {"t": 0.45, "z_tip": -1.00, "u_tip_cap": 2.0, "vol_cap": 0.01,
         "ke": 1, "n_components": 1},
        {"t": 0.49, "z_tip": -0.10, "u_tip_cap": 5.0, "vol_cap": 1e-4,
         "r_jet_z0": -1e30, "ke": 1, "n_components": 1},
        {"t": 0.50, "z_tip": -0.02, "u_tip_cap": 8.0, "vol_cap": 1e-4,
         "r_jet_z0": 0.05, "ke": 1, "n_components": 1},
        {"t": 0.51, "z_tip": +0.02, "u_tip_cap": 10.0, "vol_cap": 1e-4,
         "r_jet_z0": 0.03, "ke": 1, "n_components": 1},
        {"t": 0.52, "z_tip": +0.08, "u_tip_cap": 12.0, "vol_cap": 1e-4,
         "r_jet_z0": 0.04, "ke": 1, "n_components": 1},
        {"t": 0.53, "z_tip": +0.15, "u_tip_cap": 9.0, "vol_cap": 1e-4,
         "r_jet_z0": 0.05, "ke": 1, "n_components": 1},
    ]
    d = ladder_case(tmp_path, [(0.5, [(1, True, 0.8, -1.0, 0.1)])], rows)
    m = h.jet_metrics(h.read_jet(d / "jet.log"))
    assert abs(m["t_cross"] - 0.505) < 1e-9
    assert abs(m["v_jet0"] - 9.0) < 1e-9         # midway between 8 and 10
    assert abs(m["r_jet0"] - 0.04) < 1e-9        # midway between 0.05, 0.03
    assert abs(m["v_max"] - 12.0) < 1e-12
    assert abs(m["t_vmax"] - 0.52) < 1e-12


def test_ten_column_jet_log_parses(tmp_path):
    # enough rows for a real crossing, so r_jet0 is None because the COLUMN
    # is absent, not because t_cross was never found
    rows = [{"t": 0.45, "z_tip": -1.0, "u_tip_cap": 2.0, "vol_cap": 0.01,
             "ke": 1, "n_components": 1},
            {"t": 0.50, "z_tip": -0.02, "u_tip_cap": 8.0, "vol_cap": 1e-4,
             "ke": 1, "n_components": 1},
            {"t": 0.51, "z_tip": +0.02, "u_tip_cap": 10.0, "vol_cap": 1e-4,
             "ke": 1, "n_components": 1}]
    d = ladder_case(tmp_path, [(0.5, [(1, True, 0.8, -1.0, 0.1)])],
                    rows, with_rjet=False)
    parsed = h.read_jet(d / "jet.log")
    assert len(parsed) == 3 and "r_jet_z0" not in parsed[0]
    m = h.jet_metrics(parsed)
    assert m["t_cross"] is not None
    assert abs(m["v_jet0"] - 9.0) < 1e-9
    assert m["r_jet0"] is None  # old runs simply lack the column


def test_torn_final_line_is_tolerated(tmp_path):
    d = ladder_case(tmp_path, [(0.5, [(1, True, 0.8, -1.0, 0.1)])])
    with open(d / "components.log", "a") as fh:
        fh.write("0.505,1,1,120")  # live run cut mid-row
    frames = h.read_components(d / "components.log")
    assert len(frames) == 1


def test_rim_hover_does_not_fake_a_crossing(tmp_path):
    """Early z_tip wiggles around 0 (the meniscus rim) must not be taken as
    the jet crossing; measured on case 6301 the naive first crossing landed
    at t ~ 0.005 with u ~ -0.001."""
    rows = [
        {"t": 0.005, "z_tip": -0.0002, "u_tip_cap": 0.0, "vol_cap": 0.03,
         "ke": 1, "n_components": 1},
        {"t": 0.010, "z_tip": +0.0070, "u_tip_cap": -0.08, "vol_cap": 0.03,
         "ke": 1, "n_components": 1},
        {"t": 0.450, "z_tip": -1.0, "u_tip_cap": 2.0, "vol_cap": 0.01,
         "ke": 1, "n_components": 1},
        {"t": 0.500, "z_tip": -0.3, "u_tip_cap": 20.0, "vol_cap": 1e-4,
         "ke": 1, "n_components": 1},
        {"t": 0.505, "z_tip": +0.3, "u_tip_cap": 25.0, "vol_cap": 1e-4,
         "ke": 1, "n_components": 1},
    ]
    d = ladder_case(tmp_path, [(0.5, [(1, True, 0.8, -1.0, 0.1)])],
                    rows, with_rjet=False)
    m = h.jet_metrics(h.read_jet(d / "jet.log"))
    assert m["t_cross"] is not None and m["t_cross"] > 0.5
    assert m["v_jet0"] is not None and m["v_jet0"] > 20.0


def test_one_cell_cap_spike_is_not_vmax(tmp_path):
    """A one-cell tip cap (vol_cap below VOLCAP_MIN) must be excluded from
    both v_jet0 interpolation and v_max; measured on case 6304 a 5.9e-7 cap
    reported u = 570."""
    rows = [
        {"t": 0.45, "z_tip": -1.0, "u_tip_cap": 2.0, "vol_cap": 0.01,
         "ke": 1, "n_components": 1},
        {"t": 0.490, "z_tip": -0.56, "u_tip_cap": 570.0, "vol_cap": 5.9e-7,
         "ke": 1, "n_components": 1},
        {"t": 0.495, "z_tip": +0.26, "u_tip_cap": 164.0, "vol_cap": 2.4e-6,
         "ke": 1, "n_components": 1},
        {"t": 0.500, "z_tip": +1.01, "u_tip_cap": 136.0, "vol_cap": 2.5e-6,
         "ke": 1, "n_components": 1},
    ]
    d = ladder_case(tmp_path, [(0.5, [(1, True, 0.8, -1.0, 0.1)])],
                    rows, with_rjet=False)
    m = h.jet_metrics(h.read_jet(d / "jet.log"))
    assert m["v_jet0"] == 164.0            # nearest trustworthy frame
    assert m["v_jet0_interpolated"] is False
    assert m["v_max"] == 164.0             # the 570 spike is excluded


def test_restart_replay_rows_are_superseded(tmp_path):
    """A restarted case re-appends frames for times already logged; the
    replayed (recomputed) rows must replace the pre-crash ones."""
    d = tmp_path / "9001"
    d.mkdir()
    write_components(d / "components.log", [
        (0.500, [(1, True, 0.8, -1.0, 0.1)]),
        (0.505, [(1, True, 0.8, -1.0, 0.1), (2, False, 0.03, 1.0, 5.0)]),
        (0.510, [(1, True, 0.8, -1.0, 0.1), (2, False, 0.03, 1.1, 5.0)]),
    ])
    # crash after 0.510; restart replays from 0.505 with the recomputed
    # pass giving the drop as 0.05 and reaching only 0.505 so far
    with open(d / "components.log", "a") as fh:
        import csv as _csv
        w = _csv.writer(fh)
        w.writerow([0.505, 1, 1, 9000, vol(0.8), 0.8, -1.0, 0.1, -1.8, -0.2])
        w.writerow([0.505, 2, 0, 40, vol(0.05), 0.05, 1.0, 5.0, 0.95, 1.05])
    frames = h.read_components(d / "components.log")
    # the pre-crash 0.505 AND 0.510 frames are superseded by the replay
    assert [f[0].t for f in frames] == [0.500, 0.505]
    drops = [c for c in frames[-1] if not c.is_main]
    assert len(drops) == 1 and abs(drops[0].r_eq - 0.05) < 1e-12
