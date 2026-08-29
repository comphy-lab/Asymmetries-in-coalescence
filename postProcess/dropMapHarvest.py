#!/usr/bin/env python3
"""Harvest the drop-map campaign logs into r_d(Oh), U_d(Oh) and jet metrics.

Inputs, per case directory
--------------------------
``components.log``
    ``t,component,is_main,cells,volume,r_eq,x_mean,u_x_mean,x_min,x_max``
    written every TLOG = 0.005. The component index comes from Basilisk's
    ``tag()`` and is NOT stable between frames, so drops must be tracked by
    continuity, not by index. That tracking is this module's core job.
``jet.log``
    ``t,z_tip,u_tip_cap,r_cap_max,vol_cap[,r_jet_z0],ke,n_components,
    r_largest_detached,u_largest_detached,x_largest_detached``.
    The ``r_jet_z0`` column exists only in runs made with
    ``burstingBubbleInfiniteRr.c``; the original drop-map wrapper wrote ten
    columns. Both layouts are parsed.

Pre-registered drop selection (frozen 2026-08-29, BEFORE the ladder landed)
---------------------------------------------------------------------------
A tracked component is the *first visible forward drop* when ALL hold:

1. it is not the main component in any frame of its life;
2. its equivalent radius satisfies ``r_eq >= R_VISIBLE`` (0.021005127, the
   pinned literature threshold) on ``PERSISTENCE`` (3) consecutive frames;
3. its volume-weighted axial velocity is positive on the frame the
   persistence count completes.

The *first topological pinch* is simply the first frame at which any
non-main component exists at all. Both definitions are reported for every
case; conflating them manufactures spurious structure in r_d(Oh).

Tracking
--------
Greedy frame-to-frame matching on volume conservation and centroid
continuity: candidate pairs are scored by
``|dV|/V + |dx|/max(dx_max, r_eq)`` and accepted below a cap. A drop that
merges back into the main body simply ends its track (that is physical:
end-pinch fragments can be recaptured).

Usage
-----
    python3 dropMapHarvest.py CASEDIR [CASEDIR ...] --oh 0.01 0.015 ...
    python3 dropMapHarvest.py --ladder ladder.csv OUTDIR

Outputs one row per case in ``ladder_summary.csv`` plus per-case
``drops.csv`` and ``jet_metrics.csv``. No fitting happens here; fits live in
``fitDropMapLadder.py`` so the measured points exist independently of any
model they are later compared against.
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from dataclasses import dataclass, field
from pathlib import Path

R_VISIBLE = 0.021005127     # pinned literature threshold, never mesh-derived
PERSISTENCE = 3             # consecutive frames >= R_VISIBLE
MATCH_COST_CAP = 0.75       # max volume+centroid mismatch to link two frames
FRAME_TOL = 1e-6            # frame-time comparison tolerance


# --------------------------------------------------------------------------
# parsing
# --------------------------------------------------------------------------

@dataclass
class Component:
    t: float
    index: int
    is_main: bool
    cells: float
    volume: float
    r_eq: float
    x_mean: float
    u_x_mean: float
    x_min: float
    x_max: float


def read_components(path: Path) -> list[list[Component]]:
    """Group components.log rows into frames (one list per time level).

    The solver APPENDS, so a restarted case replays times already logged
    (the restart dump precedes the crash point). Later rows supersede
    earlier ones at the same time level: the restart recomputation wins.
    """
    by_time: dict[float, list[Component]] = {}
    order: list[float] = []
    with open(path, newline="") as fh:
        reader = csv.reader(fh)
        header = next(reader, None)
        if header is None or header[0] != "t":
            raise ValueError(f"{path}: unexpected header {header!r}")
        prev_key = None
        for row in reader:
            if len(row) != 10:
                continue  # tolerate a torn final line from a live run
            c = Component(
                t=float(row[0]), index=int(row[1]), is_main=row[2] == "1",
                cells=float(row[3]), volume=float(row[4]), r_eq=float(row[5]),
                x_mean=float(row[6]), u_x_mean=float(row[7]),
                x_min=float(row[8]), x_max=float(row[9]))
            key = round(c.t / FRAME_TOL) * FRAME_TOL
            if key != prev_key:
                if order and key <= order[-1]:
                    # time regressed: a restart is replaying from its dump
                    # point, so every previously logged frame from this time
                    # onward is the pre-crash pass and is superseded
                    order = [k for k in order if k < key]
                    by_time = {k: v for k, v in by_time.items() if k < key}
                order.append(key)
                by_time[key] = []
                prev_key = key
            elif by_time[key] and c.index <= by_time[key][-1].index:
                # same time, but component indices restart: a replay whose
                # first frame is exactly the last logged one. The recomputed
                # pass replaces the pre-crash rows wholesale.
                by_time[key] = []
            by_time[key].append(c)
    return [by_time[k] for k in order]


def read_jet(path: Path) -> list[dict]:
    """Parse jet.log, handling both the 10- and 11-column layouts.

    Same restart-replay rule as ``read_components``: jet.log is opened in
    append mode, so a restart concatenates a recomputed series onto the
    pre-crash one. On a time regression (or an exact repeat of an existing
    time) the recomputed rows supersede everything from that time onward.
    """
    rows: list[dict] = []
    with open(path, newline="") as fh:
        reader = csv.reader(fh)
        header = next(reader, None)
        if header is None:
            return rows
        names = [h.strip() for h in header]
        for row in reader:
            if len(row) != len(names):
                continue
            try:
                parsed = {k: float(v) for k, v in zip(names, row)}
            except ValueError:
                continue
            t = parsed["t"]
            while rows and rows[-1]["t"] >= t - FRAME_TOL:
                rows.pop()
            rows.append(parsed)
    return rows


# --------------------------------------------------------------------------
# tracking
# --------------------------------------------------------------------------

@dataclass
class Track:
    """One liquid component followed through time."""
    track_id: int
    history: list[Component] = field(default_factory=list)
    visible_run: int = 0            # consecutive frames >= R_VISIBLE
    confirmed_at: Component | None = None  # frame completing persistence
    ever_main: bool = False

    @property
    def born(self) -> Component:
        return self.history[0]

    @property
    def last(self) -> Component:
        return self.history[-1]

    def feed(self, c: Component) -> None:
        self.history.append(c)
        self.ever_main = self.ever_main or c.is_main
        if not c.is_main and c.r_eq >= R_VISIBLE:
            self.visible_run += 1
            if self.visible_run >= PERSISTENCE and self.confirmed_at is None:
                self.confirmed_at = c
        else:
            self.visible_run = 0


def _match_cost(prev: Component, cur: Component, dt_frames: float) -> float:
    dv = abs(cur.volume - prev.volume) / max(prev.volume, 1e-30)
    # allowed drift: a generous multiple of the drop's own size plus its
    # ballistic displacement estimated from the previous frame's velocity
    drift = abs(cur.x_mean - (prev.x_mean + prev.u_x_mean * dt_frames))
    dx = drift / max(prev.r_eq, 1e-3)
    return dv + 0.25 * dx


def track_components(frames: list[list[Component]]) -> list[Track]:
    """Greedy continuity tracking. The main component is tracked too (it is
    the reference against which 'detached' is defined) but is flagged."""
    tracks: list[Track] = []
    live: dict[int, Track] = {}
    next_id = 0
    prev_frame: list[Component] = []
    prev_assign: dict[int, int] = {}  # position in prev frame -> track id

    for frame in frames:
        dt_frames = (frame[0].t - prev_frame[0].t) if prev_frame else 0.0
        # score all pairs
        pairs = []
        for pi, pc in enumerate(prev_frame):
            for ci, cc in enumerate(frame):
                cost = _match_cost(pc, cc, dt_frames)
                if cost <= MATCH_COST_CAP:
                    pairs.append((cost, pi, ci))
        pairs.sort()
        used_prev: set[int] = set()
        used_cur: set[int] = set()
        assign: dict[int, int] = {}
        for cost, pi, ci in pairs:
            if pi in used_prev or ci in used_cur:
                continue
            used_prev.add(pi)
            used_cur.add(ci)
            tid = prev_assign[pi]
            assign[ci] = tid
        # unmatched current components are new tracks
        for ci, cc in enumerate(frame):
            if ci not in assign:
                tr = Track(track_id=next_id)
                next_id += 1
                tracks.append(tr)
                live[tr.track_id] = tr
                assign[ci] = tr.track_id
        for ci, cc in enumerate(frame):
            live[assign[ci]].feed(cc)
        prev_frame = frame
        prev_assign = {ci: tid for ci, tid in assign.items()}
    return tracks


# --------------------------------------------------------------------------
# per-case measurements
# --------------------------------------------------------------------------

def first_topological_pinch(frames: list[list[Component]]) -> dict | None:
    """First frame holding any non-main component, and that component's
    radius: the 'first pinch' definition of drop size."""
    for frame in frames:
        detached = [c for c in frame if not c.is_main]
        if detached:
            biggest = max(detached, key=lambda c: c.volume)
            return {"t_pinch_topo": frame[0].t,
                    "r_pinch_topo": biggest.r_eq,
                    "u_pinch_topo": biggest.u_x_mean,
                    "cells_pinch_topo": biggest.cells}
    return None


def first_visible_drop(tracks: list[Track]) -> dict | None:
    """Earliest confirmed track under the pre-registered selector."""
    confirmed = [tr for tr in tracks
                 if tr.confirmed_at is not None and not tr.ever_main
                 and tr.confirmed_at.u_x_mean > 0.]
    if not confirmed:
        return None
    tr = min(confirmed, key=lambda tr: tr.confirmed_at.t)
    birth = tr.born
    conf = tr.confirmed_at
    return {"t_drop_first_seen": birth.t,
            "t_drop_confirmed": conf.t,
            "r_drop_at_birth": birth.r_eq,
            "r_drop_at_confirm": conf.r_eq,
            "u_drop_at_birth": birth.u_x_mean,
            "u_drop_at_confirm": conf.u_x_mean,
            "cells_at_confirm": conf.cells,
            "track_id": tr.track_id,
            "n_confirmed_drops": len(confirmed)}


def _interp_crossing(ts, ys, level=0.0):
    """First upward crossing of `level` in the (t, y) series; linear interp."""
    for (t0, y0), (t1, y1) in zip(zip(ts, ys), list(zip(ts, ys))[1:]):
        if y0 is None or y1 is None:
            continue
        if y0 < level <= y1:
            w = (level - y0) / (y1 - y0) if y1 != y0 else 0.0
            return t0 + w * (t1 - t0)
    return None


VOLCAP_MIN = 1e-6   # tip caps smaller than this are one or two cells: their
                    # momentum/volume ratio is numerically meaningless
Z0_CLEAR = 0.05     # the tip must be this far past x = 0 before r_jet_z0
                    # measures the STEM: at the crossing instant the tip cap
                    # itself sits inside the station slab, so the minimum
                    # interface radius there is the cap curvature, not the
                    # jet radius


def jet_metrics(rows: list[dict]) -> dict:
    """v_jet(t)-derived scalars.

    - t_cross:   jet tip crosses the undisturbed surface x = 0. Before the
                 jet exists, z_tip tracks the meniscus rim hovering at
                 x ~ 0, so a naive first-crossing search fires at t ~ 0.
                 The physical crossing is the first upward crossing AFTER
                 the global minimum of z_tip (the cavity-collapse
                 excursion), which is how it is defined here.
    - v_jet0:    tip-cap velocity at t_cross (the Gordillo &
                 Blanco-Rodríguez v_jet protocol). Interpolated between the
                 bracketing frames when both have a trustworthy cap
                 (vol_cap >= VOLCAP_MIN); otherwise the first trustworthy
                 post-crossing frame is used and flagged.
    - r_jet0:    the jet STEM radius at the x = 0 station: r_jet_z0 read at
                 the first trustworthy frame with z_tip >= Z0_CLEAR after the
                 crossing, once the tip cap has left the slab. Values logged
                 before the crossing are the collapsing rim, not the jet, and
                 the value AT the crossing is the tip-cap curvature; neither
                 is used. (Column exists only on burstingBubbleInfiniteRr.c
                 runs; None otherwise.)
    - v_max:     max tip-cap velocity over trustworthy frames after
                 t_cross, with the frame time
    """
    out: dict = {"t_cross": None, "v_jet0": None, "v_jet0_interpolated": None,
                 "r_jet0": None, "t_rjet0": None, "v_max": None,
                 "t_vmax": None}
    if not rows:
        return out

    def valid(v):
        return v is not None and math.isfinite(v) and abs(v) < 1e29

    ts = [r["t"] for r in rows]
    ztip = [r["z_tip"] if valid(r.get("z_tip")) else None for r in rows]
    vtip = [r["u_tip_cap"] if valid(r.get("u_tip_cap")) else None for r in rows]
    volc = [r.get("vol_cap") for r in rows]
    trusty = [v is not None and valid(vc) and vc >= VOLCAP_MIN
              for v, vc in zip(vtip, volc)]

    # collapse excursion: global minimum of the (valid) z_tip record
    zvals = [(z, i) for i, z in enumerate(ztip) if z is not None]
    if not zvals:
        return out
    i_zmin = min(zvals)[1]

    t_cross = _interp_crossing(ts[i_zmin:], ztip[i_zmin:], 0.0)
    out["t_cross"] = t_cross

    if t_cross is not None:
        for i in range(1, len(ts)):
            if ts[i - 1] <= t_cross <= ts[i]:
                w = ((t_cross - ts[i - 1]) / (ts[i] - ts[i - 1])
                     if ts[i] != ts[i - 1] else 0.0)
                if trusty[i - 1] and trusty[i]:
                    out["v_jet0"] = vtip[i - 1] + w * (vtip[i] - vtip[i - 1])
                    out["v_jet0_interpolated"] = True
                else:
                    # first trustworthy frame at or after the crossing
                    for k in range(i, len(ts)):
                        if trusty[k]:
                            out["v_jet0"] = vtip[k]
                            out["v_jet0_interpolated"] = False
                            break
                break
        post = [(v, tt) for tt, v, ok in zip(ts, vtip, trusty)
                if ok and tt >= t_cross]
        if post:
            out["v_max"], out["t_vmax"] = max(post)
        # stem radius: first trustworthy frame after the tip clears the slab
        for i, r in enumerate(rows):
            if (ztip[i] is not None and ztip[i] >= Z0_CLEAR and trusty[i]
                    and valid(r.get("r_jet_z0"))):
                out["r_jet0"] = r["r_jet_z0"]
                out["t_rjet0"] = ts[i]
                break
    return out


def load_case(case_dir: Path):
    """Parse and track once; harvest_case and write_case_details share it."""
    frames = read_components(case_dir / "components.log")
    tracks = track_components(frames)
    jrows = read_jet(case_dir / "jet.log")
    return frames, tracks, jrows


def harvest_case(case_dir: Path, oh: float, loaded=None) -> dict:
    frames, tracks, jrows = loaded if loaded is not None else load_case(case_dir)

    row: dict = {"case": case_dir.name, "oh": oh,
                 "t_last": frames[-1][0].t if frames else None,
                 "n_frames": len(frames), "n_tracks": len(tracks)}
    topo = first_topological_pinch(frames)
    row.update(topo or {"t_pinch_topo": None, "r_pinch_topo": None,
                        "u_pinch_topo": None, "cells_pinch_topo": None})
    vis = first_visible_drop(tracks)
    row.update(vis or {"t_drop_first_seen": None, "t_drop_confirmed": None,
                       "r_drop_at_birth": None, "r_drop_at_confirm": None,
                       "u_drop_at_birth": None, "u_drop_at_confirm": None,
                       "cells_at_confirm": None, "track_id": None,
                       "n_confirmed_drops": 0})
    row.update(jet_metrics(jrows))
    return row


def write_case_details(case_dir: Path, out_dir: Path, oh: float,
                       loaded=None) -> None:
    """Per-case drops.csv (every confirmed track) and jet_metrics.csv."""
    frames, tracks, jrows = loaded if loaded is not None else load_case(case_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    with open(out_dir / "drops.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["track_id", "t_first_seen", "t_confirmed",
                    "r_eq_birth", "r_eq_confirm", "u_x_birth", "u_x_confirm",
                    "cells_confirm", "n_frames", "ever_main"])
        for tr in tracks:
            if tr.confirmed_at is None:
                continue
            w.writerow([tr.track_id, tr.born.t, tr.confirmed_at.t,
                        tr.born.r_eq, tr.confirmed_at.r_eq,
                        tr.born.u_x_mean, tr.confirmed_at.u_x_mean,
                        tr.confirmed_at.cells, len(tr.history),
                        int(tr.ever_main)])
    with open(out_dir / "jet_metrics.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        m = jet_metrics(jrows)
        w.writerow(["oh"] + list(m.keys()))
        w.writerow([oh] + list(m.values()))


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------

def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("cases", nargs="+", help="case directories")
    ap.add_argument("--oh", nargs="+", type=float, required=True,
                    help="Oh value per case directory, same order")
    ap.add_argument("--out", default="ladder_summary.csv")
    ap.add_argument("--details", default=None,
                    help="also write per-case drops.csv under this directory")
    args = ap.parse_args(argv)
    if len(args.cases) != len(args.oh):
        ap.error("need one --oh per case directory")

    rows = []
    for case, oh in zip(args.cases, args.oh):
        case_dir = Path(case)
        try:
            loaded = load_case(case_dir)
        except FileNotFoundError as exc:
            print(f"SKIP {case_dir}: {exc}", file=sys.stderr)
            continue
        rows.append(harvest_case(case_dir, oh, loaded=loaded))
        if args.details:
            write_case_details(case_dir,
                               Path(args.details) / case_dir.name, oh,
                               loaded=loaded)

    rows.sort(key=lambda r: r["oh"])
    if rows:
        with open(args.out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)
        print(f"wrote {args.out} ({len(rows)} cases)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
