r"""
# Video Post-Processing for Bubble Coalescence

Render interface geometry, strain-rate invariant, and velocity magnitude from
Basilisk snapshots of asymmetric bubble coalescence. These visualizations
support analysis of capillary-wave focusing, jetting, and droplet pinch-off by
tracking interface evolution alongside a dissipation proxy and flow speed.

## Workflow

- `postProcess/getFacet`: extract PLIC interface segments
- `postProcess/getData-generic`: sample $\log_{10}(D^2)$ and $|u|$ on a grid
- `postProcess/getCOM`: compute axial center-of-mass position and velocity

## Usage

```bash
python3 postProcess/Video-generic.py --caseToProcess simulationCases/3000 --Rr 1.0
```

## Outputs

- `simulationCases/<case>/Video/*.png`: rendered frames
- `simulationCases/<case>/<case>_COMData.csv`: COM time series
- `simulationCases/<case>/<case>.mp4`: encoded movie (if enabled)

## Author

Vatsal Sanjay
vatsal.sanjay@comphy-lab.org
CoMPhy Lab, Durham University
Last updated: Jan 2026
"""

import argparse
import csv
import math
import multiprocessing as mp
import os
import subprocess as sp
from dataclasses import dataclass
from functools import partial
from datetime import datetime
from typing import Sequence, Tuple, Optional

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.ticker import StrMethodFormatter

# Use mathtext (parallel-safe, no external LaTeX subprocess)
matplotlib.rcParams["font.family"] = "serif"
matplotlib.rcParams["mathtext.fontset"] = "cm"  # Computer Modern style


@dataclass(frozen=True)
class DomainBounds:
    """
    Symmetry-aware domain description in cylindrical coordinates.

    The code expects r in [rmin, rmax] with rmin <= 0 to leverage the axis of
    symmetry; z spans freely between zmin and zmax. Packaging these values avoids
    argument soup when plotting or generating overlays.
    """

    rmin: float
    rmax: float
    zmin: float
    zmax: float


@dataclass(frozen=True)
class RuntimeConfig:
    """
    Run-time knobs collected from CLI arguments.

    Multiprocessing workers only need a single instance of this struct, making
    later CLI additions painless. The accessors keep mirror operations (e.g.,
    computing rmin) in one place.
    """

    cpus: int
    n_snapshots: int
    grids_per_r: int
    tsnap: float
    zmin: float
    zmax: float
    rmax: float
    rr: float  # Radius ratio for coalescence
    case_dir: str
    output_dir: str
    skip_video_encode: bool
    framerate: int
    output_fps: int

    @property
    def rmin(self) -> float:
        return -self.rmax

    @property
    def bounds(self) -> DomainBounds:
        return DomainBounds(self.rmin, self.rmax, self.zmin, self.zmax)


@dataclass(frozen=True)
class PlotStyle:
    """
    Single source of truth for plot-level choices.

    Matplotlib tweaks become traceable: alter colours, fonts, or geometry here
    and every rendered frame will stay consistent without touching plotting
    logic.
    """

    figure_size: Tuple[float, float] = (19.20, 10.80)
    tick_label_size: int = 20
    zero_axis_color: str = "grey"
    axis_color: str = "black"
    line_width: float = 2.0
    interface_color: str = "#00B2FF"
    com_color: str = "red"
    com_marker_size: float = 15.0
    colorbar_width: float = 0.03
    left_colorbar_offset: float = 0.04
    right_colorbar_offset: float = 0.01


@dataclass(frozen=True)
class SnapshotInfo:
    """
    Metadata for an input snapshot and its output image.

    Storing paths and the physical time together simplifies filename logic and
    ensures logging statements stay informative.
    """

    index: int
    time: float
    source: str
    target: str


@dataclass
class FieldData:
    """
    Structured holder around the grids returned by getData-generic.

    Attributes provide convenient min/max queries so the plotting routine does
    not need to recalculate extents or worry about array shapes.
    """

    R: np.ndarray
    Z: np.ndarray
    strain_rate: np.ndarray
    velocity: np.ndarray
    nz: int

    @property
    def radial_extent(self) -> Tuple[float, float]:
        return self.R.min(), self.R.max()

    @property
    def axial_extent(self) -> Tuple[float, float]:
        return self.Z.min(), self.Z.max()


@dataclass(frozen=True)
class COMData:
    """Center of mass position and velocity at a given time."""
    time: float
    z_com: float
    u_com: float


PLOT_STYLE = PlotStyle()


def log_status(message: str, *, level: str = "INFO") -> None:
    """Print timestamped status messages for long-running CLI workflows."""
    timestamp = datetime.now().strftime("%H:%M:%S")
    print(f"[{timestamp}] [{level}] {message}", flush=True)


def parse_arguments() -> RuntimeConfig:
    """
    Parse CLI arguments and package them into an immutable runtime config.

    The returned dataclass can be shared safely across multiprocessing
    workers. Domain bounds default to values derived from `Rr`, but each bound
    can be overridden explicitly from the command line.

    #### Returns

    - `RuntimeConfig`: Parsed execution settings for the full rendering run.

    #### Example

    ```python
    config = parse_arguments()
    print(config.cpus)
    print(config.bounds)
    ```
    """
    parser = argparse.ArgumentParser(description="Generate snapshot videos for coalescence.")
    parser.add_argument("--CPUs", type=int, default=4, help="Number of CPUs to use")
    parser.add_argument(
        "--nGFS", type=int, default=4000, help="Number of restart files to process"
    )
    parser.add_argument(
        "--GridsPerR", type=int, default=256, help="Number of grids per R"
    )
    parser.add_argument("--Rr", type=float, default=1.0, help="Radius ratio (affects domain bounds)")
    parser.add_argument("--ZMIN", type=float, default=None, help="Minimum Z value (default: -2.0)")
    parser.add_argument("--ZMAX", type=float, default=None, help="Maximum Z value (default: 5*Rr-2.0)")
    parser.add_argument("--RMAX", type=float, default=None, help="Maximum R value (default: 2.5*Rr)")
    parser.add_argument("--tsnap", type=float, default=0.01, help="Time snap")
    parser.add_argument(
        "--caseToProcess",
        type=str,
        default="simulationCases/3000",
        help="Case to process",
    )
    parser.add_argument(
        "--folderToSave", type=str, default=None, help="Folder to save (default: <case>/Video)"
    )
    parser.add_argument(
        "--skip-video-encode", action="store_true",
        help="Skip ffmpeg video encoding after frame generation"
    )
    parser.add_argument(
        "--framerate", type=int, default=90,
        help="Input framerate for ffmpeg (default: 90)"
    )
    parser.add_argument(
        "--output-fps", type=int, default=30,
        help="Output video framerate (default: 30)"
    )
    args = parser.parse_args()

    # Calculate domain bounds based on Rr if not explicitly provided
    rr = args.Rr
    zmin = args.ZMIN if args.ZMIN is not None else -2.0
    zmax = args.ZMAX if args.ZMAX is not None else 5.0 * rr - 2.0
    rmax = args.RMAX if args.RMAX is not None else 2.5 * rr

    # Default output directory
    output_dir = args.folderToSave if args.folderToSave else os.path.join(args.caseToProcess, "Video")

    return RuntimeConfig(
        cpus=args.CPUs,
        n_snapshots=args.nGFS,
        grids_per_r=args.GridsPerR,
        tsnap=args.tsnap,
        zmin=zmin,
        zmax=zmax,
        rmax=rmax,
        rr=rr,
        case_dir=args.caseToProcess,
        output_dir=output_dir,
        skip_video_encode=args.skip_video_encode,
        framerate=args.framerate,
        output_fps=args.output_fps,
    )


def ensure_directory(path: str) -> None:
    """
    Create an output directory if it does not exist.

    #### Args

    - `path`: Directory to create, including any missing parents.
    """
    if not os.path.isdir(path):
        os.makedirs(path, exist_ok=True)


def run_helper(command: Sequence[str]) -> Sequence[str]:
    """
    Run a helper executable and return its stderr as decoded lines.

    The compiled helpers deliberately emit their payload to stderr, so stdout is
    ignored (it is typically empty) and we bubble up informative stderr content.

    #### Args

    - `command`: Command vector passed to `subprocess`.

    #### Returns

    - `Sequence[str]`: Decoded stderr lines produced by the helper.

    #### Raises

    - `RuntimeError`: If the helper exits with a non-zero return code.
    """
    process = sp.Popen(command, stdout=sp.PIPE, stderr=sp.PIPE)
    _, stderr = process.communicate()
    if process.returncode != 0:
        raise RuntimeError(
            f"Command {' '.join(command)} failed with code {process.returncode}:\n"
            f"{stderr.decode('utf-8')}"
        )
    return stderr.decode("utf-8").split("\n")


def get_facets(filename: str):
    """
    Collect interface facets from `getFacet` and mirror them about the axis.

    Shells out to the compiled ``getFacet`` executable, which extracts the
    volume-of-fluid (VOF) interface as a sequence of line segments. Since
    the Basilisk simulation uses axisymmetric coordinates, only the r >= 0
    half is computed. This function mirrors each segment about r=0 to create
    the full visualization.

    #### Args

    - `filename`: Path to a Basilisk snapshot file.

    #### Returns

    - `list[tuple]`: Mirrored line segments stored as
      `((r1, z1), (r2, z2))` pairs.

    #### Notes

    - `getFacet` emits coordinates as `(z, r)` pairs.
    - This function swaps them to `(r, z)` to match the plotting pipeline.
    """
    temp2 = run_helper(["postProcess/getFacet", filename])
    segs = []
    skip = False
    if len(temp2) > 1e2:
        for n1 in range(len(temp2)):
            temp3 = temp2[n1].split(" ")
            if temp3 == [""]:
                skip = False
                continue
            if not skip and n1 + 1 < len(temp2):
                temp4 = temp2[n1 + 1].split(" ")
                r1, z1 = np.array([float(temp3[1]), float(temp3[0])])
                r2, z2 = np.array([float(temp4[1]), float(temp4[0])])
                segs.append(((r1, z1), (r2, z2)))
                segs.append(((-r1, z1), (-r2, z2)))
                skip = True
    return segs


def get_com(filename: str) -> Optional[COMData]:
    """
    Extract center-of-mass data from the `getCOM` helper.

    Shells out to the compiled ``getCOM`` executable, which computes the
    volume-weighted center of mass position and velocity.

    #### Args

    - `filename`: Path to a Basilisk snapshot file.

    #### Returns

    - `COMData | None`: Time, axial position, and axial velocity when
      extraction succeeds, otherwise `None`.
    """
    try:
        temp2 = run_helper(["postProcess/getCOM", filename])
        if len(temp2) > 0 and temp2[0]:
            parts = temp2[0].split()
            if len(parts) >= 3:
                return COMData(
                    time=float(parts[0]),
                    z_com=float(parts[1]),
                    u_com=float(parts[2])
                )
    except Exception as e:
        log_status(f"getCOM failed for {filename}: {e}", level="WARN")
    return None


def get_field(filename: str, zmin: float, zmax: float, rmax: float, nr: int) -> FieldData:
    """
    Sample structured field arrays for a single snapshot.

    Shells out to the compiled ``getData-generic`` executable, which samples
    the velocity and strain-rate fields on a structured grid. Returns a
    `FieldData` object with reshaped 2D arrays, abstracting away the flattened
    helper output.

    #### Args

    - `filename`: Path to a Basilisk snapshot file.
    - `zmin`: Minimum axial coordinate in the sampling window.
    - `zmax`: Maximum axial coordinate in the sampling window.
    - `rmax`: Maximum radial coordinate on the positive branch.
    - `nr`: Number of radial samples.

    #### Returns

    - `FieldData`: Reshaped coordinate, strain-rate, and velocity arrays.

    #### Raises

    - `ValueError`: If `nr` is invalid or helper output cannot be reshaped
      into a regular `(nz, nr)` grid.
    """
    if nr <= 0:
        raise ValueError(f"nr must be positive, got {nr}")

    temp2 = run_helper(
        [
            "postProcess/getData-generic",
            filename,
            str(zmin),
            str(0),
            str(zmax),
            str(rmax),
            str(nr),
        ]
    )
    Rtemp, Ztemp, D2temp, veltemp = [], [], [], []

    for n1 in range(len(temp2)):
        temp3 = temp2[n1].split(" ")
        if temp3 == [""]:
            continue
        Ztemp.append(float(temp3[0]))
        Rtemp.append(float(temp3[1]))
        D2temp.append(float(temp3[2]))
        veltemp.append(float(temp3[3]))

    R = np.asarray(Rtemp)
    Z = np.asarray(Ztemp)
    D2 = np.asarray(D2temp)
    vel = np.asarray(veltemp)

    if len(Z) == 0:
        raise ValueError(f"No field samples were returned for {filename}")
    if len(Z) % nr != 0:
        raise ValueError(
            f"Field sample count {len(Z)} is not divisible by nr={nr} for {filename}"
        )

    nz = int(len(Z) / nr)

    log_status(f"{os.path.basename(filename)}: nz = {nz}")

    R.resize((nz, nr))
    Z.resize((nz, nr))
    D2.resize((nz, nr))
    vel.resize((nz, nr))

    return FieldData(R=R, Z=Z, strain_rate=D2, velocity=vel, nz=nz)


def build_snapshot_info(index: int, config: RuntimeConfig) -> SnapshotInfo:
    """
    Construct file paths for a given timestep index.

    #### Args

    - `index`: Timestep index used with `tsnap` to recover physical time.
    - `config`: Shared runtime configuration.

    #### Returns

    - `SnapshotInfo`: Input and output paths plus the recovered snapshot time.
    """
    time = config.tsnap * index
    # Note: coalescence uses intermediate/ directly (no results/ subdir)
    source = os.path.join(config.case_dir, "intermediate", f"snapshot-{time:.4f}")
    target = os.path.join(config.output_dir, f"{int(time * 1000):08d}.png")
    return SnapshotInfo(index=index, time=time, source=source, target=target)


def draw_domain_outline(ax, bounds: DomainBounds, style: PlotStyle) -> None:
    """
    Outline computational domain and symmetry line.

    #### Args

    - `ax`: Matplotlib axes used for plotting.
    - `bounds`: Domain extents in radial and axial directions.
    - `style`: Plotting style parameters.
    """
    ax.plot(
        [0, 0],
        [bounds.zmin, bounds.zmax],
        "-.",
        color=style.zero_axis_color,
        linewidth=style.line_width,
    )
    ax.plot(
        [bounds.rmin, bounds.rmin],
        [bounds.zmin, bounds.zmax],
        "-",
        color=style.axis_color,
        linewidth=style.line_width,
    )
    ax.plot(
        [bounds.rmin, bounds.rmax],
        [bounds.zmin, bounds.zmin],
        "-",
        color=style.axis_color,
        linewidth=style.line_width,
    )
    ax.plot(
        [bounds.rmin, bounds.rmax],
        [bounds.zmax, bounds.zmax],
        "-",
        color=style.axis_color,
        linewidth=style.line_width,
    )
    ax.plot(
        [bounds.rmax, bounds.rmax],
        [bounds.zmin, bounds.zmax],
        "-",
        color=style.axis_color,
        linewidth=style.line_width,
    )


def add_colorbar(fig, ax, mappable, *, align: str, label: str, style: PlotStyle):
    """
    Attach a vertical colorbar on the requested side of the axis.

    Using manual axes lets us keep the main axis square while still showing two
    distinct colour scales.

    #### Args

    - `fig`: Figure hosting the axes.
    - `ax`: Primary axes used to anchor the colorbar.
    - `mappable`: Image or artist represented by the color scale.
    - `align`: Either `'left'` or `'right'`.
    - `label`: Colorbar label.
    - `style`: Plot styling parameters.

    #### Returns

    - `matplotlib.colorbar.Colorbar`: The constructed colorbar object.
    """
    l, b, w, h = ax.get_position().bounds
    if align == "left":
        position = [l - style.left_colorbar_offset, b, style.colorbar_width, h]
    else:
        position = [l + w + style.right_colorbar_offset, b, style.colorbar_width, h]
    cb_ax = fig.add_axes(position)
    colorbar = plt.colorbar(mappable, cax=cb_ax, orientation="vertical")
    colorbar.set_label(label, fontsize=style.tick_label_size, labelpad=5)
    colorbar.ax.tick_params(labelsize=style.tick_label_size)
    colorbar.ax.yaxis.set_major_formatter(StrMethodFormatter("{x:,.2f}"))
    if align == "left":
        colorbar.ax.yaxis.set_ticks_position("left")
        colorbar.ax.yaxis.set_label_position("left")
    return colorbar


def plot_snapshot(
    field_data: FieldData,
    facets,
    com_data: Optional[COMData],
    bounds: DomainBounds,
    snapshot: SnapshotInfo,
    style: PlotStyle,
) -> None:
    """
    Render and persist a single snapshot figure.

    All artist construction lives here so multiprocessing workers only need to
    fetch data and call this function.

    #### Args

    - `field_data`: Structured strain-rate and velocity arrays.
    - `facets`: Interface line segments returned by `get_facets()`.
    - `com_data`: Optional center-of-mass diagnostics.
    - `bounds`: Domain extents used to keep the plot square.
    - `snapshot`: Snapshot metadata, including output path.
    - `style`: Shared plotting style.
    """
    fig, ax = plt.subplots()
    fig.set_size_inches(*style.figure_size)

    draw_domain_outline(ax, bounds, style)
    line_segments = LineCollection(
        facets, linewidths=4, colors=style.interface_color, linestyle="solid"
    )
    ax.add_collection(line_segments)

    rminp, rmaxp = field_data.radial_extent
    zminp, zmaxp = field_data.axial_extent

    cntrl1 = ax.imshow(
        field_data.strain_rate,
        cmap="hot_r",
        interpolation="Bilinear",
        origin="lower",
        extent=[-rminp, -rmaxp, zminp, zmaxp],
        vmax=2.0,
        vmin=-2.0,
    )

    cntrl2 = ax.imshow(
        field_data.velocity,
        interpolation="Bilinear",
        cmap="Purples",
        origin="lower",
        extent=[rminp, rmaxp, zminp, zmaxp],
        vmax=1.0,
        vmin=0.0,
    )

    # Add COM marker if available
    if com_data is not None:
        ax.plot(0, com_data.z_com, 'o', color=style.com_color,
                markersize=style.com_marker_size, zorder=10)

    ax.set_aspect("equal")
    ax.set_xlim(bounds.rmin, bounds.rmax)
    ax.set_ylim(bounds.zmin, bounds.zmax)
    ax.set_title(f"$t/\\tau_0$ = {snapshot.time:4.3f}", fontsize=style.tick_label_size)
    ax.axis("off")

    add_colorbar(
        fig,
        ax,
        cntrl1,
        align="left",
        label=r"$\log_{10}\left(\boldsymbol{\mathcal{D}:\mathcal{D}}\right)$",
        style=style,
    )
    add_colorbar(
        fig,
        ax,
        cntrl2,
        align="right",
        label=r"$\|\boldsymbol{u}\|$",
        style=style,
    )

    plt.savefig(snapshot.target, bbox_inches="tight")
    plt.close(fig)


def process_timestep(index: int, config: RuntimeConfig, style: PlotStyle) -> Optional[COMData]:
    """
    Worker executed for every timestep index.

    Performs availability checks, loads helper outputs, and calls plot_snapshot.

    #### Args

    - `index`: Timestep index relative to `tsnap`.
    - `config`: Shared runtime configuration.
    - `style`: Shared plotting style.

    #### Returns

    - `COMData | None`: Extracted COM diagnostics, or `None` when the
      snapshot is missing.
    """
    snapshot = build_snapshot_info(index, config)
    if not os.path.exists(snapshot.source):
        log_status(f"Missing: {os.path.basename(snapshot.source)}", level="WARN")
        return None
    if os.path.exists(snapshot.target):
        log_status(f"Exists, skipping: {os.path.basename(snapshot.target)}")
        # Still extract COM data for existing frames
        return get_com(snapshot.source)

    # Show relative path: CaseNo/intermediate/filename
    src_parts = snapshot.source.split(os.sep)
    src_rel = os.sep.join(src_parts[-3:]) if len(src_parts) >= 3 else snapshot.source
    log_status(f"Processing {src_rel}")

    try:
        facets = get_facets(snapshot.source)
        com_data = get_com(snapshot.source)
        nr = max(1, math.ceil(config.grids_per_r * config.rmax))
        field_data = get_field(
            snapshot.source, config.zmin, config.zmax, config.rmax, nr
        )
        plot_snapshot(field_data, facets, com_data, config.bounds, snapshot, style)

        # Show relative path: CaseNo/Video/filename
        tgt_parts = snapshot.target.split(os.sep)
        tgt_rel = os.sep.join(tgt_parts[-3:]) if len(tgt_parts) >= 3 else snapshot.target
        log_status(f"Saved: {tgt_rel}")

        return com_data

    except Exception as err:
        log_status(
            f"Error at {src_rel} (t={snapshot.time:.4f}): {err}", level="ERROR"
        )
        raise


def write_com_data(com_data_list: list, config: RuntimeConfig) -> None:
    """
    Write collected COM data to CSV file.

    #### Args

    - `com_data_list`: Sequence of `COMData` entries, with optional `None`s.
    - `config`: Runtime configuration providing output paths.
    """
    # Extract case number from path
    case_no = os.path.basename(config.case_dir)
    output_path = os.path.join(config.case_dir, f"{case_no}_COMData.csv")

    # Filter out None entries and sort by time
    valid_data = [d for d in com_data_list if d is not None]
    valid_data.sort(key=lambda x: x.time)

    with open(output_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['time', 'zCOM', 'uCOM'])
        for data in valid_data:
            writer.writerow([f"{data.time:.4f}", f"{data.z_com:.6g}", f"{data.u_com:.6g}"])

    log_status(f"COM data saved: {output_path} ({len(valid_data)} entries)")


def encode_video(config: RuntimeConfig) -> None:
    """
    Run ffmpeg to stitch PNG frames into an MP4 video.

    The output video is saved in the case directory with the case number
    as filename (e.g., simulationCases/3000/3000.mp4).

    #### Args

    - `config`: Runtime configuration containing paths and encoding settings.

    #### Raises

    - `RuntimeError`: If `ffmpeg` fails to encode the final video.
    """
    # Extract case number from path
    case_no = os.path.basename(config.case_dir)

    # Output path: <case_dir>/3000.mp4
    output_path = os.path.join(config.case_dir, f"{case_no}.mp4")
    input_pattern = os.path.join(config.output_dir, "*.png")

    cmd = [
        "ffmpeg", "-y",
        "-framerate", str(config.framerate),
        "-pattern_type", "glob",
        "-i", input_pattern,
        "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2",
        "-c:v", "libx264",
        "-r", str(config.output_fps),
        "-pix_fmt", "yuv420p",
        output_path
    ]

    log_status(f"Encoding video: {output_path}")
    result = sp.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        log_status(f"ffmpeg error: {result.stderr}", level="ERROR")
        raise RuntimeError(f"ffmpeg failed with code {result.returncode}")
    log_status(f"Video saved: {output_path}")


def main():
    """
    Entry point used by the CLI and documentation tooling.

    #### Workflow

    1. Parse arguments into `RuntimeConfig`.
    2. Render all requested snapshots in parallel.
    3. Write COM diagnostics to CSV.
    4. Optionally encode frames into an MP4.
    """
    config = parse_arguments()
    ensure_directory(config.output_dir)

    log_status(f"Processing case: {config.case_dir}")
    log_status(f"Domain: R=[{config.rmin:.2f},{config.rmax:.2f}], Z=[{config.zmin:.2f},{config.zmax:.2f}]")
    log_status(f"Radius ratio (Rr): {config.rr}")

    with mp.Pool(processes=config.cpus) as pool:
        worker = partial(process_timestep, config=config, style=PLOT_STYLE)
        com_data_list = pool.map(worker, range(config.n_snapshots))

    # Write COM data to CSV
    write_com_data(com_data_list, config)

    # Encode video unless skipped
    if not config.skip_video_encode:
        encode_video(config)


if __name__ == "__main__":
    main()
