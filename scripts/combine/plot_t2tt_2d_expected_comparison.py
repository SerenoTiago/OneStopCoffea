#!/usr/bin/env python3
"""Plot the Run-3 and CMS-SUS-19-006 expected T2tt limit surfaces.

Only evaluated Run-3 mass points are drawn.  The Run-2 value at each coordinate
is read from the containing MassScan2DExp TH2 bin; no interpolation is used.
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D
import matplotlib.patheffects as path_effects
from matplotlib.ticker import FixedFormatter, FixedLocator, NullLocator
import numpy as np
import uproot

from collect_split_signal_limits import cross_section_fb
from plot_t2tt_hepdata_comparison import (
    centers_to_edges,
    match_official_points,
    output_stem,
    read_figure14_surface,
)


DEFAULT_LIMITS = Path(
    "analysis_products/combine/"
    "26-10-11/run3_2024_t2tt_grid_expanded_1007bkg_automcstats_expected_137fb/limits_expected.csv"
)
DEFAULT_XSEC_YAML_DIR = Path("analyzer_resources/datasets/2024/2stop_split")
DEFAULT_ROOT = Path(
    "data/SUS-19-006/audit_2026-09-10/"
    "CMS-SUS-19-006_Figure_014-a.root"
)
DEFAULT_OUTPUT = DEFAULT_LIMITS.parent / "t2tt_run3_vs_run2_expected_2d_with_contours"
EXISTING_RATIO_LIMITS = (0.6, 0.9)
# Matches the dotted mStop - mLSP = mt reference line in the official
# CMS-SUS-19-006 Figure 14a, which is drawn at 175 GeV.  Presentation only.
TOP_MASS_GEV = 175.0
# Contours compare the absolute limit with the physical theory cross section
# from the dataset YAMLs.  The CSV xsec_fb column is only the normalization the
# Combine inputs were produced with (0.7526 fb for the legacy 1250 GeV cards),
# so it is not used as the theory curve.


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--limits-csv",
        type=Path,
        default=DEFAULT_LIMITS,
        help="Expected-limit CSV produced by collect_split_signal_limits.py.",
    )
    parser.add_argument(
        "--xsec-yaml-dir",
        type=Path,
        default=DEFAULT_XSEC_YAML_DIR,
        help="Dataset YAML directory providing the theory cross section per stop mass.",
    )
    parser.add_argument(
        "--run2-root",
        type=Path,
        default=DEFAULT_ROOT,
        help="Official CMS-SUS-19-006 Figure 14a ROOT file.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help="Output stem, or a .png/.pdf path.",
    )
    parser.add_argument(
        "--write-table",
        type=Path,
        help="Comparison CSV path (default: <output>_values.csv).",
    )
    parser.add_argument("--limit-cmap", default="rainbow")
    parser.add_argument("--ratio-cmap", default="RdBu_r")
    parser.add_argument(
        "--no-annotations",
        action="store_true",
        help="Do not annotate finite cells with their numerical values.",
    )
    return parser.parse_args()


def read_run3_points(path: Path, xsec_yaml_dir: Path) -> list[dict[str, float]]:
    points: list[dict[str, float]] = []
    theory_by_mass: dict[float, float] = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        fields = set(reader.fieldnames or [])
        required = {"mstop", "mlsp", "xsec_fb", "exp_m1", "exp", "exp_p1"}
        if not required <= fields:
            raise ValueError(f"{path} lacks columns: {sorted(required - fields)}")
        direct_pb = "exp_limit_pb" in fields
        if not direct_pb and not {"exp", "xsec_fb"} <= fields:
            raise ValueError(
                f"{path} needs exp_limit_pb, or both exp and xsec_fb"
            )
        for row in reader:
            mstop = float(row["mstop"])
            mlsp = float(row["mlsp"])
            limit_pb = (
                float(row["exp_limit_pb"])
                if direct_pb
                else float(row["exp"]) * float(row["xsec_fb"]) / 1000.0
            )
            if not all(math.isfinite(value) for value in (mstop, mlsp, limit_pb)):
                raise ValueError(f"Non-finite Run-3 row: {row}")
            if limit_pb <= 0:
                raise ValueError(f"Non-positive Run-3 expected limit: {row}")
            quantiles = {
                "r_exp_m1": float(row["exp_m1"]),
                "r_exp": float(row["exp"]),
                "r_exp_p1": float(row["exp_p1"]),
            }
            if not all(math.isfinite(value) and value > 0 for value in quantiles.values()):
                raise ValueError(f"Invalid Run-3 expected quantile row: {row}")
            reference_xsec_fb = float(row["xsec_fb"])
            if mstop not in theory_by_mass:
                theory_by_mass[mstop] = cross_section_fb(xsec_yaml_dir, int(mstop))
            theory_xsec_fb = theory_by_mass[mstop]
            quantile_limits_pb = {
                key: (
                    float(row[column])
                    if column in fields and row[column]
                    else quantiles[r_key] * reference_xsec_fb / 1000.0
                )
                for key, column, r_key in (
                    ("contour_r_exp_m1", "exp_m1_limit_pb", "r_exp_m1"),
                    ("contour_r_exp", "exp_limit_pb", "r_exp"),
                    ("contour_r_exp_p1", "exp_p1_limit_pb", "r_exp_p1"),
                )
            }
            points.append(
                {
                    "mstop": mstop,
                    "mlsp": mlsp,
                    "limit_pb": limit_pb,
                    "reference_xsec_fb": reference_xsec_fb,
                    "theory_xsec_fb": theory_xsec_fb,
                    **quantiles,
                    **{
                        key: value / (theory_xsec_fb / 1000.0)
                        for key, value in quantile_limits_pb.items()
                    },
                }
            )
    if not points:
        raise ValueError(f"No Run-3 points found in {path}")
    keys = [(point["mstop"], point["mlsp"]) for point in points]
    if len(keys) != len(set(keys)):
        raise ValueError("Duplicate Run-3 mass coordinates in the limits CSV")
    return points


def format_limit(value: float) -> str:
    return f"{value:.2g}"


def ratio_scale(values: np.ndarray) -> tuple[LogNorm, list[float], bool]:
    """Keep the established scale when safe, otherwise center it on one."""
    old_min, old_max = EXISTING_RATIO_LIMITS
    log_width = math.log(old_max / old_min)
    comfortable_min = old_min * math.exp(0.05 * log_width)
    comfortable_max = old_max / math.exp(0.05 * log_width)
    value_min = float(values.min())
    value_max = float(values.max())
    if value_min >= comfortable_min and value_max <= comfortable_max:
        return LogNorm(vmin=old_min, vmax=old_max), [0.6, 0.7, 0.8, 0.9], False

    # A logarithmically symmetric range makes equal factors above and below
    # unity occupy equal color distances. Retain the old range and add 3%
    # headroom so evaluated cells cannot sit against or beyond the colorbar.
    half_span = max(
        abs(math.log(value_min)),
        abs(math.log(value_max)),
        abs(math.log(old_min)),
        abs(math.log(old_max)),
    ) * 1.03
    vmin = math.exp(-half_span)
    vmax = math.exp(half_span)
    candidates = np.array(
        [0.02, 0.03, 0.05, 0.07, 0.1, 0.15, 0.2, 0.3, 0.5, 0.7,
         1.0, 1.5, 2.0, 3.0, 5.0, 7.0, 10.0, 15.0, 20.0, 30.0, 50.0]
    )
    ticks = candidates[(candidates >= vmin) & (candidates <= vmax)].tolist()
    if len(ticks) > 9:
        indices = np.unique(np.rint(np.linspace(0, len(ticks) - 1, 9)).astype(int))
        ticks = [ticks[index] for index in indices]
        if vmin <= 1.0 <= vmax and 1.0 not in ticks:
            ticks.append(1.0)
            ticks.sort()
    return LogNorm(vmin=vmin, vmax=vmax), ticks, True


def annotation_color(cmap: matplotlib.colors.Colormap, norm: LogNorm, value: float) -> str:
    red, green, blue, _ = cmap(norm(value))
    luminance = 0.2126 * red + 0.7152 * green + 0.0722 * blue
    return "black" if luminance > 0.52 else "white"


def draw_grid(
    ax: plt.Axes,
    x_edges: np.ndarray,
    y_edges: np.ndarray,
    x_centers: np.ndarray,
    y_centers: np.ndarray,
    values: np.ndarray,
    *,
    norm: LogNorm,
    cmap_name: str,
    title: str,
    formatter,
    annotate: bool,
):
    cmap = plt.get_cmap(cmap_name).copy()
    cmap.set_bad("white")
    mesh = ax.pcolormesh(
        x_edges,
        y_edges,
        np.ma.masked_invalid(values),
        cmap=cmap,
        norm=norm,
        shading="flat",
        edgecolors="0.78",
        linewidth=0.35,
    )
    if annotate:
        for iy, y in enumerate(y_centers):
            for ix, x in enumerate(x_centers):
                value = values[iy, ix]
                if not np.isfinite(value):
                    continue
                ax.text(
                    0.5 * (x_edges[ix] + x_edges[ix + 1]),
                    0.5 * (y_edges[iy] + y_edges[iy + 1]),
                    formatter(float(value)),
                    ha="center",
                    va="center",
                    fontsize=5.2,
                    color=annotation_color(cmap, norm, float(value)),
                )
    ax.set_title(title, fontsize=11)
    ax.set_xlabel(r"$m_{\tilde{t}}$ [GeV]")
    ax.set_xlim(x_edges[0], x_edges[-1])
    ax.set_ylim(y_edges[0], y_edges[-1])
    ax.set_xticks(x_centers)
    ax.set_xticklabels([f"{value:g}" for value in x_centers], rotation=45, ha="right")
    ax.get_xticklabels()[0].set_ha("left")
    ax.set_yticks(y_centers)
    ax.set_yticklabels([f"{value:g}" for value in y_centers])
    ax.tick_params(direction="in", top=True, right=True, labelsize=8)
    return mesh


def write_table(path: Path, rows: list[dict[str, float]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "mstop",
        "mlsp",
        "run3_expected_limit_pb",
        "run2_expected_limit_pb",
        "run3_over_run2",
        "run2_bin_center_mstop",
        "run2_bin_center_mlsp",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in sorted(rows, key=lambda item: (item["mstop"], item["mlsp"])):
            writer.writerow(
                {
                    "mstop": f"{row['mstop']:g}",
                    "mlsp": f"{row['mlsp']:g}",
                    "run3_expected_limit_pb": f"{row['run3_limit_pb']:.10g}",
                    "run2_expected_limit_pb": f"{row['run2_hepdata_limit_pb']:.10g}",
                    "run3_over_run2": f"{row['ratio']:.10g}",
                    "run2_bin_center_mstop": f"{row['run2_hepdata_mstop']:.10g}",
                    "run2_bin_center_mlsp": f"{row['run2_hepdata_mlsp']:.10g}",
                }
            )


def read_tgraph(path: Path, name: str) -> tuple[np.ndarray, np.ndarray]:
    with uproot.open(path) as root_file:
        graph = root_file[name]
        if not graph.classname.startswith("TGraph"):
            raise TypeError(f"{path}:{name} is {graph.classname}, expected TGraph")
        x = np.asarray(graph.member("fX"), dtype=float)
        y = np.asarray(graph.member("fY"), dtype=float)
    finite = np.isfinite(x) & np.isfinite(y)
    if not np.any(finite):
        raise ValueError(f"{path}:{name} has no finite coordinates")
    return x[finite], y[finite]


def add_run3_contours(
    ax: plt.Axes,
    x_centers: np.ndarray,
    y_centers: np.ndarray,
    grids: list[tuple[str, np.ndarray, str]],
    *,
    color: str = "black",
    label_prefix: str = "",
    legend_handles: list[Line2D] | None = None,
) -> list[str]:
    """Draw r=1 with piecewise-linear interpolation and no extrapolation.

    With ``legend_handles`` the handles are appended there for a combined
    legend instead of drawing this function's own legend.
    """
    absent: list[str] = []
    for label, values, linestyle in grids:
        finite = values[np.isfinite(values)]
        if not (len(finite) and float(finite.min()) <= 1.0 <= float(finite.max())):
            absent.append(label)
            continue
        contour = ax.contour(
            x_centers,
            y_centers,
            np.ma.masked_invalid(values),
            levels=[1.0],
            colors=color,
            linewidths=1.8,
            linestyles=linestyle,
            corner_mask=False,
            zorder=5,
        )
        contour.set_path_effects(
            [path_effects.Stroke(linewidth=3.2, foreground="white"), path_effects.Normal()]
        )
    handles = [
        Line2D([], [], color=color, linewidth=1.8, linestyle="-", label=f"{label_prefix}Median expected"),
        Line2D(
            [], [],
            color=color,
            linewidth=1.8,
            linestyle="--",
            label=label_prefix + (
                r"Expected $\pm 1\sigma$ ($+1\sigma$ outside grid)"
                if r"Expected $+1\sigma$" in absent
                else r"Expected $\pm 1\sigma$"
            ),
        ),
    ]
    if legend_handles is None:
        ax.legend(handles=handles, loc="upper right", fontsize=7.5, framealpha=0.9)
    else:
        legend_handles.extend(handles)
    return absent


def add_top_mass_line(ax: plt.Axes, *, label: bool = True) -> None:
    """Draw the on-shell top boundary mLSP = mStop - mTop."""
    x_min, x_max = ax.get_xlim()
    y_min, y_max = ax.get_ylim()
    line_x_min = max(x_min, y_min + TOP_MASS_GEV)
    line_x_max = min(x_max, y_max + TOP_MASS_GEV)
    if line_x_min >= line_x_max:
        return
    x = np.array([line_x_min, line_x_max])
    y = x - TOP_MASS_GEV
    ax.plot(
        x,
        y,
        color="0.35",
        linewidth=1.3,
        linestyle="--",
        zorder=4,
    )
    if not label:
        return
    label_x = line_x_min + 0.24 * (line_x_max - line_x_min)
    label_y = label_x - TOP_MASS_GEV
    start = ax.transData.transform((line_x_min, line_x_min - TOP_MASS_GEV))
    end = ax.transData.transform((line_x_max, line_x_max - TOP_MASS_GEV))
    angle = math.degrees(math.atan2(end[1] - start[1], end[0] - start[0]))
    text = ax.text(
        label_x,
        label_y,
        r"$m_{\tilde{t}}-m_{\tilde{\chi}^{0}_{1}}=m_{t}$",
        color="0.3",
        fontsize=7.5,
        rotation=angle,
        rotation_mode="anchor",
        ha="center",
        va="bottom",
        zorder=6,
    )
    text.set_path_effects(
        [path_effects.Stroke(linewidth=2.2, foreground="white"), path_effects.Normal()]
    )


def add_run2_contours(
    ax: plt.Axes,
    root_path: Path,
    *,
    color: str = "black",
    label_prefix: str = "",
    legend_handles: list[Line2D] | None = None,
) -> dict[str, tuple[int, bool]]:
    graph_specs = (
        ("ExpLim", "-"),
        ("ExpLimSdn", "--"),
        ("ExpLimSup", "--"),
    )
    statuses: dict[str, tuple[int, bool]] = {}
    x_min, x_max = ax.get_xlim()
    y_min, y_max = ax.get_ylim()
    for name, linestyle in graph_specs:
        x, y = read_tgraph(root_path, name)
        visible = bool(
            x.max() >= x_min
            and x.min() <= x_max
            and y.max() >= y_min
            and y.min() <= y_max
        )
        statuses[name] = (len(x), visible)
        line, = ax.plot(
            x,
            y,
            color=color,
            linewidth=1.8,
            linestyle=linestyle,
            zorder=5,
        )
        line.set_path_effects(
            [path_effects.Stroke(linewidth=3.2, foreground="white"), path_effects.Normal()]
        )
    handles = [
        Line2D(
            [], [], color="0.45" if not statuses["ExpLim"][1] else color,
            linewidth=1.8, linestyle="-",
            label=label_prefix + (
                "Median expected (left of range)"
                if not statuses["ExpLim"][1]
                else "Median expected"
            ),
        ),
        Line2D(
            [], [], color=color, linewidth=1.8, linestyle="--",
            label=label_prefix + (
                r"Expected $\pm 1\sigma$ ($+1\sigma$ left of range)"
                if not statuses["ExpLimSup"][1]
                else r"Expected $\pm 1\sigma$"
            ),
        ),
    ]
    if legend_handles is None:
        ax.legend(handles=handles, loc="upper right", fontsize=7.5, framealpha=0.9)
    else:
        legend_handles.extend(handles)
    return statuses


def main() -> int:
    args = parse_args()
    run3_points = read_run3_points(args.limits_csv, args.xsec_yaml_dir)
    surface = read_figure14_surface(args.run2_root, "expected")
    matched, unmatched, mode = match_official_points(run3_points, surface)
    if unmatched:
        raise ValueError(f"Run-3 points without a usable Run-2 bin: {unmatched}")
    if len(matched) != len(run3_points):
        raise ValueError(
            f"Matched {len(matched)} of {len(run3_points)} Run-3 points"
        )

    x_centers = np.array(sorted({row["mstop"] for row in matched}), dtype=float)
    y_centers = np.array(sorted({row["mlsp"] for row in matched}), dtype=float)
    x_edges = centers_to_edges(x_centers, 50.0)
    y_edges = centers_to_edges(y_centers, 50.0)
    y_edges[0] = max(0.0, y_edges[0])
    x_index = {value: index for index, value in enumerate(x_centers)}
    y_index = {value: index for index, value in enumerate(y_centers)}

    run3 = np.full((len(y_centers), len(x_centers)), np.nan)
    run2 = np.full_like(run3, np.nan)
    ratio = np.full_like(run3, np.nan)
    run3_r_m1 = np.full_like(run3, np.nan)
    run3_r_median = np.full_like(run3, np.nan)
    run3_r_p1 = np.full_like(run3, np.nan)
    for row in matched:
        ix = x_index[row["mstop"]]
        iy = y_index[row["mlsp"]]
        run3[iy, ix] = row["run3_limit_pb"]
        run2[iy, ix] = row["run2_hepdata_limit_pb"]
        ratio[iy, ix] = row["ratio"]
    for row in run3_points:
        ix = x_index[row["mstop"]]
        iy = y_index[row["mlsp"]]
        run3_r_m1[iy, ix] = row["contour_r_exp_m1"]
        run3_r_median[iy, ix] = row["contour_r_exp"]
        run3_r_p1[iy, ix] = row["contour_r_exp_p1"]

    # The absent local point must remain visibly masked even though Run 2 has a bin.
    missing_key = (1250.0, 1000.0)
    if missing_key[0] in x_index and missing_key[1] in y_index:
        ix, iy = x_index[missing_key[0]], y_index[missing_key[1]]
        if any(np.isfinite(grid[iy, ix]) for grid in (run3, run2, ratio)):
            raise ValueError("Known missing (1250,1000) point was unexpectedly filled")

    limit_values = np.concatenate((run3[np.isfinite(run3)], run2[np.isfinite(run2)]))
    limit_norm = LogNorm(vmin=float(limit_values.min()), vmax=float(limit_values.max()))

    ratio_values = ratio[np.isfinite(ratio)]
    ratio_norm, ratio_ticks, ratio_limits_changed = ratio_scale(ratio_values)

    fig = plt.figure(figsize=(26.0, 7.8), constrained_layout=True)
    grid = fig.add_gridspec(2, 3, height_ratios=[1.0, 0.055])
    axes = [fig.add_subplot(grid[0, index]) for index in range(3)]
    cax_limits = fig.add_subplot(grid[1, 0:2])
    cax_ratio = fig.add_subplot(grid[1, 2])
    for ax in axes[1:]:
        ax.sharex(axes[0])
        ax.sharey(axes[0])

    limit_mesh = draw_grid(
        axes[0], x_edges, y_edges, x_centers, y_centers, run3,
        norm=limit_norm,
        cmap_name=args.limit_cmap,
        title=(
            r"Run 3: 13.6 TeV, 137 fb$^{-1}$"
            "\nMedian expected 95% CL upper limit"
        ),
        formatter=format_limit,
        annotate=not args.no_annotations,
    )
    draw_grid(
        axes[1], x_edges, y_edges, x_centers, y_centers, run2,
        norm=limit_norm,
        cmap_name=args.limit_cmap,
        title=(
            r"Run 2: 13 TeV, 137 fb$^{-1}$"
            "\nCMS-SUS-19-006 median expected 95% CL upper limit"
        ),
        formatter=format_limit,
        annotate=not args.no_annotations,
    )
    ratio_mesh = draw_grid(
        axes[2], x_edges, y_edges, x_centers, y_centers, ratio,
        norm=ratio_norm,
        cmap_name=args.ratio_cmap,
        title="Run 3 expected/Run 2 expected",
        formatter=lambda value: f"{value:.2f}",
        annotate=not args.no_annotations,
    )
    for ax in axes:
        add_top_mass_line(ax)
    absent_run3 = add_run3_contours(
        axes[0],
        x_centers,
        y_centers,
        [
            (r"Expected $-1\sigma$", run3_r_m1, "--"),
            ("Median expected", run3_r_median, "-"),
            (r"Expected $+1\sigma$", run3_r_p1, "--"),
        ],
    )
    run2_graph_counts = add_run2_contours(axes[1], args.run2_root)
    axes[0].set_ylabel(r"$m_{\tilde{\chi}^{0}_{1}}$ [GeV]")
    for ax in axes[1:]:
        ax.tick_params(labelleft=False)

    fig.colorbar(limit_mesh, cax=cax_limits, orientation="horizontal").set_label(
        "Median expected 95% CL upper limit on cross section [pb]"
    )
    ratio_bar = fig.colorbar(ratio_mesh, cax=cax_ratio, orientation="horizontal")
    ratio_bar.set_label("Run 3 / Run 2 median expected limit")
    ratio_bar.ax.xaxis.set_major_locator(FixedLocator(ratio_ticks))
    ratio_bar.ax.xaxis.set_major_formatter(
        FixedFormatter([f"{value:g}" for value in ratio_ticks])
    )
    ratio_bar.ax.xaxis.set_minor_locator(NullLocator())

    stem, extensions = output_stem(args.output)
    stem.parent.mkdir(parents=True, exist_ok=True)
    for extension in extensions:
        fig.savefig(stem.with_suffix(f".{extension}"), dpi=220)
    plt.close(fig)

    table_path = args.write_table or stem.with_name(f"{stem.name}_values.csv")
    write_table(table_path, matched)
    print(f"Run-3 points: {len(run3_points)}")
    print(f"Run-2 matches: {len(matched)} ({mode}, no interpolation)")
    print(f"Masked grid cells: {int(np.isnan(run3).sum())}")
    print(
        f"Shared limit scale: {limit_norm.vmin:.6g} to {limit_norm.vmax:.6g} pb"
    )
    print(
        f"Evaluated ratio range: {ratio_values.min():.6g} to {ratio_values.max():.6g}"
    )
    print(f"Ratio color scale: {ratio_norm.vmin:.6g} to {ratio_norm.vmax:.6g}")
    print(f"Ratio color limits changed: {'yes' if ratio_limits_changed else 'no'}")
    print(
        "Run-3 contours: absolute exp_m1, exp, exp_p1 limits / theory cross section "
        "at r=1 (linear cell interpolation)"
    )
    renormalized = sorted(
        {
            (row["mstop"], row["reference_xsec_fb"], row["theory_xsec_fb"])
            for row in run3_points
            if not math.isclose(row["reference_xsec_fb"], row["theory_xsec_fb"])
        }
    )
    if renormalized:
        print(
            "Run-3 contours using a theory xsec that differs from the input "
            "normalization [fb]: "
            + ", ".join(
                f"mStop={mass:g}: {reference:g} -> {theory:g}"
                for mass, reference, theory in renormalized
            )
        )
    if absent_run3:
        print("Run-3 contours without an in-grid crossing: " + ", ".join(absent_run3))
    print(
        "Run-2 contour TGraphs: "
        + ", ".join(
            f"{name} ({count} points; {'visible' if visible else 'outside axes'})"
            for name, (count, visible) in run2_graph_counts.items()
        )
    )
    print(
        "Wrote plot: "
        + ", ".join(str(stem.with_suffix(f".{extension}")) for extension in extensions)
    )
    print(f"Wrote values: {table_path}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
