#!/usr/bin/env python3
"""Plot a Run-2-calibrated full-Run-3 T2tt expected-sensitivity projection.

Each Run-3 MC expected limit at the target luminosity is divided by the
point-by-point ratio R_137 = sigma95(Run 3 MC, 137/fb) / sigma95(Run 2
expected, 137/fb), where the Run-2 value is the containing CMS-SUS-19-006
Figure 14a MassScan2DExp bin.  Equivalently,
sigma95_calibrated = sigma95_Run2_137 * sigma95_Run3_L / sigma95_Run3_137.
This is a calibrated sensitivity projection, not a predicted Run-3 limit.
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
from matplotlib.ticker import FixedFormatter, FixedLocator, NullLocator
import numpy as np

from plot_t2tt_2d_expected_comparison import (
    DEFAULT_ROOT,
    DEFAULT_XSEC_YAML_DIR,
    add_run2_contours,
    add_run3_contours,
    add_top_mass_line,
    draw_grid,
    format_limit,
    ratio_scale,
    read_run3_points,
)
from plot_t2tt_hepdata_comparison import (
    centers_to_edges,
    match_official_points,
    output_stem,
    read_figure14_surface,
)


COMBINE_ROOT = Path("analysis_products/combine")
DEFAULT_REFERENCE_DIR = (
    COMBINE_ROOT / "26-10-11/run3_2024_t2tt_grid_expanded_1007bkg_automcstats_expected_137fb"
)
DEFAULT_TARGET_DIR = (
    COMBINE_ROOT / "26-10-11/run3_2024_t2tt_grid_expanded_1007bkg_automcstats_expected_283fb"
)
REFERENCE_VALUES_NAME = "t2tt_run3_vs_run2_expected_2d_cms_fig14a_style_values.csv"
# Run-2 red used for the Run-2 curves in the search-bin comparison plots.
RUN2_COLOR = "#b2182b"
QUANTILES = (("m1", "contour_r_exp_m1"), ("median", "contour_r_exp"), ("p1", "contour_r_exp_p1"))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--reference-limits",
        type=Path,
        default=DEFAULT_REFERENCE_DIR / "limits_expected.csv",
        help="Final 137/fb Run-3 expected limits defining R_137.",
    )
    parser.add_argument(
        "--reference-values",
        type=Path,
        default=DEFAULT_REFERENCE_DIR / REFERENCE_VALUES_NAME,
        help="Values CSV of the final 137/fb comparison; R_137 must reproduce it.",
    )
    parser.add_argument(
        "--target-limits",
        type=Path,
        default=DEFAULT_TARGET_DIR / "limits_expected.csv",
        help="Run-3 expected limits at the target luminosity.",
    )
    parser.add_argument("--run2-root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--xsec-yaml-dir", type=Path, default=DEFAULT_XSEC_YAML_DIR)
    parser.add_argument("--reference-luminosity-fb", type=float, default=137.0)
    parser.add_argument("--target-luminosity-fb", type=float, default=283.0)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_TARGET_DIR,
        help="Directory for the plots and values CSV.",
    )
    parser.add_argument("--limit-cmap", default="rainbow")
    # The factor is 1/R, so the comparison plot's RdBu_r is reversed here:
    # blue still means Run 3 more sensitive than Run 2.
    parser.add_argument("--factor-cmap", default="RdBu")
    parser.add_argument(
        "--no-annotations",
        action="store_true",
        help="Do not annotate finite cells with their numerical values.",
    )
    return parser.parse_args()


def keyed(points: list[dict[str, float]]) -> dict[tuple[float, float], dict[str, float]]:
    return {(point["mstop"], point["mlsp"]): point for point in points}


def check_reference_ratio(path: Path, ratio: dict[tuple[float, float], float]) -> None:
    """Require R_137 to reproduce the final 137/fb comparison values CSV."""
    with path.open(newline="") as handle:
        published = {
            (float(row["mstop"]), float(row["mlsp"])): float(row["run3_over_run2"])
            for row in csv.DictReader(handle)
        }
    if published.keys() != ratio.keys():
        raise ValueError(f"{path} covers different mass points than the reference limits")
    bad = [key for key, value in ratio.items() if not math.isclose(value, published[key], rel_tol=1e-8)]
    if bad:
        raise ValueError(f"R_137 differs from {path} at {bad[:5]}")


def quantile_limits_pb(point: dict[str, float]) -> dict[str, float]:
    theory_pb = point["theory_xsec_fb"] / 1000.0
    return {name: point[key] * theory_pb for name, key in QUANTILES}


def write_values(path: Path, rows: list[dict[str, float]]) -> None:
    fields = [
        "mstop",
        "mlsp",
        "sigma95_run3_137_pb",
        "sigma95_run2_expected_137_pb",
        "R_137",
        "calibration_factor",
        "sigma95_run3_283_pb",
        "sigma95_run3_283_m1_pb",
        "sigma95_run3_283_p1_pb",
        "sigma95_calibrated_283_pb",
        "sigma95_calibrated_283_m1_pb",
        "sigma95_calibrated_283_p1_pb",
        "theory_xsec_fb",
        "run3_283_over_theory",
        "run3_283_m1_over_theory",
        "run3_283_p1_over_theory",
        "calibrated_283_over_theory",
        "calibrated_283_m1_over_theory",
        "calibrated_283_p1_over_theory",
        "run2_bin_center_mstop",
        "run2_bin_center_mlsp",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in sorted(rows, key=lambda item: (item["mstop"], item["mlsp"])):
            writer.writerow(
                {
                    key: (f"{row[key]:g}" if key in ("mstop", "mlsp") else f"{row[key]:.10g}")
                    for key in fields
                }
            )


def to_grid(rows, key, x_index, y_index, shape) -> np.ndarray:
    grid = np.full(shape, np.nan)
    for row in rows:
        grid[y_index[row["mlsp"]], x_index[row["mstop"]]] = row[key]
    return grid


def summary(name: str, values: np.ndarray) -> str:
    return f"{name}: median {np.median(values):.4g}, range {values.min():.4g}-{values.max():.4g}"


def main() -> int:
    args = parse_args()
    lumi_ref = args.reference_luminosity_fb
    lumi = args.target_luminosity_fb

    reference = keyed(read_run3_points(args.reference_limits, args.xsec_yaml_dir))
    target = keyed(read_run3_points(args.target_limits, args.xsec_yaml_dir))
    if reference.keys() != target.keys():
        missing = sorted(set(reference) ^ set(target))
        raise ValueError(f"Reference and target grids differ at {missing[:10]}")

    surface = read_figure14_surface(args.run2_root, "expected")
    matched, unmatched, mode = match_official_points(list(reference.values()), surface)
    if unmatched or len(matched) != len(reference):
        raise ValueError(f"Run-3 points without a usable Run-2 bin: {unmatched}")
    ratio = {(row["mstop"], row["mlsp"]): row["ratio"] for row in matched}
    check_reference_ratio(args.reference_values, ratio)

    rows = []
    for row in matched:
        key = (row["mstop"], row["mlsp"])
        ref, tgt = reference[key], target[key]
        if not math.isclose(ref["theory_xsec_fb"], tgt["theory_xsec_fb"]):
            raise ValueError(f"Theory cross sections differ at {key}")
        theory_pb = tgt["theory_xsec_fb"] / 1000.0
        r_137 = row["ratio"]
        run3 = quantile_limits_pb(tgt)
        calibrated = {name: value / r_137 for name, value in run3.items()}
        rows.append(
            {
                "mstop": key[0],
                "mlsp": key[1],
                "sigma95_run3_137_pb": row["run3_limit_pb"],
                "sigma95_run2_expected_137_pb": row["run2_hepdata_limit_pb"],
                "R_137": r_137,
                "calibration_factor": 1.0 / r_137,
                "sigma95_run3_283_pb": run3["median"],
                "sigma95_run3_283_m1_pb": run3["m1"],
                "sigma95_run3_283_p1_pb": run3["p1"],
                "sigma95_calibrated_283_pb": calibrated["median"],
                "sigma95_calibrated_283_m1_pb": calibrated["m1"],
                "sigma95_calibrated_283_p1_pb": calibrated["p1"],
                "theory_xsec_fb": tgt["theory_xsec_fb"],
                "run3_283_over_theory": run3["median"] / theory_pb,
                "run3_283_m1_over_theory": run3["m1"] / theory_pb,
                "run3_283_p1_over_theory": run3["p1"] / theory_pb,
                "calibrated_283_over_theory": calibrated["median"] / theory_pb,
                "calibrated_283_m1_over_theory": calibrated["m1"] / theory_pb,
                "calibrated_283_p1_over_theory": calibrated["p1"] / theory_pb,
                "run2_bin_center_mstop": row["run2_hepdata_mstop"],
                "run2_bin_center_mlsp": row["run2_hepdata_mlsp"],
            }
        )

    x_centers = np.array(sorted({row["mstop"] for row in rows}), dtype=float)
    y_centers = np.array(sorted({row["mlsp"] for row in rows}), dtype=float)
    x_edges = centers_to_edges(x_centers, 50.0)
    y_edges = centers_to_edges(y_centers, 50.0)
    y_edges[0] = max(0.0, y_edges[0])
    x_index = {value: index for index, value in enumerate(x_centers)}
    y_index = {value: index for index, value in enumerate(y_centers)}
    shape = (len(y_centers), len(x_centers))
    grid = lambda key: to_grid(rows, key, x_index, y_index, shape)

    run3 = grid("sigma95_run3_283_pb")
    factor = grid("calibration_factor")
    calibrated = grid("sigma95_calibrated_283_pb")
    limit_values = np.concatenate((run3[np.isfinite(run3)], calibrated[np.isfinite(calibrated)]))
    limit_norm = LogNorm(vmin=float(limit_values.min()), vmax=float(limit_values.max()))
    factor_norm, factor_ticks, _ = ratio_scale(factor[np.isfinite(factor)])
    annotate = not args.no_annotations
    lumi_label = rf"{lumi:g} fb$^{{-1}}$"
    run3_title = rf"Run 3: 13.6 TeV, {lumi_label}" "\nMedian expected 95% CL upper limit (MC)"
    run3_contours = [
        (r"Expected $-1\sigma$", grid("run3_283_m1_over_theory"), "--"),
        ("Median expected", grid("run3_283_over_theory"), "-"),
        (r"Expected $+1\sigma$", grid("run3_283_p1_over_theory"), "--"),
    ]
    calibrated_contours = [
        (r"Expected $-1\sigma$", grid("calibrated_283_m1_over_theory"), "--"),
        ("Median expected", grid("calibrated_283_over_theory"), "-"),
        (r"Expected $+1\sigma$", grid("calibrated_283_p1_over_theory"), "--"),
    ]
    limit_label = "Median expected 95% CL upper limit on cross section [pb]"
    args.output_dir.mkdir(parents=True, exist_ok=True)
    written = []

    # A: the full-luminosity Run-3 MC limit surface on its own.
    fig = plt.figure(figsize=(9.6, 7.8), constrained_layout=True)
    layout = fig.add_gridspec(2, 1, height_ratios=[1.0, 0.055])
    ax = fig.add_subplot(layout[0, 0])
    mesh = draw_grid(
        ax, x_edges, y_edges, x_centers, y_centers, run3,
        norm=limit_norm, cmap_name=args.limit_cmap, title=run3_title,
        formatter=format_limit, annotate=annotate,
    )
    add_top_mass_line(ax)
    absent_single = add_run3_contours(ax, x_centers, y_centers, run3_contours)
    ax.set_ylabel(r"$m_{\tilde{\chi}^{0}_{1}}$ [GeV]")
    fig.colorbar(mesh, cax=fig.add_subplot(layout[1, 0]), orientation="horizontal").set_label(limit_label)
    stem, extensions = output_stem(args.output_dir / f"t2tt_run3_{lumi:g}fb_expected_2d")
    for extension in extensions:
        fig.savefig(stem.with_suffix(f".{extension}"), dpi=220)
        written.append(stem.with_suffix(f".{extension}"))
    plt.close(fig)

    # B: Run-3 MC limit, calibration factor, and calibrated projection.
    fig = plt.figure(figsize=(26.0, 7.8), constrained_layout=True)
    layout = fig.add_gridspec(2, 3, height_ratios=[1.0, 0.055])
    axes = [fig.add_subplot(layout[0, index]) for index in range(3)]
    caxes = [fig.add_subplot(layout[1, index]) for index in range(3)]
    for ax in axes[1:]:
        ax.sharex(axes[0])
        ax.sharey(axes[0])

    run3_mesh = draw_grid(
        axes[0], x_edges, y_edges, x_centers, y_centers, run3,
        norm=limit_norm, cmap_name=args.limit_cmap, title=run3_title,
        formatter=format_limit, annotate=annotate,
    )
    factor_mesh = draw_grid(
        axes[1], x_edges, y_edges, x_centers, y_centers, factor,
        norm=factor_norm, cmap_name=args.factor_cmap,
        title=(
            rf"Run-2 calibration factor 1/R$_{{{lumi_ref:g}}}$"
            rf"""
Run 2 expected / Run 3 expected at {lumi_ref:g} fb$^{{-1}}$"""
        ),
        formatter=lambda value: f"{value:.2f}", annotate=annotate,
    )
    calibrated_mesh = draw_grid(
        axes[2], x_edges, y_edges, x_centers, y_centers, calibrated,
        norm=limit_norm, cmap_name=args.limit_cmap,
        title=(
            rf"Run-2-calibrated full Run 3: {lumi_label}"
            "\nProjected median expected 95% CL upper limit"
        ),
        formatter=format_limit, annotate=annotate,
    )
    for ax in axes:
        add_top_mass_line(ax)
    absent_run3 = add_run3_contours(axes[0], x_centers, y_centers, run3_contours)
    handles: list[Line2D] = []
    absent_calibrated = add_run3_contours(
        axes[2], x_centers, y_centers, calibrated_contours,
        label_prefix="Calibrated: ", legend_handles=handles,
    )
    run2_graphs = add_run2_contours(
        axes[2], args.run2_root, color=RUN2_COLOR,
        label_prefix=rf"Run 2 ({lumi_ref:g} fb$^{{-1}}$): ", legend_handles=handles,
    )
    axes[2].legend(handles=handles, loc="upper right", fontsize=7.5, framealpha=0.9)
    axes[0].set_ylabel(r"$m_{\tilde{\chi}^{0}_{1}}$ [GeV]")
    for ax in axes[1:]:
        ax.tick_params(labelleft=False)

    fig.colorbar(run3_mesh, cax=caxes[0], orientation="horizontal").set_label(limit_label)
    factor_bar = fig.colorbar(factor_mesh, cax=caxes[1], orientation="horizontal")
    factor_bar.set_label(rf"Calibration factor 1/R$_{{{lumi_ref:g}}}$")
    factor_bar.ax.xaxis.set_major_locator(FixedLocator(factor_ticks))
    factor_bar.ax.xaxis.set_major_formatter(FixedFormatter([f"{value:g}" for value in factor_ticks]))
    factor_bar.ax.xaxis.set_minor_locator(NullLocator())
    fig.colorbar(calibrated_mesh, cax=caxes[2], orientation="horizontal").set_label(
        "Projected " + limit_label[0].lower() + limit_label[1:]
    )
    stem, extensions = output_stem(args.output_dir / "t2tt_full_run3_calibrated_projection")
    for extension in extensions:
        fig.savefig(stem.with_suffix(f".{extension}"), dpi=220)
        written.append(stem.with_suffix(f".{extension}"))
    plt.close(fig)

    values_path = stem.with_name(f"{stem.name}_values.csv")
    write_values(values_path, rows)

    array = lambda key: np.array([row[key] for row in rows])
    scaling = array("sigma95_run3_283_pb") / array("sigma95_run3_137_pb")
    print(f"Mass points: {len(rows)} ({mode}, no heatmap interpolation)")
    print(summary(f"Run3_{lumi:g} / Run3_{lumi_ref:g}", scaling))
    print(summary(f"Calibration factor 1/R_{lumi_ref:g}", array("calibration_factor")))
    print(summary(
        f"Calibrated_{lumi:g} / Run2_{lumi_ref:g}",
        array("sigma95_calibrated_283_pb") / array("sigma95_run2_expected_137_pb"),
    ))
    for label, absent in (("Run-3 MC", absent_run3), ("calibrated", absent_calibrated), ("single-panel", absent_single)):
        if absent:
            print(f"{label} contours without an in-grid crossing: {', '.join(absent)}")
    print(
        "Run-2 contour TGraphs: "
        + ", ".join(
            f"{name} ({count} points; {'visible' if visible else 'outside axes'})"
            for name, (count, visible) in run2_graphs.items()
        )
    )
    print("Wrote plots: " + ", ".join(str(path) for path in written))
    print(f"Wrote values: {values_path}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
