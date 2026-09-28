#!/usr/bin/env python3
"""Plot split-signal Combine limits as a mass-plane heatmap."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
import numpy as np


LIMIT_COLUMNS = ("observed", "exp", "exp_m1", "exp_p1", "exp_m2", "exp_p2")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Plot Combine signal-strength limits from collect_split_signal_limits.py. "
            "By default, r limits are converted to pb using xsec_fb / 1000."
        )
    )
    parser.add_argument("limits_csv", type=Path)
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("analysis_products/combine/split_signals/limits_heatmap"),
        help="Output path without extension, or with .png/.pdf extension.",
    )
    parser.add_argument(
        "--quantity",
        choices=LIMIT_COLUMNS,
        default="observed",
        help="Limit column to plot.",
    )
    parser.add_argument(
        "--xsec-fb",
        type=float,
        default=0.7526,
        help="Signal cross section in fb used to convert r to pb.",
    )
    parser.add_argument(
        "--plot-r",
        action="store_true",
        help="Plot raw Combine r limits instead of cross-section limits.",
    )
    parser.add_argument(
        "--strip-width",
        type=float,
        default=50.0,
        help="Width in GeV used when plotting a single stop-mass column.",
    )
    parser.add_argument(
        "--linear",
        action="store_true",
        help="Use a linear color scale instead of a log scale.",
    )
    parser.add_argument(
        "--paper_scale",
        "--paper-scale",
        action="store_true",
        help=(
            "Use the paper color scale for cross-section limits: "
            "1e-4 to 10 pb, and append '-paper_scale' to output names."
        ),
    )
    parser.add_argument(
        "--cmap",
        default="rainbow",
        help="Matplotlib colormap name.",
    )
    parser.add_argument(
        "--title",
        default=r"$pp \rightarrow \tilde{t}\tilde{t}, \tilde{t} \rightarrow t\tilde{\chi}^{0}_{1}$",
    )
    parser.add_argument(
        "--write-table",
        type=Path,
        default=None,
        help="Optional output CSV containing the plotted value for each point.",
    )
    return parser.parse_args()


def read_rows(path: Path) -> list[dict[str, float | str]]:
    rows: list[dict[str, float | str]] = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"mstop", "mlsp", *LIMIT_COLUMNS}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path} is missing required columns: {sorted(missing)}")
        for row in reader:
            parsed: dict[str, float | str] = dict(row)
            parsed["mstop"] = float(row["mstop"])
            parsed["mlsp"] = float(row["mlsp"])
            for column in LIMIT_COLUMNS:
                parsed[column] = float(row[column]) if row[column] else math.nan
            rows.append(parsed)
    if not rows:
        raise ValueError(f"No rows found in {path}")
    return rows


def centers_to_edges(centers: np.ndarray, single_width: float) -> np.ndarray:
    centers = np.asarray(sorted(float(x) for x in centers), dtype=float)
    if len(centers) == 1:
        return np.array([centers[0] - single_width / 2.0, centers[0] + single_width / 2.0])

    midpoints = (centers[1:] + centers[:-1]) / 2.0
    first = centers[0] - (midpoints[0] - centers[0])
    last = centers[-1] + (centers[-1] - midpoints[-1])
    return np.concatenate([[first], midpoints, [last]])


def output_stem(path: Path) -> tuple[Path, list[str]]:
    if path.suffix.lower() in {".png", ".pdf"}:
        return path.with_suffix(""), [path.suffix.lower().lstrip(".")]
    return path, ["png", "pdf"]


def write_augmented_table(path: Path, rows: list[dict[str, float | str]], value_key: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["mstop", "mlsp", *LIMIT_COLUMNS, value_key, "limit_log"]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


def main() -> int:
    args = parse_args()
    if args.paper_scale and args.plot_r:
        raise ValueError("--paper_scale applies to cross-section limits, not raw r limits")

    rows = read_rows(args.limits_csv)

    value_key = args.quantity if args.plot_r else f"{args.quantity}_xsec_pb"
    for row in rows:
        limit_r = float(row[args.quantity])
        row[value_key] = limit_r if args.plot_r else limit_r * args.xsec_fb / 1000.0

    x_centers = np.array(sorted({float(row["mstop"]) for row in rows}), dtype=float)
    y_centers = np.array(sorted({float(row["mlsp"]) for row in rows}), dtype=float)
    x_edges = centers_to_edges(x_centers, args.strip_width)
    y_edges = centers_to_edges(y_centers, args.strip_width)
    y_edges[0] = max(0.0, y_edges[0])

    z = np.full((len(y_centers), len(x_centers)), np.nan, dtype=float)
    x_index = {x: i for i, x in enumerate(x_centers)}
    y_index = {y: i for i, y in enumerate(y_centers)}
    for row in rows:
        z[y_index[float(row["mlsp"])]][x_index[float(row["mstop"])]] = float(row[value_key])

    finite = z[np.isfinite(z) & (z > 0)]
    if len(finite) == 0:
        raise ValueError(f"No positive finite values found for {value_key}")

    vmin = float(np.nanmin(finite))
    vmax = float(np.nanmax(finite))
    if args.paper_scale:
        vmin = 1.0e-4
        vmax = 10.0

    norm = Normalize(vmin=vmin, vmax=vmax)
    if not args.linear:
        norm = LogNorm(vmin=vmin, vmax=vmax)

    fig, ax = plt.subplots(figsize=(6.8, 6.0), constrained_layout=True)
    mesh = ax.pcolormesh(x_edges, y_edges, z, cmap=args.cmap, norm=norm, shading="flat")
    cbar = fig.colorbar(mesh, ax=ax)
    cbar.set_label(
        "95% CL upper limit on cross section [pb]" if not args.plot_r else "95% CL upper limit on signal strength r"
    )

    ax.set_xlabel(r"$m_{\tilde{t}}$ [GeV]")
    ax.set_ylabel(r"$m_{\tilde{\chi}^{0}_{1}}$ [GeV]")
    ax.set_title(args.title)
    ax.set_xlim(float(np.nanmin(x_edges)), float(np.nanmax(x_edges)))
    ax.set_ylim(float(np.nanmin(y_edges)), float(np.nanmax(y_edges)))
    ax.tick_params(direction="in", top=True, right=True)

    y_min = float(np.nanmin(y_edges))
    y_max = float(np.nanmax(y_edges))
    y_margin = 0.03 * (y_max - y_min)
    for row in rows:
        value = float(row[value_key])
        label = f"{value:.2g}" if not args.plot_r else f"{value:.2f}"
        text_y = min(max(float(row["mlsp"]), y_min + y_margin), y_max - y_margin)
        ax.text(
            float(row["mstop"]),
            text_y,
            label,
            ha="center",
            va="center",
            fontsize=8,
            color="black",
        )

    stem, extensions = output_stem(args.output)
    if args.paper_scale:
        stem = stem.with_name(f"{stem.name}-paper_scale")
    stem.parent.mkdir(parents=True, exist_ok=True)
    for extension in extensions:
        fig.savefig(stem.with_suffix(f".{extension}"), dpi=200)
    plt.close(fig)

    table_path = args.write_table or stem.with_name(f"{stem.name}_values.csv")
    write_augmented_table(table_path, rows, value_key)

    print(f"Wrote plot: {', '.join(str(stem.with_suffix(f'.{ext}')) for ext in extensions)}")
    print(f"Wrote values: {table_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
