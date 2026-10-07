#!/usr/bin/env python3
"""Diagnose variable dependence of Run-3/Run-2 174-bin background ratios."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np


KINEMATIC = (
    ("300-350", "300-600"),
    ("300-350", "600-1200"),
    ("300-350", ">=1200"),
    ("350-600", "350-600"),
    ("350-600", "600-1200"),
    ("350-600", ">=1200"),
    ("600-850", "600-1200"),
    ("600-850", ">=1200"),
    (">=850", "850-1700"),
    (">=850", ">=1700"),
)
NJETS = ("2-3", "4-5", "6-7", "8-9", ">=10")
NBJETS = ("0", "1", "2", ">=2", ">=3")
MHT_BINS = ("300-350", "350-600", "600-850", ">=850")
HT_BINS = ("300-600", "350-600", "600-1200", ">=1200", "850-1700", ">=1700")
BACKGROUNDS = ("lost_lepton", "total")


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("analysis_products/plots/174searchbinperbackground/background_validation_174_bins.csv"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("analysis_products/plots/174searchbinperbackground/ratio_pattern_diagnostics"),
    )
    parser.add_argument("--max-relative-uncertainty", type=float, default=0.5)
    parser.add_argument("--min-run2-yield", type=float, default=1.0)
    return parser.parse_args()


def load_rows(path: Path):
    with path.open() as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 174:
        raise ValueError(f"Expected 174 rows, found {len(rows)}")
    for row in rows:
        pair = (row["MHT"], row["HT"])
        if pair not in KINEMATIC:
            raise ValueError(f"Unknown kinematic interval in bin {row['bin']}: {pair}")
        row["kinematic_interval"] = KINEMATIC.index(pair) + 1
    return rows


def arrays(rows, background):
    def get(name):
        return np.asarray([float(row[name]) for row in rows])

    run3 = get(f"run3_{background}")
    run3_unc = get(f"run3_{background}_stat")
    run2 = get(f"run2_{background}")
    run2_unc = get(f"run2_{background}_total_unc")
    ratio = np.divide(run3, run2, out=np.full(174, np.nan), where=run2 > 0)
    relative = np.hypot(
        np.divide(run3_unc, run3, out=np.full(174, np.inf), where=run3 > 0),
        np.divide(run2_unc, run2, out=np.full(174, np.inf), where=run2 > 0),
    )
    return run3, run3_unc, run2, run2_unc, ratio, relative


def aggregate(rows, background, field, order):
    run3, run3_unc, run2, run2_unc, _, _ = arrays(rows, background)
    labels = np.asarray([str(row[field]) for row in rows])
    if field == "kinematic_interval":
        labels = np.asarray([str(row[field]) for row in rows])
    output = []
    for label in order:
        mask = labels == str(label)
        a, b = run3[mask].sum(), run2[mask].sum()
        au, bu = np.sqrt(np.sum(run3_unc[mask] ** 2)), np.sqrt(np.sum(run2_unc[mask] ** 2))
        ratio = a / b if b > 0 else np.nan
        rel = np.hypot(au / a, bu / b) if a > 0 and b > 0 else np.inf
        output.append((str(label), a, au, b, bu, ratio, ratio * rel))
    return output


def plot_heatmap(rows, background, reliable, output):
    _, _, _, _, ratio, _ = arrays(rows, background)
    row_groups = []
    for njet in NJETS:
        present = []
        for row in rows:
            if row["Njets"] == njet and row["Nbjets"] not in present:
                present.append(row["Nbjets"])
        row_groups.extend((njet, nbjet) for nbjet in present)
    matrix = np.full((len(row_groups), 10), np.nan)
    quality = np.zeros_like(matrix, dtype=bool)
    bins = np.zeros_like(matrix, dtype=int)
    for index, row in enumerate(rows):
        y = row_groups.index((row["Njets"], row["Nbjets"]))
        x = int(row["kinematic_interval"]) - 1
        matrix[y, x] = ratio[index]
        quality[y, x] = reliable[index]
        bins[y, x] = int(row["bin"])

    fig, ax = plt.subplots(figsize=(14, 9))
    image = ax.imshow(matrix, aspect="auto", cmap="RdBu_r", norm=LogNorm(0.25, 4.0))
    for y in range(matrix.shape[0]):
        for x in range(matrix.shape[1]):
            if not np.isfinite(matrix[y, x]):
                continue
            text = f"{matrix[y, x]:.2g}\n[{bins[y, x]}]"
            ax.text(x, y, text, ha="center", va="center", fontsize=7.2, color="black")
            if not quality[y, x]:
                ax.plot(x, y, marker="x", color="0.3", markersize=11, markeredgewidth=1.2)
    ax.set_xticks(range(10), [str(i) for i in range(1, 11)])
    ax.set_yticks(range(len(row_groups)), [f"{njet}; {nbjet}" for njet, nbjet in row_groups])
    ax.set_xlabel(r"Kinematic interval index ($H_T^{miss}$, $H_T$)")
    ax.set_ylabel(r"$N_{jet}$; $N_{b\mathrm{-jet}}$")
    ax.set_title(f"{background.replace('_', ' ').title()}: Run-3 / Run-2 ratio by search-bin definition\nNumbers in brackets are search-bin IDs; × marks low-precision bins")
    colorbar = fig.colorbar(image, ax=ax, pad=0.015, extend="both")
    colorbar.set_label("Run 3 / Run 2 (color clipped outside 0.25–4)")
    fig.tight_layout()
    fig.savefig(output.with_suffix(".png"), dpi=180)
    fig.savefig(output.with_suffix(".pdf"))
    plt.close(fig)


def plot_profiles(rows, background, output):
    specs = (
        ("kinematic_interval", tuple(range(1, 11)), "Kinematic interval"),
        ("Njets", NJETS, r"$N_{jet}$"),
        ("Nbjets", NBJETS, r"$N_{b\mathrm{-jet}}$"),
    )
    global_ratio = sum(float(row[f"run3_{background}"]) for row in rows) / sum(float(row[f"run2_{background}"]) for row in rows)
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8))
    for ax, (field, order, xlabel) in zip(axes, specs, strict=True):
        points = aggregate(rows, background, field, order)
        x = np.arange(len(points))
        y = np.asarray([point[5] for point in points])
        yerr = np.asarray([point[6] for point in points])
        ax.errorbar(x, y, yerr=yerr, fmt="o-", color="#333333", capsize=2)
        ax.axhline(global_ratio, color="#2166ac", linestyle="--", label=f"Overall = {global_ratio:.2f}")
        ax.axhline(1, color="0.5", linestyle=":")
        ax.set_xticks(x, [point[0] for point in points], rotation=35 if field != "kinematic_interval" else 0)
        ax.set_xlabel(xlabel)
        ax.set_yscale("log")
        ax.grid(axis="y", which="both", alpha=0.2)
        ax.legend(frameon=False, fontsize=9)
    axes[0].set_ylabel("Ratio of summed yields (Run 3 / Run 2)")
    fig.suptitle(f"{background.replace('_', ' ').title()}: aggregate ratio dependence", fontsize=15)
    fig.tight_layout()
    fig.savefig(output.with_suffix(".png"), dpi=180)
    fig.savefig(output.with_suffix(".pdf"))
    plt.close(fig)


def plot_ht_mht_profiles(rows, background, output):
    specs = (("MHT", MHT_BINS, r"$H_T^{miss}$ interval [GeV]"), ("HT", HT_BINS, r"$H_T$ interval [GeV]"))
    global_ratio = sum(float(row[f"run3_{background}"]) for row in rows) / sum(float(row[f"run2_{background}"]) for row in rows)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    for ax, (field, order, xlabel) in zip(axes, specs, strict=True):
        points = aggregate(rows, background, field, order)
        x = np.arange(len(points))
        y = np.asarray([point[5] for point in points])
        yerr = np.asarray([point[6] for point in points])
        ax.errorbar(x, y, yerr=yerr, fmt="o-", color="#333333", capsize=2)
        ax.axhline(global_ratio, color="#2166ac", linestyle="--", label=f"Overall = {global_ratio:.2f}")
        ax.axhline(1, color="0.5", linestyle=":")
        ax.set_xticks(x, [point[0] for point in points], rotation=30)
        ax.set_xlabel(xlabel)
        ax.set_yscale("log")
        ax.grid(axis="y", which="both", alpha=0.2)
        ax.legend(frameon=False, fontsize=9)
    axes[0].set_ylabel("Ratio of summed yields (Run 3 / Run 2)")
    fig.suptitle(f"{background.replace('_', ' ').title()}: separate HTmiss and HT dependence", fontsize=15)
    fig.tight_layout()
    fig.savefig(output.with_suffix(".png"), dpi=180)
    fig.savefig(output.with_suffix(".pdf"))
    plt.close(fig)


def weighted_additive_model(rows, background, reliable):
    _, _, _, _, ratio, relative = arrays(rows, background)
    groups = {
        "kinematic": [str(row["kinematic_interval"]) for row in rows],
        "njet": [row["Njets"] for row in rows],
        "nbjet": [row["Nbjets"] for row in rows],
    }
    mask = reliable & np.isfinite(ratio) & (ratio > 0)
    y = np.log(ratio[mask])
    weights = 1 / np.maximum(relative[mask], 0.05) ** 2

    def design(include):
        columns = [np.ones(np.count_nonzero(mask))]
        for name in include:
            values = np.asarray(groups[name])[mask]
            for level in sorted(set(values))[1:]:
                columns.append((values == level).astype(float))
        return np.column_stack(columns)

    def fit(include):
        x = design(include)
        root_w = np.sqrt(weights)
        beta = np.linalg.lstsq(x * root_w[:, None], y * root_w, rcond=None)[0]
        residual = y - x @ beta
        sse = np.sum(weights * residual**2)
        mean = np.average(y, weights=weights)
        sst = np.sum(weights * (y - mean) ** 2)
        return 1 - sse / sst, sse

    names = tuple(groups)
    full_r2, full_sse = fit(names)
    drops = {}
    for name in names:
        reduced_r2, reduced_sse = fit(tuple(item for item in names if item != name))
        drops[name] = {"delta_r2": full_r2 - reduced_r2, "delta_chi2": reduced_sse - full_sse}
    return np.count_nonzero(mask), full_r2, drops


def main():
    args = parse_args()
    rows = load_rows(args.input)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    report = []
    report.append(f"Input: {args.input}")
    report.append(f"Reliable-bin definition: Run-2 yield >= {args.min_run2_yield:g}; local and reference relative uncertainties <= {args.max_relative_uncertainty:g}")
    report.append("Kinematic intervals: " + "; ".join(f"{i}: MHT {mht}, HT {ht}" for i, (mht, ht) in enumerate(KINEMATIC, 1)))

    for background in BACKGROUNDS:
        run3, run3_unc, run2, run2_unc, ratio, relative = arrays(rows, background)
        local_rel = np.divide(run3_unc, run3, out=np.full(174, np.inf), where=run3 > 0)
        paper_rel = np.divide(run2_unc, run2, out=np.full(174, np.inf), where=run2 > 0)
        reliable = (run2 >= args.min_run2_yield) & (local_rel <= args.max_relative_uncertainty) & (paper_rel <= args.max_relative_uncertainty)
        plot_heatmap(rows, background, reliable, args.output_dir / f"{background}_ratio_heatmap")
        plot_profiles(rows, background, args.output_dir / f"{background}_aggregate_profiles")
        plot_ht_mht_profiles(rows, background, args.output_dir / f"{background}_ht_mht_profiles")
        count, r2, drops = weighted_additive_model(rows, background, reliable)
        report.append(f"\n[{background}] overall ratio={run3.sum()/run2.sum():.6g}; reliable bins={count}/174; weighted additive log-ratio R2={r2:.4f}")
        report.append("Variable importance (drop in weighted R2 when removed): " + ", ".join(f"{name}={result['delta_r2']:.4f}" for name, result in drops.items()))
        for field, order in (("kinematic_interval", range(1, 11)), ("MHT", MHT_BINS), ("HT", HT_BINS), ("Njets", NJETS), ("Nbjets", NBJETS)):
            values = aggregate(rows, background, field, order)
            report.append(field + ": " + ", ".join(f"{label}={value:.3g}±{unc:.2g}" for label, _, _, _, _, value, unc in values))
        score = np.full(174, -np.inf)
        valid = reliable & np.isfinite(ratio) & (ratio > 0)
        score[valid] = np.abs(np.log(ratio[valid])) / np.maximum(relative[valid], 0.05)
        report.append("Most significant reliable-bin discrepancies:")
        for index in np.argsort(score)[-10:][::-1]:
            if not np.isfinite(score[index]):
                continue
            row = rows[index]
            report.append(
                f"  bin {index+1}: Njet={row['Njets']}, Nb={row['Nbjets']}, interval={row['kinematic_interval']} "
                f"(MHT={row['MHT']}, HT={row['HT']}), ratio={ratio[index]:.3g}, combined rel unc={relative[index]:.3g}"
            )

    report_path = args.output_dir / "quantitative_summary.txt"
    report_path.write_text("\n".join(report) + "\n")
    print("\n".join(report))
    print(f"Wrote diagnostics to {args.output_dir}")


if __name__ == "__main__":
    main()
