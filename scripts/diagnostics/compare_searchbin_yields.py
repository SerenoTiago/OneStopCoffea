"""Compare local total-background SearchBinYield (174 bins) before/after the
njet_bin_index fix, and merge in the already-downloaded published pre-fit
background/observed values from background_comparison_174_bins.csv.

Usage:
    uv run --no-sync python3 scripts/diagnostics/compare_searchbin_yields.py \\
        --before analysis_products/old_results/26-08-30_hadsusy_1250_run2 \\
        --after analysis_products/old_results/26-09-10_hadsusy_backgrounds_jetbins_fixed \\
        --reference docs/investigations/t2tt_comparison_2026-09-10/background_comparison_174_bins.csv \\
        --output-csv docs/investigations/t2tt_comparison_2026-09-10/background_comparison_174_bins_before_after.csv \\
        --output-plot analysis_products/plots/diagnostics/searchbin_yield_before_after.png

Background dataset names are hardcoded to match hadronic_susy_backgrounds.yaml.
"""

from __future__ import annotations

import argparse
import csv
import glob as globmod
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from analyzer.core.results import loadResults, mergeAndScale
from analyzer.utils.structure_tools import globWithMeta

BACKGROUND_DATASETS = (
    "qcd_inclusive_2024",
    "tt_hadronic_2024",
    "tt_semileptonic_2024",
    "wjets_2024",
    "zjets_2024",
)


def loadTotalBackgroundSearchBinYield(result_dir, hist_name="SearchBinYield", pipeline="susy_full_validation"):
    paths = sorted(
        p
        for p in globmod.glob(os.path.join(result_dir, "*.result"))
        if any(ds in os.path.basename(p) for ds in BACKGROUND_DATASETS)
    )
    if not paths:
        raise RuntimeError(f"No background result files found in {result_dir}")
    results = loadResults(paths)
    mergeAndScale(results)
    total = None
    found_datasets = set()
    for item, meta in globWithMeta(results, ["*", "*", pipeline, hist_name]):
        found_datasets.add(meta.get("dataset_name"))
        h = item.histogram
        total = h.copy() if total is None else total + h
    missing = set(BACKGROUND_DATASETS) - found_datasets
    if missing:
        print(f"WARNING: {result_dir} is missing datasets: {sorted(missing)}")
    if total is None:
        raise RuntimeError(f"No '{hist_name}' histogram found in {result_dir}")
    axes = [a.name for a in total.axes]
    if "variation" in axes:
        total = total[{"variation": "central" if "central" in list(total.axes["variation"]) else 0}]
    return total.values()  # 174 bins, index 0 == bin 1


def loadReference(path):
    rows = {}
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows[int(row["bin"])] = row
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--before", required=True)
    ap.add_argument("--after", required=True)
    ap.add_argument("--reference", required=True)
    ap.add_argument("--output-csv", required=True)
    ap.add_argument("--output-plot", required=True)
    args = ap.parse_args()

    before_vals = loadTotalBackgroundSearchBinYield(args.before)
    after_vals = loadTotalBackgroundSearchBinYield(args.after)
    ref = loadReference(args.reference)

    n_bins = len(before_vals)
    bins = list(range(1, n_bins + 1))

    with open(args.output_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "bin",
                "Njets",
                "Nbjets",
                "MHT",
                "HT",
                "local_background_before",
                "local_background_after",
                "paper_background",
                "before_over_paper",
                "after_over_paper",
            ]
        )
        for b, before_v, after_v in zip(bins, before_vals, after_vals):
            r = ref.get(b, {})
            paper_bkg = float(r["paper_background"]) if r.get("paper_background") else float("nan")
            before_ratio = before_v / paper_bkg if paper_bkg else float("nan")
            after_ratio = after_v / paper_bkg if paper_bkg else float("nan")
            writer.writerow(
                [
                    b,
                    r.get("Njets", ""),
                    r.get("Nbjets", ""),
                    r.get("MHT", ""),
                    r.get("HT", ""),
                    before_v,
                    after_v,
                    paper_bkg,
                    before_ratio,
                    after_ratio,
                ]
            )
    print(f"Saved {args.output_csv}")

    paper_vals = np.array(
        [float(ref[b]["paper_background"]) if b in ref else np.nan for b in bins]
    )

    fig, ax = plt.subplots(figsize=(14, 5))
    ax.plot(
        bins,
        before_vals,
        drawstyle="steps-mid",
        label="Before fix",
        color="C1",
    )
    ax.plot(
        bins,
        after_vals,
        drawstyle="steps-mid",
        label="After Fix",
        color="C0",
    )
    ax.plot(
        bins,
        paper_vals,
        drawstyle="steps-mid",
        label="Paper",
        color="k",
        linewidth=1,
    )
    ax.set_yscale("log")
    ax.set_xlabel("Search bin")
    ax.set_ylabel("Total background yield")
    ax.set_title(
        "Total background SearchBinYield: before/after jet-bin fix vs published"
    )
    ax.legend()
    fig.tight_layout()
    output_dir = os.path.dirname(args.output_plot)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    fig.savefig(args.output_plot, dpi=150)
    plt.close(fig)
    print(f"Saved {args.output_plot}")


if __name__ == "__main__":
    main()
