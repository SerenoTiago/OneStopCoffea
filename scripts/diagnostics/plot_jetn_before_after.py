"""Compare JetN (and optionally other 1D histograms) between two saved-result
directories, e.g. before and after the njet_bin_index inclusive-endpoint fix.

Usage:
    uv run --no-sync python3 scripts/diagnostics/plot_jetn_before_after.py \\
        --before analysis_products/old_results/26-08-30_hadsusy_1250_run2 \\
        --after analysis_products/old_results/26-09-10_hadsusy_jetbins_fixed/signal \\
        --sample-glob 'signal_T2tt_1250_100' \\
        --hist-name JetN \\
        --output analysis_products/plots/diagnostics/jetn_signal_1250_100_before_after.png
"""

from __future__ import annotations

import argparse
import glob as globmod
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from analyzer.core.results import loadResults, mergeAndScale
from analyzer.utils.structure_tools import globWithMeta


def loadMergedHist(result_dir, sample_glob, hist_name, pipeline, restrict_basenames=None):
    # Only load files matching the sample glob: mergeAndScale collapses all
    # samples within a dataset into one, so unrelated samples must be excluded
    # before merging rather than filtered out afterward.
    paths = sorted(
        p
        for p in globmod.glob(os.path.join(result_dir, "*.result"))
        if sample_glob in os.path.basename(p)
        and (restrict_basenames is None or os.path.basename(p) in restrict_basenames)
    )
    if not paths:
        raise RuntimeError(
            f"No result files matching '{sample_glob}' found in {result_dir}"
        )
    results = loadResults(paths)
    mergeAndScale(results)
    total = None
    for item, meta in globWithMeta(results, ["*", "*", pipeline, hist_name]):
        h = item.histogram
        total = h.copy() if total is None else total + h
    if total is None:
        raise RuntimeError(
            f"No histogram '{hist_name}' found for sample matching '{sample_glob}' in {result_dir}"
        )
    return total


def toCentralCounts(h):
    # Select the "central" variation if present, sum over any remaining non-axis dims.
    axes = [a.name for a in h.axes]
    if "variation" in axes:
        h = h[{"variation": "central" if "central" in list(h.axes["variation"]) else 0}]
    return h.values(), [a.name for a in h.axes if a.name != "variation"], h.axes[-1].edges


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--before", required=True)
    ap.add_argument("--after", required=True)
    ap.add_argument("--sample-glob", required=True)
    ap.add_argument("--hist-name", default="JetN")
    ap.add_argument("--pipeline", default="susy_full_validation")
    ap.add_argument("--output", required=True)
    ap.add_argument(
        "--restrict-to-after-samples",
        action="store_true",
        help="Only sum 'before' samples whose result filename also exists in --after "
        "(for a fair comparison when --after is a partial/low-stats rerun).",
    )
    args = ap.parse_args()

    restrict = None
    if args.restrict_to_after_samples:
        restrict = {
            os.path.basename(p)
            for p in globmod.glob(os.path.join(args.after, "*.result"))
        }
        print(f"Restricting 'before' to {len(restrict)} sample(s) present in 'after'")

    before_h = loadMergedHist(
        args.before, args.sample_glob, args.hist_name, args.pipeline, restrict
    )
    after_h = loadMergedHist(args.after, args.sample_glob, args.hist_name, args.pipeline)

    before_vals, _, edges = toCentralCounts(before_h)
    after_vals, _, _ = toCentralCounts(after_h)

    centers = 0.5 * (edges[:-1] + edges[1:])
    width = edges[1] - edges[0]

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.bar(centers - width / 4, before_vals, width=width / 2, label="Before fix", color="C1")
    ax.bar(centers + width / 4, after_vals, width=width / 2, label="After fix", color="C0")
    for odd in (3, 5, 7, 9):
        ax.axvline(odd, color="gray", linestyle=":", linewidth=1)
    ax.set_xlabel("$N_{jets}$")
    ax.set_ylabel("Expected weighted yield")
    ax.set_title(f"{args.hist_name}: {args.sample_glob} ({args.pipeline})")
    ax.legend()
    fig.tight_layout()

    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    fig.savefig(args.output, dpi=150)
    print(f"Saved {args.output}")

    print("\nbin  before      after       after/before")
    for i, (b, a) in enumerate(zip(before_vals, after_vals)):
        ratio = a / b if b else float("nan")
        print(f"{i:3d}  {b:10.4f}  {a:10.4f}  {ratio:8.3f}")


if __name__ == "__main__":
    main()
