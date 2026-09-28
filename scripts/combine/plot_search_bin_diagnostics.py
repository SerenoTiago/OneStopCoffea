#!/usr/bin/env python3
"""Plot signal/background sensitivity and weighted-MC statistics by search bin."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mplhep
import numpy as np
import uproot


DEFAULT_BACKGROUNDS = (
    "qcd_inclusive_2024",
    "tt_hadronic_2024",
    "tt_semileptonic_2024",
    "wjets_2024",
    "zjets_2024",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("shapes", type=Path, help="Combine shapes ROOT file")
    parser.add_argument(
        "--signal", default="2stop_1250_2024_splitLSP", help="Signal histogram name"
    )
    parser.add_argument(
        "--background",
        action="append",
        dest="backgrounds",
        help="Background histogram name; repeat as needed",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("search_bin_diagnostics.png"),
    )
    parser.add_argument(
        "--important-fraction",
        type=float,
        default=0.90,
        help="Fraction of cumulative S^2/B used to mark important bins",
    )
    parser.add_argument("--lumi", type=float, default=137.0)
    parser.add_argument("--energy", type=float, default=13.6)
    return parser.parse_args()


def search_bin_boundaries() -> tuple[list[float], list[float]]:
    """Boundaries for the 174-bin CMS-style Njet/Nb layout."""
    major, minor = [], []
    offset = 0
    for njet_index in range(5):
        n_b_bins = 3 if njet_index == 0 else 4
        n_kinematic_bins = 10 if njet_index < 3 else 8
        for b_index in range(n_b_bins):
            offset += n_kinematic_bins
            if b_index < n_b_bins - 1:
                minor.append(offset + 0.5)
        if njet_index < 4:
            major.append(offset + 0.5)
    return major, minor


def values_and_variances(root_file, name: str) -> tuple[np.ndarray, np.ndarray]:
    obj = root_file[name]
    values = np.asarray(obj.values(flow=False), dtype=float)
    variances = obj.variances(flow=False)
    if variances is None:
        raise ValueError(f"Histogram {name!r} has no stored variances")
    return values, np.asarray(variances, dtype=float)


def main() -> None:
    args = parse_args()
    backgrounds = args.backgrounds or list(DEFAULT_BACKGROUNDS)

    with uproot.open(args.shapes) as root_file:
        signal, _ = values_and_variances(root_file, args.signal)
        background_parts = [values_and_variances(root_file, name) for name in backgrounds]

    background = np.sum([item[0] for item in background_parts], axis=0)
    background_variance = np.sum([item[1] for item in background_parts], axis=0)
    if signal.shape != background.shape:
        raise ValueError("Signal and background histograms have different binning")

    ratio = np.divide(signal, background, out=np.zeros_like(signal), where=background > 0)
    neff = np.divide(
        background**2,
        background_variance,
        out=np.zeros_like(background),
        where=background_variance > 0,
    )
    sensitivity = np.divide(
        signal**2, background, out=np.zeros_like(signal), where=background > 0
    )
    order = np.argsort(sensitivity)[::-1]
    cumulative = np.cumsum(sensitivity[order])
    important = np.zeros(signal.size, dtype=bool)
    if cumulative.size and cumulative[-1] > 0:
        important[order[cumulative <= args.important_fraction * cumulative[-1]]] = True
        important[order[0]] = True

    bins = np.arange(1, signal.size + 1)
    edges = np.arange(0.5, signal.size + 1.5)
    plt.style.use(mplhep.style.CMS)
    fig, axes = plt.subplots(
        3,
        1,
        figsize=(18, 12),
        sharex=True,
        gridspec_kw={"height_ratios": (1.25, 1, 1), "hspace": 0.08},
    )

    axes[0].stairs(background, edges, label="Total background", color="black")
    axes[0].stairs(signal, edges, label="Signal", color="tab:blue", linewidth=1.8)
    axes[0].set_yscale("log")
    axes[0].set_ylabel("Expected events")
    axes[0].legend(fontsize=15)

    axes[1].stairs(ratio, edges, color="tab:purple")
    axes[1].set_yscale("log")
    axes[1].set_ylabel(r"$S/B$")

    axes[2].stairs(neff, edges, color="tab:green", label=r"Background $N_{\rm eff}$")
    axes[2].scatter(
        bins[important & (neff < 25)],
        neff[important & (neff < 25)],
        color="tab:red",
        s=28,
        zorder=3,
        label=r"Important and $N_{\rm eff}<25$",
    )
    for threshold, label in ((100, "10%"), (25, "20%"), (10, "32%")):
        axes[2].axhline(threshold, color="gray", linestyle="--", linewidth=1)
        axes[2].text(174, threshold * 1.08, label, ha="right", va="bottom", fontsize=11)
    axes[2].set_yscale("log")
    axes[2].set_ylim(bottom=0.5)
    axes[2].set_ylabel(r"Background $N_{\rm eff}$")
    axes[2].set_xlabel("Search region bin number")
    axes[2].legend(fontsize=13, loc="upper right")

    major, minor = search_bin_boundaries()
    for axis in axes:
        axis.set_xlim(0.5, signal.size + 0.5)
        for boundary in minor:
            axis.axvline(boundary, color="gray", linestyle=":", linewidth=0.8)
        for boundary in major:
            axis.axvline(boundary, color="black", linestyle="--", linewidth=1)

    mplhep.cms.label(
        "Preliminary",
        data=False,
        lumi=args.lumi,
        com=args.energy,
        ax=axes[0],
    )
    fig.suptitle(
        "Search-bin sensitivity and weighted-MC statistical precision",
        y=0.995,
        fontsize=20,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=160, bbox_inches="tight")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
