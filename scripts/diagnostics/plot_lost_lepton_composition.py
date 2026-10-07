#!/usr/bin/env python3
"""Plot W+jets and semileptonic-ttbar composition of the 174-bin approximation."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import uproot


KINEMATIC = (
    ("300-350", "300-600"), ("300-350", "600-1200"), ("300-350", ">=1200"),
    ("350-600", "350-600"), ("350-600", "600-1200"), ("350-600", ">=1200"),
    ("600-850", "600-1200"), ("600-850", ">=1200"),
    (">=850", "850-1700"), (">=850", ">=1700"),
)
NJETS = ("2-3", "4-5", "6-7", "8-9", ">=10")


def parse_args():
    base = Path("analysis_products/combine/2026-10-06_run3_2024_projection_137fb_t2tt_1250_offline_trigger_approx")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--shapes", type=Path, default=base / "T2tt_mStop-1250_mLSP-1/shapes_bins_hadsusy.root")
    parser.add_argument("--bins", type=Path, default=Path("analysis_products/plots/174searchbinperbackground/background_validation_174_bins.csv"))
    parser.add_argument("--output-dir", type=Path, default=Path("analysis_products/plots/174searchbinperbackground/lost_lepton_composition"))
    return parser.parse_args()


def read_histograms(path):
    output = {}
    with uproot.open(path) as handle:
        for key, name in (("wjets", "wjets_2024"), ("tt_semileptonic", "tt_semileptonic_2024")):
            hist = handle[name]
            output[key] = (np.asarray(hist.values(), float), np.asarray(hist.variances(), float))
    if any(len(values) != 174 for values, _ in output.values()):
        raise ValueError("Expected 174-bin input histograms")
    return output


def sum_group(values, variances, mask):
    return values[mask].sum(), np.sqrt(variances[mask].sum())


def save(fig, stem):
    stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(stem.with_suffix(".png"), dpi=180, bbox_inches="tight")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def main():
    args = parse_args()
    with args.bins.open() as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 174:
        raise ValueError(f"Expected 174 bin definitions, got {len(rows)}")
    hist = read_histograms(args.shapes)
    w, wvar = hist["wjets"]
    tt, ttvar = hist["tt_semileptonic"]
    total = w + tt
    total_var = wvar + ttvar

    # Use physical inclusive categories: >=2 combines the paper's low-Njet >=2
    # and the high-Njet 2 and >=3 categories.
    nb_labels = ("0", "1", ">=2")
    raw_nb = np.asarray([row["Nbjets"] for row in rows])
    nb_masks = (raw_nb == "0", raw_nb == "1", np.isin(raw_nb, ("2", ">=2", ">=3")))
    by_nb = {key: [] for key in hist}
    by_nb_unc = {key: [] for key in hist}
    for key, (values, variances) in hist.items():
        for mask in nb_masks:
            value, uncertainty = sum_group(values, variances, mask)
            by_nb[key].append(value)
            by_nb_unc[key].append(uncertainty)

    x = np.arange(3)
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.8))
    axes[0].bar(x, by_nb["wjets"], color="#67a9cf", label="W+jets")
    axes[0].bar(x, by_nb["tt_semileptonic"], bottom=by_nb["wjets"], color="#ef8a62", label=r"Semileptonic $t\bar t$")
    axes[0].errorbar(x, [sum(v) for v in zip(by_nb["wjets"], by_nb["tt_semileptonic"])],
                     yerr=np.hypot(by_nb_unc["wjets"], by_nb_unc["tt_semileptonic"]), fmt="none", color="black", capsize=2)
    axes[0].set_yscale("log")
    axes[0].set_ylabel("Expected events at 137 fb$^{-1}$")
    axes[0].legend(frameon=False)
    fractions = np.divide(by_nb["tt_semileptonic"], np.add(by_nb["wjets"], by_nb["tt_semileptonic"]))
    axes[1].bar(x, fractions, color="#ef8a62", label=r"Semileptonic $t\bar t$")
    axes[1].bar(x, 1 - fractions, bottom=fractions, color="#67a9cf", label="W+jets")
    axes[1].set_ylim(0, 1)
    axes[1].set_ylabel("Fraction of local lost-lepton approximation")
    for index, fraction in enumerate(fractions):
        axes[1].text(index, fraction / 2, f"tt: {fraction:.1%}", ha="center", va="center", fontsize=10)
        axes[1].text(index, fraction + (1 - fraction) / 2, f"W: {1-fraction:.1%}", ha="center", va="center", fontsize=10)
    for ax in axes:
        ax.set_xticks(x, nb_labels)
        ax.set_xlabel(r"$N_{b\mathrm{-jet}}$")
    fig.suptitle("Run-3 direct-MC lost-lepton composition by b-tag multiplicity")
    fig.tight_layout()
    save(fig, args.output_dir / "composition_by_nb")

    # Aggregate Njet x inclusive-Nb composition.
    njet_values = np.asarray([row["Njets"] for row in rows])
    tt_fraction = np.full((len(NJETS), 3), np.nan)
    yields = np.zeros((len(NJETS), 3))
    for iy, njet in enumerate(NJETS):
        for ix, nb_mask in enumerate(nb_masks):
            mask = (njet_values == njet) & nb_mask
            yields[iy, ix] = total[mask].sum()
            if yields[iy, ix] > 0:
                tt_fraction[iy, ix] = tt[mask].sum() / yields[iy, ix]
    fig, ax = plt.subplots(figsize=(7.3, 5.2))
    image = ax.imshow(tt_fraction, vmin=0, vmax=1, cmap="RdBu_r", aspect="auto")
    for iy in range(len(NJETS)):
        for ix in range(3):
            ax.text(ix, iy, f"tt: {tt_fraction[iy, ix]:.1%}\nN={yields[iy, ix]:.3g}", ha="center", va="center")
    ax.set_xticks(range(3), nb_labels)
    ax.set_yticks(range(len(NJETS)), NJETS)
    ax.set_xlabel(r"$N_{b\mathrm{-jet}}$")
    ax.set_ylabel(r"$N_{jet}$")
    ax.set_title(r"Semileptonic $t\bar t$ fraction (and total expected yield)")
    fig.colorbar(image, ax=ax, label=r"Semileptonic $t\bar t$ fraction")
    fig.tight_layout()
    save(fig, args.output_dir / "tt_fraction_by_njet_nb")

    # Full 19-category x 10-kinematic-interval view.
    row_groups = []
    for njet in NJETS:
        for row in rows:
            group = (njet, row["Nbjets"])
            if row["Njets"] == njet and group not in row_groups:
                row_groups.append(group)
    matrix = np.full((len(row_groups), 10), np.nan)
    bin_ids = np.zeros_like(matrix, dtype=int)
    for index, row in enumerate(rows):
        iy = row_groups.index((row["Njets"], row["Nbjets"]))
        ix = KINEMATIC.index((row["MHT"], row["HT"]))
        if total[index] > 0:
            matrix[iy, ix] = tt[index] / total[index]
        bin_ids[iy, ix] = index + 1
    fig, ax = plt.subplots(figsize=(14, 8.8))
    image = ax.imshow(matrix, vmin=0, vmax=1, cmap="RdBu_r", aspect="auto")
    for iy in range(matrix.shape[0]):
        for ix in range(10):
            if np.isfinite(matrix[iy, ix]):
                ax.text(ix, iy, f"{matrix[iy, ix]:.0%}\n[{bin_ids[iy, ix]}]", ha="center", va="center", fontsize=7.3)
    ax.set_xticks(range(10), range(1, 11))
    ax.set_yticks(range(len(row_groups)), [f"{njet}; {nb}" for njet, nb in row_groups])
    ax.set_xlabel(r"Kinematic interval index ($H_T^{miss}$, $H_T$)")
    ax.set_ylabel(r"$N_{jet}$; $N_{b\mathrm{-jet}}$")
    ax.set_title(r"Semileptonic $t\bar t$ fraction of local lost-lepton approximation")
    fig.colorbar(image, ax=ax, label=r"Semileptonic $t\bar t$ fraction")
    fig.tight_layout()
    save(fig, args.output_dir / "tt_fraction_by_searchbin_group")

    lines = ["Validated 137 fb^-1 local lost-lepton composition", f"W+jets total: {w.sum():.8g}", f"Semileptonic ttbar total: {tt.sum():.8g}"]
    for label, mask in zip(nb_labels, nb_masks, strict=True):
        wsum, ttsum = w[mask].sum(), tt[mask].sum()
        lines.append(f"Nb {label}: W={wsum:.8g}, tt_semileptonic={ttsum:.8g}, total={wsum+ttsum:.8g}, tt_fraction={ttsum/(wsum+ttsum):.6f}")
    (args.output_dir / "composition_summary.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"Wrote {args.output_dir}")


if __name__ == "__main__":
    main()
