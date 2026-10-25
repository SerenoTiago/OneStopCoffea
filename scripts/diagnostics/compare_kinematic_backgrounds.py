#!/usr/bin/env python3
"""Compare Run-3 direct-MC and SUS-19-006 pre-fit backgrounds in HT, HTmiss,
Njet, and Nbjet.

Run-3 inputs are the per-sample background .result files, merged and scaled
with the analyzer's mergeAndScale (the normalization used for the datacards).
Run-3/Run-2 ratios use only coarse bins assembled from the 174 search bins, so
both predictions share identical boundaries and selections; the fine Run-3
histograms are drawn in a separate panel.  This is a developmental validation
plotter: it performs no event processing, datacard, or Combine work.
"""

from __future__ import annotations

import argparse
import csv
import glob
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FuncFormatter, MaxNLocator
import mplhep
import numpy as np
import uproot
import yaml

from compare_searchbin_backgrounds import (
    LOCAL_SAMPLES,
    TITLES,
    expected_bin_labels,
    read_paper_tables,
)

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.append(str(REPO))
from analyzer.core.results import loadResults, mergeAndScale  # noqa: E402


CATEGORIES = ("total", "lost_lepton", "z_nunu", "qcd")
RUN3_COLOR = "#2166ac"
RUN2_COLOR = "#b2182b"
RUN3_LABEL = r"$\mathbf{Run\!\!-\!3\ direct\ MC}$ (MC stat. unc.)"
RUN2_LABEL = r"$\mathbf{Run\!\!-\!2\ data\!\!-\!driven\ pre\!\!-\!fit}$ ($\mathrm{stat.}\oplus\mathrm{syst.}$ unc.)"
HTMISS_SPLIT = "$H_{T}^{miss}$"

# hist: Run-3 histogram name; display: fine-panel range (contents outside it,
# including the histogram overflow, are folded into the edge bins).
VARIABLES = {
    "ht": {"hist": "JetHT", "xlabel": r"$H_{T}$ [GeV]", "display": (300.0, 3000.0), "unit": "50 GeV"},
    "htmiss": {"hist": "HTMiss", "xlabel": r"$H_{T}^{miss}$ [GeV]", "display": (300.0, 1500.0), "unit": "50 GeV"},
    "njet": {"hist": "JetN", "xlabel": r"$N_{\rm jet}$", "display": (1.5, 13.5), "unit": "bin"},
    "nbjet": {"hist": "BJetN", "xlabel": r"$N_{b\rm-jet}$", "display": (-0.5, 6.5), "unit": "bin"},
}
# Coarse-bin boundaries drawn on the fine panels (HT boundaries depend on HTmiss).
FINE_BOUNDARIES = {
    "ht": (),
    "htmiss": (350.0, 600.0, 850.0),
    "njet": (3.5, 5.5, 7.5, 9.5),
    "nbjet": (0.5, 1.5, 2.5),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=Path("analysis_products/production_runs/2026-10-07_012944_hadronic_susy/current_backgrounds"),
        help="Directory of complete background .result files.",
    )
    parser.add_argument("--tables", type=Path, default=Path("data/SUS-19-006/audit_2026-09-10"))
    parser.add_argument("--projections", type=Path, default=Path("data/SUS-19-006/an18271_v7_prefit_projections.yaml"))
    parser.add_argument(
        "--reference-shapes",
        type=Path,
        default=Path(
            "analysis_products/combine/26-10-11/run3_2024_t2tt_grid_expanded_1007bkg_automcstats_expected_137fb/"
            "mStop-1300/T2tt_mStop-1300_mLSP-100/shapes_bins_hadsusy.root"
        ),
        help="Datacard shapes built from the same backgrounds, used to validate the normalization.",
    )
    parser.add_argument("--output-dir", type=Path, default=Path("analysis_products/plots/kinematic_background_comparison"))
    parser.add_argument("--lumi", type=float, default=137.0)
    return parser.parse_args()


# ----------------------------------------------------------------------------
# Coarse bins built from the 174 search-bin definitions
# ----------------------------------------------------------------------------

def coarse_bins(variable: str) -> list[tuple[str, np.ndarray]]:
    """Return (label, mask over the 174 search bins) for each coarse bin."""
    labels = expected_bin_labels()
    njet = np.array([label[0] for label in labels])
    nbjet = np.array([label[1] for label in labels])
    mht = np.array([label[2] for label in labels])
    ht = np.array([label[3] for label in labels])
    low_mht = mht != ">=850"
    if variable == "htmiss":
        groups = [(text, mht == text) for text in ("300-350", "350-600", "600-850", ">=850")]
    elif variable == "njet":
        groups = [(text, njet == text) for text in ("2-3", "4-5", "6-7", "8-9", ">=10")]
    elif variable == "nbjet":
        # The Njet 2-3 search bins merge Nbjet >= 2; AN Table 24 counts them in "2".
        groups = [
            ("0", nbjet == "0"),
            ("1", nbjet == "1"),
            ("2", (nbjet == "2") | ((njet == "2-3") & (nbjet == ">=2"))),
            (">=3", nbjet == ">=3"),
        ]
    elif variable == "ht":
        # HT intervals are defined jointly with HTmiss, so the two HTmiss regions
        # are kept separate instead of forcing overlapping HT intervals together.
        groups = [
            ("300-600", low_mht & np.isin(ht, ("300-600", "350-600"))),
            ("600-1200", low_mht & (ht == "600-1200")),
            (">=1200", low_mht & (ht == ">=1200")),
            ("850-1700", ~low_mht & (ht == "850-1700")),
            (">=1700", ~low_mht & (ht == ">=1700")),
        ]
    else:
        raise KeyError(variable)
    coverage = np.sum([mask for _, mask in groups], axis=0)
    if not np.all(coverage == 1):
        raise ValueError(f"{variable} coarse bins do not partition the 174 search bins")
    return groups


def tick_label(text: str) -> str:
    return text.replace(">=", "≥").replace("-", "–")


# ----------------------------------------------------------------------------
# Inputs
# ----------------------------------------------------------------------------

def central_1d(histogram):
    histogram = getattr(histogram, "histogram", histogram)
    if len(histogram.axes) == 2:
        categories = list(histogram.axes[0])
        if not categories:
            return None
        index = categories.index("central") if "central" in categories else 0
        histogram = histogram[{histogram.axes[0].name: index}]
    return histogram


def load_run3(results_dir: Path, lumi: float):
    """Merged, luminosity-scaled histograms per category, including flow bins."""
    paths = sorted(glob.glob(str(results_dir / "*.result")))
    if not paths:
        raise FileNotFoundError(f"No .result files in {results_dir}")
    results = mergeAndScale(loadResults(paths), luminosity_override=lumi)
    names = [spec["hist"] for spec in VARIABLES.values()] + ["SearchBinYield"]
    per_dataset: dict[str, dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]]] = {}
    for dataset in results.keys():
        pipeline = results[dataset]["pipelines"]["susy_full_validation"]
        per_dataset[dataset] = {}
        for name in names:
            histogram = central_1d(pipeline[name])
            edges = np.asarray(histogram.axes[0].edges, dtype=float)
            values = np.asarray(histogram.values(flow=True), dtype=float)
            variances = np.asarray(histogram.variances(flow=True), dtype=float)
            per_dataset[dataset][name] = (values, variances, edges)
    required = {sample for samples in LOCAL_SAMPLES.values() for sample in samples}
    missing = sorted(required - set(per_dataset))
    if missing:
        raise ValueError(f"Missing Run-3 datasets: {missing}")
    run3 = {}
    for key in CATEGORIES:
        run3[key] = {}
        for name in names:
            parts = [per_dataset[sample][name] for sample in LOCAL_SAMPLES[key]]
            run3[key][name] = (
                np.sum([part[0] for part in parts], axis=0),
                np.sum([part[1] for part in parts], axis=0),
                parts[0][2],
            )
    return run3, per_dataset, len(paths)


def load_an_projections(path: Path) -> dict[str, dict[str, dict[str, np.ndarray]]]:
    with path.open() as handle:
        tables = yaml.safe_load(handle)
    out = {}
    for variable in ("htmiss", "njet", "nbjet"):
        bins = tables[variable]["bins"]
        out[variable] = {"labels": [entry["label"] for entry in bins]}
        for key in CATEGORIES:
            rows = np.array([entry[key] for entry in bins], dtype=float)
            out[variable][key] = {
                "value": rows[:, 0],
                "up": np.hypot(rows[:, 1], rows[:, 3]),
                "down": np.hypot(rows[:, 2], rows[:, 4]),
            }
    return out


def project_paper(values, stat, total_unc, mask):
    """Coarse projection of the 174-bin table: stat in quadrature, syst linear.

    Treating the systematic uncertainties as fully correlated across the
    merged bins is a conservative approximation; the AN tables, which account
    for the actual correlations, are used wherever they exist.
    """
    syst = np.sqrt(np.clip(total_unc ** 2 - stat ** 2, 0.0, None))
    value = values[mask].sum()
    unc = np.hypot(np.sqrt(np.sum(stat[mask] ** 2)), syst[mask].sum())
    return value, unc


# ----------------------------------------------------------------------------
# Projections
# ----------------------------------------------------------------------------

def run3_coarse(run3_category, groups):
    values, variances, _ = run3_category["SearchBinYield"]
    inner_values, inner_variances = values[1:-1], variances[1:-1]
    return (
        np.array([inner_values[mask].sum() for _, mask in groups]),
        np.sqrt(np.array([inner_variances[mask].sum() for _, mask in groups])),
    )


def run2_coarse(variable, key, groups, an, paper):
    if variable in an:
        if [label for label, _ in groups] != an[variable]["labels"]:
            raise ValueError(f"AN {variable} bins do not match the search-bin projection")
        entry = an[variable][key]
        return entry["value"], entry["up"], entry["down"]
    values, stat, total_unc = paper
    projected = [project_paper(values[key], stat[key], total_unc[key], mask) for _, mask in groups]
    value = np.array([item[0] for item in projected])
    unc = np.array([item[1] for item in projected])
    return value, unc, unc


def edge_index(edges: np.ndarray, value: float) -> int:
    """Index of the bin edge equal to ``value`` (edges carry rounding error)."""
    matches = np.flatnonzero(np.isclose(edges, value, rtol=0.0, atol=1e-6))
    if len(matches) != 1:
        raise ValueError(f"{value} is not a bin edge of {edges}")
    return int(matches[0])


def fine_display(run3_category, hist_name, display):
    """Fine Run-3 histogram with everything outside the display range folded in."""
    values, variances, edges = run3_category[hist_name]
    low, high = display
    first = edge_index(edges, low)
    last = edge_index(edges, high)
    # values/variances carry underflow at index 0 and overflow at index -1.
    shown_edges = edges[first:last + 1]
    shown = values[first + 1:last + 1].copy()
    shown_var = variances[first + 1:last + 1].copy()
    below = values[:first + 1].sum(), variances[:first + 1].sum()
    above = values[last + 1:].sum(), variances[last + 1:].sum()
    shown[0] += below[0]
    shown_var[0] += below[1]
    shown[-1] += above[0]
    shown_var[-1] += above[1]
    return shown_edges, shown, np.sqrt(shown_var), below[0], above[0]


def ratio_with_unc(numerator, numerator_unc, denominator, denominator_up, denominator_down):
    valid = (numerator > 0) & (denominator > 0)
    ratio = np.full(len(numerator), np.nan)
    unc = np.full(len(numerator), np.nan)
    ratio[valid] = numerator[valid] / denominator[valid]
    # Symmetrize the Run-2 uncertainty conservatively, as the 174-bin comparison does.
    denominator_unc = np.maximum(denominator_up, denominator_down)
    unc[valid] = ratio[valid] * np.hypot(numerator_unc[valid] / numerator[valid], denominator_unc[valid] / denominator[valid])
    return ratio, unc


# ----------------------------------------------------------------------------
# Plotting
# ----------------------------------------------------------------------------

def ratio_ticks(low: float, high: float, max_ticks: int = 4) -> list[float]:
    """Between two and ``max_ticks`` readable ticks inside a log-scale range."""
    decades = range(int(np.floor(np.log10(low))), int(np.ceil(np.log10(high))) + 1)
    options = [(subs, 1) for subs in ((1, 1.5, 2, 3, 5, 7), (1, 2, 3, 5), (1, 2, 5), (1, 3), (1,))]
    options += [((1,), stride) for stride in (2, 3, 4, 5)]
    for subs, stride in options:
        ticks = [sub * 10.0 ** exponent for exponent in decades if exponent % stride == 0 for sub in subs]
        ticks = [tick for tick in ticks if low <= tick <= high]
        if 2 <= len(ticks) <= max_ticks:
            return ticks
    # Narrow ranges: evenly spaced plain-number ticks.
    ticks = [tick for tick in MaxNLocator(nbins=3).tick_values(low, high) if low <= tick <= high and tick > 0]
    return ticks[:max_ticks]


def ratio_tick_label(value: float, _position=None) -> str:
    """Powers of ten in scientific notation, other ticks as plain numbers."""
    exponent = np.log10(value) if value > 0 else np.nan
    if np.isfinite(exponent) and np.isclose(exponent, round(exponent)):
        return rf"$10^{{{int(round(exponent))}}}$"
    return f"{value:g}"


def positive_for_log(values):
    return np.where(values > 0, values, np.nan)


def legend(axis, loc):
    axis.legend(frameon=True, loc=loc, prop={"size": 11.5, "weight": "normal"}, facecolor="white", edgecolor="0.45", framealpha=0.82)


def draw_variable(variable, comparisons, fine, output, *, shape):
    spec = VARIABLES[variable]
    plt.style.use(mplhep.style.CMS)
    fig = plt.figure(figsize=(34, 16))
    outer = fig.add_gridspec(1, 4, wspace=0.26)
    first_fine_axis = None
    for column, key in enumerate(CATEGORIES):
        inner = outer[column].subgridspec(4, 1, height_ratios=(1.0, 0.34, 1.0, 0.55), hspace=0.06)
        ax_fine = fig.add_subplot(inner[0])
        ax_coarse = fig.add_subplot(inner[2])
        ax_ratio = fig.add_subplot(inner[3], sharex=ax_coarse)
        first_fine_axis = first_fine_axis or ax_fine
        comparison = comparisons[key]

        # Fine Run-3 distribution.
        edges, values, errors, _, _ = fine[key]
        norm = values.sum() if shape else 1.0
        centers = 0.5 * (edges[1:] + edges[:-1])
        ax_fine.stairs(positive_for_log(values / norm), edges, color=RUN3_COLOR, linewidth=2.0, label=RUN3_LABEL)
        ax_fine.errorbar(centers, values / norm, yerr=errors / norm, fmt="none", ecolor=RUN3_COLOR, elinewidth=0.9, alpha=0.65)
        ax_fine.set_yscale("log")
        for boundary in FINE_BOUNDARIES[variable]:
            ax_fine.axvline(boundary, color="0.45", linestyle=":", linewidth=1.0, zorder=0)
        ax_fine.set_xlim(edges[0], edges[-1])
        ax_fine.set_xlabel(spec["xlabel"], fontsize=17)
        ax_fine.set_ylabel(("Fraction / " if shape else "Events / ") + spec["unit"], fontsize=17)
        ax_fine.set_title(TITLES[key], fontsize=19, pad=38)
        note = "Fine Run-3 bins; outermost bins include under/overflow"
        if variable == "ht":
            note += "\nInclusive in " + HTMISS_SPLIT
        ax_fine.text(0.97, 0.95, note, transform=ax_fine.transAxes, ha="right", va="top", fontsize=10.5,
                     bbox={"facecolor": "white", "alpha": 0.8, "edgecolor": "none", "pad": 1.5})
        legend(ax_fine, "lower left")
        positive = values[values > 0] / norm
        if len(positive):
            ax_fine.set_ylim(positive.min() * 0.3, positive.max() * 30)

        # Coarse Run-3 and Run-2 comparison on identical search-region bins.
        labels = comparison["labels"]
        x = np.arange(len(labels))
        coarse_edges = np.arange(len(labels) + 1) - 0.5
        r3, r3_unc = comparison["run3"], comparison["run3_unc"]
        r2, r2_up, r2_down = comparison["run2"], comparison["run2_up"], comparison["run2_down"]
        if shape:
            r3_norm, r2_norm = r3.sum(), r2.sum()
            r3, r3_unc = r3 / r3_norm, r3_unc / r3_norm
            r2, r2_up, r2_down = r2 / r2_norm, r2_up / r2_norm, r2_down / r2_norm
        ax_coarse.stairs(positive_for_log(r3), coarse_edges, color=RUN3_COLOR, linewidth=2.2, label=RUN3_LABEL)
        ax_coarse.errorbar(x - 0.06, r3, yerr=r3_unc, fmt="none", ecolor=RUN3_COLOR, elinewidth=1.2, capsize=3)
        ax_coarse.stairs(positive_for_log(r2), coarse_edges, color=RUN2_COLOR, linewidth=2.2, label=RUN2_LABEL)
        ax_coarse.errorbar(x + 0.06, r2, yerr=[r2_down, r2_up], fmt="none", ecolor=RUN2_COLOR, elinewidth=1.2, capsize=3)
        ax_coarse.set_yscale("log")
        both = np.concatenate((r3[r3 > 0], r2[r2 > 0]))
        ax_coarse.set_ylim(both.min() * 0.25, both.max() * 60)
        ax_coarse.set_ylabel("Fraction of events" if shape else "Expected events", fontsize=17)
        ax_coarse.tick_params(labelbottom=False)
        legend(ax_coarse, "upper right")
        if variable == "ht":
            for axis in (ax_coarse, ax_ratio):
                axis.axvline(2.5, color="black", linestyle=(0, (4, 4)), linewidth=1.2)
            for center, text in ((1.0, HTMISS_SPLIT + " < 850 GeV"), (3.5, HTMISS_SPLIT + r" $\geq$ 850 GeV")):
                ax_coarse.text(center, 0.04, text, transform=ax_coarse.get_xaxis_transform(), ha="center", va="bottom", fontsize=12,
                               bbox={"facecolor": "white", "alpha": 0.8, "edgecolor": "none", "pad": 1})

        ratio, ratio_unc = ratio_with_unc(r3, r3_unc, r2, r2_up, r2_down)
        valid = np.isfinite(ratio)
        ax_ratio.axhline(1.0, color="black", linestyle="--", linewidth=1.5)
        ax_ratio.errorbar(x[valid], ratio[valid], yerr=ratio_unc[valid], fmt="o", markersize=6, color="#4d4d4d", ecolor="#777777", elinewidth=1.2)
        for index in np.flatnonzero(~valid):
            ax_ratio.text(index, 1.0, "n/a", ha="center", va="bottom", fontsize=10, color="0.4")
        ax_ratio.set_yscale("log")
        if np.any(valid):
            low = np.nanmin(ratio[valid] - np.minimum(ratio_unc[valid], 0.9 * ratio[valid]))
            high = np.nanmax(ratio[valid] + ratio_unc[valid])
            ax_ratio.set_ylim(min(low, 1.0) * 0.6, max(high, 1.0) * 1.7)
        ax_ratio.yaxis.set_major_locator(FixedLocator(ratio_ticks(*ax_ratio.get_ylim())))
        ax_ratio.yaxis.set_major_formatter(FuncFormatter(ratio_tick_label))
        ax_ratio.yaxis.set_minor_formatter(FuncFormatter(lambda value, _: ""))
        ax_ratio.set_ylabel("Run 3 / Run 2" + (" (shape)" if shape else ""), fontsize=15)
        ax_ratio.set_xticks(x)
        ax_ratio.set_xticklabels([tick_label(text) for text in labels], fontsize=14)
        ax_ratio.set_xlim(-0.5, len(labels) - 0.5)
        ax_ratio.set_xlabel(spec["xlabel"] + " (search-region bins)", fontsize=17)
        ax_ratio.grid(axis="y", which="both", alpha=0.18)

    mplhep.cms.label("Developmental", data=False, lumi=137, com=13.6, ax=first_fine_axis, fontsize=15)
    kind = "shape comparison (unit area)" if shape else "absolute expected yields"
    fig.suptitle(f"Run-3 direct MC vs Run-2 pre-fit background: {spec['xlabel']}, {kind}", y=0.97, fontsize=24)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output.with_suffix(".png"), dpi=150, bbox_inches="tight")
    fig.savefig(output.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


# ----------------------------------------------------------------------------
# Tables
# ----------------------------------------------------------------------------

def write_coarse_csv(path, variable, comparisons):
    fields = ["variable", "bin", "category", "run3", "run3_mc_stat", "run2", "run2_unc_up", "run2_unc_down",
              "run2_source", "ratio", "ratio_unc", "run3_fraction", "run2_fraction", "shape_ratio", "shape_ratio_unc"]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for key in CATEGORIES:
            c = comparisons[key]
            ratio, ratio_unc = ratio_with_unc(c["run3"], c["run3_unc"], c["run2"], c["run2_up"], c["run2_down"])
            n3, n2 = c["run3"].sum(), c["run2"].sum()
            shape, shape_unc = ratio_with_unc(c["run3"] / n3, c["run3_unc"] / n3, c["run2"] / n2, c["run2_up"] / n2, c["run2_down"] / n2)
            for index, label in enumerate(c["labels"]):
                writer.writerow({
                    "variable": variable, "bin": label, "category": key,
                    "run3": f"{c['run3'][index]:.6g}", "run3_mc_stat": f"{c['run3_unc'][index]:.6g}",
                    "run2": f"{c['run2'][index]:.6g}", "run2_unc_up": f"{c['run2_up'][index]:.6g}",
                    "run2_unc_down": f"{c['run2_down'][index]:.6g}", "run2_source": c["run2_source"],
                    "ratio": f"{ratio[index]:.6g}", "ratio_unc": f"{ratio_unc[index]:.6g}",
                    "run3_fraction": f"{c['run3'][index] / n3:.6g}", "run2_fraction": f"{c['run2'][index] / n2:.6g}",
                    "shape_ratio": f"{shape[index]:.6g}", "shape_ratio_unc": f"{shape_unc[index]:.6g}",
                })


def write_fine_csv(path, fine):
    fields = ["bin_low", "bin_high"] + [f"{key}{suffix}" for key in CATEGORIES for suffix in ("", "_mc_stat")]
    edges = fine[CATEGORIES[0]][0]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for index in range(len(edges) - 1):
            row = {"bin_low": f"{edges[index]:g}", "bin_high": f"{edges[index + 1]:g}"}
            for key in CATEGORIES:
                row[key] = f"{fine[key][1][index]:.6g}"
                row[f"{key}_mc_stat"] = f"{fine[key][2][index]:.6g}"
            writer.writerow(row)


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------

def main() -> int:
    args = parse_args()
    if not np.isclose(args.lumi, 137.0):
        raise ValueError("The Run-2 references correspond to 137 fb^-1")
    run3, per_dataset, n_files = load_run3(args.results_dir, args.lumi)
    _, paper_values, paper_stat, paper_total = read_paper_tables(args.tables)
    paper = (paper_values, paper_stat, paper_total)
    an = load_an_projections(args.projections)
    report: list[str] = [f"Run-3 inputs: {n_files} result files in {args.results_dir}, scaled to {args.lumi:g} fb^-1 with mergeAndScale"]

    # Normalization against the datacard shapes built from the same backgrounds.
    if args.reference_shapes.is_file():
        worst = 0.0
        with uproot.open(args.reference_shapes) as shapes:
            for dataset, hists in per_dataset.items():
                reference = np.asarray(shapes[dataset].values(), dtype=float)
                ours = hists["SearchBinYield"][0][1:-1]
                worst = max(worst, float(np.max(np.abs(ours - reference) / np.maximum(np.abs(reference), 1e-12))))
        report.append(f"Search-bin yields vs datacard shapes ({args.reference_shapes}): max relative difference {worst:.2e}")

    # Totals: every fine histogram (with flow) must reproduce the selected yield.
    for key in CATEGORIES:
        sb_values = run3[key]["SearchBinYield"][0]
        total = sb_values[1:-1].sum()
        flow = sb_values[0] + sb_values[-1]
        diffs = {spec["hist"]: run3[key][spec["hist"]][0].sum() - total for spec in VARIABLES.values()}
        report.append(
            f"{key}: selected yield {total:.6g}; search-bin flow {flow:.3g}; fine-histogram totals (with flow) minus selected yield: "
            + ", ".join(f"{name} {diff:+.2e}" for name, diff in diffs.items())
        )
        for spec in VARIABLES.values():
            values = run3[key][spec["hist"]][0]
            report.append(f"    {spec['hist']}: underflow {values[0]:.4g}, overflow {values[-1]:.4g}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    for variable, spec in VARIABLES.items():
        groups = coarse_bins(variable)
        comparisons, fine = {}, {}
        for key in CATEGORIES:
            r3, r3_unc = run3_coarse(run3[key], groups)
            r2, r2_up, r2_down = run2_coarse(variable, key, groups, an, paper)
            comparisons[key] = {
                "labels": [label for label, _ in groups], "run3": r3, "run3_unc": r3_unc,
                "run2": r2, "run2_up": r2_up, "run2_down": r2_down,
                "run2_source": "AN-18-271 v7 projection table" if variable in an else "sum of 174-bin HEPData table (stat quad., syst linear)",
            }
            fine[key] = fine_display(run3[key], spec["hist"], spec["display"])

        # Cross-checks of the coarse projections.
        if variable in an:
            for key in CATEGORIES:
                direct = np.array([paper_values[key][mask].sum() for _, mask in groups])
                published = an[variable][key]["value"]
                rel = np.abs(direct - published) / published
                projected_unc = np.array([project_paper(paper_values[key], paper_stat[key], paper_total[key], mask)[1] for _, mask in groups])
                published_unc = np.maximum(an[variable][key]["up"], an[variable][key]["down"])
                report.append(
                    f"{variable}/{key}: AN table vs summed 174-bin values: max |diff| {np.max(np.abs(direct - published)):.3g} events "
                    f"(max rel {rel.max():.2%}); 174-bin uncertainty prescription / AN uncertainty: "
                    + ", ".join(f"{value:.2f}" for value in projected_unc / published_unc)
                )
        values_inner = {key: run3[key][spec["hist"]] for key in CATEGORIES}
        if variable in ("htmiss", "njet", "nbjet"):
            edges = values_inner["total"][2]
            boundaries = {"htmiss": [300, 350, 600, 850], "njet": [1.5, 3.5, 5.5, 7.5, 9.5], "nbjet": [-0.5, 0.5, 1.5, 2.5]}[variable]
            for key in CATEGORIES:
                values = values_inner[key][0]
                inner = values[1:-1]
                rebinned = []
                for low, high in zip(boundaries, boundaries[1:] + [np.inf]):
                    lo_idx = edge_index(edges, low)
                    hi_idx = len(inner) if np.isinf(high) else edge_index(edges, high)
                    rebinned.append(inner[lo_idx:hi_idx].sum() + (values[-1] if np.isinf(high) else 0.0))
                diff = np.array(rebinned) - comparisons[key]["run3"]
                report.append(f"{variable}/{key}: fine histogram rebinned minus search-bin projection: " + ", ".join(f"{d:+.3g}" for d in diff))

        draw_variable(variable, comparisons, fine, args.output_dir / f"{variable}_absolute_run3_mc_vs_run2_prefit", shape=False)
        draw_variable(variable, comparisons, fine, args.output_dir / f"{variable}_shape_run3_mc_vs_run2_prefit", shape=True)
        write_coarse_csv(args.output_dir / f"{variable}_coarse_comparison.csv", variable, comparisons)
        write_fine_csv(args.output_dir / f"{variable}_run3_fine.csv", fine)

    summary = args.output_dir / "validation_summary.txt"
    summary.write_text("\n".join(report) + "\n")
    print("\n".join(report))
    print(f"Wrote plots, tables, and {summary.name} to {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
