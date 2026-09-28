"""Compare saved SR total and component backgrounds with CMS pre-fit projections.

Run from the repository root with the analysis environment active:
    PYTHONPATH=. python scripts/diagnostics/compare_background_projections.py
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
import yaml

from analyzer.core.results import loadResults, mergeAndScale
from analyzer.utils.structure_tools import globWithMeta


BACKGROUND_DATASETS = (
    "qcd_inclusive_2024",
    "tt_hadronic_2024",
    "tt_semileptonic_2024",
    "wjets_2024",
    "zjets_2024",
)
COMPONENTS = {
    "z_nunu": ("Z to neutrinos", "Z->nunu", ("zjets_2024",), "#d94747"),
    "lost_lepton": ("Lost lepton", "Lost-lepton", ("wjets_2024", "tt_semileptonic_2024"), "#3274b8"),
    "qcd": ("QCD", "QCD", ("qcd_inclusive_2024",), "#e4ba2f"),
}
TABLE_FILES = (
    "yields_and_pre_fit_bg_njets_2_3.yaml",
    "yields_and_pre_fit_bg_njets_4_5.yaml",
    "yields_and_pre_fit_bg_njets_6_7.yaml",
    "yields_and_pre_fit_bg_njets_8_9.yaml",
    "yields_and_pre_fit_bg_njets_10.yaml",
)
OBSERVABLES = {
    "njet": ("JetN", "$N_{\\mathrm{jet}}$", (2, 4, 6, 8, 10), ("2-3", "4-5", "6-7", "8-9", ">=10")),
    "nbjet": ("BJetN", "$N_{b\\mathrm{-jet}}$", (0, 1, 2), ("0", "1", ">=2")),
    "htmiss": ("HTMiss", "$H_T^{\\mathrm{miss}}$ [GeV]", (300, 350, 600, 850), ("300-350", "350-600", "600-850", ">=850")),
}


def cms_category(label: str) -> str:
    return label.replace("#", "").replace("$", "").replace("\\geq", ">=").replace("--", "-").strip()


def paper_projections(table_dir: Path) -> tuple[dict[str, np.ndarray], dict[str, dict[str, np.ndarray]]]:
    totals = {key: np.zeros(len(spec[3])) for key, spec in OBSERVABLES.items()}
    components = {
        key: {category: np.zeros(len(spec[3])) for category in COMPONENTS}
        for key, spec in OBSERVABLES.items()
    }
    n_rows = 0
    for filename in TABLE_FILES:
        with (table_dir / filename).open() as f:
            table = yaml.safe_load(f)
        columns = {entry["header"]["name"]: entry["values"] for entry in table["dependent_variables"]}
        for row in zip(*(columns[name] for name in ("Njets", "Nbjets", "MHT", "Total BG", "Z->nunu", "Lost-lepton", "QCD")), strict=True):
            njet, nbjet, mht, total, z_nunu, lost_lepton, qcd = row
            njet_label = cms_category(str(njet["value"]))
            nbjet_label = cms_category(str(nbjet["value"]))
            mht_label = cms_category(str(mht["value"]))
            # The 2-3-jet table groups >=2 b jets; this is the finest
            # inclusive Nb projection common to all five CMS tables.
            nbjet_label = nbjet_label if nbjet_label in ("0", "1") else ">=2"
            value = float(total["value"])
            category_values = {"z_nunu": z_nunu, "lost_lepton": lost_lepton, "qcd": qcd}
            for key, label in (("njet", njet_label), ("nbjet", nbjet_label), ("htmiss", mht_label)):
                labels = OBSERVABLES[key][3]
                if label not in labels:
                    raise ValueError(f"Unexpected CMS {key} category: {label}")
                index = labels.index(label)
                totals[key][index] += value
                for category, entry in category_values.items():
                    components[key][category][index] += float(entry["value"])
            n_rows += 1
    if n_rows != 174:
        raise ValueError(f"Expected 174 CMS pre-fit search bins, found {n_rows}")
    return totals, components


def project_histogram(hist, key: str) -> np.ndarray:
    _, _, lower_edges, _ = OBSERVABLES[key]
    if "variation" in [axis.name for axis in hist.axes]:
        hist = hist[{"variation": "central"}]
    if len(hist.axes) != 1:
        raise ValueError(f"Expected one-dimensional {key} histogram")
    values = np.asarray(hist.values(flow=True), dtype=float)
    centers = np.asarray(hist.axes[0].centers, dtype=float)
    if key in ("njet", "nbjet"):
        centers = np.rint(centers)
    if values[0] != 0:
        raise ValueError(f"Unexpected {key} underflow: {values[0]}")
    if np.any(values[1:-1][centers < lower_edges[0]]):
        raise ValueError(f"Unexpected {key} events below CMS projection range")
    projected = np.histogram(centers, bins=[*lower_edges, np.inf], weights=values[1:-1])[0]
    projected[-1] += values[-1]
    return projected


def local_projections(result_dir: Path, pipeline: str) -> tuple[dict[str, np.ndarray], dict[str, dict[str, np.ndarray]]]:
    paths = sorted(p for p in result_dir.glob("*.result") if p.name.split("__", 1)[0] in BACKGROUND_DATASETS)
    if not paths:
        raise ValueError(f"No background .result files in {result_dir}")
    results = loadResults([str(p) for p in paths], keep_patterns=[("*", "*", "_provenance"), ("*", "*", "pipelines", pipeline, "*")])
    mergeAndScale(results)
    totals = {}
    components = {}
    for key, (hist_name, _, _, _) in OBSERVABLES.items():
        by_dataset = {}
        found = set()
        for item, meta in globWithMeta(results, ["*", "pipelines", pipeline, hist_name]):
            dataset = meta["dataset_name"]
            found.add(dataset)
            by_dataset[dataset] = project_histogram(item.histogram, key)
        missing = set(BACKGROUND_DATASETS) - found
        if missing:
            raise ValueError(f"Missing {hist_name} datasets in {result_dir}: {sorted(missing)}")
        totals[key] = sum(by_dataset.values())
        components[key] = {
            category: sum(by_dataset[name] for name in sample_names)
            for category, (_, _, sample_names, _) in COMPONENTS.items()
        }
    return totals, components


def ratio(numerator: np.ndarray, denominator: np.ndarray) -> np.ndarray:
    return np.divide(numerator, denominator, out=np.full_like(numerator, np.nan), where=denominator > 0)


def plot_projection(key: str, local: np.ndarray, paper: np.ndarray, output: Path, normalized: bool) -> None:
    _, xlabel, _, labels = OBSERVABLES[key]
    if normalized:
        local = local / local.sum()
        paper = paper / paper.sum()
    x = np.arange(len(labels))
    fig, (upper, lower) = plt.subplots(
        2, 1, figsize=(7.3, 5.5), sharex=True,
        gridspec_kw={"height_ratios": [3, 1], "hspace": 0.08},
    )
    for position, counts, name, hatch in (
        (x - 0.19, paper, "CMS pre-fit total BG", None),
        (x + 0.19, local, "Local total BG", "///"),
    ):
        upper.bar(position, counts, 0.34, label=name, color="#686868",
                  edgecolor="0.2", linewidth=0.55, hatch=hatch)
    upper.set_ylabel("Fraction of total" if normalized else "Expected events")
    upper.set_title(f"{xlabel}: full signal region")
    upper.legend(frameon=False, loc="upper right")
    if not normalized:
        upper.set_yscale("log")
        upper.set_ylim(bottom=max(min(np.r_[paper[paper > 0], local[local > 0]]) * 0.5, 1e-3))
    lower.axhline(1, color="0.5", linewidth=1, linestyle="--")
    lower.plot(x, ratio(local, paper), linestyle="none", marker="o", markersize=5,
               color="#555555", markeredgecolor="0.2", markeredgewidth=0.4)
    lower.set_ylabel("Local / CMS")
    lower.set_xticks(x, labels)
    lower.set_xlabel(xlabel)
    lower.grid(axis="y", alpha=0.2)
    fig.subplots_adjust(left=0.14, right=0.98, bottom=0.15, top=0.91)
    fig.savefig(output, dpi=170)
    plt.close(fig)


def plot_components(
    key: str,
    local: dict[str, np.ndarray],
    paper: dict[str, np.ndarray],
    output: Path,
    normalized: bool,
) -> None:
    _, xlabel, _, labels = OBSERVABLES[key]
    x = np.arange(len(labels))
    fig, (upper, lower) = plt.subplots(
        2, 1, figsize=(8.1, 6.0), sharex=True,
        gridspec_kw={"height_ratios": [3, 1.2], "hspace": 0.08},
    )
    categories = tuple(COMPONENTS)
    offsets = {"z_nunu": -0.13, "lost_lepton": 0, "qcd": 0.13}

    if normalized:
        local_values = {name: local[name] / local[name].sum() for name in categories}
        paper_values = {name: paper[name] / paper[name].sum() for name in categories}
    else:
        local_values, paper_values = local, paper

    if normalized:
        edges = np.arange(len(labels) + 1) - 0.5
        for name in categories:
            color = COMPONENTS[name][3]
            upper.stairs(paper_values[name], edges, color=color, linewidth=2)
            upper.stairs(local_values[name], edges, color=color, linewidth=2, linestyle="--")
        upper.set_ylabel("Fraction of component")
        upper.set_ylim(bottom=0)
        color_handles = [Line2D([], [], color=COMPONENTS[name][3], linewidth=2,
                                label=COMPONENTS[name][0]) for name in categories]
        style_handles = [Line2D([], [], color="0.3", linewidth=2, label="CMS"),
                         Line2D([], [], color="0.3", linewidth=2, linestyle="--", label="Local")]
    else:
        for position, source, hatch in ((x - 0.19, paper_values, None), (x + 0.19, local_values, "///")):
            bottom = np.zeros(len(labels))
            for name in categories:
                values = source[name]
                upper.bar(position, values, 0.34, bottom=bottom, color=COMPONENTS[name][3],
                          edgecolor="0.25", linewidth=0.45, hatch=hatch)
                bottom = bottom + values
        upper.set_yscale("log")
        upper.set_ylabel("Expected events")
        color_handles = [Patch(facecolor=COMPONENTS[name][3], label=COMPONENTS[name][0])
                         for name in categories]
        style_handles = [Patch(facecolor="white", edgecolor="0.3", label="CMS"),
                         Patch(facecolor="white", edgecolor="0.3", hatch="///", label="Local")]

    title = (f"{xlabel}: {'independently normalized component shapes' if normalized else 'component yields'}"
             " in full signal region")
    upper.set_title(title)
    upper.legend(handles=color_handles + style_handles, frameon=False, ncol=2,
                 fontsize=9, loc="upper right")
    lower.axhline(1, color="0.45", linewidth=1, linestyle="--")
    for name in categories:
        values = ratio(local_values[name], paper_values[name])
        lower.plot(x + offsets[name], values, linestyle="none", marker="o", markersize=5,
                   color=COMPONENTS[name][3], markeredgecolor="0.25", markeredgewidth=0.4)
    lower.set_yscale("log")
    lower.set_ylabel("Local / CMS")
    lower.set_xticks(x, labels)
    lower.set_xlim(-0.55, len(labels) - 0.45)
    lower.set_xlabel(xlabel)
    lower.grid(axis="y", alpha=0.2)
    fig.subplots_adjust(left=0.13, right=0.98, bottom=0.13, top=0.92)
    fig.savefig(output, dpi=170)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--results",
        type=Path,
        default=Path(
            "analysis_products/old_results/26-09-10_hadsusy_backgrounds_jetbins_fixed"
        ),
    )
    parser.add_argument("--tables", type=Path, default=Path("data/SUS-19-006/audit_2026-09-10"))
    parser.add_argument("--pipeline", default="susy_full_validation")
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("analysis_products/plots/diagnostics/local_vs_cms_projections"),
    )
    args = parser.parse_args()
    paper, paper_components = paper_projections(args.tables)
    local, local_components = local_projections(args.results, args.pipeline)
    args.output.mkdir(parents=True, exist_ok=True)
    with (args.output / "projection_values.csv").open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(("observable", "bin", "local", "cms_pre_fit", "local_over_cms", "local_shape", "cms_shape", "shape_ratio"))
        for key, (_, _, _, labels) in OBSERVABLES.items():
            local_shape = local[key] / local[key].sum()
            paper_shape = paper[key] / paper[key].sum()
            for label, l, p, ls, ps in zip(labels, local[key], paper[key], local_shape, paper_shape, strict=True):
                writer.writerow((key, label, l, p, l / p if p else "", ls, ps, ls / ps if ps else ""))
            for normalized, suffix in ((False, "absolute"), (True, "shape")):
                plot_projection(key, local[key], paper[key], args.output / f"{key}_{suffix}.png", normalized)
                plot_components(key, local_components[key], paper_components[key],
                                args.output / f"{key}_components_{suffix}.png", normalized)
            print(f"{key}: local={local[key].sum():.6g}, CMS={paper[key].sum():.6g}")
    with (args.output / "component_projection_values.csv").open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(("observable", "bin", "component", "local", "cms_pre_fit", "local_over_cms",
                         "local_shape", "cms_shape", "shape_ratio"))
        for key, (_, _, _, labels) in OBSERVABLES.items():
            for name in COMPONENTS:
                local_values = local_components[key][name]
                paper_values = paper_components[key][name]
                local_shape = local_values / local_values.sum()
                paper_shape = paper_values / paper_values.sum()
                for label, l, p, ls, ps in zip(labels, local_values, paper_values, local_shape, paper_shape, strict=True):
                    writer.writerow((key, label, name, l, p, l / p if p else "", ls, ps, ls / ps if ps else ""))
    for source_name, total, by_category in (("CMS", paper, paper_components), ("Local", local, local_components)):
        remainder = total["njet"].sum() - sum(by_category["njet"][name].sum() for name in COMPONENTS)
        print(f"{source_name} total minus plotted components: {remainder:.6g} events")
    print(f"Saved twelve plots and two projection CSVs in {args.output}")


if __name__ == "__main__":
    main()
