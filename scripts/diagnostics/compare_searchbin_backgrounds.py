#!/usr/bin/env python3
"""Compare Run-3 direct-MC and SUS-19-006 pre-fit yields in all 174 bins.

This is a developmental validation plotter.  It reads an existing Combine
shapes file and the locally archived HEPData Appendix-A tables; it performs no
analysis, postprocessing, datacard, or Combine work.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mplhep
import numpy as np
import uproot
import yaml


TABLE_FILES = (
    "yields_and_pre_fit_bg_njets_2_3.yaml",
    "yields_and_pre_fit_bg_njets_4_5.yaml",
    "yields_and_pre_fit_bg_njets_6_7.yaml",
    "yields_and_pre_fit_bg_njets_8_9.yaml",
    "yields_and_pre_fit_bg_njets_10.yaml",
)
PAPER_COLUMNS = {
    "qcd": "QCD",
    "z_nunu": "Z->nunu",
    "lost_lepton": "Lost-lepton",
    "total": "Total BG",
}
LOCAL_SAMPLES = {
    "qcd": ("qcd_inclusive_2024",),
    "z_nunu": ("zjets_2024",),
    # Deliberately excludes tt_hadronic_2024.  No single-top or dileptonic-ttbar
    # sample is available in these validated shapes.
    "lost_lepton": ("tt_semileptonic_2024", "wjets_2024"),
    "total": (
        "qcd_inclusive_2024",
        "zjets_2024",
        "tt_semileptonic_2024",
        "wjets_2024",
        "tt_hadronic_2024",
    ),
}
TITLES = {
    "qcd": "QCD Background",
    "z_nunu": r"$Z\to\nu\bar{\nu}$ + Jets Background",
    "lost_lepton": "Lost-Lepton Background",
    "total": "Total Background",
}


def parse_args() -> argparse.Namespace:
    base = Path(
        "analysis_products/combine/"
        "superseded/26-10-06/run3_2024_projection_137fb_t2tt_1250_offline_trigger_approx"
    )
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--shapes",
        type=Path,
        default=base / "T2tt_mStop-1250_mLSP-1/shapes_bins_hadsusy.root",
        help="One validated 137 fb^-1 shapes file (backgrounds are identical across mass points).",
    )
    parser.add_argument(
        "--tables",
        type=Path,
        default=Path("data/SUS-19-006/audit_2026-09-10"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("analysis_products/plots/174searchbinperbackground"),
    )
    parser.add_argument("--lumi", type=float, default=137.0)
    return parser.parse_args()


def clean_label(value: object) -> str:
    return (
        str(value)
        .replace("#", "")
        .replace("$", "")
        .replace("\\geq", ">=")
        .replace("--", "-")
        .strip()
    )


def symmetric_error(entry: dict, label: str) -> float:
    for error in entry.get("errors", []):
        if error.get("label") != label:
            continue
        value = error.get("symerror", error.get("asymerror"))
        if isinstance(value, dict):
            return max(abs(float(value["minus"])), abs(float(value["plus"])))
        return abs(float(value))
    return 0.0


def expected_bin_labels() -> list[tuple[str, str, str, str]]:
    kinematic = (
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
    rows = []
    for njet_index, njet in enumerate(("2-3", "4-5", "6-7", "8-9", ">=10")):
        nbjets = ("0", "1", ">=2") if njet_index == 0 else ("0", "1", "2", ">=3")
        allowed = range(10) if njet_index < 3 else (1, 2, 4, 5, 6, 7, 8, 9)
        for nbjet in nbjets:
            for index in allowed:
                mht, ht = kinematic[index]
                rows.append((njet, nbjet, mht, ht))
    return rows


def read_paper_tables(table_dir: Path):
    labels: list[tuple[str, str, str, str]] = []
    values = {key: [] for key in PAPER_COLUMNS}
    stat = {key: [] for key in PAPER_COLUMNS}
    total_unc = {key: [] for key in PAPER_COLUMNS}
    for filename in TABLE_FILES:
        with (table_dir / filename).open() as handle:
            table = yaml.safe_load(handle)
        columns = {item["header"]["name"]: item["values"] for item in table["dependent_variables"]}
        lengths = {len(column) for column in columns.values()}
        if len(lengths) != 1:
            raise ValueError(f"Mismatched column lengths in {filename}: {sorted(lengths)}")
        for index in range(lengths.pop()):
            labels.append(tuple(clean_label(columns[name][index]["value"]) for name in ("Njets", "Nbjets", "MHT", "HT")))
            for key, column_name in PAPER_COLUMNS.items():
                entry = columns[column_name][index]
                value = float(entry["value"])
                stat_error = symmetric_error(entry, "stat")
                syst_error = symmetric_error(entry, "syst")
                values[key].append(value)
                stat[key].append(stat_error)
                total_unc[key].append(np.hypot(stat_error, syst_error))
    expected = expected_bin_labels()
    if len(labels) != 174 or labels != expected:
        mismatches = [(i + 1, actual, wanted) for i, (actual, wanted) in enumerate(zip(labels, expected)) if actual != wanted]
        raise ValueError(f"HEPData/local search-bin ordering mismatch: {mismatches[:5]}")
    return labels, *( {key: np.asarray(source[key], dtype=float) for key in PAPER_COLUMNS} for source in (values, stat, total_unc) )


def read_local_shapes(path: Path):
    required = sorted({sample for samples in LOCAL_SAMPLES.values() for sample in samples})
    pieces = {}
    with uproot.open(path) as root_file:
        missing = [name for name in required if name not in root_file]
        if missing:
            raise ValueError(f"Missing local shape histograms: {missing}")
        for name in required:
            histogram = root_file[name]
            values = np.asarray(histogram.values(flow=False), dtype=float)
            variances = histogram.variances(flow=False)
            if len(values) != 174 or variances is None:
                raise ValueError(f"{name} does not contain 174 bins with variances")
            pieces[name] = (values, np.asarray(variances, dtype=float))
    local_values, local_stat = {}, {}
    for key, samples in LOCAL_SAMPLES.items():
        local_values[key] = np.sum([pieces[name][0] for name in samples], axis=0)
        local_stat[key] = np.sqrt(np.sum([pieces[name][1] for name in samples], axis=0))
    return local_values, local_stat


def boundaries_and_labels():
    major, minor, njet_labels, nbjet_labels = [], [], [], []
    offset = 0
    for njet_index, njet in enumerate(("2–3", "4–5", "6–7", "8–9", "≥10")):
        start = offset
        nkin = 10 if njet_index < 3 else 8
        nbjets = ("0", "1", "≥2") if njet_index == 0 else ("0", "1", "2", "≥3")
        for nb_index, nbjet in enumerate(nbjets):
            nb_start = offset
            offset += nkin
            nbjet_labels.append((nb_start + (nkin + 1) / 2, nbjet))
            if nb_index < len(nbjets) - 1:
                minor.append(offset + 0.5)
        njet_labels.append((start + (offset - start + 1) / 2, njet))
        if njet_index < 4:
            major.append(offset + 0.5)
    return major, minor, njet_labels, nbjet_labels


def positive_for_log(values: np.ndarray) -> np.ndarray:
    return np.where(values > 0, values, np.nan)


def plot_category(key, local, local_stat, paper, paper_stat, paper_total_unc, output):
    bins = np.arange(1, 175)
    edges = np.arange(0.5, 175.5, 1.0)
    ratio = np.divide(local, paper, out=np.full(174, np.nan), where=paper > 0)
    ratio_rel_local = np.divide(local_stat, local, out=np.zeros(174), where=local > 0)
    ratio_rel_paper = np.divide(paper_total_unc, paper, out=np.zeros(174), where=paper > 0)
    ratio_unc = ratio * np.hypot(ratio_rel_local, ratio_rel_paper)

    positive = np.concatenate((local[local > 0], paper[paper > 0]))
    ymin = max(float(np.min(positive)) * 0.035, 1e-5)
    # Reserve ample logarithmic headroom above the largest yield so the search-
    # region category labels remain visually separate from the data.
    ymax = float(np.max(positive)) * 300.0

    plt.style.use(mplhep.style.CMS)
    fig, axes = plt.subplots(
        3, 1, figsize=(19, 12), sharex=True,
        gridspec_kw={"height_ratios": (1, 1, 0.72), "hspace": 0.06},
    )
    for axis, values, errors, color, label in (
        (axes[0], local, local_stat, "#2166ac", r"$\mathbf{Run\!\!-\!3\ direct\ MC}$ (stat. unc.)"),
        (axes[1], paper, paper_total_unc, "#b2182b", r"$\mathbf{Run\!\!-\!2\ data\!\!-\!driven\ pre\!\!-\!fit}$ ($\mathrm{stat.}\oplus\mathrm{syst.}$ unc.)"),
    ):
        axis.stairs(positive_for_log(values), edges, color=color, linewidth=2.025, label=label)
        mask = values > 0
        axis.errorbar(bins[mask], values[mask], yerr=errors[mask], fmt="none", ecolor=color, elinewidth=0.975, alpha=0.65)
        axis.set_yscale("log")
        axis.set_ylim(ymin, ymax)
        axis.set_ylabel("Expected events")
        is_paper_panel = axis is axes[1]
        axis.legend(
            frameon=True,
            loc="upper right" if is_paper_panel else "lower left",
            prop={"size": 15.3, "weight": "normal"},
            facecolor="white",
            edgecolor="0.45",
            framealpha=0.82,
        )

    valid = np.isfinite(ratio) & (ratio > 0)
    axes[2].axhline(1.0, color="black", linestyle="--", linewidth=1.5)
    axes[2].errorbar(bins[valid], ratio[valid], yerr=ratio_unc[valid], fmt="o", markersize=3.6, color="#4d4d4d", ecolor="#777777", elinewidth=0.825)
    axes[2].set_yscale("log")
    if np.any(valid):
        low = max(np.nanmin(ratio[valid] - np.minimum(ratio_unc[valid], 0.9 * ratio[valid])) * 0.075, 1e-4)
        high = np.nanmax(ratio[valid] + ratio_unc[valid]) * 13.5
        axes[2].set_ylim(low, high)
    axes[2].set_ylabel("Run 3 / Run 2")
    axes[2].set_xlabel("Search region bin number")
    axes[2].grid(axis="y", which="both", alpha=0.18)

    major, minor, njet_labels, nbjet_labels = boundaries_and_labels()
    for axis in axes:
        axis.set_xlim(0.5, 174.5)
        for boundary in minor:
            axis.axvline(boundary, color="0.45", linestyle=":", linewidth=0.75, zorder=0)
        for boundary in major:
            axis.axvline(boundary, color="black", linestyle=(0, (4, 4)), linewidth=1.0, zorder=0)
    for x, text in njet_labels:
        axes[0].text(x, 0.975, rf"$N_{{\rm jet}}={text}$", transform=axes[0].get_xaxis_transform(), ha="center", va="top", fontsize=16.5, fontweight="semibold", bbox={"facecolor": "white", "alpha": 0.8, "edgecolor": "none", "pad": 1})
    for index, (x, text) in enumerate(nbjet_labels):
        label = rf"$N_{{b\rm-jet}}={text}$" if index == 0 else text
        axes[0].text(x, 0.86, label, transform=axes[0].get_xaxis_transform(), ha="center", va="top", fontsize=13.5, fontweight="semibold", bbox={"facecolor": "white", "alpha": 0.72, "edgecolor": "none", "pad": 0.6})

    mplhep.cms.label("Developmental", data=False, lumi=137, com=13.6, ax=axes[0])
    fig.suptitle(
        TITLES[key],
        y=0.955, fontsize=22, fontweight="normal",
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output.with_suffix(".png"), dpi=170, bbox_inches="tight")
    fig.savefig(output.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    return ratio, ratio_unc


def main() -> None:
    args = parse_args()
    if not np.isclose(args.lumi, 137.0):
        raise ValueError("These validated projection shapes must be plotted at 137 fb^-1")
    labels, paper, paper_stat, paper_total_unc = read_paper_tables(args.tables)
    local, local_stat = read_local_shapes(args.shapes)

    # In HEPData the published total should be the sum of the three reported
    # data-driven components, modulo table rounding.
    component_sum = paper["qcd"] + paper["z_nunu"] + paper["lost_lepton"]
    # The table prints independently rounded component and total values (up to
    # O(10) events in its largest bins), so equality is not expected digit for
    # digit.  A larger difference would indicate a parsing/category problem.
    component_rounding_difference = component_sum - paper["total"]
    if np.max(np.abs(component_rounding_difference)) > 100:
        raise ValueError("HEPData total background is inconsistent with its components")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = args.output_dir / "background_validation_174_bins.csv"
    ratios = {}
    with csv_path.open("w", newline="") as handle:
        fields = ["bin", "Njets", "Nbjets", "MHT", "HT"]
        for key in PAPER_COLUMNS:
            fields += [f"run3_{key}", f"run3_{key}_stat", f"run2_{key}", f"run2_{key}_stat", f"run2_{key}_total_unc", f"ratio_{key}"]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        category_ratios = {}
        for key in PAPER_COLUMNS:
            stem = args.output_dir / f"{key}_absolute_run3_mc_vs_run2_prefit"
            category_ratios[key], _ = plot_category(key, local[key], local_stat[key], paper[key], paper_stat[key], paper_total_unc[key], stem)
        for index, label in enumerate(labels):
            row = {"bin": index + 1, "Njets": label[0], "Nbjets": label[1], "MHT": label[2], "HT": label[3]}
            for key in PAPER_COLUMNS:
                row.update({
                    f"run3_{key}": local[key][index], f"run3_{key}_stat": local_stat[key][index],
                    f"run2_{key}": paper[key][index], f"run2_{key}_stat": paper_stat[key][index],
                    f"run2_{key}_total_unc": paper_total_unc[key][index], f"ratio_{key}": category_ratios[key][index],
                })
            writer.writerow(row)
            ratios[index + 1] = row

    print("Verified exact 174-bin (Njet, Nbjet, MHT, HT) ordering against HEPData.")
    print(f"Largest published-total minus rounded-component difference: {np.max(np.abs(component_rounding_difference)):.3g} events")
    print(f"Read validated Run-3 shapes: {args.shapes}")
    for key in PAPER_COLUMNS:
        finite = np.array([row[f"ratio_{key}"] for row in ratios.values()], dtype=float)
        score = np.full_like(finite, -np.inf)
        valid = np.isfinite(finite) & (finite > 0)
        score[valid] = np.abs(np.log10(finite[valid]))
        order = np.argsort(score)[::-1]
        notable = [(int(i + 1), float(finite[i])) for i in order[:5] if np.isfinite(finite[i]) and finite[i] > 0]
        print(f"{key}: Run-3 sum={local[key].sum():.6g}, Run-2 sum={paper[key].sum():.6g}; largest log-ratio bins={notable}")
    print(f"Wrote plots and values to {args.output_dir}")


if __name__ == "__main__":
    main()
