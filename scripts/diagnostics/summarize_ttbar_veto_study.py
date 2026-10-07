#!/usr/bin/env python3
"""Summarize the isolated TTtoLNu2Q veto diagnostic result."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from analyzer.core.results import Histogram, UnscaledHistogram, loadResults
from analyzer.utils.structure_tools import globWithMeta


HEADLINE = (
    "baseline",
    "electron_veto_id",
    "muon_medium_no_pfiso",
    "track_pfmet",
    "all_paper_like",
)
NB_LABELS = ("0", "1", ">=2")
DECAYS = ("direct_e", "direct_mu", "tau_hadronic", "tau_leptonic", "unclassified", "ambiguous")


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--luminosity", type=float, default=137.0)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--compare-to", type=Path, help="Earlier summary directory, e.g. the 100k smoke test")
    return parser.parse_args()


def locate_result(group, name, expected_type):
    matches = []
    for item, metadata in globWithMeta(group, ["*", "*", "pipelines", "*", name]):
        if isinstance(item, expected_type):
            matches.append((item, metadata))
    if len(matches) != 1:
        raise RuntimeError(f"Expected exactly one {name}, found {len(matches)}")
    return matches[0]


def project_hist(histogram, variation, outcome="survive", decay=sum):
    view = histogram[{"veto_variation": variation, "outcome": outcome}]
    if decay is sum:
        view = view[{"w_decay": sum}]
    else:
        view = view[{"w_decay": decay}]
    return view


def values_and_variances(view):
    return np.asarray(view.values(), float), np.asarray(view.variances(), float)


def paper_lost_lepton_by_nb(path):
    with path.open() as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 174:
        raise RuntimeError(f"Expected 174 reference bins, found {len(rows)}")
    output = []
    for label in NB_LABELS:
        allowed = {label} if label != ">=2" else {"2", ">=2", ">=3"}
        output.append(sum(float(r["run2_lost_lepton"]) for r in rows if r["Nbjets"] in allowed))
    return np.asarray(output)


def main():
    args = parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    paths = sorted(args.input.glob("*.result"))
    if not paths:
        paths = sorted(args.input.rglob("*.result"))
    if len(paths) != 1:
        raise RuntimeError(f"Expected one diagnostic .result file, found {len(paths)}")
    results = loadResults([str(path) for path in paths])
    weighted, metadata = locate_result(results, "VetoDiagnosticWeightedYields", Histogram)
    raw, _ = locate_result(results, "VetoDiagnosticRawYields", UnscaledHistogram)
    tau_track_weighted, _ = locate_result(results, "TauTrackFailureReasonsWeighted", Histogram)
    tau_track_raw, _ = locate_result(results, "TauTrackFailureReasonsRaw", UnscaledHistogram)

    provenance = globWithMeta(results, ["*", "*", "_provenance"])
    if len(provenance) != 1:
        raise RuntimeError(f"Expected one provenance record, found {len(provenance)}")
    processed = provenance[0].item.chunked_events
    xsec_fb = float(metadata["x_sec"])
    scale = args.luminosity * xsec_fb / float(processed)

    rows = []
    yields = {}
    errors = {}
    for variation in weighted.histogram.axes["veto_variation"]:
        view = project_hist(weighted.histogram, variation)
        values, variances = values_and_variances(view)
        raw_values = project_hist(raw.histogram, variation).values()
        yields[variation] = values * scale
        errors[variation] = np.sqrt(variances) * scale
        for index, nb in enumerate(NB_LABELS):
            rows.append({
                "variation": variation,
                "nb": nb,
                "raw_events": int(raw_values[index]),
                "sumw": values[index],
                "sumw2": variances[index],
                "projected_137fb": yields[variation][index],
                "projected_stat_unc": errors[variation][index],
            })
    with (args.output / "veto_yields_by_nb.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    efficiency_rows = []
    baseline_view = project_hist(weighted.histogram, "baseline")
    baseline_values, baseline_variances = values_and_variances(baseline_view)
    for variation in weighted.histogram.axes["veto_variation"]:
        variation_view = project_hist(weighted.histogram, variation)
        variation_values, variation_variances = values_and_variances(variation_view)
        lost_view = project_hist(
            weighted.histogram,
            variation,
            outcome="baseline_survive_variant_fail",
        )
        _, lost_variances = values_and_variances(lost_view)
        covariance = baseline_variances - lost_variances
        for index, nb in enumerate(NB_LABELS):
            denominator = baseline_values[index]
            ratio = variation_values[index] / denominator
            ratio_variance = (
                variation_variances[index] / denominator**2
                + variation_values[index] ** 2 * baseline_variances[index] / denominator**4
                - 2.0 * variation_values[index] * covariance[index] / denominator**3
            )
            efficiency_rows.append({
                "variation": variation,
                "nb": nb,
                "yield_relative_to_baseline": ratio,
                "paired_stat_unc": np.sqrt(max(ratio_variance, 0.0)),
                "fractional_reduction": 1.0 - ratio,
            })
    with (args.output / "veto_efficiencies_relative_to_baseline.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=efficiency_rows[0].keys())
        writer.writeheader()
        writer.writerows(efficiency_rows)

    if args.compare_to is not None:
        previous_path = args.compare_to / "veto_efficiencies_relative_to_baseline.csv"
        with previous_path.open() as handle:
            previous = {
                (row["variation"], row["nb"]): row
                for row in csv.DictReader(handle)
            }
        comparison_rows = []
        for row in efficiency_rows:
            old = previous[(row["variation"], row["nb"])]
            comparison_rows.append({
                "variation": row["variation"],
                "nb": row["nb"],
                "previous_efficiency": float(old["yield_relative_to_baseline"]),
                "previous_stat_unc": float(old["paired_stat_unc"]),
                "current_efficiency": row["yield_relative_to_baseline"],
                "current_stat_unc": row["paired_stat_unc"],
            })
        with (args.output / "comparison_100k_vs_1m.csv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=comparison_rows[0].keys())
            writer.writeheader()
            writer.writerows(comparison_rows)

    fig, ax = plt.subplots(figsize=(9.2, 5.6))
    offsets = np.linspace(-0.28, 0.28, len(HEADLINE))
    for offset, variation in zip(offsets, HEADLINE):
        selected = [row for row in efficiency_rows if row["variation"] == variation]
        ax.errorbar(
            np.arange(3) + offset,
            [row["yield_relative_to_baseline"] for row in selected],
            yerr=[row["paired_stat_unc"] for row in selected],
            marker="o",
            capsize=2,
            label=variation,
        )
    ax.axhline(1, color="black", linestyle="--", linewidth=1)
    ax.set_xticks(np.arange(3), NB_LABELS)
    ax.set_xlabel(r"$N_b$ category")
    ax.set_ylabel("Yield relative to baseline")
    ax.set_title("Paired selection efficiencies")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(args.output / "veto_efficiencies_relative_to_baseline.png", dpi=180)
    fig.savefig(args.output / "veto_efficiencies_relative_to_baseline.pdf")
    plt.close(fig)

    migrations = []
    for variation in weighted.histogram.axes["veto_variation"]:
        for outcome in ("baseline_survive_variant_fail", "variant_survive_baseline_fail"):
            for decay in DECAYS:
                view = project_hist(weighted.histogram, variation, outcome=outcome, decay=decay)
                values, variances = values_and_variances(view)
                raw_values = project_hist(raw.histogram, variation, outcome=outcome, decay=decay).values()
                for index, nb in enumerate(NB_LABELS):
                    migrations.append({
                        "variation": variation,
                        "direction": outcome,
                        "nb": nb,
                        "w_decay": decay,
                        "raw_events": int(raw_values[index]),
                        "sumw": values[index],
                        "sumw2": variances[index],
                    })
    with (args.output / "baseline_migrations_by_nb_and_decay.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=migrations[0].keys())
        writer.writeheader()
        writer.writerows(migrations)

    lost_variations = HEADLINE[1:]
    fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharey=True)
    for nb_index, (ax, nb_label) in enumerate(zip(axes, NB_LABELS)):
        bottom = np.zeros(len(lost_variations))
        for decay in DECAYS:
            vals = []
            for variation in lost_variations:
                view = project_hist(
                    raw.histogram,
                    variation,
                    outcome="baseline_survive_variant_fail",
                    decay=decay,
                )
                vals.append(view.values()[nb_index])
            vals = np.asarray(vals)
            ax.bar(np.arange(len(lost_variations)), vals, bottom=bottom, label=decay)
            bottom += vals
        ax.set_title(rf"$N_b={nb_label}$")
        ax.set_xticks(np.arange(len(lost_variations)), lost_variations, rotation=55, ha="right", fontsize=8)
    axes[0].set_ylabel("Baseline survivors rejected (raw events)")
    axes[-1].legend(frameon=False, fontsize=8, bbox_to_anchor=(1.02, 1))
    fig.suptitle("Veto migrations by generator-level W decay")
    fig.tight_layout()
    fig.savefig(args.output / "baseline_rejected_migrations_by_decay.png", dpi=180)
    fig.savefig(args.output / "baseline_rejected_migrations_by_decay.pdf")
    plt.close(fig)

    track_rows = [
        row for row in migrations
        if row["variation"] == "track_pfmet"
    ]
    with (args.output / "isolated_track_bidirectional_migrations.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=track_rows[0].keys())
        writer.writeheader()
        writer.writerows(track_rows)

    tau_rows = []
    for nb_index, nb_label in enumerate(NB_LABELS):
        for reason in tau_track_raw.histogram.axes["reason"]:
            raw_value = tau_track_raw.histogram[{"nb": nb_label, "reason": reason}]
            weighted_value = tau_track_weighted.histogram[{"nb": nb_label, "reason": reason}]
            tau_rows.append({
                "nb": nb_label,
                "reason": reason,
                "raw_events": int(raw_value),
                "sumw": weighted_value.value,
                "sumw2": weighted_value.variance,
            })
    with (args.output / "hadronic_tau_isotrack_failure_reasons.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=tau_rows[0].keys())
        writer.writeheader()
        writer.writerows(tau_rows)

    reasons = list(tau_track_raw.histogram.axes["reason"])
    fig, ax = plt.subplots(figsize=(9.6, 5.8))
    bottom = np.zeros(3)
    for reason in reasons:
        values = np.asarray([
            tau_track_raw.histogram[{"nb": nb_label, "reason": reason}]
            for nb_label in NB_LABELS
        ], dtype=float)
        ax.bar(np.arange(3), values, bottom=bottom, label=reason)
        bottom += values
    ax.set_xticks(np.arange(3), NB_LABELS)
    ax.set_xlabel(r"$N_b$ category")
    ax.set_ylabel("Baseline-surviving hadronic-tau events")
    ax.set_title("Event-level IsoTrack failure reason (not gen-tau matched)")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(args.output / "hadronic_tau_isotrack_failure_reasons.png", dpi=180)
    fig.savefig(args.output / "hadronic_tau_isotrack_failure_reasons.pdf")
    plt.close(fig)

    x = np.arange(len(NB_LABELS))
    fig, ax = plt.subplots(figsize=(9.2, 5.8))
    offsets = np.linspace(-0.28, 0.28, len(HEADLINE))
    for offset, variation in zip(offsets, HEADLINE):
        ax.errorbar(x + offset, yields[variation], yerr=errors[variation], marker="o", capsize=2, label=variation)
    ax.set_yscale("log")
    ax.set_xticks(x, NB_LABELS)
    ax.set_xlabel(r"$N_b$ category")
    ax.set_ylabel(r"Projected TTtoLNu2Q yield at 137 fb$^{-1}$")
    ax.set_title("100k-event technical smoke test; not a precision prediction")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(args.output / "veto_yields_by_nb.png", dpi=180)
    fig.savefig(args.output / "veto_yields_by_nb.pdf")
    plt.close(fig)

    fig, axes = plt.subplots(1, len(NB_LABELS), figsize=(15, 5), sharey=True)
    for index, (ax, nb) in enumerate(zip(axes, NB_LABELS)):
        bottom = np.zeros(len(HEADLINE))
        for decay in DECAYS:
            vals = []
            for variation in HEADLINE:
                view = project_hist(weighted.histogram, variation, decay=decay)
                vals.append(view.values()[index])
            vals = np.asarray(vals)
            ax.bar(np.arange(len(HEADLINE)), vals, bottom=bottom, label=decay)
            bottom += vals
        ax.set_title(rf"$N_b={nb}$")
        ax.set_xticks(np.arange(len(HEADLINE)), HEADLINE, rotation=55, ha="right", fontsize=8)
    axes[0].set_ylabel("Signed weighted events in 100k sample")
    axes[-1].legend(frameon=False, fontsize=8, bbox_to_anchor=(1.02, 1))
    fig.suptitle("Generator-level leptonic-W composition of surviving events")
    fig.tight_layout()
    fig.savefig(args.output / "survivor_w_decay_composition.png", dpi=180)
    fig.savefig(args.output / "survivor_w_decay_composition.pdf")
    plt.close(fig)

    paper = paper_lost_lepton_by_nb(args.reference)
    ratio_rows = []
    fig, ax = plt.subplots(figsize=(9.2, 5.4))
    for variation in HEADLINE:
        ratio = yields[variation] / paper
        ratio_unc = errors[variation] / paper
        ax.errorbar(x, ratio, yerr=ratio_unc, marker="o", capsize=2, label=variation)
        for index, nb in enumerate(NB_LABELS):
            ratio_rows.append({"variation": variation, "nb": nb, "ttbar_over_paper_lost_lepton": ratio[index], "stat_unc": ratio_unc[index]})
    ax.axhline(1, color="black", linestyle="--", linewidth=1)
    ax.set_xticks(x, NB_LABELS)
    ax.set_xlabel(r"$N_b$ category")
    ax.set_ylabel("TTtoLNu2Q diagnostic / paper total lost-lepton")
    ax.set_title("Non-like-for-like smoke-test ratio; statistical errors only")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(args.output / "ttbar_over_paper_lost_lepton_by_nb.png", dpi=180)
    fig.savefig(args.output / "ttbar_over_paper_lost_lepton_by_nb.pdf")
    plt.close(fig)
    with (args.output / "ttbar_over_paper_lost_lepton_by_nb.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=ratio_rows[0].keys())
        writer.writeheader()
        writer.writerows(ratio_rows)

    summary = {
        "input_result": str(paths[0]),
        "processed_events": processed,
        "luminosity_fb": args.luminosity,
        "cross_section_fb": xsec_fb,
        "normalization_scale": scale,
        "warning": (
            "Technical smoke test; not statistically representative"
            if processed < 500000
            else "Finite diagnostic subsample; diagnostic extrapolations are not full-background predictions"
        ),
    }
    (args.output / "run_summary.json").write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
