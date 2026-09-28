from __future__ import annotations

from pathlib import Path
from fnmatch import fnmatch
import functools as ft

import hist
import matplotlib.pyplot as plt
import mplhep
import numpy as np
from attrs import define, field

from analyzer.modules.common.hadronic_susy import KINEMATIC_BINS, NJET_BINS
from analyzer.postprocessing.plots.common import PlotConfiguration
from analyzer.postprocessing.plots.utils import addLegend, saveFigVariants
from analyzer.postprocessing.processors import BasePostprocessor
from analyzer.postprocessing.style import StyleSet, Styler, cms_colors
from analyzer.utils.structure_tools import commonDict, dictToDot, dotFormat


def _is_data(meta):
    sample_type = meta.get("sample_type")
    value = sample_type.value if hasattr(sample_type, "value") else str(sample_type)
    return value == "Data"


def _matches_any(value, patterns):
    return any(fnmatch(value, pattern) for pattern in patterns)


def _hist_values(item):
    h = item.histogram
    if "variation" in h.axes.name:
        if "central" in h.axes["variation"]:
            h = h[{"variation": hist.loc("central")}]
        elif len(h.axes["variation"]) == 0:
            return np.zeros(174), np.zeros(174)
        else:
            h = h[{"variation": 0}]
    values = np.asarray(h.values(flow=False), dtype=float)
    variances = np.asarray(h.variances(flow=False), dtype=float)
    if values.size != 174:
        raise ValueError(f"Expected 174 search bins, got {values.size}.")
    return values, variances


def _search_bin_boundaries():
    major = []
    minor = []
    labels = []
    b_labels = []
    offset = 0
    for njet_idx, (njet_low, njet_high) in enumerate(NJET_BINS):
        start = offset
        n_b_bins = 3 if njet_idx == 0 else 4
        n_kin = len(KINEMATIC_BINS) if njet_idx < 3 else len(KINEMATIC_BINS) - 2
        for b_idx in range(n_b_bins):
            b_start = offset
            offset += n_kin
            if b_idx < n_b_bins - 1:
                minor.append(offset + 0.5)
            b_text = str(b_idx) if b_idx < n_b_bins - 1 else rf"$\geq {b_idx}$"
            b_labels.append((b_start + (n_kin + 1) / 2, b_text))
        if njet_idx < len(NJET_BINS) - 1:
            major.append(offset + 0.5)
        njet_text = (
            rf"${njet_low} \leq N_{{\mathrm{{jet}}}}$"
            if njet_high is None
            else rf"${njet_low} \leq N_{{\mathrm{{jet}}}} \leq {njet_high}$"
        )
        labels.append((start + (offset - start + 1) / 2, njet_text))
    return major, minor, labels, b_labels


def _draw_search_bin_separators(ax, ratio_ax=None):
    major, minor, njet_labels, b_labels = _search_bin_boundaries()
    axes = (ax,) if ratio_ax is None else (ax, ratio_ax)
    for axis in axes:
        for x in minor:
            axis.axvline(x, color="black", linestyle=":", linewidth=1.0, zorder=0)
        for x in major:
            axis.axvline(x, color="black", linestyle=(0, (4, 4)), linewidth=1.0, zorder=0)

    y_top = ax.get_ylim()[1]
    for x, text in njet_labels:
        ax.text(x, y_top / 2.0, text, ha="center", va="top", fontsize=12)
    ax.text(4.0, y_top / 10.0, r"$N_{\mathrm{b-jet}}$", ha="left", va="top", fontsize=10)
    for x, text in b_labels:
        ax.text(x, y_top / 35.0, text, ha="center", va="top", fontsize=10)


@define
class PaperSearchBinYield(BasePostprocessor):
    output_name: str
    signal_dataset_patterns: list[str] = field(factory=lambda: ["2stop*", "signal*"])
    ratio_ylim: tuple[float, float] = (-1.0, 3.0)
    ratio_height: float = 0.3
    style_set: str | StyleSet = field(factory=StyleSet)

    def getRunFuncs(self, group, prefix=None):
        common_meta = commonDict(group)
        output_path = dotFormat(
            self.output_name, **dict(dictToDot(common_meta)), prefix=prefix
        )
        pc = self.plot_configuration.makeFormatted(common_meta)
        yield ft.partial(
            plotPaperSearchBinYield,
            group,
            common_meta,
            output_path,
            self.style_set,
            self.signal_dataset_patterns,
            self.ratio_ylim,
            self.ratio_height,
            pc,
        )


def plotPaperSearchBinYield(
    group,
    common_meta,
    output_path,
    style_set,
    signal_dataset_patterns,
    ratio_ylim,
    ratio_height,
    plot_configuration=None,
):
    pc = plot_configuration or PlotConfiguration()
    styler = Styler(style_set)
    backgrounds = []
    signals = []
    data = []
    for item, meta in group:
        dataset_name = meta["dataset_name"]
        values, variances = _hist_values(item)
        entry = (values, variances, meta)
        if _is_data(meta):
            data.append(entry)
        elif _matches_any(dataset_name, signal_dataset_patterns):
            signals.append(entry)
        else:
            backgrounds.append(entry)

    if not backgrounds:
        raise ValueError("PaperSearchBinYield needs at least one background histogram.")

    fig, ax = plt.subplots(figsize=(14, 6))
    x = np.arange(1, 175)
    edges = np.arange(0.5, 175.5, 1.0)

    cumulative = np.zeros(174)
    total_var = np.zeros(174)
    ordered_backgrounds = sorted(backgrounds, key=lambda entry: np.sum(entry[0]))
    paper_colors = [
        cms_colors["cms-orange-light"],
        cms_colors["cms-cyan"],
        cms_colors["cms-red"],
        cms_colors["cms-purple"],
        cms_colors["cms-gray"],
    ]
    for idx, (values, variances, meta) in enumerate(ordered_backgrounds):
        top = cumulative + values
        label = meta.get("title") or meta["dataset_title"]
        color = paper_colors[idx % len(paper_colors)]
        ax.stairs(
            top,
            edges,
            baseline=cumulative,
            fill=True,
            label=label,
            color=color,
            edgecolor="black",
            linewidth=0.9,
        )
        cumulative = top
        total_var += variances

    total_bkg = cumulative
    total_unc = np.sqrt(total_var)
    positive = total_bkg > 0
    if np.any(positive):
        ax.bar(
            x[positive],
            2 * total_unc[positive],
            bottom=np.maximum(total_bkg[positive] - total_unc[positive], 1e-12),
            width=1.0,
            fill=False,
            edgecolor="gray",
            hatch="////",
            linewidth=0.0,
            label="Bkg. stat. unc.",
        )

    observed = None
    observed_var = None
    if data:
        observed = np.sum([entry[0] for entry in data], axis=0)
        observed_var = np.sum([entry[1] for entry in data], axis=0)
        ax.errorbar(
            x,
            observed,
            yerr=np.sqrt(observed_var),
            color="black",
            linestyle="none",
            marker="o",
            markersize=3,
            label="Data",
            zorder=5,
        )
    elif signals:
        observed = np.sum([entry[0] for entry in signals], axis=0)
        observed_var = np.sum([entry[1] for entry in signals], axis=0)
        for values, variances, meta in signals:
            style = styler.getStyle(meta)
            ax.errorbar(
                x,
                values,
                yerr=np.sqrt(variances),
                label=meta.get("title") or meta["dataset_title"],
                linestyle="-",
                drawstyle="steps-mid",
                marker=None,
                color=style.color,
                linewidth=1.4,
                zorder=4,
            )

    ax.set_yscale("log")
    nonzero_values = total_bkg[total_bkg > 0]
    if observed is not None:
        nonzero_values = np.concatenate([nonzero_values, observed[observed > 0]])
    ax.set_ylim(max(np.min(nonzero_values) * 0.05, 1e-3), np.max(total_bkg) * 80.0)
    ax.set_ylabel("Events")
    ax.set_xlabel("Search region bin number")
    ax.set_xlim(0.5, 174.5)
    ax.set_xticks(np.arange(20, 181, 20))
    addLegend(ax, pc)
    _draw_search_bin_separators(ax)

    all_meta = [meta for _, _, meta in backgrounds + signals + data]
    saveFigVariants(
        fig,
        ax,
        output_path,
        all_meta,
        plot_configuration=pc,
        metadata=common_meta,
    )
    plt.close(fig)
