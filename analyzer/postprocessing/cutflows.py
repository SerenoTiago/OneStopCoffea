from __future__ import annotations

from pathlib import Path
import functools as ft
from typing import Literal
from .style import StyleSet
from analyzer.utils.structure_tools import (
    commonDict,
    dictToDot,
    dotFormat,
)
from .processors import BasePostprocessor
from .plots.plots_1d import plotDictAsBars
from attrs import define, field


def _getCutflow(x, count_type="weighted"):
    if count_type == "weighted" and getattr(x, "weighted_cutflow", None) is not None:
        return x.weighted_cutflow
    return x.cutflow


def _getRelativeCutflow(x, count_type="weighted"):
    cutflow = _getCutflow(x, count_type=count_type)
    initial = cutflow["initial"]
    if initial == 0:
        return {cut: 0 for cut in cutflow}
    return {cut: value / initial for cut, value in cutflow.items()}


@define
class PlotSelectionFlow(BasePostprocessor):
    output_name: str
    style_set: str | StyleSet = field(factory=StyleSet)
    scale: Literal["log", "linear"] = "linear"
    normalize: bool = False
    count_type: Literal["weighted", "raw"] = "weighted"
    relative_to_initial: bool = False
    show_errors: bool = True

    def getRunFuncs(self, group, prefix=None):
        common_meta = commonDict(group)
        output_path = dotFormat(
            self.output_name, **dict(dictToDot(common_meta)), prefix=prefix
        )
        pc = self.plot_configuration.makeFormatted(common_meta)
        getter = _getRelativeCutflow if self.relative_to_initial else _getCutflow

        yield ft.partial(
            plotDictAsBars,
            group,
            common_meta,
            output_path,
            getter=ft.partial(getter, count_type=self.count_type),
            style_set=self.style_set,
            normalize=self.normalize,
            scale=self.scale,
            show_errors=self.show_errors,
            plot_configuration=pc,
        )


ALLOWED_COLS = Literal["count", "rel", "abs"]


@define
class CutflowTable(BasePostprocessor):
    output_name: str
    format: Literal["markdown", "csv", "latex"] = "csv"
    key: str = "{dataset_name}"
    standalone: bool = False
    highlight_rows: list[tuple[int, str]] | None = None
    count_type: Literal["weighted", "raw", "both"] = "both"
    cols: list[ALLOWED_COLS] = ["count", "rel", "abs"]

    def getRunFuncs(self, group, prefix=None):
        common_meta = commonDict(group)
        output_path = dotFormat(
            self.output_name, **dict(dictToDot(common_meta)), prefix=prefix
        )

        yield ft.partial(
            makeAndSaveCutflowTable,
            group,
            common_meta,
            output_path,
            format=self.format,
            key=self.key,
            standalone=self.standalone,
            highlight_rows=self.highlight_rows,
            count_type=self.count_type,
            cols=self.cols,
        )


def makeCutflowDf(group, key="{dataset_name}", count_type="both", cols=None):
    import pandas as pd
    import numpy as np

    cols = cols or ["count", "rel", "abs"]
    dataset_cutflows = {}
    dataset_raw_cutflows = {}
    cut_order = None
    for selection_flow, metadata in group:
        k = dotFormat(key, **dict(dictToDot(metadata)))
        dataset_cutflows[k] = _getCutflow(selection_flow, "weighted")
        dataset_raw_cutflows[k] = _getCutflow(selection_flow, "raw")
        if cut_order is None:
            cut_order = list(selection_flow.cuts)
        else:
            if cut_order != list(selection_flow.cuts):
                raise ValueError("Cutflows are not consistent across datasets.")
    all_data = {}
    for dataset_name, cutflow in dataset_cutflows.items():
        if count_type in ("weighted", "both"):
            all_data[dataset_name, "Yield"] = cutflow
        if count_type in ("raw", "both"):
            all_data[dataset_name, "Raw"] = dataset_raw_cutflows[dataset_name]

    df = pd.DataFrame(all_data)
    for dataset_name in dataset_cutflows:
        eff_source = "Yield" if (dataset_name, "Yield") in df.columns else "Raw"
        values = df.loc[:, (dataset_name, eff_source)]
        values_array = values.to_numpy(dtype=float)
        previous_array = values.shift(1).to_numpy(dtype=float)
        if "abs" in cols:
            df.loc[:, (dataset_name, "Eff. Abs.")] = np.divide(
                values_array,
                values_array[0],
                out=np.zeros_like(values_array, dtype=float),
                where=values_array[0] != 0,
            )
        if "rel" in cols:
            df.loc[:, (dataset_name, "Eff. Rel.")] = np.divide(
                values_array,
                previous_array,
                out=np.zeros_like(values_array, dtype=float),
                where=previous_array != 0,
            )
            df.loc[df.index[0], (dataset_name, "Eff. Rel.")] = 1.0
    df.sort_index(axis=1, level=[0, 1], ascending=[True, False], inplace=True)
    return df


STANDALONE_TOP = r"""\documentclass{standalone}
\usepackage{booktabs}
\usepackage[table,usenames,svgnames]{xcolor}
\begin{document}
"""

STANDALONE_BOTTOM = r"""
\end{document}
"""


def makeAndSaveCutflowTable(
    group,
    common_meta,
    output_path,
    format="csv",
    key="{dataset_name}",
    standalone=False,
    highlight_rows=None,
    count_type="both",
    cols=None,
):
    import numpy as np

    cols = cols or ["count", "rel", "abs"]

    highlight_rows = highlight_rows or []

    df = makeCutflowDf(group, key=key, count_type=count_type, cols=cols)
    output_path = Path(output_path)
    output_path.parent.mkdir(exist_ok=True, parents=True)

    s = (
        df.style.apply(
            lambda x: np.where(
                (np.arange(len(x)) % (2 * len(cols))) >= len(cols),
                "background-color: lightgray",
                "",
            ),
            axis=1,
        )
        .format("{:0.2f}", escape="latex")
        .format_index(escape="latex", axis=0)
        # .format_index(escape="latex",axis=1)
    )
    for row, color in highlight_rows:
        s = s.apply(
            lambda x: np.where(
                (np.arange(len(x)) == row), f"background-color: {color}", ""
            ),
            axis=0,
        )

    if format == "csv":
        df.to_csv(output_path)
    elif format == "markdown":
        s.to_markdown(output_path, convert_css=True)
    elif format == "latex":
        text = s.to_latex(None, convert_css=True)
        if standalone:
            text = STANDALONE_TOP + text + STANDALONE_BOTTOM

        with open(output_path, "w") as f:
            f.write(text)
