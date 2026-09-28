#!/usr/bin/env python3
"""Compare split-signal T2tt limits with CMS-SUS-19-006 HEPData."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import sys
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
import numpy as np
import yaml


LIMIT_COLUMNS = ("observed", "exp", "exp_m1", "exp_p1", "exp_m2", "exp_p2")
DEFAULT_TITLE = r"$pp \rightarrow \tilde{t}\tilde{t}, \tilde{t} \rightarrow t\tilde{\chi}^{0}_{1}$"

# Integrated luminosities in fb^-1. Edit these values to match the two inputs.
LOCAL_LUMINOSITY_FB = 137.0
HEPDATA_LUMINOSITY_FB = 137.0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Plot a three-panel comparison of local split-signal T2tt cross-section "
            "limits and the official CMS-SUS-19-006 Figure 14a HEPData table."
        )
    )
    parser.add_argument(
        "--limits-csv",
        type=Path,
        default=Path("analysis_products/combine/split_signals/limits.csv"),
        help="CSV produced by collect_split_signal_limits.py.",
    )
    parser.add_argument(
        "--hepdata-json",
        type=Path,
        default=Path("data/SUS-19-006"),
        help=(
            "HEPData table JSON/YAML file, or a directory containing one. The script "
            "prefers a file whose name/metadata matches T2tt cross-section limits."
        ),
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("analysis_products/combine/split_signals/t2tt_hepdata_comparison"),
        help="Output path without extension, or with .png/.pdf extension.",
    )
    parser.add_argument(
        "--quantity",
        choices=("observed", "expected"),
        default="observed",
        help="Use observed or expected limits for both inputs.",
    )
    parser.add_argument(
        "--xsec-fb",
        type=float,
        default=0.7526,
        help="Signal cross section in fb used to convert local Combine r limits to pb.",
    )
    parser.add_argument(
        "--strip-width",
        type=float,
        default=50.0,
        help="Width in GeV used when plotting each single stop-mass column.",
    )
    parser.add_argument(
        "--linear",
        action="store_true",
        help="Use a linear color scale for the cross-section heatmaps.",
    )
    parser.add_argument(
        "--paper-scale",
        action="store_true",
        default=True,
        help="Use the paper-like cross-section color scale, 1e-4 to 10 pb.",
    )
    parser.add_argument(
        "--auto-scale",
        action="store_false",
        dest="paper_scale",
        help="Scale the cross-section heatmaps to the finite values in this comparison.",
    )
    parser.add_argument("--cmap", default="rainbow", help="Matplotlib colormap name.")
    parser.add_argument("--title", default=DEFAULT_TITLE)
    parser.add_argument(
        "--write-table",
        type=Path,
        default=None,
        help="Optional CSV path for the matched comparison values.",
    )
    return parser.parse_args()


def read_rows(path: Path) -> list[dict[str, float | str]]:
    rows: list[dict[str, float | str]] = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"mstop", "mlsp", *LIMIT_COLUMNS}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path} is missing required columns: {sorted(missing)}")
        for row in reader:
            parsed: dict[str, float | str] = dict(row)
            parsed["mstop"] = float(row["mstop"])
            parsed["mlsp"] = float(row["mlsp"])
            for column in LIMIT_COLUMNS:
                parsed[column] = float(row[column]) if row[column] else math.nan
            rows.append(parsed)
    if not rows:
        raise ValueError(f"No rows found in {path}")
    return rows


def centers_to_edges(centers: np.ndarray, single_width: float) -> np.ndarray:
    centers = np.asarray(sorted(float(x) for x in centers), dtype=float)
    if len(centers) == 1:
        return np.array([centers[0] - single_width / 2.0, centers[0] + single_width / 2.0])

    midpoints = (centers[1:] + centers[:-1]) / 2.0
    first = centers[0] - (midpoints[0] - centers[0])
    last = centers[-1] + (centers[-1] - midpoints[-1])
    return np.concatenate([[first], midpoints, [last]])


def output_stem(path: Path) -> tuple[Path, list[str]]:
    if path.suffix.lower() in {".png", ".pdf"}:
        return path.with_suffix(""), [path.suffix.lower().lstrip(".")]
    return path, ["png", "pdf"]


def normalize(text: Any) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(text).lower())


def numeric_value(item: Any) -> float:
    if isinstance(item, dict):
        item = item.get("value")
    if item is None:
        return math.nan
    if isinstance(item, (int, float)):
        return float(item)
    text = str(item).strip()
    if not text:
        return math.nan
    match = re.search(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?", text)
    if not match:
        return math.nan
    return float(match.group(0))


def header_text(variable: dict[str, Any]) -> str:
    parts: list[str] = []
    header = variable.get("header") or {}
    if isinstance(header, dict):
        parts.extend(str(header.get(key, "")) for key in ("name", "units"))
    else:
        parts.append(str(header))
    for qualifier in variable.get("qualifiers") or []:
        parts.extend(str(qualifier.get(key, "")) for key in ("name", "value", "units"))
    return " ".join(parts)


def load_hepdata_file(path: Path) -> Any:
    with path.open() as handle:
        if path.suffix.lower() in {".yaml", ".yml"}:
            return yaml.safe_load(handle)
        return json.load(handle)


def score_hepdata_file(path: Path) -> tuple[int, str]:
    score = 0
    name = normalize(path.name)
    for token in ("t2tt", "crosssection", "upperlimits", "figure14a", "fig14a"):
        if token in name:
            score += 10
    try:
        data = load_hepdata_file(path)
    except Exception:
        return score, path.name
    text = normalize(json.dumps({key: data.get(key) for key in ("name", "location", "description") if isinstance(data, dict)}))
    for token in ("t2tt", "crosssection", "upperlimits", "figure14a", "datafromfigure14a"):
        if token in text:
            score += 20
    if isinstance(data, dict) and {"independent_variables", "dependent_variables"} <= set(data):
        score += 100
    return score, path.name


def find_hepdata_json(path: Path) -> Path:
    if path.is_file():
        return path
    if not path.exists():
        raise FileNotFoundError(f"HEPData path does not exist: {path}")
    candidates = sorted(
        candidate
        for pattern in ("*.json", "*.yaml", "*.yml")
        for candidate in path.rglob(pattern)
    )
    if not candidates:
        raise FileNotFoundError(
            f"No JSON/YAML files found under {path}. Put the HEPData Figure 14a "
            "'T2tt cross section upper limits' table JSON or YAML there, or pass --hepdata-json."
        )
    scored = sorted((score_hepdata_file(candidate), candidate) for candidate in candidates)
    best_score, best = scored[-1]
    if best_score[0] <= 0 and len(candidates) > 1:
        raise ValueError(
            "Could not identify the T2tt cross-section limit JSON automatically. "
            f"Pass --hepdata-json explicitly. Candidates: {', '.join(str(c) for c in candidates)}"
        )
    return best


def unwrap_hepdata_table(data: Any) -> dict[str, Any]:
    if isinstance(data, dict) and {"independent_variables", "dependent_variables"} <= set(data):
        return data
    if isinstance(data, dict) and "data_tables" in data:
        matches = []
        for table in data["data_tables"]:
            text = normalize(" ".join(str(table.get(key, "")) for key in ("name", "location", "description")))
            if all(token in text for token in ("t2tt", "crosssection")) and ("figure14a" in text or "datafromfigure14a" in text):
                matches.append(table)
        if matches:
            link = (matches[0].get("data") or {}).get("json")
            raise ValueError(
                "This looks like a HEPData record metadata JSON, not the table JSON with values. "
                "Download the table JSON for 'T2tt cross section upper limits' and pass that file. "
                f"Metadata points to: {link}"
            )
    raise ValueError("Unrecognized HEPData JSON structure; expected independent_variables/dependent_variables.")


def choose_independent_variables(table: dict[str, Any]) -> tuple[dict[str, Any], dict[str, Any]]:
    variables = table.get("independent_variables") or []
    if len(variables) < 2:
        raise ValueError("HEPData table must contain at least two independent variables for stop and LSP masses.")

    stop_var = None
    lsp_var = None
    for variable in variables:
        text = normalize(header_text(variable))
        if any(
            token in text
            for token in (
                "mstop",
                "mtop",
                "mtildet",
                "mtildeq",
                "squark",
                "msquark",
                "stop",
                "topquark",
            )
        ):
            stop_var = variable
        if any(
            token in text
            for token in (
                "mlsp",
                "mchi",
                "mtildechi",
                "mtildechi01",
                "neutralino",
                "chi01",
                "chitilde01",
                "lsp",
            )
        ):
            lsp_var = variable
    if stop_var is None:
        stop_var = variables[0]
    if lsp_var is None:
        lsp_var = variables[1]
    if stop_var is lsp_var:
        raise ValueError(f"Could not distinguish stop and LSP variables: {[header_text(v) for v in variables]}")
    return stop_var, lsp_var


def choose_dependent_variable(table: dict[str, Any], quantity: str) -> dict[str, Any]:
    variables = table.get("dependent_variables") or []
    if not variables:
        raise ValueError("HEPData table has no dependent variables.")

    want_observed = quantity == "observed"
    candidates = []
    for variable in variables:
        text = normalize(header_text(variable))
        if "upper" not in text and "limit" not in text and "xsec" not in text and "crosssection" not in text:
            continue
        if want_observed and "observed" in text:
            candidates.append(variable)
        if not want_observed and ("expected" in text or "median" in text):
            candidates.append(variable)
    if not candidates:
        candidates = [
            variable
            for variable in variables
            if (want_observed and "observed" in normalize(header_text(variable)))
            or (not want_observed and ("expected" in normalize(header_text(variable)) or "median" in normalize(header_text(variable))))
        ]
    if not candidates:
        raise ValueError(
            f"Could not find a {quantity} HEPData dependent variable. "
            f"Available variables: {[header_text(v) for v in variables]}"
        )
    return candidates[0]


def value_scale_to_pb(variable: dict[str, Any]) -> float:
    text = normalize(header_text(variable))
    if "fb" in text and "pb" not in text:
        return 1.0e-3
    return 1.0


def read_hepdata_points(path: Path, quantity: str) -> list[dict[str, float]]:
    data = load_hepdata_file(path)
    table = unwrap_hepdata_table(data)
    stop_var, lsp_var = choose_independent_variables(table)
    dep_var = choose_dependent_variable(table, quantity)

    stops = [numeric_value(value) for value in stop_var.get("values", [])]
    lsps = [numeric_value(value) for value in lsp_var.get("values", [])]
    limits = [numeric_value(value) for value in dep_var.get("values", [])]
    if not (len(stops) == len(lsps) == len(limits)):
        raise ValueError(
            "HEPData variable lengths differ: "
            f"stop={len(stops)}, lsp={len(lsps)}, limit={len(limits)}"
        )

    scale = value_scale_to_pb(dep_var)
    points = []
    for mstop, mlsp, limit in zip(stops, lsps, limits):
        if all(np.isfinite([mstop, mlsp, limit])):
            points.append({"mstop": float(mstop), "mlsp": float(mlsp), "limit_pb": float(limit) * scale})
    if not points:
        raise ValueError(f"No finite HEPData points found in {path}")
    print(f"Using HEPData file: {path}")
    print(f"HEPData stop variable: {header_text(stop_var)}")
    print(f"HEPData LSP variable: {header_text(lsp_var)}")
    print(f"HEPData limit variable: {header_text(dep_var)}")
    return points


def point_key(mstop: float, mlsp: float) -> tuple[int, int]:
    return (int(round(mstop)), int(round(mlsp)))


def median_spacing(values: list[float]) -> float:
    unique = np.array(sorted(set(float(value) for value in values)), dtype=float)
    if len(unique) < 2:
        return math.inf
    diffs = np.diff(unique)
    diffs = diffs[diffs > 0]
    if len(diffs) == 0:
        return math.inf
    return float(np.median(diffs))


def match_official_points(
    local_points: list[dict[str, float]],
    official_points: list[dict[str, float]],
) -> tuple[list[dict[str, float]], list[tuple[int, int]], str]:
    official_by_key = {point_key(point["mstop"], point["mlsp"]): point for point in official_points if point["limit_pb"] > 0}
    exact_matches: list[dict[str, float]] = []
    unmatched_local: list[tuple[int, int]] = []
    for local in local_points:
        key = point_key(local["mstop"], local["mlsp"])
        official = official_by_key.get(key)
        if official is None:
            unmatched_local.append(key)
            continue
        exact_matches.append(
            {
                "mstop": local["mstop"],
                "mlsp": local["mlsp"],
                "my_limit_pb": local["limit_pb"],
                "hepdata_limit_pb": official["limit_pb"],
                "ratio": local["limit_pb"] / official["limit_pb"],
                "hepdata_mstop": official["mstop"],
                "hepdata_mlsp": official["mlsp"],
            }
        )
    if exact_matches:
        return exact_matches, unmatched_local, "exact"

    stop_tol = 0.5001 * median_spacing([point["mstop"] for point in official_points])
    lsp_tol = 0.5001 * median_spacing([point["mlsp"] for point in official_points])
    positive_official = [point for point in official_points if point["limit_pb"] > 0]
    nearest_matches = []
    nearest_unmatched = []
    for local in local_points:
        candidates = [
            point
            for point in positive_official
            if abs(point["mstop"] - local["mstop"]) <= stop_tol
            and abs(point["mlsp"] - local["mlsp"]) <= lsp_tol
        ]
        if not candidates:
            nearest_unmatched.append(point_key(local["mstop"], local["mlsp"]))
            continue
        official = min(
            candidates,
            key=lambda point: (
                abs(point["mstop"] - local["mstop"]) / stop_tol if np.isfinite(stop_tol) else 0.0
            )
            ** 2
            + (
                abs(point["mlsp"] - local["mlsp"]) / lsp_tol if np.isfinite(lsp_tol) else 0.0
            )
            ** 2,
        )
        nearest_matches.append(
            {
                "mstop": local["mstop"],
                "mlsp": local["mlsp"],
                "my_limit_pb": local["limit_pb"],
                "hepdata_limit_pb": official["limit_pb"],
                "ratio": local["limit_pb"] / official["limit_pb"],
                "hepdata_mstop": official["mstop"],
                "hepdata_mlsp": official["mlsp"],
            }
        )
    return nearest_matches, nearest_unmatched, "nearest-bin"


def one_column_grid(points: list[dict[str, float]], value_key: str, stop_mass: float) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[dict[str, float]]]:
    filtered = [point for point in points if point_key(point["mstop"], 0)[0] == int(round(stop_mass))]
    y_centers = np.array(sorted({point["mlsp"] for point in filtered}), dtype=float)
    x_centers = np.array([stop_mass], dtype=float)
    z = np.full((len(y_centers), 1), np.nan, dtype=float)
    y_index = {y: i for i, y in enumerate(y_centers)}
    for point in filtered:
        z[y_index[point["mlsp"]], 0] = point[value_key]
    return x_centers, y_centers, z, filtered


def draw_limit_panel(
    fig: plt.Figure,
    ax: plt.Axes,
    points: list[dict[str, float]],
    value_key: str,
    stop_mass: float,
    norm: Normalize | LogNorm,
    cmap: str,
    strip_width: float,
    title: str,
    show_ylabel: bool,
) -> Any:
    x_centers, y_centers, z, filtered = one_column_grid(points, value_key, stop_mass)
    x_edges = centers_to_edges(x_centers, strip_width)
    y_edges = centers_to_edges(y_centers, strip_width)
    y_edges[0] = max(0.0, y_edges[0])

    mesh = ax.pcolormesh(x_edges, y_edges, z, cmap=cmap, norm=norm, shading="flat")

    ax.set_xlabel(r"$m_{\tilde{t}}$ [GeV]")
    if show_ylabel:
        ax.set_ylabel(r"$m_{\tilde{\chi}^{0}_{1}}$ [GeV]")
    ax.set_title(title)
    ax.set_xlim(float(np.nanmin(x_edges)), float(np.nanmax(x_edges)))
    ax.set_ylim(float(np.nanmin(y_edges)), float(np.nanmax(y_edges)))
    ax.tick_params(direction="in", top=True, right=True)

    y_min = float(np.nanmin(y_edges))
    y_max = float(np.nanmax(y_edges))
    y_margin = 0.03 * (y_max - y_min)
    for point in filtered:
        value = point[value_key]
        text_y = min(max(point["mlsp"], y_min + y_margin), y_max - y_margin)
        ax.text(stop_mass, text_y, f"{value:.2g}", ha="center", va="center", fontsize=7.5, color="black")
    return mesh


def write_comparison_table(path: Path, matched: list[dict[str, float]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["mstop", "mlsp", "my_limit_pb", "hepdata_limit_pb", "ratio", "hepdata_mstop", "hepdata_mlsp"]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in matched:
            writer.writerow({key: row[key] for key in fieldnames})


def main() -> int:
    args = parse_args()

    local_rows = read_rows(args.limits_csv)
    stop_masses = sorted({float(row["mstop"]) for row in local_rows})
    if len(stop_masses) != 1:
        raise ValueError(f"Expected exactly one stop mass in local limits, found {stop_masses}")
    stop_mass = stop_masses[0]

    local_column = "observed" if args.quantity == "observed" else "exp"
    local_points = []
    for row in local_rows:
        limit_pb = float(row[local_column]) * args.xsec_fb / 1000.0
        local_points.append({"mstop": float(row["mstop"]), "mlsp": float(row["mlsp"]), "limit_pb": limit_pb})

    local_lsp_values = sorted({point["mlsp"] for point in local_points})
    lsp_min, lsp_max = min(local_lsp_values), max(local_lsp_values)

    hepdata_path = find_hepdata_json(args.hepdata_json)
    hepdata_points_all = read_hepdata_points(hepdata_path, args.quantity)
    hepdata_points = hepdata_points_all
    if not hepdata_points:
        raise ValueError("No HEPData points found.")

    matched, unmatched_local, match_mode = match_official_points(local_points, hepdata_points)
    if not matched:
        raise ValueError("No local/HEPData mass points matched; interpolation is intentionally not used.")

    matched_local_points = [
        {"mstop": row["mstop"], "mlsp": row["mlsp"], "limit_pb": row["my_limit_pb"]}
        for row in matched
    ]
    matched_hepdata_points = [
        {"mstop": row["mstop"], "mlsp": row["mlsp"], "limit_pb": row["hepdata_limit_pb"]}
        for row in matched
    ]

    finite_values = np.array(
        [point["limit_pb"] for point in matched_local_points]
        + [point["limit_pb"] for point in matched_hepdata_points],
        dtype=float,
    )
    finite_values = finite_values[np.isfinite(finite_values) & (finite_values > 0)]
    if len(finite_values) == 0:
        raise ValueError("No positive finite cross-section limits found for plotting.")

    if args.paper_scale:
        vmin, vmax = 1.0e-4, 10.0
    else:
        vmin, vmax = float(np.nanmin(finite_values)), float(np.nanmax(finite_values))
    norm: Normalize | LogNorm = Normalize(vmin=vmin, vmax=vmax) if args.linear else LogNorm(vmin=vmin, vmax=vmax)

    fig = plt.figure(figsize=(11.5, 6.2), constrained_layout=True)
    grid = fig.add_gridspec(1, 4, width_ratios=[1.0, 1.0, 0.08, 1.35], wspace=0.08)
    ax_my = fig.add_subplot(grid[0, 0])
    ax_hep = fig.add_subplot(grid[0, 1], sharey=ax_my)
    cax = fig.add_subplot(grid[0, 2])
    ax_ratio = fig.add_subplot(grid[0, 3], sharey=ax_my)

    mesh = draw_limit_panel(
        fig,
        ax_my,
        matched_local_points,
        "limit_pb",
        stop_mass,
        norm,
        args.cmap,
        args.strip_width,
        f"My {args.quantity} limits ({LOCAL_LUMINOSITY_FB:g} "
        + r"fb$^{-1}$)"
        + "\n"
        + rf"$m_{{\tilde{{t}}}}={stop_mass:g}$ GeV",
        True,
    )
    draw_limit_panel(
        fig,
        ax_hep,
        matched_hepdata_points,
        "limit_pb",
        stop_mass,
        norm,
        args.cmap,
        args.strip_width,
        f"CMS-SUS-19-006 {args.quantity} ({HEPDATA_LUMINOSITY_FB:g} "
        + r"fb$^{-1}$)"
        + "\nFigure 14a HEPData",
        False,
    )
    fig.colorbar(mesh, cax=cax).set_label("95% CL upper limit on cross section [pb]")
    plt.setp(ax_hep.get_yticklabels(), visible=False)
    plt.setp(ax_ratio.get_yticklabels(), visible=False)

    ratio_x = np.array([row["ratio"] for row in matched], dtype=float)
    ratio_y = np.array([row["mlsp"] for row in matched], dtype=float)
    ax_ratio.axvline(1.0, color="0.35", linewidth=1.2, linestyle="--", zorder=1)
    ax_ratio.plot(ratio_x, ratio_y, "o", color="black", markersize=5.0, zorder=2)
    ax_ratio.set_xlabel("My limit / HEPData limit")
    ax_ratio.set_title("Matched mass-point ratio")
    ratio_finite = ratio_x[np.isfinite(ratio_x)]
    rmin = min(0.8, float(np.nanmin(ratio_finite)) * 0.92)
    rmax = max(1.2, float(np.nanmax(ratio_finite)) * 1.08)
    ax_ratio.set_xlim(rmin, rmax)
    ax_ratio.grid(True, axis="x", color="0.85", linewidth=0.8)
    ax_ratio.tick_params(direction="in", top=True, right=True)

    fig.suptitle(args.title, fontsize=14)

    stem, extensions = output_stem(args.output)
    stem.parent.mkdir(parents=True, exist_ok=True)
    for extension in extensions:
        fig.savefig(stem.with_suffix(f".{extension}"), dpi=200)
    plt.close(fig)

    table_path = args.write_table or stem.with_name(f"{stem.name}_values.csv")
    write_comparison_table(table_path, matched)

    print(f"Identified local stop mass: {stop_mass:g} GeV")
    print(f"Local LSP range: {lsp_min:g}-{lsp_max:g} GeV")
    print(f"Matched mass points: {len(matched)} ({match_mode})")
    if match_mode == "nearest-bin":
        print("Nearest HEPData bins used:")
        for row in matched:
            print(
                "  "
                f"local ({row['mstop']:g}, {row['mlsp']:g}) -> "
                f"HEPData ({row['hepdata_mstop']:g}, {row['hepdata_mlsp']:g})"
            )
    if unmatched_local:
        print("Local points without HEPData match:", ", ".join(f"({m},{l})" for m, l in unmatched_local))
    print(f"Wrote plot: {', '.join(str(stem.with_suffix(f'.{ext}')) for ext in extensions)}")
    print(f"Wrote values: {table_path}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
