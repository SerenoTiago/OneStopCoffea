#!/usr/bin/env python3
"""Collect and plot expected Combine significances for split T2tt signals."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


POINT_RE = re.compile(r"T2tt_mStop-(?P<mstop>\d+)_mLSP-(?P<mlsp>\d+)")
SIGNIFICANCE_RE = re.compile(r"^Significance:\s*(?P<value>[0-9.eE+-]+)\s*$", re.MULTILINE)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot expected Asimov significance versus LSP mass."
    )
    parser.add_argument("combine_dir", type=Path)
    parser.add_argument("-o", "--output", type=Path, required=True)
    parser.add_argument(
        "--subtitle",
        default=r"Signal trigger efficiency assumed to be 100% for $H_T^{miss}>300$ GeV",
    )
    return parser.parse_args()


def read_points(combine_dir: Path) -> list[tuple[int, int, float, Path]]:
    points = []
    for log_path in sorted(combine_dir.glob("T2tt_mStop-*_mLSP-*/*_significance.log")):
        point_match = POINT_RE.search(log_path.parent.name)
        significance_match = SIGNIFICANCE_RE.search(log_path.read_text(errors="replace"))
        if point_match is None or significance_match is None:
            continue
        points.append(
            (
                int(point_match.group("mstop")),
                int(point_match.group("mlsp")),
                float(significance_match.group("value")),
                log_path,
            )
        )
    if not points:
        raise ValueError(f"No parseable significance logs found below {combine_dir}")
    return sorted(points, key=lambda point: (point[0], point[1]))


def main() -> int:
    args = parse_args()
    points = read_points(args.combine_dir.resolve())
    output_stem = args.output.with_suffix("")
    output_stem.parent.mkdir(parents=True, exist_ok=True)

    csv_path = output_stem.with_suffix(".csv")
    with csv_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["mstop", "mlsp", "expected_significance", "significance_log"])
        writer.writerows(points)

    stop_masses = sorted({point[0] for point in points})
    fig, ax = plt.subplots(figsize=(8.0, 5.5), constrained_layout=True)
    for mstop in stop_masses:
        selected = [point for point in points if point[0] == mstop]
        ax.plot(
            [point[1] for point in selected],
            [point[2] for point in selected],
            marker="o",
            linewidth=2,
            label=rf"$m_{{\tilde{{t}}}}={mstop}$ GeV",
        )

    max_significance = max(point[2] for point in points)
    y_max = max(1.0, 1.2 * max_significance)
    for level in (1, 2, 3, 5):
        if level > y_max:
            continue
        ax.axhline(level, color="0.75", linewidth=0.8, linestyle="--", zorder=0)
    ax.set_xlabel(r"$m_{\tilde{\chi}^{0}_{1}}$ [GeV]")
    ax.set_ylabel(r"Expected Asimov significance [$\sigma$]")
    title = "Split T2tt signal significance"
    if args.subtitle:
        title += f"\n{args.subtitle}"
    ax.set_title(title, fontsize=12, pad=10)
    ax.set_ylim(0, y_max)
    ax.grid(axis="x", color="0.9", linewidth=0.8)
    ax.legend(frameon=False)

    for extension in ("png", "pdf"):
        fig.savefig(output_stem.with_suffix(f".{extension}"), dpi=180)
    plt.close(fig)

    print(f"Wrote {len(points)} points to {csv_path}")
    print(f"Wrote {output_stem.with_suffix('.png')}")
    print(f"Wrote {output_stem.with_suffix('.pdf')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
