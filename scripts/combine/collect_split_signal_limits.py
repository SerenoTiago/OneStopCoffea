#!/usr/bin/env python3
"""Collect per-point Combine AsymptoticLimits logs into a CSV table."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import sys
from pathlib import Path

import yaml


# Written into each point directory by make_split_signal_datacards.py.
NORMALIZATION_FILE = "signal_normalization.json"
POINT_RE = re.compile(r"T2tt_mStop-(?P<mstop>\d+)_mLSP-(?P<mlsp>\d+)")
OBSERVED_RE = re.compile(r"Observed Limit:\s*r\s*<\s*(?P<value>[0-9.eE+-]+)")
EXPECTED_RE = re.compile(
    r"Expected\s+(?P<quantile>2\.5|16\.0|50\.0|84\.0|97\.5)%:\s*"
    r"r\s*<\s*(?P<value>[0-9.eE+-]+)"
)

QUANTILE_TO_FIELD = {
    "2.5": "exp_m2",
    "16.0": "exp_m1",
    "50.0": "exp",
    "84.0": "exp_p1",
    "97.5": "exp_p2",
}

FIELDNAMES = [
    "mstop",
    "mlsp",
    "xsec_fb",
    "observed",
    "exp",
    "exp_m1",
    "exp_p1",
    "exp_m2",
    "exp_p2",
    "exp_limit_pb",
    "exp_m1_limit_pb",
    "exp_p1_limit_pb",
    "exp_m2_limit_pb",
    "exp_p2_limit_pb",
    "limit_log",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Parse Combine *_limits.log files for split T2tt signals."
    )
    parser.add_argument(
        "combine_dir",
        type=Path,
        help="Root directory containing T2tt_mStop-*_mLSP-* combine output directories.",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=Path,
        default=Path("analysis_products/combine/split_signals/limits.csv"),
        help="CSV output path.",
    )
    parser.add_argument(
        "--log-glob",
        default="**/*_limits.log",
        help="Glob, relative to combine_dir, used to find Combine AsymptoticLimits logs.",
    )
    parser.add_argument(
        "--expected-only",
        action="store_true",
        help="Accept blinded logs with no observed limit and validate expected quantiles only.",
    )
    parser.add_argument(
        "--xsec-yaml-dir",
        type=Path,
        default=Path("analyzer_resources/datasets/2024"),
        help="Dataset YAML directory used to add mass-dependent signal cross sections.",
    )
    return parser.parse_args()


def cross_section_fb(yaml_dir: Path, mstop: int) -> float:
    path = yaml_dir / f"2stop_{mstop}_2024_splitLSP.yaml"
    datasets = yaml.safe_load(path.read_text())
    values = {
        float(sample["x_sec"])
        for dataset in datasets
        for sample in dataset.get("samples", [])
    }
    if len(values) != 1:
        raise ValueError(f"Expected one cross section in {path}, found {sorted(values)}")
    return values.pop()


def point_from_path(path: Path) -> tuple[int, int] | None:
    for part in reversed(path.parts):
        match = POINT_RE.search(part)
        if match:
            return int(match.group("mstop")), int(match.group("mlsp"))
    return None


def parse_limit_log(path: Path) -> dict[str, str | int] | None:
    point = point_from_path(path)
    if point is None:
        return None

    row: dict[str, str | int] = {
        "mstop": point[0],
        "mlsp": point[1],
        "observed": "",
        "exp": "",
        "exp_m1": "",
        "exp_p1": "",
        "exp_m2": "",
        "exp_p2": "",
        "limit_log": str(path),
    }

    for line in path.read_text(errors="replace").splitlines():
        observed = OBSERVED_RE.search(line)
        if observed:
            row["observed"] = observed.group("value")
            continue

        expected = EXPECTED_RE.search(line)
        if expected:
            row[QUANTILE_TO_FIELD[expected.group("quantile")]] = expected.group(
                "value"
            )

    return row


def main() -> int:
    args = parse_args()
    combine_dir = args.combine_dir.resolve()
    output_path = args.output.resolve()

    if not combine_dir.is_dir():
        print(f"ERROR: combine_dir does not exist: {combine_dir}", file=sys.stderr)
        return 2

    rows = []
    skipped = 0
    for path in sorted(combine_dir.glob(args.log_glob)):
        row = parse_limit_log(path)
        if row is None:
            skipped += 1
            continue
        rows.append(row)

    rows.sort(key=lambda row: (int(row["mstop"]), int(row["mlsp"])))

    if not rows:
        print(
            f"ERROR: no parseable limit logs found under {combine_dir}",
            file=sys.stderr,
        )
        return 2

    # Combine r values are relative to the cross section stored in the signal
    # .result files when they were produced.  Converting r to pb with a
    # different YAML value would silently rescale the limits, so refuse.
    problems = []
    for row in rows:
        record_path = Path(str(row["limit_log"])).parent / NORMALIZATION_FILE
        xsec_fb = cross_section_fb(args.xsec_yaml_dir, int(row["mstop"]))
        if not record_path.is_file():
            problems.append(f"{record_path}: missing")
            continue
        used_fb = float(json.loads(record_path.read_text())["x_sec_fb"])
        if not math.isclose(used_fb, xsec_fb, rel_tol=1e-6):
            problems.append(
                f"{record_path}: inputs normalized to {used_fb:g} fb, "
                f"dataset YAML says {xsec_fb:g} fb"
            )
    if problems:
        print(
            "ERROR: signal normalization does not match the dataset YAMLs for "
            f"{len(problems)} points; regenerate their datacards from results "
            "produced with the current cross sections:",
            file=sys.stderr,
        )
        for problem in problems:
            print(f"  {problem}", file=sys.stderr)
        return 2

    for row in rows:
        xsec_fb = cross_section_fb(args.xsec_yaml_dir, int(row["mstop"]))
        row["xsec_fb"] = f"{xsec_fb:.8g}"
        for source, target in (
            ("exp", "exp_limit_pb"),
            ("exp_m1", "exp_m1_limit_pb"),
            ("exp_p1", "exp_p1_limit_pb"),
            ("exp_m2", "exp_m2_limit_pb"),
            ("exp_p2", "exp_p2_limit_pb"),
        ):
            row[target] = (
                f"{float(row[source]) * xsec_fb / 1000.0:.10g}"
                if row[source] != ""
                else ""
            )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} limit rows to {output_path}")
    if skipped:
        print(f"Skipped {skipped} logs without T2tt_mStop-*_mLSP-* in their path.")

    required_fields = ["exp", "exp_m1", "exp_p1", "exp_m2", "exp_p2"]
    if not args.expected_only:
        required_fields.append("observed")
    incomplete = [
        row
        for row in rows
        if any(row[field] == "" for field in required_fields)
    ]
    if incomplete:
        print(
            f"WARNING: {len(incomplete)} rows are missing one or more limit values.",
            file=sys.stderr,
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
