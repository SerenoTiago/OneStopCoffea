#!/usr/bin/env python3
"""Collect per-point Combine AsymptoticLimits logs into a CSV table."""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path


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
    "observed",
    "exp",
    "exp_m1",
    "exp_p1",
    "exp_m2",
    "exp_p2",
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
    return parser.parse_args()


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

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} limit rows to {output_path}")
    if skipped:
        print(f"Skipped {skipped} logs without T2tt_mStop-*_mLSP-* in their path.")

    incomplete = [
        row
        for row in rows
        if any(row[field] == "" for field in FIELDNAMES if field not in {"limit_log"})
    ]
    if incomplete:
        print(
            f"WARNING: {len(incomplete)} rows are missing one or more limit values.",
            file=sys.stderr,
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
