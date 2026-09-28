#!/usr/bin/env python3
"""Sum the signal histogram across multiple per-mass-point shapes ROOT files
into a single "signal" histogram, alongside one copy of the (identical)
background histograms, for a combined-signal 174-bin diagnostic plot.

Each per-point shapes_bins_hadsusy.root only contains that one mass point's
signal (see make_split_signal_datacards.py), so the background histograms are
identical across all points and are copied from the first input file.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import uproot

DEFAULT_BACKGROUNDS = (
    "qcd_inclusive_2024",
    "tt_hadronic_2024",
    "tt_semileptonic_2024",
    "wjets_2024",
    "zjets_2024",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "shapes", type=Path, nargs="+", help="Per-mass-point shapes_bins_hadsusy.root files"
    )
    parser.add_argument(
        "--signal-histogram",
        default="2stop_1250_2024_splitLSP",
        help="Name of the signal histogram to sum across input files",
    )
    parser.add_argument(
        "--background",
        action="append",
        dest="backgrounds",
        help="Background histogram name to copy from the first file; repeat as needed",
    )
    parser.add_argument("-o", "--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    backgrounds = args.backgrounds or list(DEFAULT_BACKGROUNDS)

    summed_signal = None
    for path in args.shapes:
        with uproot.open(path) as f:
            hist = f[args.signal_histogram].to_hist()
            summed_signal = hist.copy() if summed_signal is None else summed_signal + hist

    with uproot.open(args.shapes[0]) as f:
        available = {k.split(";")[0] for k in f.keys()}
        background_hists = {
            name: f[name].to_hist() for name in backgrounds if name in available
        }

    with uproot.recreate(args.output) as out:
        out["signal"] = summed_signal
        for name, hist in background_hists.items():
            out[name] = hist

    print(
        f"Wrote {args.output}: signal summed over {len(args.shapes)} mass points, "
        f"backgrounds {sorted(background_hists)}"
    )


if __name__ == "__main__":
    main()
