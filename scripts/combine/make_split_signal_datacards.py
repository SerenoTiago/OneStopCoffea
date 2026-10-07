#!/usr/bin/env python3
"""Create one CombineDatacard postprocessing job per split T2tt signal.

Run this from the repository after entering the normal analysis environment.
The script intentionally does not manage .venv or call uv directly.
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path


SIGNAL_RE = re.compile(
    r"2stop_(?P<dataset_stop>\d+)_\d+_splitLSP__signal_T2tt_"
    r"(?P<mstop>\d+)_(?P<mlsp>\d+)\.result$"
)

DEFAULT_BACKGROUNDS = (
    "qcd_inclusive_2024",
    "tt_hadronic_2024",
    "tt_semileptonic_2024",
    "wjets_2024",
    "zjets_2024",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Generate and optionally run per-signal CombineDatacard "
            "postprocessing configs for split T2tt signal results."
        )
    )
    parser.add_argument(
        "results_dir",
        type=Path,
        help="Directory containing background .result files and split signal .result files.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("analysis_products/combine/split_signals"),
        help="Directory where datacard directories will be written.",
    )
    parser.add_argument(
        "--config-dir",
        type=Path,
        default=Path(".application_data/generated_postprocess/split_signals"),
        help="Directory for generated per-signal postprocessing YAML files.",
    )
    parser.add_argument(
        "--background",
        action="append",
        dest="backgrounds",
        default=[],
        help=(
            "Background dataset_name to include. May be repeated. "
            "Defaults to the current hadronic SUSY 2024 backgrounds."
        ),
    )
    parser.add_argument(
        "--signal-glob",
        default="2stop_*_splitLSP__signal_T2tt_*.result",
        help="Glob used to discover split signal result files inside results_dir.",
    )
    parser.add_argument(
        "--inputs-glob",
        default="*.result",
        help="Glob of result files to pass to ./osca postprocess.",
    )
    parser.add_argument(
        "--channel",
        default="bins_hadsusy",
        help="Combine channel name used by the CombineDatacard postprocessor.",
    )
    parser.add_argument(
        "--osca",
        default="./osca",
        help="Path to the osca executable. Run inside your normal setup environment.",
    )
    parser.add_argument(
        "--parallel",
        type=int,
        default=None,
        help="Forwarded to ./osca postprocess --parallel.",
    )
    parser.add_argument(
        "--run",
        action="store_true",
        help="Actually run ./osca postprocess for each generated YAML.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Regenerate datacards even if the target datacard already exists.",
    )
    parser.add_argument(
        "--allow-missing-backgrounds",
        action="store_true",
        help="Continue even if one or more requested background datasets are absent.",
    )
    parser.add_argument(
        "--luminosity-override",
        type=float,
        help=(
            "Normalize all MC inputs to this luminosity in fb^-1 instead of "
            "the luminosity stored in their era metadata. Use this only for "
            "an explicitly labelled luminosity projection."
        ),
    )
    return parser.parse_args()


def shell_quote(value: str) -> str:
    return "'" + value.replace("'", "'\"'\"'") + "'"


def discover_signals(results_dir: Path, signal_glob: str) -> list[dict[str, object]]:
    signals = []
    for path in sorted(results_dir.glob(signal_glob)):
        match = SIGNAL_RE.match(path.name)
        if not match:
            continue
        info = {
            "path": path,
            "mstop": int(match.group("mstop")),
            "mlsp": int(match.group("mlsp")),
            "sample_name": f"signal_T2tt_{match.group('mstop')}_{match.group('mlsp')}",
        }
        signals.append(info)
    return sorted(signals, key=lambda item: (item["mstop"], item["mlsp"]))


def yaml_for_signal(
    *,
    signal: dict[str, object],
    output_dir: Path,
    channel: str,
    backgrounds: list[str],
    luminosity_override: float | None = None,
) -> str:
    mstop = int(signal["mstop"])
    mlsp = int(signal["mlsp"])
    point_dir = output_dir / f"T2tt_mStop-{mstop}_mLSP-{mlsp}"

    background_lines = "\n".join(
        f"                - {{dataset_name: \"{name}\"}}" for name in backgrounds
    )

    luminosity_line = ""
    if luminosity_override is not None:
        if luminosity_override <= 0:
            raise ValueError("luminosity_override must be positive")
        luminosity_line = f"  luminosity_override: {luminosity_override:g}\n"

    return f"""Postprocessing:
{luminosity_line}  # This override changes MC yields and plot metadata; it is not just a label.
  processors:
    - name: CombineDatacard
      channel: {channel}
      inputs:
        - "*/*/*/SearchBinYield"
      structure:
        select: {{dataset_name: "*", name: "SearchBinYield"}}
        group: {{pipeline: "*", era.name: "*", name: "*"}}
        subgroups:
          signal:
            select: {{dataset_name: "2stop_{mstop}_2024_splitLSP"}}
          observation:
            select: {{dataset_name: "re:^$"}}
          background:
            select:
              or_exprs:
{background_lines}
      output_name: "{point_dir}"
"""


def run_postprocess(
    *,
    osca: str,
    config_path: Path,
    input_paths: list[Path],
    parallel: int | None,
) -> int:
    cmd = [osca, "postprocess", str(config_path), *(str(path) for path in input_paths)]
    if parallel is not None:
        cmd.extend(["--parallel", str(parallel)])
    print("Running:", " ".join(shell_quote(part) for part in cmd), flush=True)
    return subprocess.run(cmd).returncode


def input_paths_for_signal(
    *,
    signal: dict[str, object],
    all_inputs: list[Path],
    backgrounds: list[str],
) -> list[Path]:
    signal_path = Path(signal["path"])
    background_prefixes = tuple(f"{name}__" for name in backgrounds)
    background_paths = [
        path for path in all_inputs if path.name.startswith(background_prefixes)
    ]
    return [signal_path, *background_paths]


def main() -> int:
    args = parse_args()
    results_dir = args.results_dir.resolve()
    output_dir = args.output_dir.resolve()
    config_dir = args.config_dir.resolve()
    backgrounds = args.backgrounds or list(DEFAULT_BACKGROUNDS)

    if not results_dir.is_dir():
        print(f"ERROR: results_dir does not exist: {results_dir}", file=sys.stderr)
        return 2

    input_paths = sorted(results_dir.glob(args.inputs_glob))
    if not input_paths:
        print(
            f"ERROR: no input files matched {args.inputs_glob!r} in {results_dir}",
            file=sys.stderr,
        )
        return 2

    signals = discover_signals(results_dir, args.signal_glob)
    if not signals:
        print(
            f"ERROR: no split signal result files matched {args.signal_glob!r} in {results_dir}",
            file=sys.stderr,
        )
        return 2

    output_dir.mkdir(parents=True, exist_ok=True)
    config_dir.mkdir(parents=True, exist_ok=True)

    print(f"Found {len(signals)} split signal result files.")
    print(f"Found {len(input_paths)} total result files in {results_dir}.")
    print(f"Datacard output root: {output_dir}")
    print(f"Generated config root: {config_dir}")

    failures = 0
    for signal in signals:
        mstop = int(signal["mstop"])
        mlsp = int(signal["mlsp"])
        config_path = config_dir / f"combine_T2tt_mStop-{mstop}_mLSP-{mlsp}.yaml"
        datacard_path = (
            output_dir
            / f"T2tt_mStop-{mstop}_mLSP-{mlsp}"
            / f"datacard_{args.channel}.txt"
        )

        config_text = yaml_for_signal(
            signal=signal,
            output_dir=output_dir,
            channel=args.channel,
            backgrounds=backgrounds,
            luminosity_override=args.luminosity_override,
        )
        config_path.write_text(config_text)
        current_inputs = input_paths_for_signal(
            signal=signal, all_inputs=input_paths, backgrounds=backgrounds
        )

        found_backgrounds = {
            path.name.split("__", 1)[0]
            for path in current_inputs
            if path != signal["path"] and "__" in path.name
        }
        missing_backgrounds = sorted(set(backgrounds) - found_backgrounds)
        if missing_backgrounds:
            msg = (
                "No .result files found for backgrounds: "
                + ", ".join(missing_backgrounds)
            )
            if args.allow_missing_backgrounds:
                print(f"WARNING: {msg}", file=sys.stderr)
            else:
                print(f"ERROR: {msg}", file=sys.stderr)
                print(
                    "Put the background .result files in the same results_dir, "
                    "or pass --allow-missing-backgrounds if this is intentional.",
                    file=sys.stderr,
                )
                return 2

        if datacard_path.exists() and not args.overwrite:
            print(f"Skipping existing datacard: {datacard_path}")
            continue

        if args.run:
            rc = run_postprocess(
                osca=args.osca,
                config_path=config_path,
                input_paths=current_inputs,
                parallel=args.parallel,
            )
            if rc != 0:
                failures += 1
                print(
                    f"ERROR: postprocess failed for mStop={mstop}, mLSP={mlsp} "
                    f"with exit code {rc}",
                    file=sys.stderr,
                )
        else:
            cmd = [
                args.osca,
                "postprocess",
                str(config_path),
                *(str(p) for p in current_inputs),
            ]
            if args.parallel is not None:
                cmd.extend(["--parallel", str(args.parallel)])
            print("Prepared:", " ".join(shell_quote(part) for part in cmd))

    if failures:
        print(f"Completed with {failures} failed postprocess jobs.", file=sys.stderr)
        return 1

    if not args.run:
        print("Dry run only. Add --run to execute the generated postprocessing jobs.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
