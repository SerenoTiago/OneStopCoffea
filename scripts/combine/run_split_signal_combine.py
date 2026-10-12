#!/usr/bin/env python3
"""Run combine for split T2tt datacards."""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import shutil
import subprocess
import sys
from pathlib import Path

COMBINE_IMAGE = (
    "/cvmfs/unpacked.cern.ch/gitlab-registry.cern.ch/cms-analysis/general/"
    "combine-container:latest"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run combine for each split-signal datacard directory."
    )
    parser.add_argument(
        "combine_dir",
        type=Path,
        help="Root directory produced by make_split_signal_datacards.py.",
    )
    parser.add_argument(
        "--merge",
        action="store_true",
        help="Pass --merge when using --use-generated-script.",
    )
    parser.add_argument(
        "--use-generated-script",
        action="store_true",
        help=(
            "Run each generated run_combine.sh instead of invoking the Combine "
            "container directly."
        ),
    )
    parser.add_argument(
        "--combine-image",
        default=COMBINE_IMAGE,
        help="CVMFS unpacked Combine container image to use.",
    )
    parser.add_argument(
        "--container-runtime",
        choices=("auto", "apptainer", "singularity"),
        default="auto",
        help="Container runtime used for the Combine container.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print commands without executing them.",
    )
    parser.add_argument(
        "--keep-going",
        action="store_true",
        help="Continue after a failed point instead of stopping immediately.",
    )
    parser.add_argument(
        "--significance-only",
        action="store_true",
        help="Run expected Asimov significance only; skip AsymptoticLimits.",
    )
    parser.add_argument(
        "--expected-only",
        action="store_true",
        help="Run only blinded expected AsymptoticLimits; skip significance and observed limits.",
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=1,
        help="Number of mass points to run concurrently (default: 1).",
    )
    return parser.parse_args()


def shell_quote(value: str) -> str:
    return "'" + value.replace("'", "'\"'\"'") + "'"


def combine_container_command(
    point_dir: Path,
    datacard: Path,
    method: str,
    image: str,
    runtime: str,
    expected_only: bool = False,
) -> list[str]:
    point_dir = point_dir.resolve()
    return [
        runtime,
        "exec",
        "-B",
        str(point_dir),
        "-B",
        "/cvmfs",
        "--pwd",
        str(point_dir),
        image,
        "/bin/bash",
        "-lc",
        (
            "source /cvmfs/cms.cern.ch/cmsset_default.sh && "
            'cd "$(ls -d /home/cmsusr/CMSSW_* | head -1)" && '
            "eval $(scramv1 runtime -sh) && "
            f"cd {shell_quote(str(point_dir))} && "
            f"combine -M {method} {shell_quote(datacard.name)}"
            + (" -t -1 --expectSignal=1" if method == "Significance" else "")
            + (" --run blind" if method == "AsymptoticLimits" and expected_only else "")
        ),
    ]


def run_and_log(cmd: list[str], log_path: Path, dry_run: bool) -> int:
    print("Running:", " ".join(shell_quote(part) for part in cmd))
    if dry_run:
        return 0

    with log_path.open("w") as handle:
        completed = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
        print(completed.stdout, end="")
        handle.write(completed.stdout)
    return completed.returncode


def main() -> int:
    args = parse_args()
    combine_dir = args.combine_dir.resolve()

    if not combine_dir.is_dir():
        print(f"ERROR: combine_dir does not exist: {combine_dir}", file=sys.stderr)
        return 2

    runtime = args.container_runtime
    if runtime == "auto":
        runtime = shutil.which("apptainer") or shutil.which("singularity") or ""
        if runtime:
            runtime = Path(runtime).name
    elif shutil.which(runtime) is None:
        runtime = ""

    if not args.use_generated_script and not runtime:
        print(
            "ERROR: neither 'apptainer' nor 'singularity' is available in this shell. "
            "Run this script from the host bash shell, not from inside the OSCA "
            "Singularity container, or pass --use-generated-script in an environment "
            "where combine is already installed.",
            file=sys.stderr,
        )
        return 2

    point_dirs = sorted(
        path for path in combine_dir.glob("mStop-*/T2tt_mStop-*_mLSP-*")
        if path.is_dir()
    )
    if not point_dirs:
        print(
            f"ERROR: no T2tt_mStop-*_mLSP-* directories found under {combine_dir}",
            file=sys.stderr,
        )
        return 2

    if args.jobs < 1:
        print("ERROR: --jobs must be at least 1", file=sys.stderr)
        return 2

    def run_point(point_dir: Path) -> int:
        print(f"Point: {point_dir.name}")

        if args.use_generated_script:
            script = point_dir / "run_combine.sh"
            if not script.exists():
                print(f"ERROR: missing {script}", file=sys.stderr)
                return 2

            cmd = ["bash", script.name]
            if args.merge:
                cmd.append("--merge")

            print("Running:", " ".join(shell_quote(part) for part in cmd))
            if args.dry_run:
                return 0
            completed = subprocess.run(cmd, cwd=point_dir)
            rc = completed.returncode
        else:
            datacards = sorted(point_dir.glob("datacard_*.txt"))
            if len(datacards) != 1:
                print(
                    f"ERROR: expected exactly one datacard_*.txt in {point_dir}, "
                    f"found {len(datacards)}",
                    file=sys.stderr,
                )
                return 2
            datacard = datacards[0]
            if args.expected_only:
                rc = run_and_log(
                    combine_container_command(
                        point_dir,
                        datacard,
                        "AsymptoticLimits",
                        args.combine_image,
                        runtime,
                        expected_only=True,
                    ),
                    point_dir / f"{datacard.stem}_limits.log",
                    args.dry_run,
                )
            else:
                rc = run_and_log(
                    combine_container_command(
                        point_dir, datacard, "Significance", args.combine_image, runtime
                    ),
                    point_dir / f"{datacard.stem}_significance.log",
                    args.dry_run,
                )
                if rc == 0 and not args.significance_only:
                    rc = run_and_log(
                        combine_container_command(
                            point_dir,
                            datacard,
                            "AsymptoticLimits",
                            args.combine_image,
                            runtime,
                        ),
                        point_dir / f"{datacard.stem}_limits.log",
                        args.dry_run,
                    )

        if rc != 0:
            print(
                f"ERROR: combine failed in {point_dir} with exit code {rc}",
                file=sys.stderr,
            )
        return rc

    failures = 0
    if args.jobs == 1:
        for point_dir in point_dirs:
            rc = run_point(point_dir)
            if rc != 0:
                failures += 1
                if not args.keep_going:
                    return rc
    else:
        with ThreadPoolExecutor(max_workers=args.jobs) as executor:
            futures = {executor.submit(run_point, point_dir): point_dir for point_dir in point_dirs}
            for future in as_completed(futures):
                rc = future.result()
                if rc != 0:
                    failures += 1

    if failures:
        print(f"Completed with {failures} failed combine points.", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
