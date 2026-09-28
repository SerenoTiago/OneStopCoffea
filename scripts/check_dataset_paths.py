#!/usr/bin/env python3
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path
from typing import Any
from urllib.parse import urlparse

import yaml


def load_yaml(path: Path) -> list[dict[str, Any]]:
    text = path.read_text()
    if "\t" in text:
        print(f"[WARN] {path}: replacing tab characters with spaces before YAML parsing")
        text = text.replace("\t", " ")
    data = yaml.safe_load(text)
    if data is None:
        return []
    if not isinstance(data, list):
        raise TypeError(f"{path}: expected top-level YAML list, got {type(data).__name__}")
    return data


def is_remote(path: str) -> bool:
    parsed = urlparse(path)
    return bool(parsed.scheme and parsed.scheme != "file")


def check_local_path(path: str, yaml_path: Path) -> tuple[bool, str]:
    if is_remote(path):
        return True, "remote"

    p = Path(path)
    if not p.is_absolute():
        p = yaml_path.parent / p
    if p.exists():
        return True, "exists"
    return False, "missing"


def check_open_root(path: str, tree_name: str) -> tuple[bool, str]:
    try:
        import uproot
    except Exception as exc:
        return False, f"could not import uproot: {exc}"

    try:
        with uproot.open(path) as root_file:
            if tree_name not in root_file:
                return False, f"opened, but missing tree '{tree_name}'"
        return True, f"opened tree '{tree_name}'"
    except Exception as exc:
        return False, str(exc)


def query_das_files(
    das_path: str, limit: int | None, instance: str
) -> tuple[bool, str, list[str]]:
    query = f"file dataset={das_path}"
    if instance:
        query = f"{query} instance={instance}"
    try:
        proc = subprocess.run(
            ["dasgoclient", f"--query={query}"],
            check=False,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
    except FileNotFoundError:
        return False, "dasgoclient not found in PATH", []

    if proc.returncode != 0:
        msg = proc.stderr.strip() or proc.stdout.strip() or f"exit code {proc.returncode}"
        return False, f"{msg}; attempted query: {query}", []

    files = [line.strip() for line in proc.stdout.splitlines() if line.strip()]
    if limit is not None:
        files = files[:limit]
    if not files:
        return False, "DAS query returned no files", []
    return True, f"DAS returned {len(files)} file(s)", files


def iter_samples(dataset: dict[str, Any]) -> list[dict[str, Any]]:
    samples = dataset.get("samples")
    if samples is None:
        return [dataset]
    if not isinstance(samples, list):
        raise TypeError(f"dataset {dataset.get('name', '<unknown>')}: samples is not a list")
    return samples


def check_dataset_yaml(path: Path, args: argparse.Namespace) -> tuple[int, int]:
    passed = 0
    failed = 0
    print(f"\n== {path} ==")

    for dataset in load_yaml(path):
        dataset_name = dataset.get("name", "<unknown dataset>")
        for sample in iter_samples(dataset):
            sample_name = sample.get("name", dataset_name)
            tree_name = sample.get("tree_name", dataset.get("tree_name", "Events"))
            label = f"{dataset_name} / {sample_name}"

            files = sample.get("files")
            das_path = sample.get("das_path")

            if files is not None:
                if not isinstance(files, list):
                    print(f"[FAIL] {label}: files is not a list")
                    failed += 1
                    continue

                for file_path in files:
                    ok, msg = check_local_path(str(file_path), path)
                    status = "OK" if ok else "FAIL"
                    print(f"[{status}] {label}: {file_path} ({msg})")
                    passed += int(ok)
                    failed += int(not ok)

                    if ok and args.open_root:
                        ok_open, msg_open = check_open_root(str(file_path), tree_name)
                        status_open = "OK" if ok_open else "FAIL"
                        print(f"  [{status_open}] open: {msg_open}")
                        passed += int(ok_open)
                        failed += int(not ok_open)

            elif das_path is not None:
                if not args.check_das:
                    print(f"[SKIP] {label}: DAS path present, not queried: {das_path}")
                    continue

                ok, msg, das_files = query_das_files(
                    str(das_path), args.das_limit, args.das_instance
                )
                status = "OK" if ok else "FAIL"
                print(f"[{status}] {label}: {das_path} ({msg})")
                passed += int(ok)
                failed += int(not ok)

                if ok and args.open_root and das_files:
                    redirector = args.redirector.rstrip("/")
                    for das_file in das_files:
                        root_url = f"{redirector}/{das_file}"
                        ok_open, msg_open = check_open_root(root_url, tree_name)
                        status_open = "OK" if ok_open else "FAIL"
                        print(f"  [{status_open}] open: {root_url} ({msg_open})")
                        passed += int(ok_open)
                        failed += int(not ok_open)
            else:
                print(f"[FAIL] {label}: no 'files' or 'das_path' entry")
                failed += 1

    return passed, failed


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Check file path entries in OneStopCoffea dataset YAML files."
    )
    parser.add_argument("yamls", nargs="+", type=Path, help="Dataset YAML file(s).")
    parser.add_argument(
        "--check-das",
        action="store_true",
        help="Query das_path entries with dasgoclient.",
    )
    parser.add_argument(
        "--das-limit",
        type=int,
        default=1,
        help="Maximum DAS files to query/open per sample when --check-das is set.",
    )
    parser.add_argument(
        "--das-instance",
        default="",
        help="Optional DAS instance to include in the query, e.g. prod/phys03.",
    )
    parser.add_argument(
        "--open-root",
        action="store_true",
        help="Try opening ROOT files with uproot and check the tree exists.",
    )
    parser.add_argument(
        "--redirector",
        default="root://cms-xrd-global.cern.ch/",
        help="Redirector prepended to DAS /store paths when --open-root is used.",
    )
    args = parser.parse_args()

    total_passed = 0
    total_failed = 0
    for yaml_path in args.yamls:
        passed, failed = check_dataset_yaml(yaml_path, args)
        total_passed += passed
        total_failed += failed

    print(f"\nSummary: {total_passed} passed, {total_failed} failed")
    return 1 if total_failed else 0


if __name__ == "__main__":
    sys.exit(main())
