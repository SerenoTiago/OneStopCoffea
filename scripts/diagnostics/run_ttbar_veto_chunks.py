#!/usr/bin/env python3
"""Run the isolated veto diagnostic over distinct, fixed-size input chunks."""

from __future__ import annotations

import argparse
import copy
import time
from pathlib import Path

from analyzer.core.analysis import loadAnalysis
from analyzer.core.event_collection import FileChunk
from analyzer.core.executors.finalizers import basicFinalizer
from analyzer.core.running import getRepos, getTasks
from analyzer.utils.querying import Pattern, PatternMode


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--chunks", type=int, default=10)
    parser.add_argument("--events-per-chunk", type=int, default=100000)
    return parser.parse_args()


def main():
    args = parse_args()
    analysis = loadAnalysis(str(args.config))
    datasets, eras = getRepos(analysis.extra_dataset_paths, analysis.extra_era_paths)
    tasks = getTasks(
        datasets,
        eras,
        analysis.event_collections,
        location_priorities=analysis.location_priorities,
        filter_dataset=Pattern("tt_semileptonic_2024", PatternMode.LITERAL),
        filter_sample=Pattern("TTtoLNu2Q", PatternMode.LITERAL),
    )
    if len(tasks) != 1:
        raise RuntimeError(f"Expected one TTtoLNu2Q task, found {len(tasks)}")
    task = tasks[0]
    task.file_set.updateFromCache()
    candidates = [
        info for info in task.file_set.files.values()
        if info.nevents is not None and info.nevents >= args.events_per_chunk
    ]
    if len(candidates) < args.chunks:
        raise RuntimeError(
            f"Only {len(candidates)} cached distinct files can supply "
            f"{args.events_per_chunk} events; need {args.chunks}"
        )
    candidates.sort(key=lambda item: item.file_path)
    # Evenly span the deterministic file list instead of taking its prefix.
    indices = [round(i * (len(candidates) - 1) / (args.chunks - 1)) for i in range(args.chunks)]
    chosen = [candidates[index] for index in indices]
    if len({item.file_path for item in chosen}) != args.chunks:
        raise RuntimeError("File selection produced duplicate input files")

    analysis.analyzer.initModules(task.metadata)
    merged = None
    started = time.monotonic()
    completed = []
    for index, info in enumerate(chosen):
        chunk = FileChunk(
            file_path=info.file_path,
            event_start=0,
            event_stop=args.events_per_chunk,
            tree_name=info.tree_name,
            schema_name=info.schema_name,
            file_nevents=info.nevents,
        )
        print(f"[{index + 1}/{args.chunks}] {info.file_path}", flush=True)
        result = analysis.analyzer.run(chunk, task.metadata, task.pipelines)
        result.finalize(basicFinalizer)
        merged = result if merged is None else merged + result
        completed.append(info.file_path)

    if merged["tt_semileptonic_2024"]["TTtoLNu2Q"]["_provenance"].chunked_events != args.chunks * args.events_per_chunk:
        raise RuntimeError("Merged provenance does not match requested event count")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    if args.output.exists():
        raise FileExistsError(args.output)
    args.output.write_bytes(merged.toBytes())
    elapsed = time.monotonic() - started
    manifest = args.output.with_suffix(".files.txt")
    manifest.write_text("\n".join(completed) + "\n")
    print(f"Wrote {args.output}")
    print(f"Processed {args.chunks * args.events_per_chunk} events in {elapsed:.1f} seconds")


if __name__ == "__main__":
    main()
