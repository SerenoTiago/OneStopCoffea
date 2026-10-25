# Silent data loss in the Dask executor (open framework issue)

Date: 2026-10-15.  Status: **documented only; no code changed.**
Found while running the LO `zjets_2024` production
(`production_runs/2026-10-07_012944_hadronic_susy/logs/zjets_lo_2024.log`).

The Dask executor (`analyzer/core/executors/dask_exec.py`) can drop input
events without failing the run.  The process still exits 0, so the exit code
alone cannot be trusted; every output must be checked against the expected
event count (provenance `chunked_events` vs DAS `n_events`).

## 1. A failed file-info (`prep`) task drops the whole file

- `run()` catches the exception from a prep future and only logs it
  (`dask_exec.py:360-364`, `logger.warning(e)`); the warning does not name the
  file.
- With `result is None` the file's `nevents` is never filled in
  (`dask_exec.py:369-370`).
- `FileSet.toChunked` keeps only files whose `nevents` is known
  (`event_collection.py:433-439`), so the file is never chunked or processed.
- Observed: one XRootD `Operation expired` at 00:06:47 for an HT100-200 file.

## 2. A failed chunk is dropped from its sample

- `getAnalyzerRunFunc.inner` turns any exception in a chunk into
  `DaskRunResult(None, [DaskRunException(chunk, e)])`
  (`dask_exec.py:201-202`).
- When the merged sample result is non-empty it is yielded and saved
  (`dask_exec.py:396-397`); the accumulated `maybe_exceptions` are discarded
  without a log message.  The saved provenance correctly omits the failed
  chunks, so this is detectable afterwards, but nothing reports it.

## 3. A failed merge or completion discards the entire sample

- If a merge raises (here `ResultIntegrityError: Overlapping Provenance` from
  `ResultGroup.checkOk`, `results.py:191-199`), the exception is logged as a
  plain warning (`dask_exec.py:363`) and the completion future's `result` is
  `None`.
- `dask_exec.py:391-393` then `continue`s: no output file and no further
  message.  The `Result was None` warning (`dask_exec.py:399-401`) only covers
  the different case where the completion succeeded but held no result.
- Observed: the complete `ZNuNu4Jets_HT100to200` sample (1787 files,
  ~332 M events) was lost at 00:30:16 after ~57 min.  The overlap involved
  chunk (188632, 282948) of `.../90000/9a138ff3-e550-4f76-a45c-d8d13b314cb0.root`.
  It coincided with heavy adaptive-scaling churn (679 Condor workers removed by
  dask-jobqueue, many "lost already computed task(s)" reschedules), but the
  exact mechanism producing a duplicated chunk was not established.

## Recovery of HT100-200 (same day)

- Rerun with `configurations/production/hadronic_susy/backgrounds_zjets_rerun.py`
  (cluster pinned at 250 workers, `min_workers == max_workers`): 0 workers
  removed, no overlap failure, but 6 files were dropped by failed prep tasks
  (path 1) and one chunk by a failed chunk (path 2): 390,462 events short.
- Two `osca patch` rounds (`buildMissingFileset`) recovered exactly the
  missing files/chunks; the first round itself dropped one file to another
  XRootD timeout.  `osca merge` of the three results has 332,233,649 processed
  events = DAS, 1787/1787 files, no overlap.
- Inference (not proven): `iaddMany` mutates its first input in place, so a
  merge recomputed after a worker loss can re-add already-merged inputs; this
  would explain the overlap failure under adaptive churn.

## Suggested direction (not implemented)

- Fail or at least summarise loudly: collect prep/chunk/merge failures and
  report them per sample at the end, with a nonzero exit code if any sample is
  incomplete.
- Name the file in prep-failure messages.
- Investigate how a chunk can appear twice in one merge under worker churn
  before relaxing the overlap check; the check itself correctly prevents
  double counting.
