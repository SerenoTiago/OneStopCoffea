"""Rerun configuration for single Z->nunu samples (2026-10-15).

Identical to backgrounds.yaml (loaded verbatim, so the selection, triggers,
location priorities and pipelines are unchanged).  It only adds an executor
equal to the premade lpc-dask-condor-4G-100000 except that min_workers equals
max_workers, which pins the Dask cluster size.

Reason: the 2026-10-14 zjets_lo_2024 production lost ZNuNu4Jets_HT100to200
to an "Overlapping Provenance" merge failure during heavy adaptive-scaling
churn (679 workers removed by dask-jobqueue).  See
docs/investigations/dask_exec_silent_data_loss_2026-10-15.md.

Usage:
  ./osca run configurations/production/hadronic_susy/backgrounds_zjets_rerun.py \
      <output> -e lpc-dask-condor-4G-100000-fixed \
      --filter-dataset 'zjets*2024*' --filter-sample ZNuNu4Jets_HT100to200
"""

from pathlib import Path

from analyzer.utils.yamlload import loadTemplateYaml

ANALYSIS = loadTemplateYaml(Path(__file__).with_name("backgrounds.yaml"))
ANALYSIS["extra_executors"] = {
    "lpc-dask-condor-4G-100000-fixed": {
        "executor_name": "LPCCondorDask",
        "container": "/cvmfs/unpacked.cern.ch/registry.hub.docker.com/coffeateam/coffea-dask-almalinux9:2025.10.2-py3.12",
        "chunk_size": 100000,
        "min_workers": 250,
        "max_workers": 250,
        "worker_memory": "4GB",
    },
}
