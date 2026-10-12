# Run-3 2024 T2tt grid-extension audit

Date: 2026-10-07.  No physics selection, datacard, or Combine job was run.

## Source of truth and the 1200/1250 mismatch

The strongest mass identifier available in these NanoAOD files is the event-level
`Events/GenModel_T2tt_<mStop>_<mLSP>` flag.  The same labels are independently
repeated in `Runs/genEventCount_T2tt_<mStop>_<mLSP>`.  A 25-event check per mass
point also found last-copy `GenPart` stops (`abs(pdgId)==1000006`) and neutralinos
(`abs(pdgId)==1000022`) at exactly the masses encoded by each flag.

The DAS dataset whose primary name says `Par-mStop-1200` actually contains only
`GenModel_T2tt_1250_*` points.  Thus the existing split files are correctly named
1250 at generator level; the DAS primary-dataset name is misleading.  The same
offset occurs systematically in the higher-mass samples listed below, so primary
dataset names must not be used as mass truth.

The 1250 DAS dataset has two files and 488,217 total events.  The earlier local
split used only the 413,008-event file.  Its nonempty split shards have the exact
selected counts and correct generator masses, but the source file has unreadable
baskets near entry 170,669 and its 1250/1000 output has no `Events` tree.  The
second 75,209-event file was downloaded, split, and validated for all twelve LSP
flags.  The configured eleven usable points now include both source-file shards;
1250/1000 remains excluded because the first-file shard cannot be recovered.

## Run-3 2024 inventory verified from ROOT content

`Files/events` are DAS totals.  LSP lists are unions of the `GenModel_T2tt_*`
branches across every file in each dataset.  Representative files for every row,
and every file where branch coverage differed, were opened with ROOT.
The per-file LFNs, event counts, sizes, and DBS checksums are recorded in
[`run3_2024_files.tsv`](run3_2024_files.tsv).

| DAS name says mStop | ROOT mStop | LSP masses in ROOT (GeV) | Files / events |
|---:|---:|---|---:|
| 1200 | 1250 | 1, 100-1100 by 100 | 2 / 488,217 |
| 1250 | 1300 | 1, 100-1100 by 100 | 3 / 512,229 |
| 1300 | 1350 | 1, 100-1200 by 100 | 2 / 491,181 |
| 1350 | 1400 | 1, 100-1200 by 100 | 2 / 401,505 |
| 1400 | 1450 | 1, 100-1300 by 100 | 2 / 478,264 |
| 1450 | 1500 | 1, 100-1300 by 100 | 3 / 478,261 |
| 1500 | 1600 | 1, 100-1400 by 100 | 2 / 446,441 |
| 1600 | 1700 | 1, 100-1500 by 100 | 3 / 539,465 |
| 1700 | 1800 | 1, 100-1600 by 100 | 2 / 506,695 |
| 1800 | 1900 | 1, 100-1700 by 100 | 2 / 528,098 |
| 1900 | 2000 | 1, 100-1800 by 100 | 3 / 565,753 |

There is no generic RunIII2024Summer24NanoAODv15 dataset yielding mStop=1550,
1650, 1750, 1850, or 1950 in the generator flags.  For mStop=1600 the canonical
choice is the 2-file `madgraphMLM-pythia8` v1 dataset (446,441 events, prep ID
`NPS-RunIII2024Summer24NanoAODv15-00725`).  It is part of the contiguous newer
NPS-00719--00728 campaign used for every other mass in this grid.  The overlapping
5-file `madgraph-pythia8` v2 dataset (565,793 events, prep ID
`SUS-RunIII2024Summer24NanoAODv15-00058`) uses an older, different production
configuration.  It was neither downloaded nor merged because there is no evidence
that the two productions are statistically disjoint or intended to be combined.

## Run-3 coordinates on the Run-2 published surface

The official Figure 14a file stores the expected surface as the 128 x 160 TH2D
`MassScan2DExp`, with x edges [150, 2000) GeV and y edges [1, 1150) GeV.  The
published HEPData coordinates are the raster's bin centers, not a list of nominal
simulated mass points.  Nominal masses must therefore be looked up in their
containing bin (ROOT `FindBin` semantics), without requiring equality to a center
and without interpolation.

All 120 coordinates formed from the ten requested stop masses and the twelve LSP
masses fall in finite, positive expected-limit bins.  The usable local grid has
119 points after excluding the known-bad (1250,1000) signal shard.  The previous
zero-match conclusion was solely an artifact of comparing nominal masses to bin
centers.

## Local files and validation

Unsplit files are under `analysis_inputs/signals/raw/2024/mStop-XXXX/`.  Split
shards are under
`analysis_inputs/signals/split/2024/mStop-XXXX/<source UUID>/`.  The source UUID
layer is retained deliberately: shards from separate input files are not merged,
and provenance remains unambiguous.  Only mLSP=1 and 100--1100 GeV by 100 were
split for this grid.

For the newly split 75,209-event file, every output has one readable `Events`
tree; each output count exactly equals the sum of its input `GenModel` flag; and
25 checked events per output have exact last-copy stop/LSP generator masses.
Legacy and failed artifacts are retained separately under `legacy-preaudit/`,
`failed-regeneration/`, and `regeneration-audit/` for provenance.

For mStop=1300--1900, 21 source files produced 237 source-preserving shards.
All 237 have readable, nonempty `Events` trees.  Per-branch input counts total
3,664,601 selected events and exactly equal the summed output counts.  The first
event in every shard has the requested last-copy stop and neutralino masses in
`GenPart`.  The actual file set equals the set predicted by the input branches,
with no duplicate or extra shards.  Including the 22 usable canonical 1250
shards, 259 shards were reopened successfully in the final validation pass.

## Later Combine implications

The five Run-3 background result sets are signal independent in the current
workflow and can be reused for every signal point, provided the selection,
luminosity override, and background systematic model remain identical.

To group future outputs as
`mStop-XXXX/T2tt_mStop-XXXX_mLSP-YYY/`, change the point path in
`make_split_signal_datacards.py`, then make `run_split_signal_combine.py`, the
per-point loop in `run_full_pipeline.sh`, and
`plot_split_signal_significances.py` discover nested directories.  The limit
collector already uses a recursive log glob and extracts the point from any path
component, so it is already compatible.

The expected eventual Combine workload is 119 signal points.  No datacard or
Combine job was run while constructing or validating this file-level grid.
