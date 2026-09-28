# Audit of the T2tt equal-luminosity comparison

Date: 2026-09-10. Scope: investigate comparison fairness, trace the 1250/100 GeV signal point, and compare local and published search-bin yields. No analysis code, datacards, or plots were changed; no plot-producing code was written. The CSV files here are diagnostic tables, not corrected physics results.

## Findings that change the interpretation

The existing plot is not a like-for-like reproduction of CMS-SUS-19-006. Its local input is a background-only artificial dataset with no nuisance parameters; its published input is an observed limit. The local samples are configured and saved at 13.6 TeV, whereas CMS used 13 TeV. Both use 137 fb^-1, but equal luminosity does not resolve these differences.

More importantly, a confirmed jet-count selection bug excludes events containing exactly 3, 5, 7, or 9 selected jets. The saved histograms demonstrate that this affected the inputs behind the plot. The signal efficiency and background distribution must be revisited before interpreting a limit ratio as a physics-performance difference.

The background comparison also revises the working hypothesis: local total background is *higher* than the paper's pre-fit central value in 167 of 174 nominally corresponding bins, and lower in only seven. QCD is a major contributor. This comparison is diagnostic: the local jet selections are currently wrong, and the collision energies and reconstruction differ.

## 1. Is the limit comparison fair?

### Expected versus observed

All 11 local cards use `kmax 0`, contain no nuisance-parameter rows, and have no `autoMCStats` directive. Their `data_obs` histogram equals the sum of the five nominal background templates bin by bin. Thus, their reported “observed” limits describe artificial background-only data, not detector observations. Counting fluctuations still enter the likelihood; template statistical and systematic uncertainties are not modeled.

The existing HEPData file already contains **both observed and expected cross-section limits**:

`data/SUS-19-006/HEPData-ins1749379-v1-T2tt_cross_section_upper_limits.yaml`.

I extracted the expected values and independently checked them against the official Figure 14a ROOT histogram `MassScan2DExp`. At all 11 chosen positions they agree to within 0.0026%, consistent with rounding in the YAML. The local CSV's `exp` field is the median expected Combine result.

| Neutralino mass [GeV] | Local expected limit [pb] | CMS expected map value [pb] | Local/CMS |
|---:|---:|---:|---:|
| 1 | 0.00255177 | 0.0016372 | 1.559 |
| 100 | 0.00259880 | 0.0016856 | 1.542 |
| 200 | 0.00270469 | 0.0017587 | 1.538 |
| 300 | 0.00290458 | 0.0018774 | 1.547 |
| 400 | 0.00339844 | 0.0021262 | 1.598 |
| 500 | 0.00419808 | 0.0025792 | 1.628 |
| 600 | 0.00543287 | 0.0034070 | 1.595 |
| 700 | 0.00717318 | 0.0047780 | 1.501 |
| 800 | 0.0106305 | 0.0076859 | 1.383 |
| 900 | 0.0165102 | 0.013498 | 1.223 |
| 1100 | 0.0536228 | 0.053465 | 1.003 |

This removes much of the apparent high-neutralino-mass disagreement in the original observed comparison. It does not fix the selections, energy difference, or missing uncertainty model. These are comparisons of the available expected outputs, not claims of equivalent likelihoods.

### Mass coordinates

The HEPData values represent the published map's grid. There are no exact local mass-coordinate matches in the original comparison. At local (1250,100) GeV, the selected published limit-map coordinate is (1255.7,97.947) GeV. The official ROOT cell containing (1250,100) has edges [1248.4375,1262.890625] and [94.35625,101.5375] GeV, respectively, and reproduces the selected HEPData value.

This verifies that the unusual coordinates are consistent with the published ROOT map, rather than evidence of an arbitrary downloaded table. It does not establish the original simulated mass-grid point or remove interpolation already present in the published map. A future plot should say “published map sampled at local masses,” not “exact matched mass points.”

### Luminosity, collision energy, and reference cross section

Saved metadata for all 28 audited input samples (one signal plus 27 background samples) contain luminosity 137.0 fb^-1 and energy 13.6 TeV. This is not inferred just from today's YAML or the plot title. Background dataset names/DAS paths identify Run 3 production. The era named `run2` retains 2024 reconstruction/correction metadata and a UParT b-tag working point. The paper uses 2016–2018, 13 TeV collisions and a DeepCSV medium b-tag working point.

The local signal normalization and limit conversion both use 0.7526 fb. The conversion to pb divides by 1000 correctly. The physical provenance/perturbative order of that reference cross section was not established. However, simply changing a consistently used reference cross section rescales the signal template and its fitted multiplier inversely; it is not by itself a remedy for a cross-section-limit discrepancy. Signal acceptance and background predictions must also be valid for the intended energy and model.

The published analysis and selections are documented in the repository's `data/references/SUS-19-006/1908.04722v2.pdf` (Sections 4–6 and Appendix A) and on the [CMS publication page](https://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-19-006/index.html).

## 2. Trace of the 1250/100 point

### Signal source and normalization

Source: `analysis_inputs/signals/split/2024/624134b4-1901-4471-b93d-76e7d308946e/signal_T2tt_1250_100.root`.

The Events tree contains 36,246 entries, matching both configured and processed event counts. All `genWeight` values are +1. A spot check of the first 100 events finds two last-copy stops with mass 1250 GeV and two last-copy neutralinos with mass 100 GeV in every event. This is a model/mass sanity check, not a complete generator validation. The split file has no Runs tree; full original production provenance was not reconstructed.

The saved result is now archived at `analysis_products/old_results/26-08-30_hadsusy_1250_run2/2stop_1250_2024_splitLSP__signal_T2tt_1250_100.result`.

| Quantity | Value |
|---|---:|
| Processed generated events | 36,246 |
| Selected events in SearchBinYield | 9,541 |
| Raw selection fraction | 26.3229% |
| Reference produced events, 137 × 0.7526 | 103.1062 |
| Normalization factor per simulated event | 0.002844622855 |
| Selected expected signal yield | 27.140546659 |
| Median expected Combine multiplier | 3.4531 |
| Expected cross-section limit | 0.00259880306 pb |
| Artificial-data “observed” multiplier | 3.4626 |
| Cross-section value in the original plot | 0.00260595276 pb |

The scaling code is `analyzer/core/results.py::mergeAndScale`: luminosity × cross section / processed events. I independently applied that formula to each raw sample histogram, then summed samples by process. This reproduced **every bin exactly** in all six corresponding ROOT signal/background templates. Template integrals agree with datacard rates to floating-point precision. No SearchBinYield underflow or overflow is present in the 28 audited raw inputs.

Thus, for this point, there is no evidence of an additional luminosity scaling, an fb/pb conversion error, or dropped template bins between saved results and the card. This does not validate the upstream selections, physical background cross sections, or the representativeness of partially processed samples.

As a separate statistical sanity check, evaluating the fixed-background Poisson Asimov expression `q(r) = 2 sum[r*s - b*log(1+r*s/b)]` and solving `q = 1.959964^2` gives r = 3.461574. This is within 0.25% of the saved Combine median expected value 3.4531. The check uses the same asymptotic approximation and is not an exact low-count coverage test or a replacement Combine run, but it provides no indication of an order-one error in the final limit calculation.

### Confirmed jet-count bug

In `analyzer/modules/common/hadronic_susy.py`, `NJET_BINS` stores `(2,3), (4,5), (6,7), (8,9), (10,None)`. But `njet_bin_index` passes them to `_in_range`, which uses `value < high`. Consequently:

| Intended category | Current accepted jet counts |
|---|---|
| 2–3 | 2 only |
| 4–5 | 4 only |
| 6–7 | 6 only |
| 8–9 | 8 only |
| >=10 | >=10 |

Direct evaluation confirms that 3, 5, 7, and 9 map to invalid category -1. All saved selected JetN histograms in the audited signal/background inputs have zero entries at those four counts. The selected signal JetN counts are 831 at 2 jets, 4,110 at 4, 3,380 at 6, and 1,025 at 8, with the remaining events at >=10 jets.

The signal cutflow has 18,664 events after the isolated-track veto but only 9,541 after `search_bin_valid`: a loss of 9,123 (48.88%). **Not all of this loss can be assigned to the bug from the saved one-dimensional histograms**, because the final validity cut also enforces legitimate kinematic exclusions, including the dropped intervals at >=8 jets. Exact recovery requires event-level reevaluation.

The complete saved signal cutflow is in `normalization_audit.json`. In order: initial 36,246; electron veto 31,615; muon veto 26,988; photon veto 26,793; jet multiplicity 25,972; HT 25,845; MHT 22,003; MHT<HT 21,745; angular cuts 18,908; isolated-track veto 18,664; valid search bin 9,541.

The downloaded CMS T2tt efficiency map gives 0.3927 at its nearest coordinate (1251.5,97.947) GeV, compared with local 0.263229. The local acceptance is about 0.670 of that published value. Different map coordinates, collision energy, reconstruction, and corrections prevent treating this as an exact isolated bug measurement. The published light-LSP T2tt cutflow is at (950,100), not (1250,100), so it is unsuitable for a direct numerical cut-by-cut comparison at our chosen point.

### Other upstream differences requiring follow-up

The local isolated-track transverse mass uses jet-based HTMiss and its azimuth; the paper defines it with PF missing transverse momentum. The active configuration also uses Run 3 b tagging and first-pass lepton definitions. It does not explicitly schedule trigger, event-noise-filter, or event-weight correction modules in `susy_full_validation`; the presence of correction filenames in era metadata does not prove those corrections were applied.

All 28 raw selected histograms have sum of weights equal to sum of squared weights and integer bin contents, consistent with unit-weight filling. Histogram construction only applies weights if a `Weights` collection exists. For the signal, constant generator weights make count normalization appropriate; generator-weight handling and sample cross sections for the backgrounds remain to be audited. The current five background processes also omit explicit dileptonic ttbar, single-top, and other smaller samples; the paper's lost-lepton estimate is not identical to a sum of the local process labels.

## 3. Background comparison in all 174 bins

I downloaded all five HEPData tables of **pre-fit** backgrounds and observations. The first table's 30 total-background and observed values agree exactly with the CSV already in the repository. All 174 published Njet, Nb, HT, and MHT labels agree with the *intended* local layout, including omitted kinematic intervals. The actual Njet acceptance differs because of the bug above.

The complete comparison, with local signal, artificial observation, process yields, weighted-MC standard deviation, effective MC count, and separate published asymmetric statistical/systematic uncertainties, is in `background_comparison_174_bins.csv`.

| Nominal Njet range | Local total background | Published pre-fit total |
|---|---:|---:|
| 2–3 | 5,876,295.98 | 309,777.40 |
| 4–5 | 1,446,637.94 | 91,707.30 |
| 6–7 | 66,014.83 | 16,810.40 |
| 8–9 | 5,007.16 | 2,213.30 |
| >=10 | 416.60 | 249.13 |
| All | 7,394,372.50 | 420,757.53 |

Totals are sums of tabulated central values, with no assertion about aggregate uncertainty. Correlations would be required to construct that uncertainty. QCD totals are 7,103,894 locally versus 8,330.81 in the paper. Large inclusive totals alone do not measure their influence on the limit.

Examples in bins carrying substantial local signal:

| Bin | Local signal | Local background | Published background | Published observed | Local/published background |
|---:|---:|---:|---:|---:|---:|
| 59 | 1.9002 | 87.8875 | 13.3 | 15 | 6.608 |
| 60 | 0.4551 | 9.8654 | 1.5 | 1 | 6.577 |
| 99 | 1.4166 | 31.0855 | 4.8 | 9 | 6.476 |
| 100 | 0.4182 | 5.9824 | 1.7 | 1 | 3.519 |

These local backgrounds include QCD contributions of 65.0806, 7.2350, 17.0952, and 3.0830 events, respectively. This makes inflated background a plausible contributor to weaker sensitivity in addition to lost signal. Its quantitative contribution has not been isolated with reruns.

Only bins 26, 134, 144, 145, 146, 151, and 158 have smaller local total backgrounds. Published central backgrounds are zero in bins 149, 157, 165, and 173 but have nonzero upper uncertainties; ratios are deliberately left undefined there. A zero published central value must not be treated as perfectly known zero background.

### Sparse QCD tails

The low-pT QCD samples produce enormous extrapolations from very few passing MC events:

| QCD sample | Processed events | Passing events | Expected selected yield |
|---|---:|---:|---:|
| pT 15–20 | 95,768,884 | 2 | 2,534,036.00 |
| pT 20–30 | 66,835,110 | 1 | 852,110.51 |
| pT 30–50 | 98,951,656 | 7 | 1,088,366.83 |

For example, bin 10 has 155,581.43 local background events versus 50.7 in the paper, and local effective MC count about 1.001. Effective count here means `(sum weights)^2 / sum weights squared`, which indicates the statistical precision of the weighted estimate, not the number of real collision events. The current Combine model treats the resulting large background expectation as fixed. These rare passing QCD events deserve event-level inspection for reconstruction, filtering, pileup, and normalization issues; their origin cannot be diagnosed from these histograms alone. They should not simply be deleted or capped.

## Proposed plots — specification only

### A. Search-bin background and signal audit

Use the 174 published bin numbers on the horizontal axis with boundaries at 30/31, 70/71, 110/111, and 142/143, and lighter subdivisions for b-jet categories. Label the comparison explicitly as local 13.6 TeV MC at 137 fb^-1 versus published 13 TeV pre-fit estimates at 137 fb^-1. Any display of the current templates must prominently note the missing odd jet multiplicities.

The first panel would show the local total background and local QCD separately, the published total-background central values with their per-bin asymmetric uncertainty, and the published real observations. A logarithmic vertical axis is needed. The local signal at its reference normalization could be a clearly labeled unstacked line. Do not plot local `data_obs` as real observations: it duplicates local background.

The second panel would show local total background / published total background and, separately, local QCD / published QCD. Zero denominators would be marked as undefined, with their original yields available in the table. These ratios compare estimates, not statistical significances. Published per-bin total uncertainty may be visualized by combining its tabulated statistical and systematic components in quadrature on each side, explicitly documenting that display convention; no inter-bin independence should be implied.

The third panel would show local MC relative precision or effective count, alongside signal yield, to distinguish discrepancies in irrelevant bins from potentially consequential ones. MC errors shown here are diagnostic and must be labeled as absent from the current Combine likelihood. Avoid using S/sqrt(B) as an exact limit or significance, particularly in sparse bins.

### B. Revised limit comparison

Use neutralino mass horizontally and cross-section limit in pb vertically at local stop mass 1250 GeV. Show local median expected and published median expected map values, with a ratio panel. This is clearer than a wide colored strip for one stop mass. The local expected bands can be shown if clearly identified; published expected bands should be added only from appropriate cross-section-band data, not inferred from exclusion-contour coordinates.

Record the published map coordinate used for each point and state that the published map is sampled, not an exact common simulated grid. An observed CMS curve can appear separately for context, but should not define the headline sensitivity ratio. With current inputs this plot would show ratios about 1.54 near 100 GeV, peaking around 1.63 near 500 GeV, then approaching 1.00 near 1100 GeV. It must not be presented as a validated reproduction until the selection and modeling differences are resolved.

### C. Jet-count validation after correction

Before interpreting a corrected limit, compare selected JetN distributions before and after the integer-range fix, and measure the restored signal/background yield in each search bin. The validation should specifically demonstrate nonzero acceptance for 3, 5, 7, and 9 jets and preservation of the intended kinematic exclusions. No corrected distribution can be recovered reliably by multiplying existing bin yields by a common factor.

## Data obtained and remaining limits

The download succeeded without requiring a manual transfer. Source files are in `data/SUS-19-006/audit_2026-09-10/`, with URLs and SHA-256 checksums in `sources.json`. These include five pre-fit-yield tables, the T2tt acceptance/efficiency map, two squark cutflow tables, and the official Figure 14a ROOT file.

For manual access, use the [HEPData record](https://www.hepdata.net/record/ins1749379) and select the named tables, then download YAML/CSV/JSON. The [CMS publication page](https://cms-results.web.cern.ch/cms-results/public-results/publications/SUS-19-006/index.html) also provides the Figure 14a ROOT file. The existing local T2tt limits YAML already supplies expected and observed maps. No additional user download is needed for the comparisons in this report.

A full reproduction would require an appropriate signal template per search bin at the relevant mass/energy, the complete background likelihood/correlations and control-region treatment, and validated event reconstruction/weights. The acceptance map and published example cutflows do not supply that full model. The current audit does not claim to have validated every dataset cross section, reconstructed the historical executable configuration, diagnosed the rare QCD events, or measured the separate effects of the identified problems on the limit.

The next concrete action should be to correct and validate jet-count acceptance, inspect the high-weight QCD survivors, and rerun the affected analysis inputs. Uncertainty modeling should follow validation of those central yields. Neither the published yields nor a change of labels can repair the current event-selection loss.

## Audit artifacts

- `normalization_audit.json`: per-sample processed counts, metadata, normalization, flows, selected counts, and signal cutflow; exact reconstruction differences for all six templates.
- `background_comparison_174_bins.csv`: all 174 nominally aligned search bins, process yields, real/pseudo observations, and uncertainty components.
- `expected_limit_comparison.csv`: all 11 local expected limits, published expected values, ROOT cross-checks, and sampling coordinates.
- `input_sha256.json`: checksums of inspected source/configuration files, saved results, cards, and ROOT templates, so later edits can be distinguished from this audit.
- `card_integrity_checks.json`: rate/integral agreement, shape lengths, and artificial-data checks across all 11 cards; all process-bin contents are finite and nonnegative.
