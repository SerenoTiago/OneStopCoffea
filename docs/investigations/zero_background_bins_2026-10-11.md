# Zero/negative-background search bins (open caveat)

Date: 2026-10-11.  Status: **open, for advisor discussion.**  No floor, bin
merging, clipping, smoothing, or Run-2-derived background has been applied.

## Affected bins

With the complete current-selection backgrounds
(`analysis_products/production_runs/2026-10-07_012944_hadronic_susy/background_results`),
two of the 174 search bins have zero or negative total predicted background
at 137 fb^-1:

| Bin | Njet | Nb | HTmiss [GeV] | HT [GeV] | Total MC background |
|----:|-----|----|--------------|----------|---------------------|
| 171 | >=10 | >=3 | 600-850 | 600-1200 | 0 exactly |
| 174 | >=10 | >=3 | >=850 | >=1700 | -0.0036 |

## Cause

Sparse MC with signed (+-1) generator weights; not a missing sample, dropped
bin, or histogram bug.

- Bin 171: only two Z->nunu (PTNuNu600) MC events survive, one with weight +1
  and one with -1, so the yield cancels to exactly zero (sumw = 0, sumw2 = 2).
- Bin 174: one W+jets (PTLNu600toInf) event with weight -1 (-0.0048 events)
  and five Z->nunu events (+0.0012 events).
- All other background samples have no MC events in either bin.
- The older 09-28 background set (previous lepton vetoes) had 0.117 and 0.235
  events here from semileptonic ttbar; the current veto removes them.

The datacard therefore carries a negative W+jets nominal in bin 174 and a
slightly negative Asimov `data_obs` (the blinded expected limits do not use
`data_obs`).  Combine sets the negative nominal to zero when building the
workspace and reports no warning.

## Comparison with CMS-SUS-19-006 (Run 2, pre-fit)

The published total background is 0.1 events in both bins (lost-lepton 0,
Z->nunu 0.05, QCD 0.04/0.02) with upward statistical uncertainties of +3.0
and +1.3 events from zero-event lost-lepton control regions.  The published
central value is also about zero; the difference is that it carries a large
uncertainty.

## Treatment in the current datacards

The datacards now use the standard `* autoMCStats 0`.  In these two bins it
creates only per-process Gaussian MC-statistical parameters whose sizes come
from the stored sumw2: background 0 +- 0.0017 events (bin 171) and
0.0012 +- 0.0055 events (bin 174).  These are far below one event, so both
bins remain effectively background-free in the fit.

## Measured impact

Without MC-statistical nuisances, a closed-form reproduction of Combine's
median expected limit (agreeing with Combine to 0.19% median, 0.73% max over
all 198 points) gives:

| Variation | Median change in r_exp | Largest change |
|-----------|-----------------------:|---------------:|
| Drop bins 171 and 174 | +1.5% | +6.5% (mStop 1050, mLSP 800) |
| Use 0.1 events (Run-2 central) in each | +1.3% | +4.2% (mStop 1050, mLSP 800) |

So treating these bins as exactly background-free strengthens the expected
limits by about 1-2% typically and up to about 6% for some compressed points.
The signal in the two bins is a median 0.06% (max 0.24%) of the total.
Because autoMCStats assigns them sub-event uncertainties, this estimate is
expected to remain approximately valid with autoMCStats; it was not
re-measured with Combine.

For reference, enabling autoMCStats over the whole grid weakens the median
expected limit by +6.6% (range +5.1% to +15.9%, largest at mStop 600); that
change comes from the other low-statistics bins, not from these two.

## Questions for discussion

1. Should zero-MC bins receive an uncertainty reflecting the zero-event
   Poisson upper bound (as the Run-2 lost-lepton estimate effectively does)?
2. Should the Njet >= 10, Nb >= 3 high-HTmiss bins be merged or dropped for
   the Run-3 projection?
3. Should negative per-process nominals be clipped explicitly in the datacard
   rather than relying on Combine's internal handling?
