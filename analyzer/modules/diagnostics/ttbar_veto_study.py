"""Isolated diagnostics for the hadronic-SUSY lepton and track vetoes.

This module intentionally does not filter events.  It evaluates several veto
definitions on the same input chunk and stores compact weighted/unweighted
histograms, including paired baseline-versus-variation migration counts.
"""

from __future__ import annotations

import awkward as ak
import hist
import numpy as np
from attrs import define

from analyzer.core.analysis_modules import AnalyzerModule
from analyzer.core.columns import Column
from analyzer.core.results import Histogram, UnscaledHistogram
from analyzer.modules.common.hadronic_susy import delta_phi


VARIATIONS = (
    "baseline",
    "electron_veto_id",
    "muon_medium_no_pfiso",
    "track_pfmet",
    "all_paper_like",
    "muon_tight_no_pfiso",
)

STAGES = (
    "input",
    "met_mht_trigger",
    "zero_electron",
    "zero_muon",
    "zero_photon",
    "njets",
    "ht_gt_300",
    "htmiss_gt_300",
    "htmiss_lt_ht",
    "htmiss_dphi",
    "zero_iso_track",
    "search_bin_valid",
)

DECAYS = (
    "direct_e",
    "direct_mu",
    "tau_hadronic",
    "tau_leptonic",
    "unclassified",
    "ambiguous",
)

OUTCOMES = (
    "survive",
    "baseline_survive_variant_fail",
    "variant_survive_baseline_fail",
)

TAU_TRACK_REASONS = (
    "no_reconstructed_isotrack",
    "no_basic_candidate",
    "fails_track_quality",
    "fails_isolation",
    "mt_ge_100",
    "other",
)


def _eventwise_take(values, indices):
    """Take a jagged value at a jagged, event-local index."""
    safe = ak.where(indices >= 0, indices, 0)
    counts = ak.to_numpy(ak.num(values, axis=1))
    offsets = np.concatenate(([0], np.cumsum(counts[:-1])))
    global_indices = ak.flatten(safe + offsets[:, None])
    selected = ak.flatten(values)[global_indices]
    return ak.unflatten(selected, counts)


def _ancestor_flags(genparts, target_ids, max_depth=30):
    """Return per-particle flags for ancestors with requested absolute IDs."""
    parent_index = genparts.genPartIdxMother
    active_index = parent_index
    flags = {pid: ak.zeros_like(genparts.pdgId, dtype=bool) for pid in target_ids}
    for _ in range(max_depth):
        valid = active_index >= 0
        ancestor_pdg = ak.fill_none(_eventwise_take(genparts.pdgId, active_index), 0)
        for pid in target_ids:
            flags[pid] = flags[pid] | (valid & (abs(ancestor_pdg) == pid))
        next_index = ak.fill_none(
            _eventwise_take(genparts.genPartIdxMother, active_index), -1
        )
        active_index = ak.where(valid, next_index, -1)
    return flags


def _first_distinct_parent_pdg(genparts, max_depth=20):
    """Find the first parent whose absolute PDG ID differs from the child."""
    child_abs_id = abs(genparts.pdgId)
    parent_index = genparts.genPartIdxMother
    for _ in range(max_depth):
        valid = parent_index >= 0
        parent_pdg = ak.fill_none(_eventwise_take(genparts.pdgId, parent_index), 0)
        same_copy = valid & (abs(parent_pdg) == child_abs_id)
        next_index = ak.fill_none(
            _eventwise_take(genparts.genPartIdxMother, parent_index), -1
        )
        parent_index = ak.where(same_copy, next_index, parent_index)
    return ak.fill_none(_eventwise_take(genparts.pdgId, parent_index), 0)


def _has_marked_ancestor(genparts, marked, max_depth=30):
    """Test whether each particle descends from any marked particle."""
    ancestor_index = genparts.genPartIdxMother
    found = ak.zeros_like(genparts.pdgId, dtype=bool)
    for _ in range(max_depth):
        valid = ancestor_index >= 0
        found = found | (
            valid
            & ak.fill_none(_eventwise_take(marked, ancestor_index), False)
        )
        next_index = ak.fill_none(
            _eventwise_take(genparts.genPartIdxMother, ancestor_index), -1
        )
        ancestor_index = ak.where(valid, next_index, -1)
    return found


def _classify_w_decay(genparts):
    """Classify the leptonic W decay using generator ancestry."""
    abs_id = abs(genparts.pdgId)
    distinct_parent_abs_id = abs(_first_distinct_parent_pdg(genparts))
    direct_w_e = (abs_id == 11) & (distinct_parent_abs_id == 24)
    direct_w_mu = (abs_id == 13) & (distinct_parent_abs_id == 24)
    direct_w_tau = (abs_id == 15) & (distinct_parent_abs_id == 24)

    direct_e = ak.any(direct_w_e, axis=1)
    direct_mu = ak.any(direct_w_mu, axis=1)
    w_tau = ak.any(direct_w_tau, axis=1)
    from_direct_w_tau = _has_marked_ancestor(genparts, direct_w_tau)
    tau_lepton = ak.any(
        ((abs_id == 11) | (abs_id == 13)) & from_direct_w_tau, axis=1
    )
    tau_leptonic = w_tau & tau_lepton
    tau_hadronic = w_tau & ~tau_lepton

    category_count = (
        ak.values_astype(direct_e, np.int8)
        + ak.values_astype(direct_mu, np.int8)
        + ak.values_astype(tau_leptonic, np.int8)
        + ak.values_astype(tau_hadronic, np.int8)
    )
    labels = np.full(len(category_count), "unclassified", dtype="U16")
    labels[ak.to_numpy(tau_hadronic)] = "tau_hadronic"
    labels[ak.to_numpy(tau_leptonic)] = "tau_leptonic"
    labels[ak.to_numpy(direct_mu)] = "direct_mu"
    labels[ak.to_numpy(direct_e)] = "direct_e"
    labels[ak.to_numpy(category_count > 1)] = "ambiguous"
    return labels


def _total_weight(columns, template):
    if "Weights" not in columns.fields or not columns["Weights"].fields:
        return ak.ones_like(template, dtype=float)
    weights = columns["Weights"]
    fields = iter(weights.fields)
    result = weights[next(fields)]
    for name in fields:
        result = result * weights[name]
    return result


def _nb_category(nbjets):
    values = ak.to_numpy(nbjets)
    labels = np.full(len(values), ">=2", dtype="U3")
    labels[values == 0] = "0"
    labels[values == 1] = "1"
    return labels


def _make_yield_hist(storage):
    return hist.Hist(
        hist.axis.StrCategory(VARIATIONS, name="veto_variation"),
        hist.axis.StrCategory(("0", "1", ">=2"), name="nb"),
        hist.axis.StrCategory(DECAYS, name="w_decay"),
        hist.axis.StrCategory(OUTCOMES, name="outcome"),
        storage=storage,
    )


def _make_cutflow_hist(storage):
    return hist.Hist(
        hist.axis.StrCategory(VARIATIONS, name="veto_variation"),
        hist.axis.StrCategory(STAGES, name="stage"),
        hist.axis.StrCategory(("0", "1", ">=2"), name="nb"),
        hist.axis.StrCategory(DECAYS, name="w_decay"),
        storage=storage,
    )


def _make_tau_track_hist(storage):
    return hist.Hist(
        hist.axis.StrCategory(("0", "1", ">=2"), name="nb"),
        hist.axis.StrCategory(TAU_TRACK_REASONS, name="reason"),
        storage=storage,
    )


def _tau_track_failure_reason(tracks, htmiss, htmiss_phi):
    """Classify why an event has no vetoing isolated track.

    This is an event-level candidate audit, not a reconstructed-to-generator
    tau match.  The hierarchy reports the latest selection stage reached by
    any NanoAOD IsoTrack in the event.
    """
    is_lepton = (abs(tracks.pdgId) == 11) | (abs(tracks.pdgId) == 13)
    min_pt = ak.where(is_lepton, 5.0, 10.0)
    max_iso = ak.where(is_lepton, 0.2, 0.1)
    basic = (
        (tracks.pt > min_pt)
        & (abs(tracks.eta) < 2.4)
        & (tracks.fromPV >= 2)
        & (tracks.charge != 0)
    )
    quality = basic & tracks.isHighPurityTrack & (abs(tracks.dxy) < 0.2) & (abs(tracks.dz) < 0.1)
    isolated = quality & (tracks.pfRelIso03_chg < max_iso)
    mt = np.sqrt(
        2.0
        * tracks.pt
        * htmiss[:, None]
        * (1.0 - np.cos(delta_phi(tracks.phi, htmiss_phi[:, None])))
    )
    mt_pass = isolated & (mt < 100.0)

    has_raw = ak.num(tracks, axis=1) > 0
    has_basic = ak.any(basic, axis=1)
    has_quality = ak.any(quality, axis=1)
    has_isolated = ak.any(isolated, axis=1)
    has_mt_pass = ak.any(mt_pass, axis=1)
    labels = np.full(len(has_raw), "other", dtype="U28")
    labels[~ak.to_numpy(has_raw)] = "no_reconstructed_isotrack"
    labels[ak.to_numpy(has_raw & ~has_basic)] = "no_basic_candidate"
    labels[ak.to_numpy(has_basic & ~has_quality)] = "fails_track_quality"
    labels[ak.to_numpy(has_quality & ~has_isolated)] = "fails_isolation"
    labels[ak.to_numpy(has_isolated & ~has_mt_pass)] = "mt_ge_100"
    return labels


@define
class TTbarVetoDiagnostic(AnalyzerModule):
    """Evaluate correlated veto variations without changing event selection."""

    def inputs(self, metadata):
        return [
            Column("GenPart"),
            Column("NBJet"),
            Column("Weights"),
            Column("tightElectron"),
            Column("vetoElectron"),
            Column("tightMuon"),
            Column("mediumMuonNoPFIso"),
            Column("tightMuonNoPFIso"),
            Column("goodIsoTrack"),
            Column("goodIsoTrackPFMET"),
            Column("IsoTrack"),
            Column("HTMiss"),
            Column("HTMissPhi"),
            *[Column(("Selection", name)) for name in STAGES[1:] if name not in (
                "zero_electron", "zero_muon", "zero_iso_track"
            )],
        ]

    def outputs(self, metadata):
        return []

    def run(self, columns, params):
        template = columns[Column(("Selection", "met_mht_trigger"))]
        weights = _total_weight(columns, template)
        decay = _classify_w_decay(columns["GenPart"])
        nb = _nb_category(columns["NBJet"])

        electron_masks = {
            "baseline": ak.num(columns["tightElectron"], axis=1) == 0,
            "paper": ak.num(columns["vetoElectron"], axis=1) == 0,
        }
        muon_masks = {
            "baseline": ak.num(columns["tightMuon"], axis=1) == 0,
            "medium": ak.num(columns["mediumMuonNoPFIso"], axis=1) == 0,
            "tight_no_pfiso": ak.num(columns["tightMuonNoPFIso"], axis=1) == 0,
        }
        track_masks = {
            "baseline": ak.num(columns["goodIsoTrack"], axis=1) == 0,
            "pfmet": ak.num(columns["goodIsoTrackPFMET"], axis=1) == 0,
        }
        veto_masks = {
            "baseline": (electron_masks["baseline"], muon_masks["baseline"], track_masks["baseline"]),
            "electron_veto_id": (electron_masks["paper"], muon_masks["baseline"], track_masks["baseline"]),
            "muon_medium_no_pfiso": (electron_masks["baseline"], muon_masks["medium"], track_masks["baseline"]),
            "track_pfmet": (electron_masks["baseline"], muon_masks["baseline"], track_masks["pfmet"]),
            "all_paper_like": (electron_masks["paper"], muon_masks["medium"], track_masks["pfmet"]),
            "muon_tight_no_pfiso": (electron_masks["baseline"], muon_masks["tight_no_pfiso"], track_masks["baseline"]),
        }

        selection_masks = {
            name: columns[Column(("Selection", name))]
            for name in STAGES[1:]
            if name not in ("zero_electron", "zero_muon", "zero_iso_track")
        }

        weighted_yields = _make_yield_hist(hist.storage.Weight())
        raw_yields = _make_yield_hist(hist.storage.Double())
        weighted_cutflow = _make_cutflow_hist(hist.storage.Weight())
        raw_cutflow = _make_cutflow_hist(hist.storage.Double())
        weighted_tau_track = _make_tau_track_hist(hist.storage.Weight())
        raw_tau_track = _make_tau_track_hist(hist.storage.Double())

        final_masks = {}
        for variation, (electron_mask, muon_mask, track_mask) in veto_masks.items():
            stage_masks = {
                **selection_masks,
                "zero_electron": electron_mask,
                "zero_muon": muon_mask,
                "zero_iso_track": track_mask,
            }
            cumulative = ak.ones_like(template, dtype=bool)
            for stage in STAGES:
                if stage != "input":
                    cumulative = cumulative & stage_masks[stage]
                weighted_cutflow.fill(
                    veto_variation=variation,
                    stage=stage,
                    nb=nb[cumulative],
                    w_decay=decay[cumulative],
                    weight=weights[cumulative],
                )
                raw_cutflow.fill(
                    veto_variation=variation,
                    stage=stage,
                    nb=nb[cumulative],
                    w_decay=decay[cumulative],
                )
            final_masks[variation] = cumulative

        baseline = final_masks["baseline"]
        for variation, mask in final_masks.items():
            outcome_masks = {
                "survive": mask,
                "baseline_survive_variant_fail": baseline & ~mask,
                "variant_survive_baseline_fail": mask & ~baseline,
            }
            for outcome, outcome_mask in outcome_masks.items():
                weighted_yields.fill(
                    veto_variation=variation,
                    nb=nb[outcome_mask],
                    w_decay=decay[outcome_mask],
                    outcome=outcome,
                    weight=weights[outcome_mask],
                )
                raw_yields.fill(
                    veto_variation=variation,
                    nb=nb[outcome_mask],
                    w_decay=decay[outcome_mask],
                    outcome=outcome,
                )

        # Diagnose baseline-surviving hadronic-tau events.  These reasons refer
        # to the reconstructed IsoTrack candidates in the event; NanoAOD does
        # not provide a direct IsoTrack-to-generator-tau association here.
        tau_reason = _tau_track_failure_reason(
            columns["IsoTrack"], columns["HTMiss"], columns["HTMissPhi"]
        )
        tau_hadronic_survivor = baseline & (decay == "tau_hadronic")
        weighted_tau_track.fill(
            nb=nb[tau_hadronic_survivor],
            reason=tau_reason[tau_hadronic_survivor],
            weight=weights[tau_hadronic_survivor],
        )
        raw_tau_track.fill(
            nb=nb[tau_hadronic_survivor],
            reason=tau_reason[tau_hadronic_survivor],
        )

        return columns, [
            Histogram("VetoDiagnosticWeightedYields", axes=[], histogram=weighted_yields),
            UnscaledHistogram("VetoDiagnosticRawYields", axes=[], histogram=raw_yields),
            Histogram("VetoDiagnosticWeightedCutflow", axes=[], histogram=weighted_cutflow),
            UnscaledHistogram("VetoDiagnosticRawCutflow", axes=[], histogram=raw_cutflow),
            Histogram("TauTrackFailureReasonsWeighted", axes=[], histogram=weighted_tau_track),
            UnscaledHistogram("TauTrackFailureReasonsRaw", axes=[], histogram=raw_tau_track),
        ]
