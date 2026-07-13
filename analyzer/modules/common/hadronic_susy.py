from __future__ import annotations

import logging

import awkward as ak
import numpy as np
from attrs import define, field

from analyzer.core.analysis_modules import AnalyzerModule
from analyzer.core.columns import Column, addSelection


logger = logging.getLogger("analyzer.modules")

INVALID_SEARCH_BIN = -1

KINEMATIC_BINS = (
    (300.0, 350.0, 300.0, 600.0),
    (300.0, 350.0, 600.0, 1200.0),
    (300.0, 350.0, 1200.0, None),
    (350.0, 600.0, 350.0, 600.0),
    (350.0, 600.0, 600.0, 1200.0),
    (350.0, 600.0, 1200.0, None),
    (600.0, 850.0, 600.0, 1200.0),
    (600.0, 850.0, 1200.0, None),
    (850.0, None, 850.0, 1700.0),
    (850.0, None, 1700.0, None),
)

NJET_BINS = (
    (2, 3),
    (4, 5),
    (6, 7),
    (8, 9),
    (10, None),
)


def _in_range(values, low, high):
    passed = values >= low
    if high is not None:
        passed = passed & (values < high)
    return passed


def angle_to_minus_pi_pi(angle):
    return (angle + np.pi) % (2 * np.pi) - np.pi


def htmiss_components(jet_pt, jet_phi):
    px = ak.sum(jet_pt * np.cos(jet_phi), axis=1)
    py = ak.sum(jet_pt * np.sin(jet_phi), axis=1)
    return -px, -py


def htmiss_pt_phi(jets):
    htmiss_x, htmiss_y = htmiss_components(jets.pt, jets.phi)
    htmiss = np.hypot(htmiss_x, htmiss_y)
    phi = np.arctan2(htmiss_y, htmiss_x)
    return htmiss, phi


def delta_phi(phi_1, phi_2):
    return np.abs(angle_to_minus_pi_pi(phi_1 - phi_2))


def leading_jet_delta_phi(htmiss_phi, jets, max_jets=4):
    jets = jets[ak.argsort(jets.pt, axis=1, ascending=False)]
    padded = ak.pad_none(jets.phi, max_jets, axis=1)
    ret = []
    for idx in range(max_jets):
        dphi = delta_phi(htmiss_phi, padded[:, idx])
        ret.append(ak.fill_none(dphi, np.nan))
    return ret


def njet_bin_index(njet):
    ret = np.full(len(njet), INVALID_SEARCH_BIN, dtype=np.int64)
    for idx, (low, high) in enumerate(NJET_BINS):
        ret = np.where(_in_range(njet, low, high), idx, ret)
    return ret


def nbjet_bin_index(njet, nbjet):
    ret = np.full(len(njet), INVALID_SEARCH_BIN, dtype=np.int64)

    low_njet = (njet >= 2) & (njet <= 3)
    ret = np.where(low_njet & (nbjet == 0), 0, ret)
    ret = np.where(low_njet & (nbjet == 1), 1, ret)
    ret = np.where(low_njet & (nbjet >= 2), 2, ret)

    high_njet = njet >= 4
    ret = np.where(high_njet & (nbjet == 0), 0, ret)
    ret = np.where(high_njet & (nbjet == 1), 1, ret)
    ret = np.where(high_njet & (nbjet == 2), 2, ret)
    ret = np.where(high_njet & (nbjet >= 3), 3, ret)
    return ret


def n_nbjet_bins_for_njet_index(njet_idx):
    return np.where(njet_idx == 0, 3, 4)


def kinematic_bin_index(htmiss, ht):
    ret = np.full(len(ht), INVALID_SEARCH_BIN, dtype=np.int64)
    for idx, (mht_low, mht_high, ht_low, ht_high) in enumerate(KINEMATIC_BINS):
        passed = _in_range(htmiss, mht_low, mht_high) & _in_range(ht, ht_low, ht_high)
        ret = np.where(passed, idx, ret)
    return ret


def search_bin_id(njet, nbjet, ht, htmiss):
    """
    Assign CMS-style hadronic SUSY search-region bins.

    Valid search bins are numbered 1 through 174. Events outside the published
    bin definitions are assigned INVALID_SEARCH_BIN (-1). The binning includes
    Njet, Nb, HT, HTmiss, the Njet >= 8 dropped kinematic bins, and the strict
    HTmiss < HT requirement. Lepton/photon/isolated-track vetoes and dPhi cuts
    are event preselection cuts and are intentionally not encoded here.
    """
    njet = np.asarray(njet)
    nbjet = np.asarray(nbjet)
    ht = np.asarray(ht)
    htmiss = np.asarray(htmiss)

    njet_idx = njet_bin_index(njet)
    nbjet_idx = nbjet_bin_index(njet, nbjet)
    kin_idx = kinematic_bin_index(htmiss, ht)

    valid = (
        (njet_idx >= 0)
        & (nbjet_idx >= 0)
        & (kin_idx >= 0)
        & (htmiss < ht)
        & ~((njet >= 8) & ((kin_idx == 0) | (kin_idx == 3)))
    )

    ret = np.full(len(njet), INVALID_SEARCH_BIN, dtype=np.int64)
    offset = 0
    for this_njet_idx in range(len(NJET_BINS)):
        n_b_bins = 3 if this_njet_idx == 0 else 4
        allowed_kin = [idx for idx in range(len(KINEMATIC_BINS))]
        if this_njet_idx >= 3:
            allowed_kin = [idx for idx in allowed_kin if idx not in (0, 3)]
        for this_nbjet_idx in range(n_b_bins):
            for this_kin_idx in allowed_kin:
                offset += 1
                in_bin = (
                    valid
                    & (njet_idx == this_njet_idx)
                    & (nbjet_idx == this_nbjet_idx)
                    & (kin_idx == this_kin_idx)
                )
                ret = np.where(in_bin, offset, ret)
    return ret


def allowed_search_bins():
    bins = []
    for njet in (2, 4, 6, 8, 10):
        nbjets = (0, 1, 2) if njet <= 3 else (0, 1, 2, 3)
        for nbjet in nbjets:
            for kin_idx, (mht_low, mht_high, ht_low, ht_high) in enumerate(
                KINEMATIC_BINS
            ):
                if njet >= 8 and kin_idx in (0, 3):
                    continue
                htmiss = mht_low + 1.0
                ht = ht_low + 1.0
                if htmiss >= ht:
                    ht = htmiss + 1.0
                bins.append(search_bin_id([njet], [nbjet], [ht], [htmiss])[0])
    return sorted(set(int(x) for x in bins if x != INVALID_SEARCH_BIN))


def get_metadata_btag_info(metadata, discriminator_col=None, wp_name=None):
    btag_meta = metadata["era"].get("btag_scale_factors", {})
    tagger = discriminator_col or btag_meta.get("tagger")
    correction_name = btag_meta.get("correction_name", "deepJet_wp_values")
    working_point = wp_name or "M"
    if tagger is None:
        tagger = "btagDeepFlavB"
    return tagger, correction_name, working_point


@define
class HTMiss(AnalyzerModule):
    input_col: Column
    output_col: Column = field(factory=lambda: Column("HTMiss"))
    output_phi_col: Column = field(factory=lambda: Column("HTMissPhi"))

    def run(self, columns, params):
        jets = columns[self.input_col]
        htmiss, phi = htmiss_pt_phi(jets)
        columns[self.output_col] = htmiss
        columns[self.output_phi_col] = phi
        return columns, []

    def inputs(self, metadata):
        return [self.input_col]

    def outputs(self, metadata):
        return [self.output_col, self.output_phi_col]


@define
class HTMissDeltaPhi(AnalyzerModule):
    jet_col: Column
    htmiss_phi_col: Column
    output_prefix: str = "DeltaPhiHTMissJet"
    max_jets: int = 4

    def run(self, columns, params):
        jets = columns[self.jet_col]
        htmiss_phi = columns[self.htmiss_phi_col]
        for idx, dphi in enumerate(
            leading_jet_delta_phi(htmiss_phi, jets, self.max_jets)
        ):
            columns[Column(f"{self.output_prefix}{idx + 1}")] = dphi
        return columns, []

    def inputs(self, metadata):
        return [self.jet_col, self.htmiss_phi_col]

    def outputs(self, metadata):
        return [
            Column(f"{self.output_prefix}{idx + 1}") for idx in range(self.max_jets)
        ]


@define
class BJetSelector(AnalyzerModule):
    input_col: Column
    output_col: Column
    discriminator_col: str | None = None
    working_point: str = "M"

    __corrections: dict = field(factory=dict)

    def run(self, columns, params):
        tagger, correction_name, wp_name = get_metadata_btag_info(
            columns.metadata, self.discriminator_col, self.working_point
        )
        jets = columns[self.input_col]
        threshold = self.getWP(columns.metadata, correction_name, wp_name)
        columns[self.output_col] = jets[jets[tagger] > threshold]
        return columns, []

    def getWP(self, metadata, correction_name, wp_name):
        file_path = metadata["era"]["btag_scale_factors"]["file"]
        key = (file_path, correction_name, wp_name)
        if key in self.__corrections:
            return self.__corrections[key]
        import correctionlib

        cset = correctionlib.CorrectionSet.from_file(file_path)
        threshold = cset[correction_name].evaluate(wp_name)
        self.__corrections[key] = threshold
        return threshold

    def inputs(self, metadata):
        tagger, _, _ = get_metadata_btag_info(
            metadata, self.discriminator_col, self.working_point
        )
        return [self.input_col, self.input_col + tagger]

    def outputs(self, metadata):
        return [self.output_col]


@define
class SearchBin(AnalyzerModule):
    """
    Add a 1-based hadronic SUSY search-bin column.

    Valid bins are 1..174 and invalid/unassigned events are -1. The assignment
    follows search_bin_id and therefore applies only the bin-definition parts of
    the selection, including strict HTmiss < HT. Other preselection cuts must be
    represented by separate selection modules before histogramming.
    """

    njet_col: Column
    nbjet_col: Column
    ht_col: Column
    htmiss_col: Column
    output_col: Column = field(factory=lambda: Column("SearchBin"))

    def run(self, columns, params):
        columns[self.output_col] = search_bin_id(
            ak.to_numpy(columns[self.njet_col]),
            ak.to_numpy(columns[self.nbjet_col]),
            ak.to_numpy(columns[self.ht_col]),
            ak.to_numpy(columns[self.htmiss_col]),
        )
        return columns, []

    def inputs(self, metadata):
        return [self.njet_col, self.nbjet_col, self.ht_col, self.htmiss_col]

    def outputs(self, metadata):
        return [self.output_col]


@define
class ScalarFilter(AnalyzerModule):
    selection_name: str
    input_col: Column
    min_value: float | None = None
    max_value: float | None = None

    def run(self, columns, params):
        values = columns[self.input_col]
        sel = ak.ones_like(values, dtype=bool)
        if self.min_value is not None:
            sel = sel & (values > self.min_value)
        if self.max_value is not None:
            sel = sel & (values < self.max_value)
        addSelection(columns, self.selection_name, sel)
        return columns, []

    def inputs(self, metadata):
        return [self.input_col]

    def outputs(self, metadata):
        return [Column(("Selection", self.selection_name))]


@define
class HTMissLessThanHT(AnalyzerModule):
    ht_col: Column
    htmiss_col: Column
    selection_name: str = "htmiss_lt_ht"

    def run(self, columns, params):
        addSelection(
            columns, self.selection_name, columns[self.htmiss_col] < columns[self.ht_col]
        )
        return columns, []

    def inputs(self, metadata):
        return [self.ht_col, self.htmiss_col]

    def outputs(self, metadata):
        return [Column(("Selection", self.selection_name))]


@define
class HTMissDeltaPhiSelection(AnalyzerModule):
    njet_col: Column
    dphi_prefix: str = "DeltaPhiHTMissJet"
    selection_name: str = "htmiss_dphi"

    def run(self, columns, params):
        njet = columns[self.njet_col]
        dphi1 = columns[Column(f"{self.dphi_prefix}1")]
        dphi2 = columns[Column(f"{self.dphi_prefix}2")]
        dphi3 = columns[Column(f"{self.dphi_prefix}3")]
        dphi4 = columns[Column(f"{self.dphi_prefix}4")]
        sel = (dphi1 > 0.5) & (dphi2 > 0.5)
        sel = sel & ((njet < 3) | (dphi3 > 0.3))
        sel = sel & ((njet < 4) | (dphi4 > 0.3))
        addSelection(columns, self.selection_name, sel)
        return columns, []

    def inputs(self, metadata):
        return [self.njet_col] + [
            Column(f"{self.dphi_prefix}{idx}") for idx in range(1, 5)
        ]

    def outputs(self, metadata):
        return [Column(("Selection", self.selection_name))]


@define
class PassAllSelection(AnalyzerModule):
    selection_name: str
    warning: str | None = None

    def run(self, columns, params):
        if self.warning:
            logger.warning(self.warning)
        addSelection(
            columns,
            self.selection_name,
            ak.ones_like(columns[Column("event")], dtype=bool),
        )
        return columns, []

    def inputs(self, metadata):
        return [Column("event")]

    def outputs(self, metadata):
        return [Column(("Selection", self.selection_name))]
