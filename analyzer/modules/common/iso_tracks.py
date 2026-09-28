from __future__ import annotations

import awkward as ak
import numpy as np
from attrs import define, field

from analyzer.core.adl import ADLBlock, ADLStatement
from analyzer.core.analysis_modules import AnalyzerModule
from analyzer.core.columns import Column
from analyzer.modules.common.hadronic_susy import delta_phi


def _abs_pdg_id_in(values, pdg_ids):
    selected = ak.zeros_like(values, dtype=bool)
    abs_values = abs(values)
    for pdg_id in pdg_ids:
        selected = selected | (abs_values == pdg_id)
    return selected


@define
class IsoTrackMaker(AnalyzerModule):
    """
    Select isolated tracks for the hadronic SUSY isolated-track veto.

    The thresholds follow the CMS SUSY prescription: tracks must satisfy
    |eta| < 2.4, mT < 100 GeV, and charged isolation in dR=0.3. PF electron
    and muon tracks use pT > 5 GeV and isolation < 0.2; other tracks use
    pT > 10 GeV and isolation < 0.1.
    """

    input_col: Column
    output_col: Column
    htmiss_col: Column = Column("HTMiss")
    htmiss_phi_col: Column = Column("HTMissPhi")
    max_abs_eta: float = 2.4
    max_mt: float = 100.0
    min_lepton_track_pt: float = 5.0
    min_hadron_track_pt: float = 10.0
    max_lepton_track_iso: float = 0.2
    max_hadron_track_iso: float = 0.1
    min_from_pv: int = 2
    max_abs_dxy: float | None = 0.2
    max_abs_dz: float | None = 0.1
    require_high_purity: bool = True
    require_pf_candidate: bool = False
    lepton_pdg_ids: tuple[int, ...] = (11, 13)

    __corrections: dict = field(factory=dict)

    def run(self, columns, params):
        tracks = columns[self.input_col]
        htmiss = columns[self.htmiss_col][:, None]
        htmiss_phi = columns[self.htmiss_phi_col][:, None]

        is_lepton_track = _abs_pdg_id_in(tracks.pdgId, self.lepton_pdg_ids)
        min_pt = ak.where(
            is_lepton_track, self.min_lepton_track_pt, self.min_hadron_track_pt
        )
        max_iso = ak.where(
            is_lepton_track, self.max_lepton_track_iso, self.max_hadron_track_iso
        )
        mt = np.sqrt(
            2.0
            * tracks.pt
            * htmiss
            * (1.0 - np.cos(delta_phi(tracks.phi, htmiss_phi)))
        )

        passed = (
            (tracks.pt > min_pt)
            & (abs(tracks.eta) < self.max_abs_eta)
            & (tracks.pfRelIso03_chg < max_iso)
            & (mt < self.max_mt)
            & (tracks.fromPV >= self.min_from_pv)
            & (tracks.charge != 0)
        )
        if self.require_high_purity:
            passed = passed & tracks.isHighPurityTrack
        if self.require_pf_candidate:
            passed = passed & tracks.isPFcand
        if self.max_abs_dxy is not None:
            passed = passed & (abs(tracks.dxy) < self.max_abs_dxy)
        if self.max_abs_dz is not None:
            passed = passed & (abs(tracks.dz) < self.max_abs_dz)

        columns[self.output_col] = tracks[passed]
        return columns, []

    def inputs(self, metadata):
        return [self.input_col, self.htmiss_col, self.htmiss_phi_col]

    def outputs(self, metadata):
        return [self.output_col]

    def adlExport(self, metadata):
        statements = [
            ADLStatement("take", self.input_col.adl_name),
            ADLStatement("select", f"abs(eta) < {self.max_abs_eta}"),
            ADLStatement("select", f"fromPV >= {self.min_from_pv}"),
            ADLStatement("select", "charge != 0"),
            ADLStatement(
                "select",
                f"pfRelIso03_chg < {self.max_lepton_track_iso} for PF e/mu tracks, otherwise < {self.max_hadron_track_iso}",
            ),
            ADLStatement(
                "select",
                f"pt > {self.min_lepton_track_pt} for PF e/mu tracks, otherwise > {self.min_hadron_track_pt}",
            ),
            ADLStatement("select", f"mT(track, HTMiss) < {self.max_mt}"),
        ]
        if self.require_high_purity:
            statements.append(ADLStatement("select", "isHighPurityTrack == True"))
        if self.require_pf_candidate:
            statements.append(ADLStatement("select", "isPFcand == True"))
        if self.max_abs_dxy is not None:
            statements.append(ADLStatement("select", f"abs(dxy) < {self.max_abs_dxy}"))
        if self.max_abs_dz is not None:
            statements.append(ADLStatement("select", f"abs(dz) < {self.max_abs_dz}"))

        return [
            ADLBlock(
                block_type="object",
                name=self.output_col.adl_name,
                statements=statements,
            )
        ]
