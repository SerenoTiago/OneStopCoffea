import enum

import awkward as ak
from analyzer.core.adl import ADLBlock, ADLStatement
from analyzer.core.analysis_modules import AnalyzerModule
from analyzer.core.columns import Column
from attrs import define, field


class CutBasedWPs(str, enum.Enum):
    fail = "fail"
    veto = "veto"
    loose = "loose"
    medium = "medium"
    tight = "tight"


cut_bitmap_mapping = dict(fail=0, veto=0b001, loose=0b001, medium=0b010, tight=0b100)
cut_based_mapping = dict(fail=0, veto=1, loose=1, medium=2, tight=3)


def pass_cut_based_bitmap(bitmap, working_point: CutBasedWPs):
    required_bit = cut_bitmap_mapping[working_point]
    if required_bit == 0:
        return ak.ones_like(bitmap, dtype=bool)
    return (bitmap & required_bit) != 0


def pass_cut_based_id(cut_based, working_point: CutBasedWPs):
    return cut_based >= cut_based_mapping[working_point]


def pass_photon_id(photons, working_point: CutBasedWPs):
    if "cutBasedBitmap" in photons.fields:
        return pass_cut_based_bitmap(photons.cutBasedBitmap, working_point)
    if "cutBased" in photons.fields:
        return pass_cut_based_id(photons.cutBased, working_point)
    raise AttributeError(
        "Photon collection has neither 'cutBasedBitmap' nor 'cutBased' fields."
    )


@define
class PhotonMaker(AnalyzerModule):
    """
    Select photons based on kinematics, VID cut-based ID, and optional isolation.

    Parameters
    ----------
    input_col : Column
        Column containing the input photon collection.
    output_col : Column
        Column where the selected photons will be stored.
    working_point : CutBasedWPs
        Cut-based ID working point (fail, veto, loose, medium, tight).
    min_pt : float, optional
        Minimum transverse momentum in GeV, by default 100.
    max_abs_eta : float, optional
        Maximum absolute pseudorapidity, by default 2.5.
    require_electron_veto : bool, optional
        Whether to require the NanoAOD electronVeto flag, by default True.
    max_pf_rel_iso03_all : float or None, optional
        Optional maximum PF relative isolation in a dR=0.3 cone.
    max_pf_rel_iso03_chg : float or None, optional
        Optional maximum charged PF relative isolation in a dR=0.3 cone.
    """

    input_col: Column
    output_col: Column
    working_point: CutBasedWPs
    min_pt: float = 100.0
    max_abs_eta: float = 2.5
    require_electron_veto: bool = True
    max_pf_rel_iso03_all: float | None = None
    max_pf_rel_iso03_chg: float | None = None

    __corrections: dict = field(factory=dict)

    def run(self, columns, params):
        photons = columns[self.input_col]
        passed = (photons.pt > self.min_pt) & (abs(photons.eta) < self.max_abs_eta)
        passed = passed & pass_photon_id(photons, self.working_point)

        if self.require_electron_veto:
            passed = passed & photons.electronVeto
        if self.max_pf_rel_iso03_all is not None:
            passed = passed & (photons.pfRelIso03_all < self.max_pf_rel_iso03_all)
        if self.max_pf_rel_iso03_chg is not None:
            passed = passed & (photons.pfRelIso03_chg < self.max_pf_rel_iso03_chg)

        columns[self.output_col] = photons[passed]
        return columns, []

    def inputs(self, metadata):
        return [self.input_col]

    def outputs(self, metadata):
        return [self.output_col]

    def adlExport(self, metadata):
        statements = [
            ADLStatement("take", self.input_col.adl_name),
            ADLStatement("select", f"pt > {self.min_pt}"),
            ADLStatement("select", f"abs(eta) < {self.max_abs_eta}"),
            ADLStatement(
                "select",
                f"cutBasedBitmap passes {self.working_point} or cutBased >= {cut_based_mapping[self.working_point]}",
            ),
        ]
        if self.require_electron_veto:
            statements.append(ADLStatement("select", "electronVeto == True"))
        if self.max_pf_rel_iso03_all is not None:
            statements.append(
                ADLStatement(
                    "select", f"pfRelIso03_all < {self.max_pf_rel_iso03_all}"
                )
            )
        if self.max_pf_rel_iso03_chg is not None:
            statements.append(
                ADLStatement(
                    "select", f"pfRelIso03_chg < {self.max_pf_rel_iso03_chg}"
                )
            )

        return [
            ADLBlock(
                block_type="object",
                name=self.output_col.adl_name,
                statements=statements,
            )
        ]
