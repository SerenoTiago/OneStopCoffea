import awkward as ak
import numpy as np

from analyzer.modules.common.hadronic_susy import (
    INVALID_SEARCH_BIN,
    allowed_search_bins,
    delta_phi,
    htmiss_pt_phi,
    leading_jet_delta_phi,
    njet_bin_index,
    search_bin_id,
)


def test_htmiss_pt_phi():
    jets = ak.Array(
        [
            {"pt": [100.0, 100.0], "phi": [0.0, np.pi / 2]},
            {"pt": [50.0, 50.0], "phi": [0.0, np.pi]},
        ]
    )

    htmiss, phi = htmiss_pt_phi(jets)

    assert np.allclose(ak.to_numpy(htmiss), [np.sqrt(20000.0), 0.0], atol=1e-10)
    assert np.isclose(float(phi[0]), -3 * np.pi / 4)


def test_delta_phi_wraps_to_pi():
    assert np.isclose(delta_phi(0.1, 2 * np.pi - 0.1), 0.2)
    assert np.isclose(delta_phi(-np.pi + 0.1, np.pi - 0.1), 0.2)


def test_leading_jet_delta_phi_pads_missing_jets():
    jets = ak.Array(
        [
            {"pt": [1.0, 1.0, 1.0], "phi": [0.1, 1.1, 2.1]},
            {"pt": [1.0], "phi": [0.5]},
        ]
    )

    dphis = leading_jet_delta_phi(ak.Array([0.0, 0.0]), jets, max_jets=4)

    assert len(dphis) == 4
    assert np.allclose(ak.to_numpy(dphis[0]), [0.1, 0.5])
    assert np.isnan(float(dphis[1][1]))
    assert np.isnan(float(dphis[3][0]))


def test_leading_jet_delta_phi_sorts_by_descending_pt():
    jets = ak.Array(
        [
            {
                "pt": [20.0, 100.0, 50.0],
                "phi": [0.2, 1.0, 2.0],
            }
        ]
    )

    dphis = leading_jet_delta_phi(ak.Array([0.0]), jets, max_jets=3)

    assert np.isclose(float(dphis[0][0]), 1.0)
    assert np.isclose(float(dphis[1][0]), 2.0)
    assert np.isclose(float(dphis[2][0]), 0.2)


def test_search_bin_mapping_has_174_allowed_bins():
    bins = allowed_search_bins()

    assert len(bins) == 174
    assert bins == list(range(1, 175))


def test_search_bin_assignment_boundaries_and_invalid_regions():
    bins = search_bin_id(
        njet=[2, 2, 4, 8, 8, 10, 1],
        nbjet=[0, 2, 3, 0, 0, 3, 0],
        ht=[400.0, 1300.0, 1800.0, 400.0, 901.0, 1800.0, 900.0],
        htmiss=[325.0, 325.0, 900.0, 325.0, 900.0, 900.0, 500.0],
    )

    assert bins[0] == 1
    assert bins[1] == 23
    assert bins[2] == 70
    assert bins[3] == INVALID_SEARCH_BIN
    assert bins[4] != INVALID_SEARCH_BIN
    assert bins[5] == 174
    assert bins[6] == INVALID_SEARCH_BIN


def test_search_bin_rejects_htmiss_greater_than_ht():
    bins = search_bin_id(
        njet=[4, 4],
        nbjet=[0, 0],
        ht=[500.0, 500.0],
        htmiss=[700.0, 500.0],
    )

    assert bins[0] == INVALID_SEARCH_BIN
    assert bins[1] == INVALID_SEARCH_BIN


def test_jet_categories_include_both_integer_endpoints():
    assert njet_bin_index(np.arange(0, 13)).tolist() == [
        -1, -1, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 4
    ]
    # Published bin IDs for Nb=0, 350<=MHT<600 and 600<=HT<1200.
    assert search_bin_id(
        njet=list(range(2, 12)),
        nbjet=[0] * 10,
        ht=[800.0] * 10,
        htmiss=[400.0] * 10,
    ).tolist() == [5, 5, 35, 35, 75, 75, 113, 113, 145, 145]


def test_restored_odd_jet_counts_preserve_kinematic_exclusions():
    # The low-HT intervals at >=8 jets remain excluded for even and odd counts.
    assert search_bin_id(
        njet=[7, 8, 9, 10, 7, 8, 9, 10],
        nbjet=[0] * 8,
        ht=[400.0] * 4 + [500.0] * 4,
        htmiss=[325.0] * 4 + [400.0] * 4,
    ).tolist() == [71, -1, -1, -1, 74, -1, -1, -1]
    assert search_bin_id(
        njet=[3, 5, 7, 9],
        nbjet=[0] * 4,
        ht=[500.0, 500.0, 500.0, 500.0],
        htmiss=[500.0, 600.0, 500.0, 600.0],
    ).tolist() == [-1] * 4


def test_continuous_kinematic_edges_remain_half_open():
    assert search_bin_id(
        njet=[3] * 6,
        nbjet=[0] * 6,
        ht=[800.0, 800.0, 800.0, 800.0, 1000.0, 1700.0],
        htmiss=[299.0, 300.0, 350.0, 600.0, 850.0, 850.0],
    ).tolist() == [-1, 2, 5, 7, 9, 10]
