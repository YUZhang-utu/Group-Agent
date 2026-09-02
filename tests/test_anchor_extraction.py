import numpy as np

from aidd_agent.anchor_extraction import atomic_sasa, donor_hydrogen_angle, protein_hbond_roles


def test_protein_backbone_and_sidechain_hbond_roles():
    assert protein_hbond_roles("ALA", "N") == (True, False)
    assert protein_hbond_roles("PRO", "N") == (False, False)
    assert protein_hbond_roles("ALA", "O") == (False, True)
    assert protein_hbond_roles("LYS", "NZ") == (True, False)
    assert protein_hbond_roles("ASP", "OD1") == (False, True)


def test_atomic_sasa_decreases_with_nearby_occluder():
    center = np.asarray([[0., 0., 0.]])
    isolated = atomic_sasa(center, ["O"], center, ["O"], target_occluder_indices=[0])
    occluded = atomic_sasa(center, ["O"], np.asarray([[0., 0., 0.], [0., 0., 3.]]),
                           ["O", "C"], target_occluder_indices=[0])
    assert 0 < occluded[0] < isolated[0]


def test_donor_hydrogen_acceptor_linear_angle():
    assert donor_hydrogen_angle(np.array([0., 0., 0.]), np.array([1., 0., 0.]),
                                np.array([2.8, 0., 0.])) == 180.0
