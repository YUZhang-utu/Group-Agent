import json

import numpy as np

import aidd_agent.anchor_extraction as anchor_extraction
from aidd_agent.anchor_extraction import atomic_sasa, donor_hydrogen_angle, protein_hbond_roles
from aidd_agent.screening_rerank import AnchorDefinition, QueryManifest


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


def test_extract_and_write_query_manifest_is_atomic(tmp_path, monkeypatch):
    expected = QueryManifest(
        "1X8B:824:A:901", "direct-hbond-atomcenter-v1",
        {"anchor": 1.5, "ordinary": 1.0, "solvent_exposed": 0.5},
        (AnchorDefinition("A0:HBA:0", (0,), "HBA", (1.0, 2.0, 3.0),
                          (), {"interaction": "test", "distance_angstrom": 2.8,
                               "angle_degrees": None, "burial": {}}, 1.5),))
    monkeypatch.setattr(anchor_extraction, "extract_query_manifest",
                        lambda *args, **kwargs: expected)
    output = tmp_path / "nested" / "query_manifest.json"

    result = anchor_extraction.extract_and_write_query_manifest(
        tmp_path / "1X8B.cif", tmp_path / "824.cif", "824",
        "1X8B:824:A:901", output)

    assert json.loads(output.read_text(encoding="utf-8")) == result
    assert result["query_id"] == "1X8B:824:A:901"
    assert len(result["anchors"]) == 1
    assert not output.with_name(output.name + ".partial").exists()
