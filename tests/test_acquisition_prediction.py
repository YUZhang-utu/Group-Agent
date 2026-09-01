from __future__ import annotations

import json

import pytest

from aidd_agent.acquisition import acquire_alphafold_db_mmcif, acquire_rcsb_mmcif
from aidd_agent.prediction import (
    inspect_alphafold3_output, load_model_profile, prediction_command, write_alphafold3_input,
    write_prediction_input,
)


def test_acquire_rcsb_mmcif_writes_scoped_content_and_metadata(tmp_path):
    payload = b"data_1ABC\n_entry.id 1ABC\n#\n" + b" " * 40
    result = acquire_rcsb_mmcif(
        "1abc", tmp_path,
        fetch=lambda url: (payload, "chemical/x-mmcif"),
    )
    assert result["pdb_id"] == "1ABC"
    assert (tmp_path / "inputs" / "structures" / "1ABC.cif").read_bytes() == payload
    assert json.loads((tmp_path / "inputs" / "structures" / "1ABC.metadata.json").read_text())["sha256"] == result["sha256"]


def test_acquire_rejects_bad_id_and_non_cif(tmp_path):
    with pytest.raises(ValueError):
        acquire_rcsb_mmcif("../../x", tmp_path)
    with pytest.raises(ValueError, match="not the requested"):
        acquire_rcsb_mmcif("1ABC", tmp_path, fetch=lambda url: (b"<html>not found</html>", "text/html"))


def test_alphafold3_input_and_profile(tmp_path):
    output = tmp_path / "runs" / "af3.json"
    write_alphafold3_input("wee1", ["ACDEFG"], output, tmp_path, seeds=(7, 11))
    payload = json.loads(output.read_text())
    assert payload["dialect"] == "alphafold3"
    assert payload["modelSeeds"] == [7, 11]
    profile_path = tmp_path / "profile.json"
    profile_path.write_text(json.dumps({
        "backend": "alphafold3", "python": "python", "runner": "run.py",
        "model_parameters": "/models", "databases": "/db", "license_acknowledged": True,
    }))
    command = prediction_command(load_model_profile(profile_path), output, tmp_path / "result")
    assert command[:2] == ["python", "run.py"]
    assert "--model_dir=/models" in command


def test_alphafold3_requires_license_acknowledgement(tmp_path):
    profile = tmp_path / "profile.json"
    profile.write_text(json.dumps({
        "backend": "alphafold3", "python": "python", "runner": "run.py",
        "model_parameters": "/models", "databases": "/db", "license_acknowledged": False,
    }))
    with pytest.raises(ValueError, match="terms"):
        load_model_profile(profile)


def test_acquire_alphafold_db_model(tmp_path):
    payload = b"data_AF-P30291-F1-model_v4\n_entry.id AF-P30291-F1\n#\n" + b" " * 40
    result = acquire_alphafold_db_mmcif(
        "p30291", tmp_path, fetch=lambda url: (payload, "chemical/x-mmcif"))
    assert result["uniprot_id"] == "P30291"
    assert result["source_type"] == "predicted"
    assert (tmp_path / "inputs" / "structures" /
            "AF-P30291-F1-model_v4.cif").read_bytes() == payload
    with pytest.raises(ValueError, match="UniProt"):
        acquire_alphafold_db_mmcif("../../x", tmp_path)


def test_boltz2_input_profile_and_command(tmp_path):
    input_path = tmp_path / "runs" / "boltz" / "wee1.yaml"
    write_prediction_input("boltz2", "wee1", ["ACDEFG"], input_path, tmp_path)
    assert "id: A" in input_path.read_text()
    profile_path = tmp_path / "boltz.json"
    profile_path.write_text(json.dumps({
        "backend": "boltz2", "executable": "/env/bin/boltz",
        "cache": "/models/boltz2", "use_msa_server": False,
        "use_potentials": True,
    }))
    command = prediction_command(
        load_model_profile(profile_path), input_path, tmp_path / "output")
    assert command[:2] == ["/env/bin/boltz", "predict"]
    assert command[-1] == "--use_potentials"
    assert "--use_msa_server" not in command


def test_chai1_input_profile_and_command(tmp_path):
    input_path = tmp_path / "runs" / "chai" / "wee1.fasta"
    write_prediction_input("chai1", "wee1", ["ACDEFG"], input_path, tmp_path)
    assert input_path.read_text().startswith(">protein|wee1_A\nACDEFG")
    profile_path = tmp_path / "chai.json"
    profile_path.write_text(json.dumps({
        "backend": "chai1", "executable": "/env/bin/chai-lab",
        "use_msa_server": True, "use_templates_server": True,
    }))
    command = prediction_command(
        load_model_profile(profile_path), input_path, tmp_path / "output")
    assert command == ["/env/bin/chai-lab", "fold", str(input_path),
                       str(tmp_path / "output"), "--use-msa-server",
                       "--use-templates-server"]


def test_inspect_alphafold3_output_builds_hashed_manifest(tmp_path):
    output = tmp_path / "output" / "wee1"
    output.mkdir(parents=True)
    structure = output / "wee1_model.cif"
    structure.write_text("data_wee1\n#\n", encoding="utf-8")
    (output / "wee1_ranking_scores.csv").write_text(
        "seed,sample,ranking_score\n42,0,0.81\n42,1,0.86\n", encoding="utf-8")
    (output / "wee1_summary_confidences.json").write_text(
        json.dumps({"chain_iptm": [0.9]}), encoding="utf-8")
    input_path = tmp_path / "wee1.json"
    input_path.write_text("{}", encoding="utf-8")
    manifest = inspect_alphafold3_output(
        output.parent, construct_name="WEE1_299_569", model_version="commit-123",
        chain_ids=["A"], input_path=input_path,
        runtime={"container": "alphafold3-local-rebuilt.sif"})
    assert manifest["backend"] == "alphafold3"
    assert manifest["ranking_score"] == 0.86
    assert manifest["confidence"] == {"chain_iptm": [0.9]}
    assert len(manifest["structure_sha256"]) == 64
    assert len(manifest["provenance"]["input_sha256"]) == 64
