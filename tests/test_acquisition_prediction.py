from __future__ import annotations

import json

import pytest

from aidd_agent.acquisition import acquire_rcsb_mmcif
from aidd_agent.prediction import load_model_profile, prediction_command, write_alphafold3_input


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
