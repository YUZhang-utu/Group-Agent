import json
from pathlib import Path
from types import SimpleNamespace

from aidd_agent.reduce_pocket import (
    _reduce_atom_name, load_reduce_profile, run_reduce, validate_reduce_run,
)


def test_reduce_profile_and_execution_manifest(tmp_path: Path):
    profile = tmp_path / "reduce.json"
    profile.write_text(json.dumps({"executable": "/opt/reduce", "het_dictionary": None,
                                   "build_arguments": ["-build"]}))
    pocket = tmp_path / "pocket.pdb"; pocket.write_text("ATOM\n")
    calls = []

    def runner(command, **kwargs):
        calls.append((command, kwargs))
        return SimpleNamespace(returncode=0, stdout="ATOM H\n", stderr="USER  MOD Single : A 10 HIS\n")

    manifest = run_reduce(profile, pocket, tmp_path / "out", runner=runner)
    assert calls[0][0] == ["/opt/reduce", "-build", str(pocket.resolve())]
    assert manifest["scope"].endswith("within_5A")
    assert Path(manifest["output"]["path"]).read_text() == "ATOM H\n"


def test_reduce_profile_rejects_non_flag_arguments(tmp_path: Path):
    profile = tmp_path / "bad.json"
    profile.write_text(json.dumps({"executable": "reduce", "build_arguments": ["input.pdb"]}))
    try:
        load_reduce_profile(profile)
    except ValueError as exc:
        assert "flags" in str(exc)
    else:
        raise AssertionError("unsafe argument was accepted")


def test_reduce_validation_requires_ligand_and_protein_hydrogens(tmp_path: Path):
    pocket = tmp_path / "pocket.pdb"
    pocket.write_text(
        "HETATM    1  O1  QT9 A 601       0.000   0.000   0.000  1.00 10.00           O\n"
        "ATOM      2  N   ALA A 100       2.800   0.000   0.000  1.00 10.00           N\n")
    out = tmp_path / "out"

    def runner(command, **kwargs):
        return SimpleNamespace(returncode=0, stdout=pocket.read_text() +
            "HETATM    3  H1  QT9 A 601       0.900   0.000   0.000  1.00 10.00           H\n"
            "ATOM      4  H   ALA A 100       1.900   0.000   0.000  1.00 10.00           H\n",
            stderr="WARNING: residues 303 and 305 appear unbonded\nUSER  MOD Flip A 101 ASN\n")

    profile = tmp_path / "reduce.json"
    profile.write_text(json.dumps({"executable": "reduce", "het_dictionary": None,
                                   "build_arguments": ["-build"]}))
    run_reduce(profile, pocket, out, runner=runner)
    report = validate_reduce_run(out, "QT9", "A", "601")
    assert report["accepted"] and report["angle_ready"]
    assert report["hydrogen_counts"]["ligand_after"] == 1
    assert report["het_dictionary_mode"] == "reduce_default_lookup"
    assert report["fatal_warnings"] == []
    assert len(report["expected_fragment_warnings"]) == 1


def test_reduce_dictionary_atom_name_format():
    assert _reduce_atom_name("C1") == " C1 "
    assert _reduce_atom_name("H123") == "H123"
