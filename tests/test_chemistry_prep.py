from pathlib import Path

import pytest

from aidd_agent import chemistry_prep


def test_enumerate_retained_ligand_instances_preserves_identity(monkeypatch, tmp_path: Path):
    data = {
        "_atom_site.group_PDB": ["HETATM", "HETATM", "HETATM", "ATOM"],
        "_atom_site.auth_comp_id": ["LIG", "LIG", "GOL", "LIG"],
        "_atom_site.auth_asym_id": ["A"] * 4,
        "_atom_site.auth_seq_id": ["501", "501", "502", "20"],
        "_atom_site.pdbx_PDB_ins_code": ["?", "?", "?", "?"],
        "_atom_site.label_alt_id": ["A", "A", ".", "."],
        "_atom_site.pdbx_PDB_model_num": ["1"] * 4,
        "_atom_site.auth_atom_id": ["C1", "N1", "C1", "CA"],
        "_atom_site.type_symbol": ["C", "N", "C", "C"],
        "_atom_site.Cartn_x": ["1", "2", "3", "4"],
        "_atom_site.Cartn_y": ["0"] * 4,
        "_atom_site.Cartn_z": ["0"] * 4,
    }
    monkeypatch.setattr(chemistry_prep, "_mmcif_dict", lambda path: data)
    result = chemistry_prep.enumerate_ligand_instances(tmp_path / "x.cif", ["LIG"])
    assert len(result) == 1
    assert result[0]["model"] == "1"
    assert result[0]["chain_id"] == "A"
    assert result[0]["residue_number"] == "501"
    assert result[0]["insertion_code"] == ""
    assert result[0]["altloc"] == "A"
    assert [atom["atom_id"] for atom in result[0]["atoms"]] == ["C1", "N1"]


def test_fetch_ccd_validates_and_caches(tmp_path: Path):
    class Response:
        def __enter__(self): return self
        def __exit__(self, *args): pass
        def read(self): return b"data_LIG\n_chem_comp.id LIG\n"
    calls = []
    def opener(url, timeout):
        calls.append((url, timeout)); return Response()
    first = chemistry_prep.fetch_ccd("lig", tmp_path, opener)
    second = chemistry_prep.fetch_ccd("LIG", tmp_path, opener)
    assert first == second
    assert first.read_bytes().startswith(b"data_LIG")
    assert len(calls) == 1


def test_fetch_ccd_rejects_wrong_component(tmp_path: Path):
    class Response:
        def __enter__(self): return self
        def __exit__(self, *args): pass
        def read(self): return b"data_OTHER\n"
    with pytest.raises(ValueError, match="not CCD LIG"):
        chemistry_prep.fetch_ccd("LIG", tmp_path, lambda *args, **kwargs: Response())


def test_dependency_failures_are_explicit(monkeypatch):
    monkeypatch.setattr(chemistry_prep, "_mmcif_dict", lambda path: {})
    assert chemistry_prep.enumerate_ligand_instances(Path("missing"), ["LIG"]) == []


def test_ccd_explicit_kekule_order_wins_over_aromatic_flag():
    class BondType:
        SINGLE = "single"
        DOUBLE = "double"
        TRIPLE = "triple"
        AROMATIC = "aromatic"

    class FakeChem:
        pass

    FakeChem.BondType = BondType
    assert chemistry_prep._ccd_bond_type(FakeChem, "SING", "Y") == "single"
    assert chemistry_prep._ccd_bond_type(FakeChem, "DOUB", "Y") == "double"
    assert chemistry_prep._ccd_bond_type(FakeChem, "AROM", "Y") == "aromatic"
    assert chemistry_prep._ccd_bond_type(FakeChem, "unknown", "Y") == "aromatic"
