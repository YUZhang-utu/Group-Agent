import json
import hashlib
from pathlib import Path

import numpy as np

import aidd_agent.chemical_companion as companion_module

from aidd_agent.chemical_companion import (
    ATOM_DTYPE, BOND_DTYPE, CHEM_META_DTYPE, FEATURE_DIRECTION_DTYPE,
    TORSION_DTYPE, ChemicalCompanionReader, _artifact_identity_lookup,
    build_chemical_companion_shard,
)
from aidd_agent.chemical_geometry import DIRECTION_AXIAL, DIRECTION_SIGNED
from aidd_agent.conformer_artifacts import META_DTYPE


def _write_v1_shard(path: Path, source_path: Path, source_sha256: str) -> None:
    path.mkdir()
    np.zeros(1, dtype=META_DTYPE).tofile(path / "meta.bin")
    np.asarray([b"C1"], dtype="S16").tofile(path / "conformer_ids.bin")
    np.asarray([b"M1"], dtype="S16").tofile(path / "molecule_ids.bin")
    (path / "manifest.json").write_text(json.dumps({
        "format": "aidd-conformer-artifact-shard",
        "source_path": str(source_path),
        "source_sha256": source_sha256,
        "conformers": 1,
    }), encoding="utf-8")


def test_missing_relocated_source_does_not_create_partial(tmp_path: Path):
    shard = tmp_path / "split_0001"
    missing = tmp_path / "missing.mol2"
    _write_v1_shard(shard, missing, "0" * 64)

    try:
        build_chemical_companion_shard(
            tmp_path / "registry.sqlite3", "LIB", shard, tmp_path / "output",
            source_path=missing)
    except FileNotFoundError as error:
        assert str(missing) in str(error)
    else:
        raise AssertionError("missing source was accepted")
    assert not (tmp_path / "output" / ".split_0001.partial").exists()


def test_source_hash_mismatch_does_not_create_partial(tmp_path: Path):
    shard = tmp_path / "split_0001"
    source = tmp_path / "split_0001.mol2"
    source.write_text("not the registered source", encoding="utf-8")
    _write_v1_shard(shard, source, "0" * 64)

    try:
        build_chemical_companion_shard(
            tmp_path / "registry.sqlite3", "LIB", shard, tmp_path / "output")
    except ValueError as error:
        assert "source SHA-256 mismatch" in str(error)
    else:
        raise AssertionError("mismatched source was accepted")
    assert not (tmp_path / "output" / ".split_0001.partial").exists()


def test_artifact_identity_lookup_uses_stable_row_order():
    meta = np.zeros(2, dtype=META_DTYPE)
    meta["global_id"] = [17, 18]
    conformers = np.asarray([b"C17", b"C18"], dtype="S16")
    molecules = np.asarray([b"M8", b"M8"], dtype="S16")

    assert _artifact_identity_lookup(meta, conformers, molecules) == {
        0: (17, "C17", "M8", None),
        1: (18, "C18", "M8", None),
    }


def test_registry_free_shard_build_checks_artifact_row_count(
        tmp_path: Path, monkeypatch):
    source = tmp_path / "split_0001.mol2"
    source.write_text("locked source bytes", encoding="utf-8")
    source_sha256 = hashlib.sha256(source.read_bytes()).hexdigest()
    shard = tmp_path / "split_0001"
    _write_v1_shard(shard, source, source_sha256)
    meta = np.memmap(shard / "meta.bin", dtype=META_DTYPE, mode="r+")
    meta[0]["global_id"] = 7
    meta.flush()

    record = {
        "global_id": 7, "conformer_id": "C1", "molecule_id": "M1",
        "atoms": np.empty(0, dtype=ATOM_DTYPE),
        "bonds": np.empty(0, dtype=BOND_DTYPE),
        "features": np.empty(0, dtype=FEATURE_DIRECTION_DTYPE),
        "feature_members": np.empty(0, dtype="<u2"),
        "torsions": np.empty(0, dtype=TORSION_DTYPE),
        "torsion_members": np.empty(0, dtype="<u2"),
    }

    class FakePool:
        def __enter__(self):
            return self

        def __exit__(self, *_args):
            return False

        def imap(self, _function, _tasks, chunksize):
            assert chunksize == 8
            return iter(((record,),))

    class FakeContext:
        def Pool(self, *_args, **_kwargs):
            return FakePool()

    monkeypatch.setattr(companion_module.mp, "get_context", lambda _name: FakeContext())
    monkeypatch.setattr(companion_module, "_grouped_tasks", lambda *_args: iter((object(),)))

    manifest = build_chemical_companion_shard(
        None, "LIB", shard, tmp_path / "output", workers=1)

    assert manifest["conformers"] == 1
    assert manifest["identity_source"] == "artifact-v1"
    assert (tmp_path / "output" / "split_0001" / "manifest.json").is_file()
    assert not (tmp_path / "output" / ".split_0001.partial").exists()


def test_companion_reader_resolves_identity_features_and_torsions(tmp_path: Path):
    shard = tmp_path / "shard-0"; shard.mkdir()
    atoms = np.asarray([(6, 0, 1, 1), (7, 1, 0, 0), (8, -1, 0, 0)],
                       dtype=ATOM_DTYPE)
    bonds = np.asarray([(0, 1, 1, 0, 0), (1, 2, 2, 0, 0)], dtype=BOND_DTYPE)
    features = np.asarray([
        ((1, 0, 0), DIRECTION_SIGNED, 0, 1),
        ((0, 0, 1), DIRECTION_AXIAL, 1, 2),
    ], dtype=FEATURE_DIRECTION_DTYPE)
    members = np.asarray([1, 0, 2], dtype="<u2")
    torsions = np.asarray([(0, 1, 2, 0, 0, 1)], dtype=TORSION_DTYPE)
    moving = np.asarray([2], dtype="<u2")
    meta = np.asarray([(4, 0, 0, 0, 0, 0, 0, 3, 2, 2, 3, 1, 1)],
                      dtype=CHEM_META_DTYPE)
    arrays = {
        "chem-meta.bin": meta, "atoms.bin": atoms, "bonds.bin": bonds,
        "feature-directions.bin": features, "feature-members.bin": members,
        "torsions.bin": torsions, "torsion-members.bin": moving,
        "conformer_ids.bin": np.asarray([b"C4"], dtype="S16"),
        "molecule_ids.bin": np.asarray([b"M2"], dtype="S16"),
    }
    for name, values in arrays.items():
        values.tofile(shard / name)
    (shard / "manifest.json").write_text(json.dumps({
        "format": "aidd-chemical-companion-shard", "version": 1,
        "global_id_start": 4, "conformers": 1,
    }), encoding="utf-8")
    catalog = tmp_path / "catalog.json"
    catalog.write_text(json.dumps({
        "format": "aidd-chemical-companion-catalog", "version": 1,
        "shards": [{"name": "shard-0", "path": str(shard),
                    "global_id_start": 4, "conformers": 1}],
    }), encoding="utf-8")
    record = ChemicalCompanionReader(catalog).get(4)
    assert (record.molecule_id, record.conformer_id) == ("M2", "C4")
    assert record.atomic_numbers.tolist() == [6, 7, 8]
    assert record.formal_charges.tolist() == [0, 1, -1]
    assert [values.tolist() for values in record.feature_members] == [[1], [0, 2]]
    assert record.feature_kinds.tolist() == [DIRECTION_SIGNED, DIRECTION_AXIAL]
    assert len(record.torsions) == 1
    assert record.torsions[0].moving_atoms == (2,)
