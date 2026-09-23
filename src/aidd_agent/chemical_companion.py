from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import multiprocessing as mp
import os
from pathlib import Path
import sqlite3
import time
from typing import Iterator

import numpy as np

from .chemical_geometry import (
    DIRECTION_AXIAL, DIRECTION_NONE, DIRECTION_SIGNED, TerminalTorsion, unit_vector,
)
from .conformer_artifacts import FEATURE_TYPES, META_DTYPE
from .mol2 import iter_mol2_records, load_rdkit_mol2


FORMAT = "aidd-chemical-companion-catalog"
SHARD_FORMAT = "aidd-chemical-companion-shard"
VERSION = 1

CHEM_META_DTYPE = np.dtype([
    ("global_id", "<i8"),
    ("atom_offset", "<u8"), ("bond_offset", "<u8"),
    ("feature_offset", "<u8"), ("member_offset", "<u8"),
    ("torsion_offset", "<u8"), ("torsion_member_offset", "<u8"),
    ("atoms", "<u2"), ("bonds", "<u2"), ("features", "<u2"),
    ("feature_members", "<u2"), ("torsions", "<u2"),
    ("torsion_members", "<u2"),
])
ATOM_DTYPE = np.dtype([
    ("atomic_number", "u1"), ("formal_charge", "i1"), ("flags", "u1"),
    ("chiral_tag", "u1"),
])
BOND_DTYPE = np.dtype([
    ("begin", "<u2"), ("end", "<u2"), ("order", "u1"), ("flags", "u1"),
    ("stereo", "u1"),
])
FEATURE_DIRECTION_DTYPE = np.dtype([
    ("direction", "<f2", (3,)), ("kind", "u1"),
    ("member_offset", "<u8"), ("member_count", "<u2"),
])
TORSION_DTYPE = np.dtype([
    ("atom_a", "<u2"), ("atom_b", "<u2"),
    ("atom_c", "<u2"), ("atom_d", "<u2"),
    ("moving_offset", "<u8"), ("moving_count", "<u2"),
])

ATOM_FLAG_AROMATIC = 1
BOND_FLAG_AROMATIC = 1
BOND_FLAG_RING = 2

_WORKER_FACTORY = None
_WORKER_MAX_MOVING = 12


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _decode(value) -> str:
    return bytes(value).split(b"\0", 1)[0].decode("utf-8")


def _worker_init(max_moving_atoms: int) -> None:
    from rdkit import RDConfig
    from rdkit.Chem import ChemicalFeatures
    global _WORKER_FACTORY, _WORKER_MAX_MOVING
    _WORKER_FACTORY = ChemicalFeatures.BuildFeatureFactory(
        str(Path(RDConfig.RDDataDir) / "BaseFeatures.fdef"))
    _WORKER_MAX_MOVING = max_moving_atoms


def _bond_order(bond) -> int:
    from rdkit import Chem
    mapping = {Chem.BondType.SINGLE: 1, Chem.BondType.DOUBLE: 2,
               Chem.BondType.TRIPLE: 3, Chem.BondType.AROMATIC: 4}
    value = mapping.get(bond.GetBondType())
    if value is None:
        raise ValueError(f"unsupported bond type {bond.GetBondType()}")
    return value


def _feature_direction(mol, family: str, atom_ids: tuple[int, ...]) -> tuple[np.ndarray, int]:
    conformer = mol.GetConformer()
    positions = lambda ids: np.asarray([
        list(conformer.GetAtomPosition(int(index))) for index in ids], dtype=np.float64)
    heavy_ids = tuple(index for index in atom_ids
                      if mol.GetAtomWithIdx(index).GetAtomicNum() > 1)
    if family == "Aromatic" and len(heavy_ids) >= 3:
        centered = positions(heavy_ids) - positions(heavy_ids).mean(axis=0)
        _, _, vh = np.linalg.svd(centered, full_matrices=False)
        vector = unit_vector(vh[-1])
        return ((np.zeros(3), DIRECTION_NONE) if vector is None
                else (vector, DIRECTION_AXIAL))
    if family not in {"Donor", "Acceptor"} or len(heavy_ids) != 1:
        return np.zeros(3), DIRECTION_NONE
    index = heavy_ids[0]
    atom = mol.GetAtomWithIdx(index)
    center = positions((index,))[0]
    hydrogens = [neighbor.GetIdx() for neighbor in atom.GetNeighbors()
                 if neighbor.GetAtomicNum() == 1]
    heavy_neighbors = [neighbor.GetIdx() for neighbor in atom.GetNeighbors()
                       if neighbor.GetAtomicNum() > 1]
    if family == "Donor" and hydrogens:
        vector = unit_vector(positions(hydrogens).mean(axis=0) - center)
    elif heavy_neighbors:
        vector = unit_vector(center - positions(heavy_neighbors).mean(axis=0))
    else:
        vector = None
    return ((np.zeros(3), DIRECTION_NONE) if vector is None
            else (vector, DIRECTION_SIGNED))


def _terminal_torsions(mol, heavy_map: dict[int, int]) -> list[tuple]:
    from rdkit.Chem import Lipinski
    adjacency = {local: set() for local in heavy_map.values()}
    for bond in mol.GetBonds():
        left, right = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if left in heavy_map and right in heavy_map:
            a, b = heavy_map[left], heavy_map[right]
            adjacency[a].add(b); adjacency[b].add(a)

    def component(start: int, blocked: tuple[int, int]) -> set[int]:
        found, stack = set(), [start]
        while stack:
            node = stack.pop()
            if node in found:
                continue
            found.add(node)
            stack.extend(neighbor for neighbor in adjacency[node]
                         if (node, neighbor) != blocked and
                         (neighbor, node) != blocked and neighbor not in found)
        return found

    seen, result = set(), []
    for match in mol.GetSubstructMatches(Lipinski.RotatableBondSmarts):
        original_b, original_c = int(match[0]), int(match[1])
        if original_b not in heavy_map or original_c not in heavy_map:
            continue
        bond = mol.GetBondBetweenAtoms(original_b, original_c)
        if bond is None or bond.IsInRing():
            continue
        b, c = heavy_map[original_b], heavy_map[original_c]
        side_c = component(c, (b, c)); side_b = component(b, (b, c))
        if len(side_b) < len(side_c):
            b, c, side_c = c, b, side_b
        if len(side_c) > _WORKER_MAX_MOVING:
            continue
        fixed_neighbors = sorted(adjacency[b] - {c})
        moving_neighbors = sorted(adjacency[c] - {b})
        if not fixed_neighbors or not moving_neighbors:
            continue
        key = tuple(sorted((b, c)))
        if key in seen:
            continue
        seen.add(key)
        result.append((fixed_neighbors[0], b, c, moving_neighbors[0],
                       tuple(sorted(side_c))))
    return sorted(result, key=lambda row: (len(row[4]), row[1], row[2], row[4]))


def _process_group(task):
    output = []
    for global_id, conformer_id, molecule_id, raw_text in task:
        mol, sanitization = load_rdkit_mol2(raw_text, conformer_id)
        from .budget_export import hydrogen_roundtrip
        hydrogen_audit = hydrogen_roundtrip(mol)
        heavy = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1]
        heavy_map = {original: local for local, original in enumerate(heavy)}
        atoms = np.empty(len(heavy), dtype=ATOM_DTYPE)
        for local, original in enumerate(heavy):
            atom = mol.GetAtomWithIdx(original)
            charge = int(atom.GetFormalCharge())
            if charge < -128 or charge > 127:
                raise ValueError(f"formal charge outside int8 at {conformer_id}")
            atoms[local] = (atom.GetAtomicNum(), charge,
                            ATOM_FLAG_AROMATIC if atom.GetIsAromatic() else 0,
                            int(atom.GetChiralTag()))
        bonds = []
        for bond in mol.GetBonds():
            left, right = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            if left not in heavy_map or right not in heavy_map:
                continue
            flags = ((BOND_FLAG_AROMATIC if bond.GetIsAromatic() else 0)
                     | (BOND_FLAG_RING if bond.IsInRing() else 0))
            bonds.append((heavy_map[left], heavy_map[right], _bond_order(bond), flags,
                          int(bond.GetStereo())))
        bond_array = np.asarray(bonds, dtype=BOND_DTYPE)
        feature_rows, members = [], []
        features = [feature for feature in _WORKER_FACTORY.GetFeaturesForMol(mol)
                    if feature.GetFamily() in FEATURE_TYPES]
        for feature in features:
            atom_ids = tuple(int(value) for value in feature.GetAtomIds())
            local_members = tuple(heavy_map[value] for value in atom_ids if value in heavy_map)
            direction, kind = _feature_direction(mol, feature.GetFamily(), atom_ids)
            feature_rows.append((direction.astype(np.float16), kind,
                                 len(members), len(local_members)))
            members.extend(local_members)
        feature_array = np.asarray(feature_rows, dtype=FEATURE_DIRECTION_DTYPE)
        torsion_rows, torsion_members = [], []
        for atom_a, atom_b, atom_c, atom_d, moving in _terminal_torsions(mol, heavy_map):
            torsion_rows.append((atom_a, atom_b, atom_c, atom_d,
                                 len(torsion_members), len(moving)))
            torsion_members.extend(moving)
        output.append({
            "global_id": global_id, "conformer_id": conformer_id,
            "molecule_id": molecule_id, "atoms": atoms, "bonds": bond_array,
            "features": feature_array,
            "feature_members": np.asarray(members, dtype="<u2"),
            "torsions": np.asarray(torsion_rows, dtype=TORSION_DTYPE),
            "torsion_members": np.asarray(torsion_members, dtype="<u2"),
            "sanitization": sanitization,
            "hydrogen_audit": hydrogen_audit,
        })
    return output


def _grouped_tasks(source_path: Path, lookup: dict[int, tuple]) -> Iterator[tuple]:
    current = None
    records = []
    for record in iter_mol2_records(source_path):
        identity = lookup.get(record.record_index)
        if identity is None:
            continue
        global_id, conformer_id, molecule_id, expected_sha256 = identity
        if expected_sha256 is not None and record.content_sha256 != expected_sha256:
            raise ValueError(
                f"source record checksum mismatch at record {record.record_index}")
        if current is not None and molecule_id != current:
            yield tuple(records); records = []
        current = molecule_id
        records.append((global_id, conformer_id, molecule_id, record.raw_text))
    if records:
        yield tuple(records)


def _v1_identity(v1_shard: Path) -> tuple[dict, np.memmap, np.memmap, np.memmap]:
    manifest = json.loads((v1_shard / "manifest.json").read_text(encoding="utf-8"))
    if manifest.get("format") != "aidd-conformer-artifact-shard":
        raise ValueError(f"not an artifact v1 shard: {v1_shard}")
    meta = np.memmap(v1_shard / "meta.bin", dtype=META_DTYPE, mode="r")
    conformers = np.memmap(v1_shard / "conformer_ids.bin", dtype="S16", mode="r")
    molecules = np.memmap(v1_shard / "molecule_ids.bin", dtype="S16", mode="r")
    if not len(meta) == len(conformers) == len(molecules) == int(manifest["conformers"]):
        raise ValueError("artifact v1 identity lengths disagree")
    return manifest, meta, conformers, molecules


def _artifact_identity_lookup(
        meta: np.ndarray, conformers: np.ndarray,
        molecules: np.ndarray) -> dict[int, tuple]:
    return {
        index: (int(row["global_id"]), _decode(conformers[index]),
                _decode(molecules[index]), None)
        for index, row in enumerate(meta)
    }


def build_chemical_companion_shard(
        db_path: Path | None, library_id: str, v1_shard: Path, output_root: Path, *,
        workers: int = 16, max_moving_atoms: int = 12,
        source_path: Path | None = None) -> dict:
    if workers <= 0 or max_moving_atoms <= 0:
        raise ValueError("workers and max_moving_atoms must be positive")
    v1_shard = v1_shard.resolve(); output_root = output_root.resolve()
    v1_manifest, v1_meta, v1_conformers, v1_molecules = _v1_identity(v1_shard)
    shard_name = v1_shard.name
    final = output_root / shard_name
    if final.is_dir():
        manifest = json.loads((final / "manifest.json").read_text(encoding="utf-8"))
        if manifest.get("format") != SHARD_FORMAT:
            raise ValueError(f"existing output is not a chemical companion: {final}")
        return manifest
    registered_source = str(v1_manifest["source_path"])
    source = (Path(source_path).resolve() if source_path is not None
              else Path(registered_source))
    if not source.is_file():
        raise FileNotFoundError(source)
    partial = output_root / f".{shard_name}.partial"
    if partial.exists():
        raise ValueError(f"incomplete companion exists; inspect before retry: {partial}")
    expected_source_sha256 = str(v1_manifest["source_sha256"])
    actual_source_sha256 = _sha256(source)
    if actual_source_sha256 != expected_source_sha256:
        raise ValueError(
            f"source SHA-256 mismatch for {shard_name}: expected "
            f"{expected_source_sha256}, got {actual_source_sha256}")
    lookup = _artifact_identity_lookup(v1_meta, v1_conformers, v1_molecules)
    source_validation = (
        "whole-file artifact SHA-256 plus ordered artifact-v1 identity and shape matched")
    if db_path is not None:
        with sqlite3.connect(db_path) as connection:
            connection.row_factory = sqlite3.Row
            rows = connection.execute(
                """SELECT c.id conformer_id,c.molecule_id,c.source_record_index,
                          c.content_sha256
                   FROM conformer c JOIN molecule m ON m.id=c.molecule_id
                   WHERE m.library_id=? AND c.source_path=?
                   ORDER BY c.source_record_index""",
                (library_id, registered_source)).fetchall()
        if len(rows) != len(v1_meta):
            raise ValueError(
                f"registry/artifact conformer count mismatch for {shard_name}: "
                f"registry={len(rows)}, artifact={len(v1_meta)}, db={db_path}, "
                f"registered_source={registered_source!r}; omit --db to use the "
                "immutable artifact-v1 identity on a relocated workstation")
        lookup = {}
        for index, row in enumerate(rows):
            gid = int(v1_meta[index]["global_id"])
            conformer_id = str(row["conformer_id"])
            molecule_id = str(row["molecule_id"])
            if (_decode(v1_conformers[index]) != conformer_id
                    or _decode(v1_molecules[index]) != molecule_id):
                raise ValueError(f"registry/artifact identity mismatch at global ID {gid}")
            lookup[int(row["source_record_index"])] = (
                gid, conformer_id, molecule_id, str(row["content_sha256"]))
        source_validation = (
            "whole-file artifact SHA-256, ordered artifact-v1 identity and shape, "
            "and all per-record registry SHA-256 values matched")

    partial.mkdir(parents=True)
    names = ("chem-meta.bin", "atoms.bin", "bonds.bin", "feature-directions.bin",
             "feature-members.bin", "torsions.bin", "torsion-members.bin",
             "conformer_ids.bin", "molecule_ids.bin", "hydrogen-audit.jsonl", "total-h.bin")
    streams = {name: (partial / name).open("wb") for name in names}
    offsets = {"atom": 0, "bond": 0, "feature": 0, "member": 0,
               "torsion": 0, "torsion_member": 0}
    written = 0; started = time.perf_counter()
    sanitization_counts = {"strict": 0, "aromatic_no_kekulize": 0}
    context = mp.get_context("spawn")
    try:
        with context.Pool(workers, initializer=_worker_init,
                          initargs=(max_moving_atoms,)) as pool:
            for group in pool.imap(_process_group,
                                   _grouped_tasks(source, lookup), chunksize=8):
                for record in group:
                    expected = v1_meta[written]
                    if (record["global_id"] != int(expected["global_id"])
                            or len(record["atoms"]) != int(expected["heavy_atoms"])
                            or len(record["features"]) != int(expected["features"])):
                        raise ValueError(
                            f"companion/artifact shape mismatch at row {written}")
                    record["atoms"].tofile(streams["atoms.bin"])
                    audit=record['hydrogen_audit']
                    streams['hydrogen-audit.jsonl'].write((json.dumps(dict(global_id=record['global_id'],**audit))+'\n').encode())
                    np.asarray([a['total_h'] for a in audit['atoms']],dtype='<u2').tofile(streams['total-h.bin'])
                    record["bonds"].tofile(streams["bonds.bin"])
                    record["features"].tofile(streams["feature-directions.bin"])
                    record["feature_members"].tofile(streams["feature-members.bin"])
                    record["torsions"].tofile(streams["torsions.bin"])
                    record["torsion_members"].tofile(streams["torsion-members.bin"])
                    sanitization_counts[record["sanitization"]] += 1
                    meta = np.asarray([(
                        record["global_id"], offsets["atom"], offsets["bond"],
                        offsets["feature"], offsets["member"], offsets["torsion"],
                        offsets["torsion_member"], len(record["atoms"]),
                        len(record["bonds"]), len(record["features"]),
                        len(record["feature_members"]), len(record["torsions"]),
                        len(record["torsion_members"])),
                    ], dtype=CHEM_META_DTYPE)
                    meta.tofile(streams["chem-meta.bin"])
                    np.asarray([record["conformer_id"].encode()], dtype="S16").tofile(
                        streams["conformer_ids.bin"])
                    np.asarray([record["molecule_id"].encode()], dtype="S16").tofile(
                        streams["molecule_ids.bin"])
                    for key, payload in (("atom", record["atoms"]),
                                         ("bond", record["bonds"]),
                                         ("feature", record["features"]),
                                         ("member", record["feature_members"]),
                                         ("torsion", record["torsions"]),
                                         ("torsion_member", record["torsion_members"])):
                        offsets[key] += len(payload)
                    written += 1
    finally:
        for stream in streams.values():
            stream.close()
    expected_rows = len(v1_meta)
    if written != expected_rows:
        raise ValueError(f"expected {expected_rows} companion rows, wrote {written}")
    files = {name: {"bytes": (partial / name).stat().st_size,
                    "sha256": _sha256(partial / name)} for name in names}
    manifest = {
        "format": SHARD_FORMAT, "version": VERSION,
        "library_id": library_id, "name": shard_name,
        "artifact_v1_shard": str(v1_shard),
        "artifact_v1_manifest_sha256": _sha256(v1_shard / "manifest.json"),
        "source_path": str(source.resolve()),
        "source_sha256": actual_source_sha256,
        "source_validation": source_validation,
        "identity_source": ("artifact-v1" if db_path is None
                            else "artifact-v1-plus-registry"),
        "global_id_start": int(v1_meta[0]["global_id"]),
        "conformers": written, "counts": offsets,
        "workers": workers, "max_moving_atoms": max_moving_atoms,
        "rdkit_sanitization_counts": sanitization_counts,
        "wall_seconds": time.perf_counter() - started,
        "direction_model": "explicit-H-or-local-heavy-geometry-v1",
        "torsion_model": "rdkit-strict-nonring-smallest-side-v1",
        "files": files,
    }
    (partial / "manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8")
    os.replace(partial, final)
    return manifest


def build_chemical_companion_catalog(
        db_path: Path | None, library_id: str, artifact_catalog: Path,
        output_root: Path, *, workers: int = 16,
        max_moving_atoms: int = 12,
        source_overrides: dict[str, Path] | None = None) -> dict:
    artifact_catalog = artifact_catalog.resolve(); output_root = output_root.resolve()
    source_catalog = json.loads(artifact_catalog.read_text(encoding="utf-8"))
    if source_catalog.get("format") != "aidd-conformer-artifact-catalog":
        raise ValueError("not an artifact v1 catalog")
    if source_catalog.get("library_id") != library_id:
        raise ValueError("artifact catalog library does not match --library")
    output_root.mkdir(parents=True, exist_ok=True)
    source_overrides = source_overrides or {}
    unknown = sorted(set(source_overrides) - {str(row["name"])
                                             for row in source_catalog["shards"]})
    if unknown:
        raise ValueError("source override names are not artifact shards: " + ", ".join(unknown))
    shards = []
    for row in source_catalog["shards"]:
        v1_shard = Path(row["path"])
        if not v1_shard.is_dir():
            v1_shard = artifact_catalog.parent / row["name"]
        manifest = build_chemical_companion_shard(
            db_path, library_id, v1_shard, output_root,
            workers=workers, max_moving_atoms=max_moving_atoms,
            source_path=source_overrides.get(str(row["name"])))
        directory = output_root / v1_shard.name
        shards.append({
            "name": v1_shard.name, "path": str(directory),
            "manifest_sha256": _sha256(directory / "manifest.json"),
            "global_id_start": manifest["global_id_start"],
            "conformers": manifest["conformers"],
        })
    document = {
        "format": FORMAT, "version": VERSION, "library_id": library_id,
        "artifact_v1_catalog": str(artifact_catalog),
        "artifact_v1_catalog_sha256": _sha256(artifact_catalog),
        "shards": shards,
    }
    temporary = output_root / "catalog.json.partial"
    temporary.write_text(json.dumps(document, indent=2), encoding="utf-8")
    os.replace(temporary, output_root / "catalog.json")
    return document


@dataclass(frozen=True)
class ChemicalConformer:
    global_id: int
    conformer_id: str
    molecule_id: str
    atomic_numbers: np.ndarray
    formal_charges: np.ndarray
    atom_flags: np.ndarray
    bonds: np.ndarray
    feature_directions: np.ndarray
    feature_kinds: np.ndarray
    feature_members: tuple[np.ndarray, ...]
    torsions: tuple[TerminalTorsion, ...]
    total_hydrogens: np.ndarray | None = None


class ChemicalCompanionReader:
    def __init__(self, catalog_path: Path):
        self.catalog_path = catalog_path.resolve()
        catalog = json.loads(self.catalog_path.read_text(encoding="utf-8"))
        if catalog.get("format") != FORMAT or catalog.get("version") != VERSION:
            raise ValueError("not a supported chemical companion catalog")
        self._shards = []
        mapped = lambda path, dtype: (np.empty(0, dtype=dtype)
                                     if path.stat().st_size == 0
                                     else np.memmap(path, dtype=dtype, mode="r"))
        for row in catalog["shards"]:
            directory = Path(row["path"])
            if not directory.is_dir():
                directory = self.catalog_path.parent / row["name"]
            manifest = json.loads((directory / "manifest.json").read_text(encoding="utf-8"))
            count, start = int(manifest["conformers"]), int(manifest["global_id_start"])
            arrays = {
                "meta": mapped(directory / "chem-meta.bin", CHEM_META_DTYPE),
                "atoms": mapped(directory / "atoms.bin", ATOM_DTYPE),
                "bonds": mapped(directory / "bonds.bin", BOND_DTYPE),
                "features": mapped(directory / "feature-directions.bin", FEATURE_DIRECTION_DTYPE),
                "members": mapped(directory / "feature-members.bin", "<u2"),
                "torsions": mapped(directory / "torsions.bin", TORSION_DTYPE),
                "torsion_members": mapped(directory / "torsion-members.bin", "<u2"),
                "conformers": mapped(directory / "conformer_ids.bin", "S16"),
                "molecules": mapped(directory / "molecule_ids.bin", "S16"),
            }
            if 'total-h.bin' in manifest.get('files',{}):
                arrays['hydrogens']=mapped(directory/'total-h.bin','<u2')
                if len(arrays['hydrogens'])!=len(arrays['atoms']):raise ValueError('Hydrogen metadata length mismatch')
            if not len(arrays["meta"]) == len(arrays["conformers"]) == len(arrays["molecules"]) == count:
                raise ValueError(f"chemical companion row count mismatch: {directory}")
            self._shards.append((start, start + count, arrays))
        self._shards.sort(key=lambda item: item[0])

    def get(self, global_id: int) -> ChemicalConformer:
        gid = int(global_id)
        for start, stop, arrays in self._shards:
            if not start <= gid < stop:
                continue
            index = gid - start; row = arrays["meta"][index]
            if int(row["global_id"]) != gid:
                raise ValueError(f"chemical companion global ID mismatch: {gid}")
            sliced = lambda name, offset, count: arrays[name][int(offset):int(offset)+int(count)]
            atoms = sliced("atoms", row["atom_offset"], row["atoms"])
            bonds = sliced("bonds", row["bond_offset"], row["bonds"])
            features = sliced("features", row["feature_offset"], row["features"])
            members = sliced("members", row["member_offset"], row["feature_members"])
            torsions = sliced("torsions", row["torsion_offset"], row["torsions"])
            torsion_members = sliced(
                "torsion_members", row["torsion_member_offset"], row["torsion_members"])
            feature_members = tuple(np.asarray(
                members[int(feature["member_offset"]):
                        int(feature["member_offset"])+int(feature["member_count"])],
                dtype=np.int64) for feature in features)
            definitions = tuple(TerminalTorsion(
                int(value["atom_a"]), int(value["atom_b"]),
                int(value["atom_c"]), int(value["atom_d"]),
                tuple(map(int, torsion_members[
                    int(value["moving_offset"]):
                    int(value["moving_offset"])+int(value["moving_count"])])))
                for value in torsions)
            return ChemicalConformer(
                gid, _decode(arrays["conformers"][index]),
                _decode(arrays["molecules"][index]),
                np.asarray(atoms["atomic_number"], dtype=np.uint8),
                np.asarray(atoms["formal_charge"], dtype=np.int8),
                np.asarray(atoms["flags"], dtype=np.uint8), np.asarray(bonds),
                np.asarray(features["direction"], dtype=np.float64),
                np.asarray(features["kind"], dtype=np.uint8),
                feature_members, definitions,
                np.asarray(sliced('hydrogens',row['atom_offset'],row['atoms'])) if 'hydrogens' in arrays else None)
        raise IndexError(f"global ID outside chemical companion: {gid}")
