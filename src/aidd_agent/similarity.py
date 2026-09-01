from __future__ import annotations

from dataclasses import dataclass
import gzip
import json
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np


class ChemistryDependencyError(RuntimeError):
    pass


def _rdkit():
    try:
        from rdkit import Chem, DataStructs
        from rdkit.Chem import AllChem, rdMolDescriptors, rdShapeHelpers
    except ImportError as exc:
        raise ChemistryDependencyError(
            "RDKit is required for ligand similarity. Install the optional chemistry environment."
        ) from exc
    return Chem, DataStructs, AllChem, rdMolDescriptors, rdShapeHelpers


@dataclass(frozen=True)
class SimilarityHit:
    molecule_id: str
    score: float


@dataclass(frozen=True)
class HierarchicalSimilarityHit:
    molecule_id: str
    conformer_id: str | None
    rank: int
    morgan_tanimoto: float
    usrcat_similarity: float | None
    combined_score: float
    stages_completed: tuple[str, ...]


def search_morgan_2d(query_smiles: str, records: Iterable[tuple[str, str]],
                     radius: int = 2, n_bits: int = 2048,
                     limit: int = 100) -> list[SimilarityHit]:
    Chem, DataStructs, AllChem, _, _ = _rdkit()
    query = Chem.MolFromSmiles(query_smiles)
    if query is None:
        raise ValueError("Invalid query SMILES")
    query_fp = AllChem.GetMorganGenerator(radius=radius, fpSize=n_bits).GetFingerprint(query)
    hits = []
    for molecule_id, smiles in records:
        molecule = Chem.MolFromSmiles(smiles)
        if molecule is None:
            continue
        fp = AllChem.GetMorganGenerator(radius=radius, fpSize=n_bits).GetFingerprint(molecule)
        hits.append(SimilarityHit(molecule_id, float(DataStructs.TanimotoSimilarity(query_fp, fp))))
    return sorted(hits, key=lambda hit: (-hit.score, hit.molecule_id))[:limit]


def usrcat_descriptor(mol) -> Sequence[float]:
    _, _, _, rdMolDescriptors, _ = _rdkit()
    if mol.GetNumConformers() == 0:
        raise ValueError("A 3D conformer is required")
    return tuple(float(x) for x in rdMolDescriptors.GetUSRCAT(mol))


MORGAN_INDEX_VERSION = 1
USRCAT_INDEX_VERSION = 1


def build_morgan_index(records: Iterable[tuple[str, str]], output: Path,
                       radius: int = 2, n_bits: int = 2048) -> dict:
    Chem, _, AllChem, _, _ = _rdkit()
    generator = AllChem.GetMorganGenerator(radius=radius, fpSize=n_bits)
    rows = []
    for molecule_id, smiles in records:
        molecule = Chem.MolFromSmiles(smiles)
        if molecule is not None:
            rows.append([molecule_id, smiles, generator.GetFingerprint(molecule).ToBitString()])
    payload = {"format": "aidd-morgan", "version": MORGAN_INDEX_VERSION,
               "radius": radius, "n_bits": n_bits, "records": rows}
    output.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(output, "wt", encoding="utf-8", newline="\n") as stream:
        json.dump(payload, stream, separators=(",", ":"))
    return {"path": str(output), "records": len(rows), "radius": radius, "n_bits": n_bits}


def query_morgan_index(query_smiles: str, index_path: Path,
                       limit: int = 100) -> list[SimilarityHit]:
    Chem, DataStructs, AllChem, _, _ = _rdkit()
    with gzip.open(index_path, "rt", encoding="utf-8") as stream:
        payload = json.load(stream)
    if payload.get("format") != "aidd-morgan" or payload.get("version") != MORGAN_INDEX_VERSION:
        raise ValueError("Unsupported Morgan index format or version")
    molecule = Chem.MolFromSmiles(query_smiles)
    if molecule is None:
        raise ValueError("Invalid query SMILES")
    query = AllChem.GetMorganGenerator(radius=payload["radius"], fpSize=payload["n_bits"]).GetFingerprint(molecule)
    hits = [SimilarityHit(row[0], float(DataStructs.TanimotoSimilarity(
        query, DataStructs.CreateFromBitString(row[2])))) for row in payload["records"]]
    return sorted(hits, key=lambda hit: (-hit.score, hit.molecule_id))[:limit]


def save_usrcat_index(records: Iterable[tuple[str, Sequence[float]]], output: Path) -> dict:
    rows = list(records)
    if not rows:
        raise ValueError("USRCAT index cannot be empty")
    descriptors = np.asarray([row[1] for row in rows], dtype=np.float32)
    if descriptors.ndim != 2 or descriptors.shape[1] != 60:
        raise ValueError("Every USRCAT descriptor must contain 60 values")
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output, format="aidd-usrcat", version=USRCAT_INDEX_VERSION,
                        molecule_ids=np.asarray([row[0] for row in rows]), descriptors=descriptors)
    return {"path": str(output), "records": len(rows), "dimensions": 60}


def query_usrcat_index(query_descriptor: Sequence[float], index_path: Path,
                       limit: int = 100) -> list[SimilarityHit]:
    query = np.asarray(query_descriptor, dtype=np.float32)
    if query.shape != (60,):
        raise ValueError("USRCAT query descriptor must contain 60 values")
    with np.load(index_path, allow_pickle=False) as index:
        if str(index["format"]) != "aidd-usrcat" or int(index["version"]) != USRCAT_INDEX_VERSION:
            raise ValueError("Unsupported USRCAT index format or version")
        ids, descriptors = index["molecule_ids"], index["descriptors"]
        distances = np.linalg.norm(descriptors - query, axis=1)
    order = np.argsort(distances, kind="stable")[:limit]
    return [SimilarityHit(str(ids[i]), float(1.0 / (1.0 + distances[i]))) for i in order]


def hierarchical_similarity_search(
    query_smiles: str,
    query_descriptor: Sequence[float] | None,
    morgan_index_path: Path,
    usrcat_index_path: Path | None = None,
    *,
    conformer_to_molecule: dict[str, str] | None = None,
    two_d_pool: int = 1000,
    limit: int = 100,
    two_d_weight: float = 0.4,
    three_d_weight: float = 0.6,
) -> list[HierarchicalSimilarityHit]:
    if two_d_pool < 1 or limit < 1 or limit > two_d_pool:
        raise ValueError("Similarity limits must satisfy 1 <= limit <= two_d_pool")
    if two_d_weight < 0 or three_d_weight < 0 or two_d_weight + three_d_weight <= 0:
        raise ValueError("Similarity weights must be non-negative with a positive sum")
    two_d = query_morgan_index(query_smiles, morgan_index_path, limit=two_d_pool)
    two_d_by_id = {hit.molecule_id: hit.score for hit in two_d}
    best_shape: dict[str, tuple[str, float]] = {}
    if query_descriptor is not None and usrcat_index_path is not None:
        shape_hits = query_usrcat_index(
            query_descriptor, usrcat_index_path, limit=2_147_483_647)
        mapping = conformer_to_molecule or {}
        for hit in shape_hits:
            molecule_id = mapping.get(hit.molecule_id, hit.molecule_id)
            if molecule_id not in two_d_by_id:
                continue
            previous = best_shape.get(molecule_id)
            if previous is None or hit.score > previous[1] or (
                    hit.score == previous[1] and hit.molecule_id < previous[0]):
                best_shape[molecule_id] = (hit.molecule_id, hit.score)
    weight_sum = two_d_weight + three_d_weight
    rows = []
    for molecule_id, score_2d in two_d_by_id.items():
        shape = best_shape.get(molecule_id)
        if shape:
            combined = (two_d_weight * score_2d + three_d_weight * shape[1]) / weight_sum
            rows.append((molecule_id, shape[0], score_2d, shape[1], combined,
                         ("morgan2d", "usrcat3d")))
        else:
            rows.append((molecule_id, None, score_2d, None, score_2d, ("morgan2d",)))
    rows.sort(key=lambda row: (-row[4], -row[2], row[0], row[1] or ""))
    return [HierarchicalSimilarityHit(
        molecule_id=row[0], conformer_id=row[1], rank=index,
        morgan_tanimoto=row[2], usrcat_similarity=row[3], combined_score=row[4],
        stages_completed=row[5]) for index, row in enumerate(rows[:limit], start=1)]
