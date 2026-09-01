from __future__ import annotations

from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import re
import sqlite3
from typing import Any, Sequence

from .campaign import _campaign_for_user
from .registry import stable_id, utc_now
from .similarity import hierarchical_similarity_search, search_morgan_2d


_SAFE_CHAIN = re.compile(r"^[A-Za-z0-9]{1,4}$")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def register_campaign_ligands(
    connection: sqlite3.Connection, user_id: str, campaign_id: str,
    records: Sequence[dict[str, Any]],
) -> list[str]:
    _campaign_for_user(connection, campaign_id, user_id, write=True)
    if not records:
        raise ValueError("At least one ligand record is required")
    inserted: list[str] = []
    for item in records:
        pdb_id = str(item.get("pdb_id", "")).strip().upper()
        ccd_id = str(item.get("ccd_id", "")).strip().upper()
        chain_id = str(item.get("chain_id", "")).strip()
        residue_number = str(item.get("residue_number", "")).strip()
        altloc = str(item.get("altloc", "")).strip()
        smiles = str(item.get("standardized_smiles", "")).strip()
        if not pdb_id or not ccd_id or not residue_number or not smiles:
            raise ValueError("pdb_id, ccd_id, residue_number, and standardized_smiles are required")
        if not _SAFE_CHAIN.fullmatch(chain_id):
            raise ValueError(f"Invalid ligand chain ID: {chain_id}")
        candidate = connection.execute(
            "SELECT * FROM structure_candidate WHERE campaign_id=? AND pdb_id=?",
            (campaign_id, pdb_id),
        ).fetchone()
        if not candidate:
            raise KeyError(f"PDB candidate is not registered in this Campaign: {pdb_id}")
        retained = set(json.loads(candidate["ligand_ids_json"]))
        if ccd_id not in retained:
            raise ValueError(
                f"Component {ccd_id} is not a retained ligand for {pdb_id}; "
                "filtered components cannot be registered")
        descriptor = item.get("usrcat")
        if descriptor is not None and (
                not isinstance(descriptor, list) or len(descriptor) != 60 or
                not all(isinstance(value, (int, float)) for value in descriptor)):
            raise ValueError("USRCAT descriptor must contain 60 numeric values")
        source_path = item.get("source_path")
        source_sha256 = item.get("source_sha256")
        if source_path:
            path = Path(str(source_path)).resolve()
            if not path.is_file():
                raise FileNotFoundError(path)
            actual = _sha256(path)
            if source_sha256 and str(source_sha256).lower() != actual:
                raise ValueError(f"Ligand source checksum mismatch: {path}")
            source_path, source_sha256 = str(path), actual
        ligand_id = stable_id("LIG")
        cursor = connection.execute(
            """INSERT OR IGNORE INTO campaign_ligand(
               id, campaign_id, structure_candidate_id, ccd_id, chain_id,
               residue_number, altloc, standardized_smiles, source_path,
               source_sha256, usrcat_json, metadata_json, created_at)
               VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
            (ligand_id, campaign_id, candidate["id"], ccd_id, chain_id,
             residue_number, altloc, smiles, source_path, source_sha256,
             json.dumps(descriptor) if descriptor is not None else None,
             json.dumps(item.get("metadata", {}), sort_keys=True), utc_now()),
        )
        if cursor.rowcount:
            inserted.append(ligand_id)
    return inserted


def campaign_ligand_status(connection: sqlite3.Connection, user_id: str,
                           campaign_id: str) -> dict[str, Any]:
    _campaign_for_user(connection, campaign_id, user_id, write=False)
    rows = connection.execute(
        """SELECT cl.*, sc.pdb_id FROM campaign_ligand cl
           JOIN structure_candidate sc ON sc.id=cl.structure_candidate_id
           WHERE cl.campaign_id=? ORDER BY sc.pdb_id, cl.ccd_id, cl.chain_id,
           cl.residue_number, cl.altloc""", (campaign_id,),
    ).fetchall()
    lock = connection.execute(
        "SELECT * FROM query_ligand_lock WHERE campaign_id=?", (campaign_id,),
    ).fetchone()
    ligands = []
    for row in rows:
        value = dict(row)
        value["usrcat"] = json.loads(value.pop("usrcat_json")) \
            if value["usrcat_json"] is not None else None
        value["metadata"] = json.loads(value.pop("metadata_json"))
        ligands.append(value)
    return {
        "campaign_id": campaign_id,
        "ligands": ligands,
        "query_lock": ({
            "ligand_id": lock["ligand_id"],
            "snapshot": json.loads(lock["snapshot_json"]),
            "rationale": lock["rationale"], "selected_by": lock["selected_by"],
            "created_at": lock["created_at"],
        } if lock else None),
    }


def compare_campaign_ligands(connection: sqlite3.Connection, user_id: str,
                             campaign_id: str) -> dict[str, Any]:
    status = campaign_ligand_status(connection, user_id, campaign_id)
    ligands = status["ligands"]
    pairs = []
    for left_index, left in enumerate(ligands):
        for right in ligands[left_index + 1:]:
            score_2d = search_morgan_2d(
                left["standardized_smiles"],
                [(right["id"], right["standardized_smiles"])], limit=1)[0].score
            score_3d = None
            if left["usrcat"] is not None and right["usrcat"] is not None:
                distance = sum((a - b) ** 2 for a, b in zip(left["usrcat"], right["usrcat"])) ** 0.5
                score_3d = 1.0 / (1.0 + distance)
            pairs.append({"left_ligand_id": left["id"],
                          "right_ligand_id": right["id"],
                          "morgan_tanimoto": score_2d,
                          "usrcat_similarity": score_3d})
    return {"campaign_id": campaign_id,
            "methods": {"morgan2d": {"radius": 2, "n_bits": 2048},
                        "usrcat3d": {"dimensions": 60, "missing": "null"}},
            "ligands": ligands, "pairs": pairs}


def select_query_ligand(connection: sqlite3.Connection, user_id: str,
                        campaign_id: str, ligand_id: str, rationale: str) -> None:
    _campaign_for_user(connection, campaign_id, user_id, write=True)
    if not rationale.strip():
        raise ValueError("Query ligand selection rationale is required")
    ligand = connection.execute(
        "SELECT * FROM campaign_ligand WHERE id=? AND campaign_id=?",
        (ligand_id, campaign_id),
    ).fetchone()
    if not ligand:
        raise KeyError(f"Unknown Campaign ligand: {ligand_id}")
    existing = connection.execute(
        "SELECT ligand_id FROM query_ligand_lock WHERE campaign_id=?", (campaign_id,),
    ).fetchone()
    if existing:
        raise ValueError("Campaign query ligand is already locked")
    snapshot = dict(ligand)
    for key in ("metadata_json", "usrcat_json"):
        snapshot[key.removesuffix("_json")] = (
            json.loads(snapshot.pop(key)) if snapshot[key] is not None else None)
    connection.execute(
        """INSERT INTO query_ligand_lock(
           campaign_id, ligand_id, snapshot_json, rationale, selected_by, created_at)
           VALUES (?, ?, ?, ?, ?, ?)""",
        (campaign_id, ligand_id, json.dumps(snapshot, sort_keys=True),
         rationale.strip(), user_id, utc_now()),
    )


def run_hierarchical_library_search(
    connection: sqlite3.Connection, user_id: str, campaign_id: str,
    library_id: str, morgan_index_path: Path, *,
    usrcat_index_path: Path | None = None,
    conformer_to_molecule: dict[str, str] | None = None,
    two_d_pool: int = 1000, limit: int = 100,
) -> dict[str, Any]:
    _campaign_for_user(connection, campaign_id, user_id, write=True)
    library = connection.execute("SELECT * FROM library WHERE id=?", (library_id,)).fetchone()
    if not library:
        raise KeyError(f"Unknown library: {library_id}")
    lock = connection.execute(
        """SELECT q.snapshot_json, q.ligand_id, cl.standardized_smiles, cl.usrcat_json
           FROM query_ligand_lock q JOIN campaign_ligand cl ON cl.id=q.ligand_id
           WHERE q.campaign_id=?""", (campaign_id,),
    ).fetchone()
    if not lock:
        raise ValueError("A locked Campaign query ligand is required")
    morgan_path = morgan_index_path.resolve()
    if not morgan_path.is_file():
        raise FileNotFoundError(morgan_path)
    shape_path = usrcat_index_path.resolve() if usrcat_index_path else None
    if shape_path is not None and not shape_path.is_file():
        raise FileNotFoundError(shape_path)
    descriptor = json.loads(lock["usrcat_json"]) if lock["usrcat_json"] else None
    hits = hierarchical_similarity_search(
        lock["standardized_smiles"], descriptor, morgan_path, shape_path,
        conformer_to_molecule=conformer_to_molecule, two_d_pool=two_d_pool,
        limit=limit)
    parameters = {"two_d_pool": two_d_pool, "limit": limit,
                  "morgan_sha256": _sha256(morgan_path),
                  "usrcat_sha256": _sha256(shape_path) if shape_path else None,
                  "ranking": "weighted_0.4_morgan_0.6_usrcat; 2d score when 3d missing"}
    result_rows = [asdict(hit) for hit in hits]
    search_id = stable_id("LSS")
    connection.execute(
        """INSERT INTO ligand_similarity_search(
           id, campaign_id, library_id, query_ligand_id, morgan_index_path,
           usrcat_index_path, parameters_json, query_snapshot_json,
           results_json, created_by, created_at)
           VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
        (search_id, campaign_id, library_id, lock["ligand_id"], str(morgan_path),
         str(shape_path) if shape_path else None, json.dumps(parameters, sort_keys=True),
         lock["snapshot_json"], json.dumps(result_rows, sort_keys=True), user_id, utc_now()),
    )
    return {"search_id": search_id, "campaign_id": campaign_id,
            "library_id": library_id, "parameters": parameters, "hits": result_rows}
