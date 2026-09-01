from __future__ import annotations

import json
from pathlib import Path
import re
import sqlite3
from typing import Any, Mapping, Sequence

import numpy as np

from .campaign import _campaign_for_user
from .project_context import ensure_within
from .registry import stable_id, utc_now

SAFE_TOKEN = re.compile(r"^[A-Za-z0-9_.-]+$")


def create_comparison_set(connection: sqlite3.Connection, user_id: str,
                          campaign_id: str, name: str, candidate_ids: Sequence[str],
                          reference_candidate_id: str, pocket_residues: Sequence[int],
                          chains: Mapping[str, str], rationale: str) -> str:
    campaign = _campaign_for_user(connection, campaign_id, user_id)
    unique = list(dict.fromkeys(candidate_ids))
    if len(unique) < 2:
        raise ValueError("A comparison set requires at least two structures")
    if reference_candidate_id not in unique:
        raise ValueError("Reference structure must be a comparison member")
    if not rationale.strip() or not pocket_residues:
        raise ValueError("Rationale and pocket residues are required")
    placeholders = ",".join("?" for _ in unique)
    rows = connection.execute(
        f"SELECT id, pdb_id, ligand_ids_json FROM structure_candidate WHERE campaign_id=? AND id IN ({placeholders})",
        (campaign_id, *unique),
    ).fetchall()
    if len(rows) != len(unique):
        raise ValueError("All comparison members must belong to the Campaign")
    set_id = stable_id("CMP")
    connection.execute(
        """INSERT INTO structure_comparison_set VALUES (?, ?, ?, ?, ?, ?, ?)""",
        (set_id, campaign_id, name.strip(), reference_candidate_id,
         json.dumps(sorted(set(int(x) for x in pocket_residues))), rationale.strip(), utc_now()),
    )
    by_id = {row["id"]: row for row in rows}
    for candidate_id in unique:
        chain = chains.get(candidate_id, "A")
        if not SAFE_TOKEN.fullmatch(chain):
            raise ValueError(f"Invalid chain ID: {chain}")
        connection.execute(
            "INSERT INTO structure_comparison_member VALUES (?, ?, ?, ?)",
            (set_id, candidate_id, chain, by_id[candidate_id]["ligand_ids_json"]),
        )
    return set_id


def create_campaign_comparison_set(
    connection: sqlite3.Connection, user_id: str, campaign_id: str, name: str,
    pdb_ids: Sequence[str], reference_pdb_id: str, pocket_residues: Sequence[int],
    chains: Mapping[str, str], rationale: str,
) -> str:
    """Create a comparison set using public PDB IDs instead of registry IDs."""
    normalized = list(dict.fromkeys(str(pdb_id).upper() for pdb_id in pdb_ids))
    reference = reference_pdb_id.upper()
    if len(normalized) < 2:
        raise ValueError("A comparison set requires at least two PDB structures")
    placeholders = ",".join("?" for _ in normalized)
    rows = connection.execute(
        f"SELECT id, pdb_id FROM structure_candidate WHERE campaign_id=? "
        f"AND pdb_id IN ({placeholders})",
        (campaign_id, *normalized),
    ).fetchall()
    by_pdb = {row["pdb_id"]: row["id"] for row in rows}
    missing = [pdb_id for pdb_id in normalized if pdb_id not in by_pdb]
    if missing:
        raise ValueError(
            "PDB candidates do not belong to the Campaign: " + ", ".join(missing)
        )
    if reference not in by_pdb:
        raise ValueError("Reference PDB must be a comparison member")
    normalized_chains = {by_pdb[pdb_id]: chains.get(pdb_id, "A") for pdb_id in normalized}
    return create_comparison_set(
        connection, user_id, campaign_id, name,
        [by_pdb[pdb_id] for pdb_id in normalized], by_pdb[reference],
        pocket_residues, normalized_chains, rationale,
    )


def comparison_set_status(connection: sqlite3.Connection, user_id: str,
                          comparison_set_id: str) -> dict[str, Any]:
    row = connection.execute(
        "SELECT * FROM structure_comparison_set WHERE id=?", (comparison_set_id,),
    ).fetchone()
    if not row:
        raise KeyError(f"Unknown comparison set: {comparison_set_id}")
    _campaign_for_user(connection, row["campaign_id"], user_id, write=False)
    members = connection.execute(
        """SELECT sc.id AS candidate_id, sc.pdb_id, scm.chain_id,
                  scm.ligand_ids_json,
                  CASE WHEN sc.id = cs.reference_candidate_id THEN 1 ELSE 0 END AS is_reference
           FROM structure_comparison_member scm
           JOIN structure_candidate sc ON sc.id = scm.candidate_id
           JOIN structure_comparison_set cs ON cs.id = scm.comparison_set_id
           WHERE scm.comparison_set_id=? ORDER BY is_reference DESC, sc.pdb_id""",
        (comparison_set_id,),
    ).fetchall()
    return {
        "comparison_set_id": row["id"], "campaign_id": row["campaign_id"],
        "name": row["name"], "rationale": row["rationale"],
        "pocket_residues": json.loads(row["pocket_residues_json"]),
        "members": [
            {"candidate_id": member["candidate_id"], "pdb_id": member["pdb_id"],
             "chain_id": member["chain_id"], "is_reference": bool(member["is_reference"]),
             "ligand_ids": json.loads(member["ligand_ids_json"])}
            for member in members
        ],
    }


def generate_comparison_pymol_review(
    connection: sqlite3.Connection, user_id: str, comparison_set_id: str,
    project_root: Path, output_path: Path | None = None,
) -> dict[str, Any]:
    status = comparison_set_status(connection, user_id, comparison_set_id)
    root = project_root.resolve()
    structures = []
    missing = []
    reference_id = ""
    for member in status["members"]:
        pdb_id = member["pdb_id"]
        path = root / "inputs" / "structures" / f"{pdb_id}.cif"
        if not path.is_file():
            missing.append(str(path))
        structure_id = f"PDB_{pdb_id}"
        if member["is_reference"]:
            reference_id = structure_id
        structures.append({"structure_id": structure_id, "path": str(path),
                           "chain_id": member["chain_id"]})
    if missing:
        raise FileNotFoundError("Missing downloaded structure files: " + ", ".join(missing))
    destination = output_path or (
        root / "target" / "reviews" / comparison_set_id / "review.pml"
    )
    script = generate_pymol_review(
        structures, reference_id, status["pocket_residues"], destination, root,
    )
    return status | {"pymol_script": str(script)}


def kabsch_align(mobile: np.ndarray, reference: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    mobile = np.asarray(mobile, dtype=float)
    reference = np.asarray(reference, dtype=float)
    if mobile.shape != reference.shape or mobile.ndim != 2 or mobile.shape[1] != 3:
        raise ValueError("Coordinate arrays must have matching shape (N, 3)")
    if len(mobile) < 3:
        raise ValueError("At least three matched coordinates are required")
    mobile_center, reference_center = mobile.mean(0), reference.mean(0)
    covariance = (mobile - mobile_center).T @ (reference - reference_center)
    u, _, vt = np.linalg.svd(covariance)
    correction = np.eye(3)
    correction[-1, -1] = np.sign(np.linalg.det(u @ vt))
    rotation = u @ correction @ vt
    translation = reference_center - mobile_center @ rotation
    aligned = mobile @ rotation + translation
    rmsd = float(np.sqrt(np.mean(np.sum((aligned - reference) ** 2, axis=1))))
    return rotation, translation, rmsd


def compare_pocket_coordinates(
    structures: Mapping[str, Mapping[int, Sequence[float]]], reference_id: str,
    pocket_residues: Sequence[int],
) -> dict[str, Any]:
    if reference_id not in structures:
        raise KeyError(reference_id)
    reference = structures[reference_id]
    results = {}
    for structure_id, coordinates in structures.items():
        matched = [r for r in pocket_residues if r in reference and r in coordinates]
        missing = [r for r in pocket_residues if r not in coordinates]
        if len(matched) < 3:
            results[structure_id] = {"status": "insufficient_overlap", "matched_residues": matched,
                                     "missing_residues": missing}
            continue
        mobile = np.array([coordinates[r] for r in matched], dtype=float)
        target = np.array([reference[r] for r in matched], dtype=float)
        rotation, translation, rmsd = kabsch_align(mobile, target)
        aligned = mobile @ rotation + translation
        deviations = np.linalg.norm(aligned - target, axis=1)
        results[structure_id] = {
            "status": "ok", "matched_residues": matched, "missing_residues": missing,
            "coverage": len(matched) / len(pocket_residues), "pocket_rmsd_angstrom": rmsd,
            "per_residue_deviation_angstrom": {
                str(residue): float(value) for residue, value in zip(matched, deviations)
            },
            "rotation": rotation.tolist(), "translation": translation.tolist(),
        }
    return {"reference_id": reference_id, "pocket_residues": list(pocket_residues),
            "structures": results}


def generate_pymol_review(structures: Sequence[dict[str, Any]], reference_id: str,
                          pocket_residues: Sequence[int], output_path: Path,
                          project_root: Path) -> Path:
    output = ensure_within(output_path, project_root)
    if not SAFE_TOKEN.fullmatch(reference_id):
        raise ValueError("Invalid reference ID")
    names = {item["structure_id"] for item in structures}
    if reference_id not in names:
        raise ValueError("Reference ID is not present")
    residues = "+".join(str(int(x)) for x in sorted(set(pocket_residues)))
    lines = ["reinitialize", "set retain_order, 1", "bg_color white"]
    for item in structures:
        name = item["structure_id"]
        chain = item.get("chain_id", "A")
        if not SAFE_TOKEN.fullmatch(name) or not SAFE_TOKEN.fullmatch(chain):
            raise ValueError("Unsafe structure or chain identifier")
        path = ensure_within(Path(item["path"]), project_root)
        lines += [f'load "{path.as_posix()}", {name}', f"hide everything, {name}",
                  f"show cartoon, {name} and chain {chain}",
                  f"select {name}_pocket, {name} and chain {chain} and resi {residues}",
                  f"show sticks, {name}_pocket",
                  f"select {name}_ligands, {name} and organic", f"show sticks, {name}_ligands"]
        if name != reference_id:
            lines.append(f"align {name} and chain {chain}, {reference_id}")
    lines += ["group comparison_structures, " + " ".join(sorted(names)),
              "orient comparison_structures", "zoom comparison_structures", ""]
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(lines), encoding="utf-8", newline="\n")
    return output


def create_receptor_ensemble(
    connection: sqlite3.Connection, user_id: str, campaign_id: str, name: str,
    reference_pdb_id: str, pocket_residues: Sequence[int],
    chains: Mapping[str, str], rationale: str,
    exclusions: Mapping[str, str] | None = None,
) -> str:
    campaign = _campaign_for_user(connection, campaign_id, user_id)
    if campaign["state"] != "structures_review":
        raise ValueError("Receptor ensembles require structures_review state")
    if not rationale.strip():
        raise ValueError("A scientific ensemble rationale is required")
    experimental = connection.execute(
        "SELECT id, pdb_id FROM structure_candidate WHERE campaign_id=? ORDER BY pdb_id",
        (campaign_id,),
    ).fetchall()
    predicted = connection.execute(
        "SELECT id, construct_name FROM predicted_structure_candidate WHERE campaign_id=? ORDER BY created_at",
        (campaign_id,),
    ).fetchall()
    reference = next(
        (row for row in experimental if row["pdb_id"] == reference_pdb_id.upper()), None)
    if not reference:
        raise ValueError(f"Reference PDB is not a Campaign candidate: {reference_pdb_id.upper()}")
    exclusions = {key.upper(): value.strip() for key, value in (exclusions or {}).items()}
    known_pdbs = {row["pdb_id"] for row in experimental}
    unknown_exclusions = sorted(set(exclusions) - known_pdbs)
    if unknown_exclusions:
        raise ValueError("Excluded PDBs are not Campaign candidates: " + ", ".join(unknown_exclusions))
    if any(not reason for reason in exclusions.values()):
        raise ValueError("Every ensemble exclusion requires a rationale")
    members = [("experimental", row["id"], row["pdb_id"])
               for row in experimental if row["pdb_id"] not in exclusions]
    members += [("predicted", row["id"], row["id"], row["construct_name"])
                for row in predicted]
    if len(members) < 2:
        raise ValueError("A receptor ensemble requires at least two candidates")
    ensemble_id = stable_id("ENS")
    connection.execute(
        "INSERT INTO receptor_ensemble VALUES (?, ?, ?, ?, ?, ?, ?)",
        (ensemble_id, campaign_id, name.strip(), reference["id"],
         json.dumps(sorted(set(int(value) for value in pocket_residues))),
         rationale.strip(), utc_now()),
    )
    for row in experimental:
        if row["pdb_id"] in exclusions:
            connection.execute(
                "INSERT INTO receptor_ensemble_exclusion VALUES (?, ?, 'experimental', ?)",
                (ensemble_id, row["id"], exclusions[row["pdb_id"]]),
            )
    normalized_members = [(*item, "") if len(item) == 3 else item for item in members]
    for kind, candidate_id, public_id, construct_name in normalized_members:
        chain = chains.get(public_id, chains.get(candidate_id, "A"))
        if not SAFE_TOKEN.fullmatch(chain):
            raise ValueError(f"Invalid chain ID: {chain}")
        match = re.search(r"_(\d+)_(\d+)$", construct_name) if kind == "predicted" else None
        residue_offset = int(match.group(1)) - 1 if match else 0
        connection.execute(
            "INSERT INTO receptor_ensemble_member VALUES (?, ?, ?, ?, ?)",
            (ensemble_id, kind, candidate_id, chain, residue_offset),
        )
    return ensemble_id


def receptor_ensemble_status(
    connection: sqlite3.Connection, user_id: str, ensemble_id: str,
) -> dict[str, Any]:
    ensemble = connection.execute(
        "SELECT * FROM receptor_ensemble WHERE id=?", (ensemble_id,)
    ).fetchone()
    if not ensemble:
        raise KeyError(f"Unknown receptor ensemble: {ensemble_id}")
    _campaign_for_user(connection, ensemble["campaign_id"], user_id, write=False)
    members = connection.execute(
        "SELECT * FROM receptor_ensemble_member WHERE ensemble_id=? ORDER BY candidate_kind, candidate_id",
        (ensemble_id,),
    ).fetchall()
    excluded_rows = connection.execute(
        """SELECT ree.candidate_id, ree.candidate_kind, ree.reason, sc.pdb_id
           FROM receptor_ensemble_exclusion ree
           LEFT JOIN structure_candidate sc ON sc.id=ree.candidate_id
           WHERE ree.ensemble_id=? ORDER BY sc.pdb_id""",
        (ensemble_id,),
    ).fetchall()
    resolved = []
    for member in members:
        if member["candidate_kind"] == "experimental":
            row = connection.execute(
                "SELECT pdb_id, ligand_ids_json FROM structure_candidate WHERE id=?",
                (member["candidate_id"],),
            ).fetchone()
            resolved.append({"candidate_kind": "experimental",
                             "candidate_id": member["candidate_id"],
                             "display_id": row["pdb_id"], "pdb_id": row["pdb_id"],
                             "chain_id": member["chain_id"],
                             "residue_offset": member["residue_offset"],
                             "ligand_ids": json.loads(row["ligand_ids_json"]),
                             "is_reference": member["candidate_id"] == ensemble["reference_candidate_id"]})
        else:
            row = connection.execute(
                "SELECT backend, construct_name, structure_path FROM predicted_structure_candidate WHERE id=?",
                (member["candidate_id"],),
            ).fetchone()
            resolved.append({"candidate_kind": "predicted",
                             "candidate_id": member["candidate_id"],
                             "display_id": member["candidate_id"],
                             "backend": row["backend"], "construct_name": row["construct_name"],
                             "structure_path": row["structure_path"],
                             "chain_id": member["chain_id"], "ligand_ids": [],
                             "residue_offset": member["residue_offset"],
                             "is_reference": False})
    return {"ensemble_id": ensemble_id, "campaign_id": ensemble["campaign_id"],
            "name": ensemble["name"], "rationale": ensemble["rationale"],
            "pocket_residues": json.loads(ensemble["pocket_residues_json"]),
            "excluded_candidates": [dict(item) for item in excluded_rows],
            "members": sorted(resolved, key=lambda item: (not item["is_reference"], item["display_id"]))}


def generate_receptor_ensemble_pymol_review(
    connection: sqlite3.Connection, user_id: str, ensemble_id: str,
    project_root: Path, output_path: Path,
) -> dict[str, Any]:
    status = receptor_ensemble_status(connection, user_id, ensemble_id)
    root = project_root.resolve()
    output = ensure_within(output_path, root)
    missing: list[str] = []
    objects: list[dict[str, Any]] = []
    reference_object = ""
    for index, member in enumerate(status["members"], start=1):
        if member["candidate_kind"] == "experimental":
            path = root / "inputs" / "structures" / f"{member['pdb_id']}.cif"
            object_name = f"PDB_{member['pdb_id']}"
        else:
            path = Path(member["structure_path"])
            object_name = f"PRED_{index}_{member['candidate_id'].replace('-', '_')}"
        if not path.is_file():
            missing.append(str(path))
        if member["is_reference"]:
            reference_object = object_name
        objects.append(member | {"path": path, "object_name": object_name})
    if missing:
        raise FileNotFoundError("Missing ensemble structure files: " + ", ".join(missing))
    lines = ["reinitialize", "bg_color white", "set retain_order, 1"]
    for item in objects:
        name, chain = item["object_name"], item["chain_id"]
        lines.extend([f'load "{item["path"].as_posix()}", {name}',
                      f"hide everything, {name}",
                      f"show cartoon, {name} and chain {chain}"])
        if name != reference_object:
            lines.append(f"cealign {reference_object} and chain {objects[0]['chain_id']} and polymer.protein, {name} and chain {chain} and polymer.protein")
        member_residues = [value - item["residue_offset"]
                           for value in status["pocket_residues"]
                           if value - item["residue_offset"] > 0]
        residues = "+".join(str(value) for value in member_residues)
        if residues:
            lines.extend([f"select {name}_canonical_pocket, {name} and chain {chain} and resi {residues}",
                          f"show sticks, {name}_canonical_pocket"])
        ligand_ids = item["ligand_ids"]
        if ligand_ids:
            ligand_selection = "+".join(ligand_ids)
            lines.extend([f"select {name}_ligands, {name} and resn {ligand_selection} and not polymer",
                          f"show sticks, {name}_ligands",
                          f"select {name}_ligand_pocket, byres (({name} and chain {chain} and polymer.protein) within 6 of {name}_ligands)",
                          f"show sticks, {name}_ligand_pocket"])
    lines.extend(["group receptor_ensemble, " + " ".join(item["object_name"] for item in objects),
                  "orient receptor_ensemble", "zoom receptor_ensemble", ""])
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(lines), encoding="utf-8", newline="\n")
    return status | {"pymol_script": str(output)}
