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
