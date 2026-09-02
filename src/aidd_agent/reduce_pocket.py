from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import subprocess
from typing import Callable

import numpy as np


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_reduce_profile(path: Path) -> dict:
    profile = json.loads(path.read_text(encoding="utf-8"))
    allowed = {"executable", "het_dictionary", "build_arguments"}
    if set(profile) - allowed or not profile.get("executable"):
        raise ValueError("Invalid Reduce profile")
    arguments = profile.get("build_arguments", ["-build"])
    if not isinstance(arguments, list) or not all(isinstance(x, str) and x.startswith("-") for x in arguments):
        raise ValueError("Reduce build_arguments must be a list of flags")
    return {"executable": str(profile["executable"]),
            "het_dictionary": profile.get("het_dictionary"),
            "build_arguments": arguments}


def export_query_pocket_pdb(mmcif: Path, output: Path, ccd_id: str,
                            chain_id: str, residue_number: str,
                            radius_angstrom: float = 5.0) -> dict:
    """Export the query ligand and complete protein residues contacting it."""
    from Bio.PDB import MMCIFParser, PDBIO, Select

    structure = MMCIFParser(QUIET=True).get_structure(mmcif.stem, str(mmcif))
    model = next(structure.get_models())
    ligand_residues = [residue for chain in model for residue in chain
                       if chain.id == chain_id and residue.get_resname().strip() == ccd_id
                       and str(residue.id[1]) == str(residue_number)]
    if len(ligand_residues) != 1:
        raise ValueError(f"Expected one {ccd_id} {chain_id}:{residue_number}, found {len(ligand_residues)}")
    ligand = ligand_residues[0]
    ligand_xyz = np.asarray([atom.coord for atom in ligand if atom.element != "H"])
    pocket = set()
    for chain in model:
        for residue in chain:
            if residue is ligand or residue.id[0].strip():
                continue
            heavy = np.asarray([atom.coord for atom in residue if atom.element != "H"])
            if len(heavy) and np.min(np.linalg.norm(heavy[:, None] - ligand_xyz[None, :], axis=2)) <= radius_angstrom:
                pocket.add((chain.id, residue.id))

    class PocketSelect(Select):
        def accept_residue(self, residue):
            chain = residue.get_parent()
            return int(residue is ligand or (chain.id, residue.id) in pocket)

    output.parent.mkdir(parents=True, exist_ok=True)
    io = PDBIO(); io.set_structure(structure); io.save(str(output), PocketSelect())
    return {"path": str(output.resolve()), "sha256": _sha256(output), "ccd_id": ccd_id,
            "chain_id": chain_id, "residue_number": str(residue_number),
            "radius_angstrom": radius_angstrom, "protein_residues": len(pocket)}


def run_reduce(profile_path: Path, pocket_pdb: Path, output_dir: Path,
               runner: Callable = subprocess.run) -> dict:
    """Run Reduce with an argument array and preserve all provenance/output."""
    profile = load_reduce_profile(profile_path)
    output_dir.mkdir(parents=True, exist_ok=True)
    environment = os.environ.copy()
    if profile["het_dictionary"]:
        environment["REDUCE_HET_DICT"] = str(Path(profile["het_dictionary"]).resolve())
    command = [profile["executable"], *profile["build_arguments"], str(pocket_pdb.resolve())]
    completed = runner(command, capture_output=True, text=True, check=False, env=environment)
    if completed.returncode != 0:
        raise RuntimeError(f"Reduce failed ({completed.returncode}): {completed.stderr.strip()}")
    hydrogenated = output_dir / "query-pocket.reduce.pdb"
    hydrogenated.write_text(completed.stdout, encoding="utf-8")
    log = output_dir / "reduce.stderr.log"
    log.write_text(completed.stderr, encoding="utf-8")
    manifest = {"format": "aidd-reduce-query-pocket", "version": 1,
                "scope": "query_ligand_and_complete_protein_residues_within_5A",
                "input": {"path": str(pocket_pdb.resolve()), "sha256": _sha256(pocket_pdb)},
                "command": command, "het_dictionary": profile["het_dictionary"],
                "returncode": int(completed.returncode),
                "output": {"path": str(hydrogenated.resolve()), "sha256": _sha256(hydrogenated)},
                "stderr_log": {"path": str(log.resolve()), "sha256": _sha256(log)}}
    manifest_path = output_dir / "reduce_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return manifest
