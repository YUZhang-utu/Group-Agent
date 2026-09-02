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


def _reduce_atom_name(name: str) -> str:
    name = name.strip()
    if not 1 <= len(name) <= 4:
        raise ValueError(f"Reduce cannot represent CCD atom name {name!r}")
    if len(name) == 4:
        return name
    return f" {name:<3}"


def build_reduce_het_dictionary(ccd_path: Path, ccd_id: str, output: Path) -> dict:
    """Create a minimal query-specific Reduce dictionary from authoritative CCD."""
    import gemmi

    ccd_id = ccd_id.upper()
    block = gemmi.cif.read_file(str(ccd_path)).sole_block()
    atom_rows = block.find(["_chem_comp_atom.atom_id", "_chem_comp_atom.type_symbol"])
    bond_rows = block.find(["_chem_comp_bond.atom_id_1", "_chem_comp_bond.atom_id_2"])
    atoms = [(str(row[0]), str(row[1]).upper()) for row in atom_rows]
    if not atoms or not bond_rows:
        raise ValueError(f"CCD {ccd_id} lacks atoms or bonds")
    neighbors = {name: [] for name, _ in atoms}
    for row in bond_rows:
        left, right = str(row[0]), str(row[1])
        if left in neighbors and right in neighbors:
            neighbors[left].append(right); neighbors[right].append(left)
    hydrogen_count = sum(element == "H" for _, element in atoms)
    if hydrogen_count == 0:
        raise ValueError(f"CCD {ccd_id} contains no hydrogen definitions")
    lines = [f"RESIDUE {ccd_id:>3}{len(atoms):7d}"]
    for name, _ in atoms:
        connected = neighbors[name]
        formatted = " ".join(_reduce_atom_name(other) for other in connected)
        lines.append(f"CONECT {_reduce_atom_name(name)}{len(connected):5d} {formatted}".rstrip())
    lines.append("END")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(lines) + "\n", encoding="ascii")
    return {"format": "aidd-reduce-het-dictionary", "version": 1,
            "ccd_id": ccd_id, "ccd_path": str(ccd_path.resolve()),
            "ccd_sha256": _sha256(ccd_path), "path": str(output.resolve()),
            "sha256": _sha256(output), "atoms": len(atoms),
            "hydrogens": hydrogen_count, "bonds": len(bond_rows)}


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
    dictionary_path = None
    if profile["het_dictionary"]:
        dictionary_path = Path(profile["het_dictionary"]).resolve()
        if not dictionary_path.is_file():
            raise FileNotFoundError(f"Reduce HET dictionary does not exist: {dictionary_path}")
        environment["REDUCE_HET_DICT"] = str(dictionary_path)
    dictionary_arguments = ["-DB", str(dictionary_path)] if dictionary_path else []
    command = [profile["executable"], *dictionary_arguments,
               *profile["build_arguments"], str(pocket_pdb.resolve())]
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
                "command": command, "het_dictionary": (str(dictionary_path) if dictionary_path else None),
                "dictionary_cli_applied": bool(dictionary_path),
                "returncode": int(completed.returncode),
                "output": {"path": str(hydrogenated.resolve()), "sha256": _sha256(hydrogenated)},
                "stderr_log": {"path": str(log.resolve()), "sha256": _sha256(log)}}
    manifest_path = output_dir / "reduce_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return manifest


def _pdb_atoms(path: Path) -> list[dict]:
    atoms = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line[:6].strip() not in {"ATOM", "HETATM"}:
            continue
        element = line[76:78].strip().upper() if len(line) >= 78 else ""
        atom_name = line[12:16].strip()
        if not element:
            element = atom_name.lstrip("0123456789")[:1].upper()
        atoms.append({"record": line[:6].strip(), "atom_name": atom_name,
                      "resname": line[17:20].strip(), "chain": line[21:22].strip(),
                      "residue": line[22:26].strip(), "element": element})
    return atoms


def validate_reduce_run(output_dir: Path, ccd_id: str, chain_id: str,
                        residue_number: str) -> dict:
    """Verify that Reduce actually hydrogenated the query ligand and pocket."""
    manifest_path = output_dir / "reduce_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    input_path = Path(manifest["input"]["path"])
    output_path = Path(manifest["output"]["path"])
    log_path = Path(manifest["stderr_log"]["path"])
    hash_checks = {"input": _sha256(input_path) == manifest["input"]["sha256"],
                   "output": _sha256(output_path) == manifest["output"]["sha256"],
                   "stderr_log": _sha256(log_path) == manifest["stderr_log"]["sha256"]}
    before, after = _pdb_atoms(input_path), _pdb_atoms(output_path)

    def ligand(atom: dict) -> bool:
        return (atom["resname"] == ccd_id and atom["chain"] == chain_id
                and atom["residue"] == str(residue_number))

    before_ligand_h = sum(a["element"] == "H" and ligand(a) for a in before)
    after_ligand_h = sum(a["element"] == "H" and ligand(a) for a in after)
    before_protein_h = sum(a["element"] == "H" and not ligand(a) for a in before)
    after_protein_h = sum(a["element"] == "H" and not ligand(a) for a in after)
    log_lines = log_path.read_text(encoding="utf-8").splitlines()
    warnings = [line for line in log_lines
                if any(token in line.lower() for token in
                       ("warning", "unknown", "not found", "no connectivity"))
                or ("het" in line.lower() and any(token in line.lower() for token in
                                                   ("error", "missing", "failed")))]
    fragment_warnings = [line for line in warnings if "unbonded" in line.lower()]
    fatal_warnings = [line for line in warnings if line not in fragment_warnings]
    flips = [line for line in log_lines if "flip" in line.lower() or "user  mod" in line.lower()]
    accepted = (all(hash_checks.values()) and after_ligand_h > before_ligand_h
                and after_protein_h > before_protein_h and not fatal_warnings)
    report = {"format": "aidd-reduce-query-pocket-validation", "version": 1,
              "query": {"ccd_id": ccd_id, "chain_id": chain_id,
                        "residue_number": str(residue_number)},
              "hash_checks": hash_checks, "hydrogen_counts": {
                  "ligand_before": before_ligand_h, "ligand_after": after_ligand_h,
                  "protein_before": before_protein_h, "protein_after": after_protein_h},
              "flip_records": flips, "warnings": warnings,
              "expected_fragment_warnings": fragment_warnings,
              "fatal_warnings": fatal_warnings,
              "angle_ready": accepted, "accepted": accepted,
              "het_dictionary_mode": ("explicit_override" if manifest.get("het_dictionary")
                                      else "reduce_default_lookup")}
    (output_dir / "reduce_validation.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8")
    return report
