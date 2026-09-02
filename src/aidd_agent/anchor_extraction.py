from __future__ import annotations

from collections import defaultdict
import json
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np

from .chemistry_prep import _ccd_molecule, _column, _mmcif_dict, enumerate_ligand_instances
from .screening_rerank import AnchorDefinition, QueryManifest
from .reduce_pocket import _pdb_atoms


SIDECHAIN_DONORS = {
    ("ARG", "NE"), ("ARG", "NH1"), ("ARG", "NH2"), ("ASN", "ND2"),
    ("GLN", "NE2"), ("HIS", "ND1"), ("HIS", "NE2"), ("LYS", "NZ"),
    ("SER", "OG"), ("THR", "OG1"), ("TRP", "NE1"), ("TYR", "OH"),
    ("CYS", "SG"),
}
SIDECHAIN_ACCEPTORS = {
    ("ASP", "OD1"), ("ASP", "OD2"), ("GLU", "OE1"), ("GLU", "OE2"),
    ("ASN", "OD1"), ("GLN", "OE1"), ("HIS", "ND1"), ("HIS", "NE2"),
    ("SER", "OG"), ("THR", "OG1"), ("TYR", "OH"), ("CYS", "SG"),
}

VDW_RADII = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "P": 1.80,
             "F": 1.47, "CL": 1.75, "BR": 1.85, "I": 1.98}


def protein_hbond_roles(residue: str, atom_name: str) -> tuple[bool, bool]:
    residue, atom_name = residue.upper(), atom_name.upper()
    donor = atom_name == "N" and residue != "PRO"
    acceptor = atom_name in {"O", "OXT"}
    return donor or (residue, atom_name) in SIDECHAIN_DONORS, \
        acceptor or (residue, atom_name) in SIDECHAIN_ACCEPTORS


def _sphere_points(count: int = 960) -> np.ndarray:
    indices = np.arange(count, dtype=np.float64) + 0.5
    z = 1.0 - 2.0 * indices / count
    radius = np.sqrt(1.0 - z * z)
    angle = np.pi * (1.0 + np.sqrt(5.0)) * indices
    return np.column_stack((radius * np.cos(angle), radius * np.sin(angle), z))


def atomic_sasa(target_centers: np.ndarray, target_elements: Iterable[str],
                occluder_centers: np.ndarray, occluder_elements: Iterable[str],
                *, probe_radius: float = 1.4, sphere_points: int = 960,
                target_occluder_indices: Sequence[int] | None = None) -> np.ndarray:
    """Deterministic Shrake--Rupley SASA for selected atoms in one environment."""
    target_centers = np.asarray(target_centers, dtype=np.float64)
    occluder_centers = np.asarray(occluder_centers, dtype=np.float64)
    target_elements, occluder_elements = list(target_elements), list(occluder_elements)
    if len(target_centers) != len(target_elements) or len(occluder_centers) != len(occluder_elements):
        raise ValueError("coordinate and element counts differ")
    sphere = _sphere_points(sphere_points)
    occluder_radii = np.asarray([VDW_RADII.get(e.upper(), 1.70) + probe_radius
                                 for e in occluder_elements])
    own = list(target_occluder_indices) if target_occluder_indices is not None else [-1] * len(target_centers)
    areas = []
    for center, element, own_index in zip(target_centers, target_elements, own):
        radius = VDW_RADII.get(element.upper(), 1.70) + probe_radius
        surface = center + sphere * radius
        exposed = np.ones(sphere_points, dtype=bool)
        for index, (other, other_radius) in enumerate(zip(occluder_centers, occluder_radii)):
            if index == own_index:
                continue
            exposed &= np.einsum("ij,ij->i", surface - other, surface - other) >= other_radius ** 2
            if not exposed.any():
                break
        areas.append(4.0 * np.pi * radius * radius * float(exposed.mean()))
    return np.asarray(areas)


def _protein_atoms(path: Path, model: str, excluded: tuple[str, str, str]) -> list[dict]:
    data = _mmcif_dict(path)
    groups = data.get("_atom_site.group_PDB", [])
    if isinstance(groups, str):
        groups = [groups]
    n = len(groups)
    columns = {name: _column(data, name, n) for name in (
        "_atom_site.auth_comp_id", "_atom_site.auth_asym_id", "_atom_site.auth_seq_id",
        "_atom_site.auth_atom_id", "_atom_site.type_symbol", "_atom_site.Cartn_x",
        "_atom_site.Cartn_y", "_atom_site.Cartn_z", "_atom_site.pdbx_PDB_model_num")}
    result = []
    for i, group in enumerate(groups):
        if columns["_atom_site.pdbx_PDB_model_num"][i] != model:
            continue
        identity = (columns["_atom_site.auth_comp_id"][i], columns["_atom_site.auth_asym_id"][i],
                    columns["_atom_site.auth_seq_id"][i])
        if identity == excluded or group.upper() != "ATOM":
            continue
        element = columns["_atom_site.type_symbol"][i].upper()
        if element == "H":
            continue
        residue, atom_name = identity[0], columns["_atom_site.auth_atom_id"][i]
        donor, acceptor = protein_hbond_roles(residue, atom_name)
        result.append({"residue": residue, "chain": identity[1], "residue_number": identity[2],
                       "atom_name": atom_name, "element": element, "donor": donor,
                       "acceptor": acceptor, "xyz": np.asarray([
                           float(columns["_atom_site.Cartn_x"][i]),
                           float(columns["_atom_site.Cartn_y"][i]),
                           float(columns["_atom_site.Cartn_z"][i])])})
    return result


def extract_query_manifest(mmcif: Path, ccd: Path, ccd_id: str,
                           query_id: str, max_distance: float = 3.5) -> QueryManifest:
    """Extract observed direct H-bond anchors; no result-admission rule is created."""
    from rdkit import RDConfig
    from rdkit.Chem import ChemicalFeatures

    instances = enumerate_ligand_instances(mmcif, [ccd_id])
    if len(instances) != 1:
        raise ValueError(f"Expected one {ccd_id} instance, found {len(instances)}")
    instance = instances[0]
    ligand = _ccd_molecule(ccd, instance["atoms"])
    factory = ChemicalFeatures.BuildFeatureFactory(str(Path(RDConfig.RDDataDir) / "BaseFeatures.fdef"))
    features = [f for f in factory.GetFeaturesForMol(ligand) if f.GetFamily() in {"Donor", "Acceptor"}]
    protein = _protein_atoms(mmcif, instance["model"],
                             (ccd_id, instance["chain_id"], instance["residue_number"]))
    protein_xyz = np.asarray([atom["xyz"] for atom in protein])
    ligand_conf = ligand.GetConformer()
    ligand_xyz = np.asarray([list(ligand_conf.GetAtomPosition(i)) for i in range(ligand.GetNumAtoms())])
    ligand_elements = [atom.GetSymbol() for atom in ligand.GetAtoms()]
    protein_elements = [atom["element"] for atom in protein]
    isolated_sasa = atomic_sasa(ligand_xyz, ligand_elements, ligand_xyz, ligand_elements,
                                target_occluder_indices=range(len(ligand_xyz)))
    complex_sasa = atomic_sasa(ligand_xyz, ligand_elements,
                               np.vstack((ligand_xyz, protein_xyz)),
                               ligand_elements + protein_elements,
                               target_occluder_indices=range(len(ligand_xyz)))
    ligand_to_protein = np.linalg.norm(ligand_xyz[:, None, :] - protein_xyz[None, :, :], axis=2)
    pocket_residues = {(protein[j]["chain"], protein[j]["residue_number"], protein[j]["residue"])
                       for j in np.where(np.any(ligand_to_protein <= 5.0, axis=0))[0]}
    anchors = []
    for feature in features:
        family = feature.GetFamily()
        ligand_role = "HBD" if family == "Donor" else "HBA"
        center = np.asarray(feature.GetPos())
        partners = [atom for atom in protein
                    if (atom["acceptor"] if ligand_role == "HBD" else atom["donor"])
                    and np.linalg.norm(atom["xyz"] - center) <= max_distance]
        if not partners:
            continue
        partner = min(partners, key=lambda atom: np.linalg.norm(atom["xyz"] - center))
        distance = float(np.linalg.norm(partner["xyz"] - center))
        burial_count = int(np.sum(np.linalg.norm(protein_xyz - center, axis=1) <= 6.0))
        atom_indices = tuple(int(i) for i in feature.GetAtomIds())
        iso = float(isolated_sasa[list(atom_indices)].sum())
        bound = float(complex_sasa[list(atom_indices)].sum())
        relative_burial = 0.0 if iso <= 1e-12 else max(0.0, min(1.0, 1.0 - bound / iso))
        anchor_id = f"A{len(anchors)}:{ligand_role}:" + "-".join(map(str, atom_indices))
        anchors.append(AnchorDefinition(
            anchor_id, atom_indices, ligand_role, tuple(map(float, center)), (),
            {"interaction": "direct_hydrogen_bond_geometry",
             "distance_angstrom": distance, "angle_degrees": None,
             "burial": {"method": "shrake_rupley_atomic_v1", "probe_radius_angstrom": 1.4,
                        "sphere_points": 960, "isolated_sasa_angstrom2": iso,
                        "complex_sasa_angstrom2": bound,
                        "relative_burial": relative_burial,
                        "protein_heavy_atoms_within_6A": burial_count,
                        "pocket_residues_within_5A": len(pocket_residues)},
             "protein_partner": {key: partner[key] for key in
                                 ("residue", "chain", "residue_number", "atom_name")},
             "angle_status": "not_evaluated_no_explicit_hydrogen"}, 1.5))
    return QueryManifest(query_id, "direct-hbond-atomcenter-v1",
                         {"anchor": 1.5, "ordinary": 1.0, "solvent_exposed": 0.5},
                         tuple(anchors))


def donor_hydrogen_angle(donor: np.ndarray, hydrogen: np.ndarray,
                         acceptor: np.ndarray) -> float:
    left, right = donor - hydrogen, acceptor - hydrogen
    cosine = np.dot(left, right) / (np.linalg.norm(left) * np.linalg.norm(right))
    return float(np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0))))


def enrich_manifest_with_reduce(query_manifest: Path, ccd_path: Path,
                                hydrogenated_pdb: Path, output: Path) -> dict:
    """Add observed D-H-A angles and protein-side projection points."""
    import gemmi

    document = json.loads(query_manifest.read_text(encoding="utf-8"))
    atoms = _pdb_atoms(hydrogenated_pdb)
    ccd = gemmi.cif.read_file(str(ccd_path)).sole_block()
    ccd_atoms = [(str(row[0]), str(row[1]).upper()) for row in
                 ccd.find(["_chem_comp_atom.atom_id", "_chem_comp_atom.type_symbol"])]
    heavy_names = [name for name, element in ccd_atoms if element != "H"]

    def match(resname, chain, residue, atom_name):
        found = [atom for atom in atoms if atom["resname"] == resname
                 and atom["chain"] == chain and atom["residue"] == str(residue)
                 and atom["atom_name"] == atom_name]
        if len(found) != 1:
            raise ValueError(f"Expected one PDB atom {resname} {chain}:{residue}:{atom_name}, found {len(found)}")
        return found[0]

    for anchor in document["anchors"]:
        ligand_indices = anchor["ligand_atom_indices"]
        if len(ligand_indices) != 1:
            anchor["evidence"]["angle_status"] = "not_evaluated_multi_atom_feature"
            continue
        ligand_name = heavy_names[int(ligand_indices[0])]
        query_parts = document["query_id"].split(":")
        ligand = match(query_parts[1], query_parts[2], query_parts[3], ligand_name)
        partner_def = anchor["evidence"]["protein_partner"]
        partner = match(partner_def["residue"], partner_def["chain"],
                        partner_def["residue_number"], partner_def["atom_name"])
        donor, acceptor = (ligand, partner) if anchor["feature_type"] == "HBD" else (partner, ligand)
        donor_xyz, acceptor_xyz = np.asarray(donor["xyz"]), np.asarray(acceptor["xyz"])
        hydrogens = [atom for atom in atoms if atom["element"] == "H"
                     and atom["resname"] == donor["resname"] and atom["chain"] == donor["chain"]
                     and atom["residue"] == donor["residue"]
                     and np.linalg.norm(np.asarray(atom["xyz"]) - donor_xyz) <= 1.35]
        if not hydrogens:
            anchor["evidence"]["angle_status"] = "not_evaluated_no_bonded_donor_hydrogen"
            continue
        evaluated = [(donor_hydrogen_angle(donor_xyz, np.asarray(h["xyz"]), acceptor_xyz), h)
                     for h in hydrogens]
        angle, hydrogen = max(evaluated, key=lambda item: item[0])
        anchor["evidence"].update({
            "angle_degrees": angle, "angle_status": "evaluated_reduce_hydrogen",
            "donor_atom": donor["atom_name"], "donor_hydrogen": hydrogen["atom_name"],
            "hydrogen_acceptor_distance_angstrom": float(np.linalg.norm(
                np.asarray(hydrogen["xyz"]) - acceptor_xyz)),
            "projection_method": "observed_protein_partner_heavy_atom",
        })
        anchor["projected_points"] = [list(map(float, partner["xyz"]))]
    document["feature_definition_version"] = "direct-hbond-reduce-angle-projection-v1"
    document["reduce_source"] = str(hydrogenated_pdb.resolve())
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(document, indent=2), encoding="utf-8")
    return document
