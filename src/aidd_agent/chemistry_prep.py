from __future__ import annotations

from collections import defaultdict
import hashlib
import json
from pathlib import Path
import sqlite3
from typing import Any, Callable, Iterable
from urllib.request import urlopen

from .campaign import _campaign_for_user
from .mol2 import iter_mol2_records
from .similarity import build_morgan_index, save_usrcat_index, usrcat_descriptor


STANDARDIZATION_VERSION = "rdkit-cleanup-fragment-parent-uncharge-v1"
CCD_URL = "https://files.rcsb.org/ligands/download/{ccd_id}.cif"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _mmcif_dict(path: Path) -> dict[str, list[str]]:
    try:
        from Bio.PDB.MMCIF2Dict import MMCIF2Dict
    except ImportError as exc:
        raise RuntimeError("Biopython is required for mmCIF ligand extraction") from exc
    return MMCIF2Dict(str(path))


def _column(data: dict, name: str, length: int, default: str = "") -> list[str]:
    value = data.get(name, [default] * length)
    if isinstance(value, str):
        value = [value]
    if len(value) != length:
        raise ValueError(f"Malformed mmCIF column {name}")
    return [default if str(x) in {".", "?"} else str(x) for x in value]


def enumerate_ligand_instances(path: Path, retained_ccd_ids: Iterable[str]) -> list[dict[str, Any]]:
    """Enumerate coordinate-bearing retained ligands without inferring bonds."""
    retained = {str(x).upper() for x in retained_ccd_ids}
    data = _mmcif_dict(path)
    groups = data.get("_atom_site.group_PDB", [])
    if isinstance(groups, str):
        groups = [groups]
    n = len(groups)
    comp = _column(data, "_atom_site.auth_comp_id", n)
    chain = _column(data, "_atom_site.auth_asym_id", n)
    seq = _column(data, "_atom_site.auth_seq_id", n)
    insertion = _column(data, "_atom_site.pdbx_PDB_ins_code", n)
    altloc = _column(data, "_atom_site.label_alt_id", n)
    model = _column(data, "_atom_site.pdbx_PDB_model_num", n, "1")
    atom = _column(data, "_atom_site.auth_atom_id", n)
    element = _column(data, "_atom_site.type_symbol", n)
    xs = _column(data, "_atom_site.Cartn_x", n)
    ys = _column(data, "_atom_site.Cartn_y", n)
    zs = _column(data, "_atom_site.Cartn_z", n)
    buckets: dict[tuple[str, ...], list[dict[str, Any]]] = defaultdict(list)
    for i in range(n):
        ccd = comp[i].upper()
        if groups[i].upper() != "HETATM" or ccd not in retained:
            continue
        key = (model[i] or "1", chain[i], seq[i], insertion[i], altloc[i], ccd)
        buckets[key].append({"atom_id": atom[i], "element": element[i],
                             "x": float(xs[i]), "y": float(ys[i]), "z": float(zs[i])})
    return [{"model": key[0], "chain_id": key[1], "residue_number": key[2],
             "insertion_code": key[3], "altloc": key[4], "ccd_id": key[5],
             "atoms": atoms} for key, atoms in sorted(buckets.items())]


def fetch_ccd(ccd_id: str, cache_dir: Path,
              opener: Callable[..., Any] = urlopen) -> Path:
    ccd = ccd_id.strip().upper()
    if not ccd or not ccd.isalnum():
        raise ValueError(f"Invalid CCD ID: {ccd_id}")
    cache_dir.mkdir(parents=True, exist_ok=True)
    destination = cache_dir / f"{ccd}.cif"
    if destination.is_file():
        return destination
    with opener(CCD_URL.format(ccd_id=ccd), timeout=60) as response:
        payload = response.read()
    if not payload or f"data_{ccd}".lower().encode() not in payload[:512].lower():
        raise ValueError(f"RCSB response is not CCD {ccd}")
    temporary = destination.with_suffix(".cif.part")
    temporary.write_bytes(payload)
    temporary.replace(destination)
    return destination


def _rdkit_tools():
    try:
        from rdkit import Chem, rdBase
        from rdkit.Chem.MolStandardize import rdMolStandardize
    except ImportError as exc:
        raise RuntimeError("RDKit is required for chemistry preparation") from exc
    return Chem, rdBase, rdMolStandardize


def standardize_parent(molecule):
    Chem, _, standardize = _rdkit_tools()
    cleaned = standardize.Cleanup(molecule)
    parent = standardize.FragmentParent(cleaned)
    parent = standardize.Uncharger().uncharge(parent)
    Chem.SanitizeMol(parent)
    return parent, Chem.MolToSmiles(parent, canonical=True, isomericSmiles=True)


def _ccd_molecule(ccd_path: Path, coordinates: list[dict[str, Any]]):
    try:
        import gemmi
    except ImportError as exc:
        raise RuntimeError("Gemmi is required for authoritative CCD parsing") from exc
    Chem, _, _ = _rdkit_tools()
    block = gemmi.cif.read_file(str(ccd_path)).sole_block()
    atom_rows = block.find(["_chem_comp_atom.atom_id", "_chem_comp_atom.type_symbol",
                            "_chem_comp_atom.charge"])
    bond_rows = block.find(["_chem_comp_bond.atom_id_1", "_chem_comp_bond.atom_id_2",
                            "_chem_comp_bond.value_order", "_chem_comp_bond.pdbx_aromatic_flag"])
    atoms = [(str(r[0]), str(r[1]), str(r[2])) for r in atom_rows]
    bonds = [(str(r[0]), str(r[1]), str(r[2]), str(r[3])) for r in bond_rows]
    if not atoms or not bonds:
        raise ValueError(f"CCD {ccd_path.stem} has no authoritative atoms/bonds")
    by_name = {row["atom_id"]: row for row in coordinates}
    missing_heavy = sorted(name for name, element, _ in atoms
                           if element.upper() != "H" and name not in by_name)
    if missing_heavy:
        raise ValueError(
            "Crystal coordinates are missing CCD heavy atoms: " + ",".join(missing_heavy))
    present = [(name, element, charge) for name, element, charge in atoms if name in by_name]
    if not present:
        raise ValueError("No crystal atom names match the CCD definition")
    rw = Chem.RWMol()
    indices = {}
    conformer = Chem.Conformer(len(present))
    conformer.Set3D(True)
    for name, element, charge in present:
        atom = Chem.Atom(element)
        atom.SetProp("_CCDAtomName", name)
        if charge not in {"", ".", "?"}:
            atom.SetFormalCharge(int(charge))
        index = rw.AddAtom(atom)
        indices[name] = index
        xyz = by_name[name]
        conformer.SetAtomPosition(index, (xyz["x"], xyz["y"], xyz["z"]))
    orders = {"SING": Chem.BondType.SINGLE, "DOUB": Chem.BondType.DOUBLE,
              "TRIP": Chem.BondType.TRIPLE, "AROM": Chem.BondType.AROMATIC}
    for left, right, order, aromatic in bonds:
        if left in indices and right in indices:
            bond_type = Chem.BondType.AROMATIC if aromatic.upper() == "Y" else orders.get(order.upper())
            if bond_type is None:
                raise ValueError(f"Unsupported CCD bond order: {order}")
            rw.AddBond(indices[left], indices[right], bond_type)
    mol = rw.GetMol()
    mol.AddConformer(conformer)
    Chem.SanitizeMol(mol)
    return mol


def prepare_campaign_ligands(connection: sqlite3.Connection, user_id: str,
                             campaign_id: str, structures_dir: Path,
                             output_dir: Path, ccd_cache: Path) -> dict[str, Any]:
    _campaign_for_user(connection, campaign_id, user_id, write=False)
    candidates = connection.execute(
        "SELECT pdb_id, ligand_ids_json FROM structure_candidate WHERE campaign_id=? ORDER BY pdb_id",
        (campaign_id,)).fetchall()
    output_dir.mkdir(parents=True, exist_ok=True)
    records, errors = [], []
    Chem, rdBase, _ = _rdkit_tools()
    writer_dir = output_dir / "sdf"
    writer_dir.mkdir(exist_ok=True)
    for candidate in candidates:
        retained = json.loads(candidate["ligand_ids_json"])
        if not retained:
            continue
        source = structures_dir / f"{candidate['pdb_id']}.cif"
        if not source.is_file():
            errors.append({"pdb_id": candidate["pdb_id"], "error": "missing_mmcif", "path": str(source)})
            continue
        for instance in enumerate_ligand_instances(source, retained):
            identity = "_".join([candidate["pdb_id"], instance["ccd_id"], instance["model"],
                                 instance["chain_id"], instance["residue_number"],
                                 instance["insertion_code"] or "-", instance["altloc"] or "-"])
            try:
                ccd_path = fetch_ccd(instance["ccd_id"], ccd_cache)
                mol = _ccd_molecule(ccd_path, instance["atoms"])
                parent, smiles = standardize_parent(mol)
                descriptor = list(usrcat_descriptor(parent))
                sdf_path = writer_dir / f"{identity}.sdf"
                writer = Chem.SDWriter(str(sdf_path)); writer.write(parent); writer.close()
                records.append({"pdb_id": candidate["pdb_id"], **{k: v for k, v in instance.items() if k != "atoms"},
                                "standardized_smiles": smiles, "usrcat": descriptor,
                                "source_path": str(sdf_path.resolve()), "source_sha256": _sha256(sdf_path),
                                "mmcif_sha256": _sha256(source), "ccd_sha256": _sha256(ccd_path),
                                "standardization_version": STANDARDIZATION_VERSION,
                                "metadata": {"model": instance["model"],
                                             "insertion_code": instance["insertion_code"],
                                             "mmcif_path": str(source.resolve()),
                                             "mmcif_sha256": _sha256(source),
                                             "ccd_path": str(ccd_path.resolve()),
                                             "ccd_sha256": _sha256(ccd_path),
                                             "standardization_version": STANDARDIZATION_VERSION}})
            except Exception as exc:
                errors.append({"identity": identity, "error": str(exc)})
    manifest = {"format": "aidd-campaign-ligands", "version": 1,
                "standardization_version": STANDARDIZATION_VERSION,
                "rdkit_version": rdBase.rdkitVersion, "ligands": records, "errors": errors}
    path = output_dir / "campaign-ligands.json"
    path.write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    manifest["manifest_path"] = str(path.resolve()); manifest["manifest_sha256"] = _sha256(path)
    return manifest


def build_library_indices(connection: sqlite3.Connection, library_id: str,
                          output_dir: Path) -> dict[str, Any]:
    rows = connection.execute(
        """SELECT m.id molecule_id, c.id conformer_id, c.source_path, c.source_record_index
           FROM molecule m JOIN conformer c ON c.molecule_id=m.id
           WHERE m.library_id=? ORDER BY m.id, c.conformer_index, c.id""", (library_id,)).fetchall()
    if not rows:
        raise KeyError(f"Unknown or empty library: {library_id}")
    Chem, rdBase, _ = _rdkit_tools()
    parents: dict[str, str] = {}
    shapes, mapping, errors = [], {}, []
    rows_by_source: dict[str, dict[int, Any]] = defaultdict(dict)
    for row in rows:
        rows_by_source[str(row["source_path"])][row["source_record_index"]] = row
    for source_path in sorted(rows_by_source):
        pending = rows_by_source[source_path]
        seen: set[int] = set()
        for record in iter_mol2_records(Path(source_path)):
            row = pending.get(record.record_index)
            if row is None:
                continue
            seen.add(record.record_index)
            try:
                mol = Chem.MolFromMol2Block(record.raw_text, sanitize=True, removeHs=False)
                if mol is None:
                    raise ValueError("RDKit rejected MOL2 record")
                parent, smiles = standardize_parent(mol)
                previous = parents.setdefault(row["molecule_id"], smiles)
                if previous != smiles:
                    raise ValueError("Conformers disagree on standardized parent SMILES")
                shapes.append((row["conformer_id"], usrcat_descriptor(parent)))
                mapping[row["conformer_id"]] = row["molecule_id"]
            except Exception as exc:
                errors.append({"molecule_id": row["molecule_id"],
                               "conformer_id": row["conformer_id"], "error": str(exc)})
        for index in sorted(set(pending) - seen):
            row = pending[index]
            errors.append({"molecule_id": row["molecule_id"],
                           "conformer_id": row["conformer_id"],
                           "error": "Registered MOL2 source record is missing"})
    if errors:
        raise ValueError(f"Library preparation failed for {len(errors)} conformers; first: {errors[0]}")
    output_dir.mkdir(parents=True, exist_ok=True)
    morgan_path, usrcat_path = output_dir / "morgan.json.gz", output_dir / "usrcat.npz"
    morgan = build_morgan_index(sorted(parents.items()), morgan_path)
    usrcat = save_usrcat_index(shapes, usrcat_path)
    map_path = output_dir / "conformer-to-molecule.json"
    map_path.write_text(json.dumps(mapping, sort_keys=True, separators=(",", ":")), encoding="utf-8")
    manifest = {"format": "aidd-library-indices", "version": 1, "library_id": library_id,
                "standardization_version": STANDARDIZATION_VERSION, "rdkit_version": rdBase.rdkitVersion,
                "molecule_count": len(parents), "conformer_count": len(shapes),
                "morgan": {**morgan, "sha256": _sha256(morgan_path)},
                "usrcat": {**usrcat, "sha256": _sha256(usrcat_path)},
                "conformer_map": {"path": str(map_path), "sha256": _sha256(map_path)}}
    manifest_path = output_dir / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    return {**manifest, "manifest_path": str(manifest_path), "manifest_sha256": _sha256(manifest_path)}
