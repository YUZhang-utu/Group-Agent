from __future__ import annotations

import argparse
import json
from pathlib import Path
import sqlite3
import time

import numpy as np

from .chemistry_prep import _ccd_molecule, enumerate_ligand_instances, standardize_parent
from .mol2 import iter_mol2_records


def run(db_path: Path, library_id: str, vectors_path: Path, conformer_ids_path: Path,
        molecule_ids_path: Path, mmcif: Path, ccd: Path, ccd_id: str,
        output_dir: Path) -> dict:
    from rdkit import Chem, DataStructs, RDLogger
    from rdkit.Chem import AllChem, rdMolDescriptors
    RDLogger.DisableLog("rdApp.warning")
    output_dir.mkdir(parents=True, exist_ok=True)
    instances = enumerate_ligand_instances(mmcif, [ccd_id])
    query, query_smiles = standardize_parent(_ccd_molecule(ccd, instances[0]["atoms"]))
    generator = AllChem.GetMorganGenerator(radius=2, fpSize=2048)
    query_fp = generator.GetFingerprint(query)
    with sqlite3.connect(db_path) as connection:
        connection.row_factory = sqlite3.Row
        rows = connection.execute(
            """SELECT c.source_path,c.source_record_index,c.molecule_id,c.id conformer_id
               FROM conformer c JOIN molecule m ON m.id=c.molecule_id
               WHERE m.library_id=? ORDER BY c.source_path,c.source_record_index""",
            (library_id,)).fetchall()
    lookup = {(r["source_path"], r["source_record_index"]): r for r in rows}
    seen, molecule_ids, morgan_scores, invalid = set(), [], [], []
    started = time.perf_counter()
    for source in sorted({r["source_path"] for r in rows}):
        for record in iter_mol2_records(Path(source)):
            row = lookup.get((source, record.record_index))
            if row is None or row["molecule_id"] in seen:
                continue
            seen.add(row["molecule_id"])
            mol = Chem.MolFromMol2Block(record.raw_text, sanitize=True, removeHs=False)
            if mol is None:
                invalid.append(row["molecule_id"]); continue
            fp = generator.GetFingerprint(mol)
            molecule_ids.append(row["molecule_id"])
            morgan_scores.append(DataStructs.TanimotoSimilarity(query_fp, fp))
    morgan_seconds = time.perf_counter() - started
    vectors = np.load(vectors_path, mmap_mode="r")
    conformer_ids = np.load(conformer_ids_path, mmap_mode="r").astype("U16")
    vector_molecules = np.load(molecule_ids_path, mmap_mode="r").astype("U16")
    query_usr = np.asarray(rdMolDescriptors.GetUSRCAT(query), dtype=np.float32)
    distances = np.linalg.norm(np.asarray(vectors) - query_usr, axis=1)
    best: dict[str, tuple[float, str]] = {}
    for index, mol_id in enumerate(vector_molecules):
        candidate = (float(distances[index]), str(conformer_ids[index]))
        if mol_id not in best or candidate < best[mol_id]: best[str(mol_id)] = candidate
    morgan_order = np.argsort(-np.asarray(morgan_scores), kind="stable")
    usrcat_order = sorted(best, key=lambda x: (best[x][0], x))
    top_morgan = [molecule_ids[i] for i in morgan_order[:1000]]
    top_usrcat = usrcat_order[:1000]
    counts = {}
    for mol_id in vector_molecules: counts[str(mol_id)] = counts.get(str(mol_id), 0) + 1
    report = {"format": "aidd-real-baselines", "version": 1,
              "query": {"pdb_id": mmcif.stem.upper(), "ccd_id": ccd_id,
                        "smiles": query_smiles},
              "library": {"molecules": len(best), "conformers": len(vectors),
                          "conformer_count_distribution": {
                              str(k): list(counts.values()).count(k) for k in sorted(set(counts.values()))}},
              "morgan": {"radius": 2, "bits": 2048, "seconds": morgan_seconds,
                         "valid": len(molecule_ids), "invalid": len(invalid)},
              "top1000_overlap_molecules": len(set(top_morgan) & set(top_usrcat)),
              "top_morgan": [{"molecule_id": molecule_ids[i], "score": morgan_scores[i]}
                             for i in morgan_order[:100]],
              "top_usrcat": [{"molecule_id": x, "conformer_id": best[x][1],
                              "distance": best[x][0], "similarity": 1/(1+best[x][0])}
                             for x in usrcat_order[:100]]}
    (output_dir / "real_baselines.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    return report


def main() -> int:
    p=argparse.ArgumentParser()
    for name in ("db","vectors","conformer-ids","molecule-ids","mmcif","ccd","output-dir"):
        p.add_argument("--"+name, type=Path, required=True)
    p.add_argument("--library", required=True); p.add_argument("--ccd-id", required=True)
    a=p.parse_args(); result=run(a.db,a.library,a.vectors,a.conformer_ids,a.molecule_ids,
                                 a.mmcif,a.ccd,a.ccd_id,a.output_dir)
    print(json.dumps({k:v for k,v in result.items() if not k.startswith("top_")},indent=2)); return 0

if __name__ == "__main__": raise SystemExit(main())
