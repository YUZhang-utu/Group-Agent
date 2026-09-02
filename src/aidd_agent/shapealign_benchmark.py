from __future__ import annotations

import argparse
import json
from pathlib import Path
import sqlite3
import time

import numpy as np

from .chemistry_prep import _ccd_molecule, enumerate_ligand_instances, standardize_parent
from .mol2 import iter_mol2_records
from .similarity import usrcat_descriptor


def run(db: Path, vectors_path: Path, conformer_ids_path: Path, mmcif: Path,
        ccd: Path, ccd_id: str, output: Path, candidate_count: int) -> dict:
    from rdkit import Chem, RDLogger
    from rdkit.Chem import rdShapeAlign
    RDLogger.DisableLog("rdApp.warning")
    instance = enumerate_ligand_instances(mmcif, [ccd_id])[0]
    query, smiles = standardize_parent(_ccd_molecule(ccd, instance["atoms"]))
    vectors = np.load(vectors_path, mmap_mode="r")
    ids = np.load(conformer_ids_path, mmap_mode="r").astype("U16")
    query_usr = np.asarray(usrcat_descriptor(query), dtype=np.float32)
    distance = np.linalg.norm(np.asarray(vectors) - query_usr, axis=1)
    selected_indices = np.argpartition(distance, candidate_count - 1)[:candidate_count]
    selected = {str(ids[i]): float(distance[i]) for i in selected_indices}
    with sqlite3.connect(db) as connection:
        connection.row_factory = sqlite3.Row
        placeholders = ",".join("?" * len(selected))
        rows = connection.execute(
            f"SELECT id,molecule_id,source_path,source_record_index FROM conformer WHERE id IN ({placeholders})",
            tuple(selected)).fetchall()
    targets = {(r["source_path"], r["source_record_index"]): r for r in rows}
    results, invalid = [], []
    started = time.perf_counter()
    for source in sorted({r["source_path"] for r in rows}):
        for record in iter_mol2_records(Path(source)):
            row = targets.get((source, record.record_index))
            if row is None: continue
            mol = Chem.MolFromMol2Block(record.raw_text, sanitize=True, removeHs=False)
            if mol is None:
                invalid.append(row["id"]); continue
            shape_probe = Chem.Mol(mol)
            t0=time.perf_counter(); shape, _ = rdShapeAlign.AlignMol(
                query, shape_probe, useColors=False, opt_param=1.0)
            shape_seconds=time.perf_counter()-t0
            color_probe = Chem.Mol(mol)
            t0=time.perf_counter(); joint_shape, color = rdShapeAlign.AlignMol(
                query, color_probe, useColors=True, opt_param=0.5)
            joint_seconds=time.perf_counter()-t0
            results.append({"molecule_id": row["molecule_id"], "conformer_id": row["id"],
                            "usrcat_distance": selected[row["id"]],
                            "shape_only_score": float(shape),
                            "joint_shape_score": float(joint_shape),
                            "joint_color_score": float(color),
                            "rdkit_joint_sum": float(joint_shape+color),
                            "shape_seconds": shape_seconds, "joint_seconds": joint_seconds})
    results.sort(key=lambda x: (-x["rdkit_joint_sum"], x["conformer_id"]))
    report={"format":"aidd-rdkit-shapealign-baseline","version":1,
            "scope":"RDKit shape-align baseline; not anchor-weighted/Tversky production score",
            "query":{"pdb_id":mmcif.stem.upper(),"ccd_id":ccd_id,"smiles":smiles},
            "candidate_source":"exact USRCAT nearest conformers",
            "candidate_count":candidate_count,"valid":len(results),"invalid":len(invalid),
            "wall_seconds":time.perf_counter()-started,
            "shape_compute_seconds":sum(x["shape_seconds"] for x in results),
            "joint_compute_seconds":sum(x["joint_seconds"] for x in results),
            "top_results":results[:100]}
    output.write_text(json.dumps(report,indent=2),encoding="utf-8"); return report


def main():
    p=argparse.ArgumentParser()
    for x in ("db","vectors","conformer-ids","mmcif","ccd","output"):
        p.add_argument("--"+x,type=Path,required=True)
    p.add_argument("--ccd-id",required=True);p.add_argument("--candidate-count",type=int,default=1000)
    a=p.parse_args(); r=run(a.db,a.vectors,a.conformer_ids,a.mmcif,a.ccd,a.ccd_id,a.output,a.candidate_count)
    print(json.dumps({k:v for k,v in r.items() if k!="top_results"},indent=2)); return 0

if __name__=="__main__": raise SystemExit(main())
