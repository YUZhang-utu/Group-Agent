from __future__ import annotations

import argparse
import json
from pathlib import Path
import sqlite3
import time

import numpy as np
from numpy.lib.format import open_memmap

from .mol2 import iter_mol2_records


def build_real_usrcat(db_path: Path, library_id: str, output_dir: Path) -> dict:
    from rdkit import Chem, RDLogger, rdBase
    from rdkit.Chem import rdMolDescriptors

    RDLogger.DisableLog("rdApp.warning")
    output_dir.mkdir(parents=True, exist_ok=True)
    checkpoint_path = output_dir / "checkpoint.json"
    with sqlite3.connect(db_path) as connection:
        connection.row_factory = sqlite3.Row
        rows = connection.execute(
            """SELECT c.id conformer_id, c.molecule_id, c.source_path,
                      c.source_record_index
               FROM conformer c JOIN molecule m ON m.id=c.molecule_id
               WHERE m.library_id=?
               ORDER BY c.source_path, c.source_record_index""", (library_id,)).fetchall()
    if not rows:
        raise ValueError(f"No conformers for library {library_id}")
    total = len(rows)
    mode = "r+" if checkpoint_path.is_file() else "w+"
    vectors = open_memmap(output_dir / "usrcat.f32.npy", mode=mode,
                          dtype=np.float32, shape=(total, 60))
    conformer_ids = open_memmap(output_dir / "conformer_ids.s16.npy", mode=mode,
                                dtype="S16", shape=(total,))
    molecule_ids = open_memmap(output_dir / "molecule_ids.s16.npy", mode=mode,
                               dtype="S16", shape=(total,))
    checkpoint = json.loads(checkpoint_path.read_text()) if checkpoint_path.is_file() else {
        "format": "aidd-real-usrcat-checkpoint", "version": 1, "completed_sources": [],
        "valid": 0, "invalid": 0, "errors": []}
    completed = set(checkpoint["completed_sources"])
    by_source: dict[str, list[tuple[int, sqlite3.Row]]] = {}
    for global_index, row in enumerate(rows):
        by_source.setdefault(row["source_path"], []).append((global_index, row))
    started = time.perf_counter()
    for source_path, source_rows in by_source.items():
        if source_path in completed:
            continue
        wanted = {row["source_record_index"]: (global_index, row)
                  for global_index, row in source_rows}
        seen = set()
        shard_started = time.perf_counter()
        for record in iter_mol2_records(Path(source_path)):
            target = wanted.get(record.record_index)
            if target is None:
                continue
            global_index, row = target; seen.add(record.record_index)
            molecule = Chem.MolFromMol2Block(record.raw_text, sanitize=True, removeHs=False)
            if molecule is None or molecule.GetNumConformers() != 1:
                checkpoint["invalid"] += 1
                checkpoint["errors"].append({"conformer_id": row["conformer_id"],
                                             "error": "RDKit rejected conformer"})
                vectors[global_index] = np.nan
            else:
                descriptor = rdMolDescriptors.GetUSRCAT(molecule)
                vectors[global_index] = descriptor
                checkpoint["valid"] += 1
            conformer_ids[global_index] = row["conformer_id"].encode("ascii")
            molecule_ids[global_index] = row["molecule_id"].encode("ascii")
        missing = sorted(set(wanted) - seen)
        if missing:
            raise ValueError(f"Missing {len(missing)} registered records in {source_path}")
        vectors.flush(); conformer_ids.flush(); molecule_ids.flush()
        checkpoint["completed_sources"].append(source_path)
        checkpoint["last_shard_seconds"] = time.perf_counter() - shard_started
        checkpoint["rdkit_version"] = rdBase.rdkitVersion
        checkpoint_path.write_text(json.dumps(checkpoint, indent=2), encoding="utf-8")
    checkpoint["total"] = total
    checkpoint["run_seconds"] = time.perf_counter() - started
    checkpoint_path.write_text(json.dumps(checkpoint, indent=2), encoding="utf-8")
    return checkpoint


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--db", type=Path, required=True)
    parser.add_argument("--library", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(build_real_usrcat(args.db, args.library, args.output_dir), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
