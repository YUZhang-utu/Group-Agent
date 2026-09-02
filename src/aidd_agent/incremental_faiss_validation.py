from __future__ import annotations

import argparse
import json
from pathlib import Path
import time

import numpy as np

from .chemistry_prep import _ccd_molecule, enumerate_ligand_instances, standardize_parent
from .similarity import usrcat_descriptor


def validate_incremental(index_dir: Path, catalog_path: Path, mmcif: Path,
                         ccd: Path, ccd_id: str, output: Path) -> dict:
    """Validate stable IDs and query recall after append-only FAISS additions."""
    import faiss

    catalog = json.loads(catalog_path.read_text(encoding="utf-8"))
    transform = np.load(index_dir / "transform.npz")
    mean, std = transform["mean"], transform["std"]
    vectors = []
    for shard in catalog["shards"]:
        raw = np.memmap(Path(shard["path"]) / "usrcat.f32.bin", dtype="<f4", mode="r").reshape(-1, 60)
        vectors.append(np.asarray(raw))
    database = np.ascontiguousarray((np.concatenate(vectors) - mean) / std, dtype=np.float32)

    instances = enumerate_ligand_instances(mmcif, [ccd_id])
    if len(instances) != 1:
        raise ValueError(f"Expected one {ccd_id} instance, found {len(instances)}")
    parent, smiles = standardize_parent(_ccd_molecule(ccd, instances[0]["atoms"]))
    raw_query = np.asarray(usrcat_descriptor(parent), dtype=np.float32)
    query = np.ascontiguousarray(((raw_query - mean) / std)[None], dtype=np.float32)
    exact_distances = np.einsum("ij,ij->i", database - query[0], database - query[0])
    stages = []
    for stage_no, shard in enumerate(catalog["shards"], start=1):
        index = faiss.read_index(str(index_dir / f"index_after_{stage_no}_shards.faiss"))
        expected = int(shard["global_id_start"] + shard["conformers"])
        if int(index.ntotal) != expected:
            raise AssertionError(f"Stage {stage_no}: ntotal {index.ntotal} != {expected}")
        index.nprobe = min(256, index.nlist)
        k = min(100_000, int(index.ntotal))
        started = time.perf_counter()
        _, found = index.search(query, k)
        elapsed = time.perf_counter() - started
        ids = found[0]
        ids = ids[ids >= 0]
        if len(ids) and (int(ids.min()) < 0 or int(ids.max()) >= expected):
            raise AssertionError(f"Stage {stage_no}: returned an ID outside the stable range")
        found_set = set(map(int, ids))
        stage_truth = np.argsort(exact_distances[:expected], kind="stable")[:10_000]
        recalls = {
            f"recall_at_{n}": float(len(set(map(int, stage_truth[:n])) & found_set) / n)
            for n in (100, 1000, 10000)
        }
        stages.append({"stage": stage_no, "ntotal": int(index.ntotal),
                       "expected_ntotal": expected, "query_seconds": float(elapsed),
                       "returned": int(len(ids)), "minimum_id": int(ids.min()),
                       "maximum_id": int(ids.max()), **recalls})

    final_ids = set(map(int, faiss.read_index(str(index_dir / "index_after_2_shards.faiss"))
                        .search(database[[0, len(vectors[0])]], 20)[1].ravel()))
    stable_sentinels = {0, len(vectors[0])}
    report = {"format": "aidd-incremental-faiss-validation", "version": 1,
              "query": {"pdb_id": mmcif.stem.upper(), "ccd_id": ccd_id,
                        "smiles": smiles}, "stages": stages,
              "sentinel_ids": sorted(stable_sentinels),
              "sentinel_self_retrieval": stable_sentinels.issubset(final_ids)}
    output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    return report


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--index-dir", type=Path, required=True)
    parser.add_argument("--catalog", type=Path, required=True)
    parser.add_argument("--mmcif", type=Path, required=True)
    parser.add_argument("--ccd", type=Path, required=True)
    parser.add_argument("--ccd-id", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(validate_incremental(args.index_dir, args.catalog, args.mmcif,
                                          args.ccd, args.ccd_id, args.output), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
