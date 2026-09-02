from __future__ import annotations

import argparse
import json
from pathlib import Path
import time

import numpy as np

from .chemistry_prep import _ccd_molecule, enumerate_ligand_instances, standardize_parent
from .similarity import usrcat_descriptor


def run(vectors_path: Path, mmcif: Path, ccd: Path, ccd_id: str,
        output: Path, seed: int = 20260902) -> dict:
    import faiss
    raw = np.load(vectors_path, mmap_mode="r")
    count, d = raw.shape
    rng = np.random.default_rng(seed)
    train_size = min(count, max(30_000, count // 100))
    train_ids = rng.choice(count, size=train_size, replace=False)
    training_raw = np.asarray(raw[train_ids], dtype=np.float32)
    mean, std = training_raw.mean(0), training_raw.std(0)
    std[std < 1e-6] = 1
    database = np.ascontiguousarray((np.asarray(raw) - mean) / std, dtype=np.float32)
    training = np.ascontiguousarray((training_raw - mean) / std, dtype=np.float32)
    instances = enumerate_ligand_instances(mmcif, [ccd_id])
    if len(instances) != 1:
        raise ValueError(f"Expected one {ccd_id} instance, found {len(instances)}")
    molecule = _ccd_molecule(ccd, instances[0]["atoms"])
    parent, smiles = standardize_parent(molecule)
    query_raw = np.asarray(usrcat_descriptor(parent), dtype=np.float32)
    query = np.ascontiguousarray(((query_raw - mean) / std)[None, :], dtype=np.float32)
    exact_started = time.perf_counter()
    distances = np.einsum("ij,ij->i", database - query[0], database - query[0])
    exact_order = np.argsort(distances, kind="stable")[:10_000]
    exact_seconds = time.perf_counter() - exact_started
    configurations = [(60, 20, "d60_m20"), (60, 30, "d60_m30"),
                      (64, 32, "d64_m32_zero_pad")]
    results = []
    for dim, m, name in configurations:
        db = np.pad(database, ((0, 0), (0, 4))) if dim == 64 else database
        tr = np.pad(training, ((0, 0), (0, 4))) if dim == 64 else training
        q = np.pad(query, ((0, 0), (0, 4))) if dim == 64 else query
        index = faiss.IndexIVFPQ(faiss.IndexFlatL2(dim), dim, 512, m, 8)
        t0 = time.perf_counter(); index.train(tr); train_seconds = time.perf_counter() - t0
        t0 = time.perf_counter(); index.add(db); add_seconds = time.perf_counter() - t0
        candidate_k = min(100_000, count)
        for nprobe in (4, 8, 16, 32, 64, 128, 256):
            index.nprobe = nprobe
            t0 = time.perf_counter(); _, found = index.search(q, candidate_k)
            query_seconds = time.perf_counter() - t0
            recalls = {}
            for budget in (10_000, candidate_k):
                found_set = set(found[0, :budget].tolist())
                for truth_k in (100, 1000, 10000):
                    recalls[f"truth{truth_k}_in_candidates{budget}"] = (
                        len(set(exact_order[:truth_k].tolist()) & found_set) / truth_k)
            results.append({"configuration": name, "m": m, "nprobe": nprobe,
                            "train_seconds": train_seconds, "add_seconds": add_seconds,
                            "candidate_k": candidate_k,
                            "query_seconds": query_seconds, **recalls})
    report = {"format": "aidd-real-faiss-query-benchmark", "version": 1,
              "library_conformers": count, "query": {"pdb_id": mmcif.stem.upper(),
              "ccd_id": ccd_id, "smiles": smiles, "heavy_atoms": parent.GetNumHeavyAtoms(),
              "model": instances[0]["model"], "chain_id": instances[0]["chain_id"],
              "residue_number": instances[0]["residue_number"]},
              "parameters": {"train_size": train_size, "nlist": 512,
                             "standardization": "per-dimension-zscore"},
              "exact_seconds": exact_seconds, "results": results}
    output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    return report


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--vectors", type=Path, required=True)
    p.add_argument("--mmcif", type=Path, required=True)
    p.add_argument("--ccd", type=Path, required=True)
    p.add_argument("--ccd-id", required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    report = run(a.vectors, a.mmcif, a.ccd, a.ccd_id, a.output)
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
