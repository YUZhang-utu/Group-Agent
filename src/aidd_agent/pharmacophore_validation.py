from __future__ import annotations

import json
from pathlib import Path
import time
from typing import Iterable

import numpy as np

from .chemistry_prep import _ccd_molecule, enumerate_ligand_instances, standardize_parent
from .pharmacophore_index import (
    build_pharmacophore_index, compile_pharmacophore_query,
    load_external_l1_ids, search_pharmacophore_index,
)
from .similarity import usrcat_descriptor


def generate_faiss_l1_ids(index_dir: Path, artifact_catalog_path: Path,
                          mmcif: Path, ccd: Path, ccd_id: str,
                          output: Path, *, k: int = 100_000,
                          nprobe: int = 256) -> dict:
    """Generate the baseline L1 ID set from an existing incremental FAISS index."""
    import faiss

    if k <= 0 or nprobe <= 0:
        raise ValueError("FAISS k and nprobe must be positive")
    catalog = json.loads(artifact_catalog_path.read_text(encoding="utf-8"))
    manifest = json.loads((index_dir / "manifest.json").read_text(encoding="utf-8"))
    if int(manifest["final_ntotal"]) != int(catalog["conformers"]):
        raise ValueError("FAISS and artifact catalog conformer counts differ")
    final_record = manifest["additions"][-1]
    configured = Path(final_record["index_path"])
    final_index = configured if configured.exists() else index_dir / configured.name
    if not final_index.exists():
        fallback = index_dir / f"index_after_{len(catalog['shards'])}_shards.faiss"
        if not fallback.exists():
            raise FileNotFoundError("final incremental FAISS index was not found")
        final_index = fallback
    transform = np.load(index_dir / "transform.npz")
    mean, std = transform["mean"], transform["std"]
    if mean.shape != (60,) or std.shape != (60,) or np.any(std <= 0):
        raise ValueError("invalid FAISS USRCAT transform")

    instances = enumerate_ligand_instances(mmcif, [ccd_id])
    if len(instances) != 1:
        raise ValueError(f"Expected one {ccd_id} instance, found {len(instances)}")
    parent, smiles = standardize_parent(_ccd_molecule(ccd, instances[0]["atoms"]))
    raw = np.asarray(usrcat_descriptor(parent), dtype=np.float32)
    query = np.ascontiguousarray(((raw - mean) / std)[None], dtype=np.float32)
    index = faiss.read_index(str(final_index))
    index.nprobe = min(int(nprobe), int(index.nlist))
    started = time.perf_counter()
    _, ids = index.search(query, min(int(k), int(index.ntotal)))
    seconds = time.perf_counter() - started
    ids = ids[0]; ids = np.unique(ids[ids >= 0].astype(np.int64))
    if len(ids) and (int(ids.min()) < 0 or int(ids.max()) >= int(catalog["conformers"])):
        raise ValueError("FAISS returned an ID outside the artifact catalog")
    output.parent.mkdir(parents=True, exist_ok=True)
    np.save(output, ids)
    return {
        "index": str(final_index.resolve()), "index_ntotal": int(index.ntotal),
        "requested_k": int(k), "returned": int(len(ids)),
        "nprobe": int(index.nprobe), "query_seconds": float(seconds),
        "ccd_id": ccd_id, "standardized_smiles": smiles,
        "output": str(output.resolve()),
    }


def validate_pharmacophore_results(index_catalog_path: Path, query_plan_path: Path,
                                    result_path: Path, result_manifest_path: Path,
                                    external_l1_ids: Iterable[int] = ()) -> dict:
    """Apply structural, nesting, range, and recall-safety acceptance checks."""
    index = json.loads(index_catalog_path.read_text(encoding="utf-8"))
    query = json.loads(query_plan_path.read_text(encoding="utf-8"))
    manifest = json.loads(result_manifest_path.read_text(encoding="utf-8"))
    arrays = np.load(result_path)
    ids = np.asarray(arrays["global_ids"], dtype=np.int64)
    flags = np.asarray(arrays["source_flags"], dtype=np.uint8)
    tiers = np.asarray(arrays["highest_tier"], dtype=np.uint8)
    matched = np.asarray(arrays["matched_pair_counts"], dtype=np.uint16)
    external = np.unique(np.fromiter(external_l1_ids, dtype=np.int64))

    checks = {
        "query_has_invariant_pairs": len(query["pairs"]) > 0,
        "array_lengths_agree": len(ids) == len(flags) == len(tiers) == len(matched),
        "ids_sorted_unique": bool(len(ids) < 2 or np.all(ids[1:] > ids[:-1])),
        "ids_in_catalog_range": bool(
            not len(ids) or (ids.min() >= 0 and ids.max() < int(index["conformers"]))),
        "tier_values_valid": bool(np.all(tiers <= 3)),
        "source_flags_valid": bool(np.all((flags >= 1) & (flags <= 3))),
        "matched_shape_valid": matched.shape == (len(ids), len(query["profiles"])),
        "external_l1_preserved": bool(set(external.tolist()) <= set(ids.tolist())),
    }
    if len(external):
        positions = np.searchsorted(ids, external)
        positions_valid = bool(np.all(positions < len(ids)))
        if positions_valid:
            checks["external_l1_flagged"] = bool(np.all((flags[positions] & 2) != 0))
        else:
            checks["external_l1_flagged"] = False
    else:
        checks["external_l1_flagged"] = True

    profile_columns = {row["name"]: column for column, row in enumerate(query["profiles"])}
    tier_thresholds = {"loose": 1, "balanced": 2, "strict": 3}
    calculated_counts = {}
    for name, tier_value in tier_thresholds.items():
        if name in profile_columns:
            calculated_counts[name] = int(np.sum(tiers >= tier_value))
    checks["manifest_profile_counts_agree"] = calculated_counts == {
        name: int(manifest["profile_counts"][name]) for name in calculated_counts}
    checks["tiers_nested"] = (
        calculated_counts.get("strict", 0) <= calculated_counts.get("balanced", 0)
        <= calculated_counts.get("loose", 0))
    checks["manifest_union_count_agrees"] = int(manifest["union_count"]) == len(ids)
    return {
        "accepted": all(checks.values()), "checks": checks,
        "counts": {"union": int(len(ids)), "external_l1": int(len(external)),
                   **calculated_counts},
    }


def run_pharmacophore_validation(artifact_catalog: Path,
                                 pharmacophore_index_dir: Path,
                                 query_manifest: Path, output_dir: Path, *,
                                 external_l1_path: Path | None = None,
                                 faiss_index_dir: Path | None = None,
                                 mmcif: Path | None = None, ccd: Path | None = None,
                                 ccd_id: str | None = None, faiss_k: int = 100_000,
                                 nprobe: int = 256, bin_width: float = 0.5,
                                 max_distance: float = 20.0) -> dict:
    """One-command build/reuse, query, tier search, and acceptance validation."""
    output_dir.mkdir(parents=True, exist_ok=True)
    timings = {}
    started = time.perf_counter()
    index = build_pharmacophore_index(
        artifact_catalog, pharmacophore_index_dir,
        bin_width=bin_width, max_distance=max_distance)
    timings["build_or_reuse_index_seconds"] = time.perf_counter() - started

    query_plan_path = output_dir / "pharmacophore-query-v1.json"
    started = time.perf_counter()
    query = compile_pharmacophore_query(query_manifest, query_plan_path)
    timings["compile_query_seconds"] = time.perf_counter() - started

    faiss_report = None
    if faiss_index_dir is not None:
        required = {"mmcif": mmcif, "ccd": ccd, "ccd_id": ccd_id}
        missing = [name for name, value in required.items() if value is None]
        if missing:
            raise ValueError("FAISS query requires " + ", ".join(missing))
        external_l1_path = output_dir / "faiss-global-ids.npy"
        faiss_report = generate_faiss_l1_ids(
            faiss_index_dir, artifact_catalog, mmcif, ccd, ccd_id,
            external_l1_path, k=faiss_k, nprobe=nprobe)
    external = load_external_l1_ids(external_l1_path)

    result_path = output_dir / "pharmacophore-hits.npz"
    started = time.perf_counter()
    result = search_pharmacophore_index(
        pharmacophore_index_dir / "catalog.json", query_plan_path,
        result_path, external_l1_ids=external)
    timings["tier_search_seconds"] = time.perf_counter() - started
    validation = validate_pharmacophore_results(
        pharmacophore_index_dir / "catalog.json", query_plan_path, result_path,
        result_path.with_suffix(".manifest.json"), external)
    report = {
        "format": "aidd-pharmacophore-end-to-end-validation", "version": 1,
        "accepted": validation["accepted"], "timings": timings,
        "index": {"library_id": index["library_id"],
                  "conformers": index["conformers"],
                  "shards": len(index["shards"])},
        "query": {"query_id": query.get("query_id"),
                  "supported_anchors": query["supported_anchors"],
                  "pairs": len(query["pairs"]), "query_hash": query["query_hash"]},
        "faiss": faiss_report, "result": result, "validation": validation,
    }
    report_path = output_dir / "validation-report.json"
    report_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    return report
