from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np


PLAN_FORMAT = "aidd-multi-cocrystal-query-plan"
RESULT_FORMAT = "aidd-multi-cocrystal-molecule-aggregation"
DEFAULT_OBJECTIVE = "atomcentered_anchored_joint"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _atomic_json(path: Path, document: Mapping) -> None:
    partial = path.with_name(path.name + ".partial")
    with partial.open("w", encoding="utf-8", newline="\n") as stream:
        json.dump(document, stream, indent=2, sort_keys=True, ensure_ascii=False)
        stream.write("\n"); stream.flush(); os.fsync(stream.fileno())
    os.replace(partial, path)


def _atomic_jsonl(path: Path, rows: Sequence[Mapping]) -> None:
    partial = path.with_name(path.name + ".partial")
    with partial.open("w", encoding="utf-8", newline="\n") as stream:
        for row in rows:
            stream.write(json.dumps(row, sort_keys=True, ensure_ascii=False,
                                    separators=(",", ":")) + "\n")
        stream.flush(); os.fsync(stream.fileno())
    os.replace(partial, path)


def _resolve(plan_path: Path, value: str) -> Path:
    path = Path(value)
    return (path if path.is_absolute() else plan_path.parent / path).resolve()


def _text(value) -> str:
    return str(value).strip()


def _load_query_result(plan_path: Path, query: Mapping,
                       primary_objective: str,
                       top_conformers: int) -> tuple[list[dict], dict]:
    for key in ("query_id", "receptor_id", "site_id", "result"):
        if not _text(query.get(key, "")):
            raise ValueError(f"query record requires {key}")
    result_path = _resolve(plan_path, str(query["result"]))
    if not result_path.is_file():
        raise FileNotFoundError(result_path)
    with np.load(result_path, allow_pickle=False) as archive:
        required = {"global_ids", "molecule_ids", "conformer_ids", "objective_names"}
        if not required.issubset(archive.files):
            raise ValueError(f"detailed Gaussian result is missing identity arrays: {result_path}")
        gids = np.asarray(archive["global_ids"], dtype=np.int64)
        molecule_ids = archive["molecule_ids"].astype(str)
        conformer_ids = archive["conformer_ids"].astype(str)
        objectives = tuple(archive["objective_names"].astype(str))
        if primary_objective not in objectives:
            raise ValueError(f"primary objective {primary_objective} absent from {result_path}")
        if not (gids.ndim == molecule_ids.ndim == conformer_ids.ndim == 1
                and len(gids) == len(molecule_ids) == len(conformer_ids)):
            raise ValueError(f"identity arrays have inconsistent shapes: {result_path}")
        if len(np.unique(gids)) != len(gids) or np.any(gids < 0):
            raise ValueError(f"global IDs must be unique and non-negative: {result_path}")
        if any(not value for value in molecule_ids) or any(not value for value in conformer_ids):
            raise ValueError(f"molecule and conformer IDs must be non-empty: {result_path}")
        scores, transforms = {}, {}
        for objective in objectives:
            score_key, transform_key = f"{objective}__objective", f"{objective}__transform"
            if score_key not in archive or transform_key not in archive:
                raise ValueError(f"detailed score/transform absent for {objective}: {result_path}")
            scores[objective] = np.asarray(archive[score_key], dtype=np.float64)
            transforms[objective] = np.asarray(archive[transform_key], dtype=np.float64)
            if (scores[objective].shape != (len(gids),)
                    or transforms[objective].shape != (len(gids), 16)
                    or not np.all(np.isfinite(scores[objective]))
                    or not np.all(np.isfinite(transforms[objective]))):
                raise ValueError(f"invalid score/transform arrays for {objective}: {result_path}")
        seed_ids = (archive["best_seed_ids"].astype(str)
                    if "best_seed_ids" in archive else None)
        if seed_ids is not None and seed_ids.shape != (len(gids), len(objectives)):
            raise ValueError(f"best seed array has inconsistent shape: {result_path}")

    by_molecule: dict[str, list[int]] = {}
    for index, molecule_id in enumerate(molecule_ids):
        by_molecule.setdefault(str(molecule_id), []).append(index)
    records = []
    for molecule_id in sorted(by_molecule):
        indices = by_molecule[molecule_id]
        objective_best = {}
        for column, objective in enumerate(objectives):
            ordered = sorted(indices, key=lambda i: (-float(scores[objective][i]), int(gids[i])))
            best = ordered[0]
            objective_best[objective] = {
                "global_id": int(gids[best]),
                "conformer_id": str(conformer_ids[best]),
                "score": float(scores[objective][best]),
                "seed_id": str(seed_ids[best, column]) if seed_ids is not None else None,
                "transform": transforms[objective][best].tolist(),
            }
        primary_order = sorted(
            indices, key=lambda i: (-float(scores[primary_objective][i]), int(gids[i])))
        alternatives = [{
            "global_id": int(gids[i]), "conformer_id": str(conformer_ids[i]),
            "score": float(scores[primary_objective][i]),
            "transform": transforms[primary_objective][i].tolist(),
        } for i in primary_order[:top_conformers]]
        records.append({
            "query_id": str(query["query_id"]),
            "receptor_id": str(query["receptor_id"]),
            "site_id": str(query["site_id"]),
            "molecule_id": molecule_id,
            "primary_objective": primary_objective,
            "objective_best": objective_best,
            "primary_conformers": alternatives,
            "conformers_observed": len(indices),
            "source_result": str(result_path),
        })
    records.sort(key=lambda row: (
        -row["objective_best"][primary_objective]["score"], row["molecule_id"]))
    denominator = max(1, len(records) - 1)
    for zero_rank, record in enumerate(records):
        record["query_molecule_rank"] = zero_rank + 1
        record["query_rank_percentile"] = (
            1.0 if len(records) == 1 else 1.0 - zero_rank / denominator)
    source = {
        "query_id": str(query["query_id"]), "receptor_id": str(query["receptor_id"]),
        "site_id": str(query["site_id"]), "result": str(result_path),
        "result_sha256": _sha256(result_path), "molecules": len(records),
        "conformers": len(gids), "objectives": list(objectives),
        "metadata": query.get("metadata", {}),
    }
    if query.get("result_manifest"):
        manifest_path = _resolve(plan_path, str(query["result_manifest"]))
        result_manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        declared_result_hash = (result_manifest.get("result_sha256")
                                or result_manifest.get("result_npz_sha256"))
        if declared_result_hash and declared_result_hash != source["result_sha256"]:
            raise ValueError(f"result manifest checksum mismatch: {manifest_path}")
        source["result_manifest"] = str(manifest_path)
        source["result_manifest_sha256"] = _sha256(manifest_path)
    return records, source


def aggregate_multi_cocrystal_results(
        plan_path: Path, output_dir: Path, *,
        primary_objective: str = DEFAULT_OBJECTIVE,
        top_conformers_per_query: int = 3,
        per_query_quota: int = 1000,
        consensus_quota: int = 5000,
        global_limit_per_site: int = 20_000,
        rrf_k: int = 60) -> dict:
    """Aggregate independent detailed 3D searches into docking-ready tasks."""
    positive = (top_conformers_per_query, per_query_quota, consensus_quota,
                global_limit_per_site, rrf_k)
    if any(value <= 0 for value in positive):
        raise ValueError("aggregation quotas, Top-M, limits, and rrf_k must be positive")
    plan_path = plan_path.resolve(); output_dir = output_dir.resolve()
    plan = json.loads(plan_path.read_text(encoding="utf-8"))
    if plan.get("format") != PLAN_FORMAT or int(plan.get("version", 0)) != 1:
        raise ValueError("not a supported multi-cocrystal query plan")
    if not _text(plan.get("library_id", "")):
        raise ValueError("plan requires library_id")
    queries = plan.get("queries")
    if not isinstance(queries, list) or not queries:
        raise ValueError("plan requires at least one query")
    query_ids = [_text(row.get("query_id", "")) for row in queries]
    if len(query_ids) != len(set(query_ids)):
        raise ValueError("query_id values must be unique")
    output_dir.mkdir(parents=True, exist_ok=True)
    existing = output_dir / "manifest.json"
    configuration = {
        "primary_objective": primary_objective,
        "top_conformers_per_query": top_conformers_per_query,
        "per_query_quota": per_query_quota, "consensus_quota": consensus_quota,
        "global_limit_per_site": global_limit_per_site, "rrf_k": rrf_k,
    }
    config_hash = hashlib.sha256(json.dumps(
        {"plan_sha256": _sha256(plan_path), **configuration},
        sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    if existing.is_file():
        old = json.loads(existing.read_text(encoding="utf-8"))
        if old.get("format") != RESULT_FORMAT or old.get("config_hash") != config_hash:
            raise ValueError("existing aggregation directory has different inputs or parameters")

    evidence, sources = [], []
    for query in queries:
        records, source = _load_query_result(
            plan_path, query, primary_objective, top_conformers_per_query)
        evidence.extend(records); sources.append(source)
    evidence.sort(key=lambda row: (
        row["site_id"], row["query_id"], row["query_molecule_rank"],
        row["molecule_id"]))
    queries_by_site: dict[str, list[str]] = {}
    for source in sources:
        queries_by_site.setdefault(source["site_id"], []).append(source["query_id"])
    for site_id, site_queries in queries_by_site.items():
        if global_limit_per_site < len(site_queries) * per_query_quota:
            raise ValueError(
                f"global limit for {site_id} cannot guarantee all per-query quotas")

    grouped: dict[tuple[str, str], list[dict]] = {}
    for row in evidence:
        grouped.setdefault((row["site_id"], row["molecule_id"]), []).append(row)
    summaries = []
    for (site_id, molecule_id), rows in grouped.items():
        query_support = sorted({row["query_id"] for row in rows})
        receptor_support = sorted({row["receptor_id"] for row in rows})
        summaries.append({
            "site_id": site_id, "molecule_id": molecule_id,
            "query_support_count": len(query_support),
            "receptor_support_count": len(receptor_support),
            "query_ids": query_support, "receptor_ids": receptor_support,
            "best_query_rank": min(row["query_molecule_rank"] for row in rows),
            "best_query_rank_percentile": max(row["query_rank_percentile"] for row in rows),
            "rrf_score": sum(1.0 / (rrf_k + row["query_molecule_rank"]) for row in rows),
        })
    summaries.sort(key=lambda row: (
        row["site_id"], -row["rrf_score"], row["best_query_rank"], row["molecule_id"]))
    site_summaries: dict[str, list[dict]] = {}
    for row in summaries:
        site_summaries.setdefault(row["site_id"], []).append(row)
    for rows in site_summaries.values():
        for rank, row in enumerate(rows, start=1):
            row["site_consensus_rank"] = rank

    admissions, admission_lookup = [], {}
    for site_id in sorted(queries_by_site):
        reasons: dict[str, list[str]] = {}
        order: list[str] = []
        for query_id in sorted(queries_by_site[site_id]):
            ranked = [row for row in evidence
                      if row["site_id"] == site_id and row["query_id"] == query_id]
            for row in ranked[:per_query_quota]:
                molecule_id = row["molecule_id"]
                if molecule_id not in reasons:
                    reasons[molecule_id] = []; order.append(molecule_id)
                reasons[molecule_id].append(f"protected_query:{query_id}")
        for row in site_summaries[site_id][:consensus_quota]:
            molecule_id = row["molecule_id"]
            if molecule_id not in reasons:
                if len(order) >= global_limit_per_site:
                    break
                reasons[molecule_id] = []; order.append(molecule_id)
            reasons[molecule_id].append("site_consensus")
        summary_lookup = {row["molecule_id"]: row for row in site_summaries[site_id]}
        for site_admission_rank, molecule_id in enumerate(order, start=1):
            summary = summary_lookup[molecule_id]
            row = {**summary, "site_admission_rank": site_admission_rank,
                   "admission_reasons": reasons[molecule_id]}
            admissions.append(row); admission_lookup[(site_id, molecule_id)] = row

    tasks = []
    for evidence_row in evidence:
        key = (evidence_row["site_id"], evidence_row["molecule_id"])
        if key not in admission_lookup:
            continue
        pose = evidence_row["objective_best"][primary_objective]
        task_key = "|".join((evidence_row["site_id"], evidence_row["molecule_id"],
                             evidence_row["query_id"], evidence_row["receptor_id"]))
        tasks.append({
            "docking_task_id": "DT-" + hashlib.sha256(task_key.encode()).hexdigest()[:16].upper(),
            "library_id": str(plan["library_id"]),
            "site_id": evidence_row["site_id"],
            "molecule_id": evidence_row["molecule_id"],
            "query_id": evidence_row["query_id"],
            "receptor_id": evidence_row["receptor_id"],
            "global_id": pose["global_id"], "conformer_id": pose["conformer_id"],
            "objective": primary_objective, "search_score": pose["score"],
            "query_molecule_rank": evidence_row["query_molecule_rank"],
            "candidate_to_query_transform": pose["transform"],
            "source_result": evidence_row["source_result"],
            "admission_reasons": admission_lookup[key]["admission_reasons"],
        })
    tasks.sort(key=lambda row: (
        row["site_id"], admission_lookup[(row["site_id"], row["molecule_id"])]["site_admission_rank"],
        row["query_id"], row["receptor_id"]))

    paths = {
        "query_molecule_evidence": output_dir / "query-molecule-evidence.jsonl",
        "molecule_summary": output_dir / "molecule-summary.jsonl",
        "docking_admission": output_dir / "docking-admission.jsonl",
        "docking_tasks": output_dir / "docking-tasks.jsonl",
    }
    payloads = {
        "query_molecule_evidence": evidence, "molecule_summary": summaries,
        "docking_admission": admissions, "docking_tasks": tasks,
    }
    for name, path in paths.items():
        _atomic_jsonl(path, payloads[name])
    manifest = {
        "format": RESULT_FORMAT, "version": 1, "status": "complete",
        "config_hash": config_hash, "configuration": configuration,
        "plan": str(plan_path), "plan_sha256": _sha256(plan_path),
        "library_id": str(plan["library_id"]), "sources": sources,
        "counts": {
            "queries": len(queries), "sites": len(queries_by_site),
            "query_molecule_evidence": len(evidence),
            "unique_site_molecules": len(summaries),
            "admitted_site_molecules": len(admissions),
            "docking_tasks": len(tasks),
        },
        "outputs": {name: {"path": str(path), "sha256": _sha256(path)}
                    for name, path in paths.items()},
        "invariants": {
            "cross_query_policy": "site-scoped union",
            "raw_scores_averaged": False,
            "query_protected_quotas": True,
            "unique_admission_keys": len(admission_lookup) == len(admissions),
            "all_tasks_reference_admission": all(
                (row["site_id"], row["molecule_id"]) in admission_lookup for row in tasks),
        },
    }
    _atomic_json(existing, manifest)
    return manifest
