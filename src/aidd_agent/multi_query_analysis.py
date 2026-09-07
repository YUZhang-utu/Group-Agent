from __future__ import annotations

from collections import Counter
import hashlib
from itertools import combinations
import json
import os
from pathlib import Path
from typing import Iterable, Mapping

import numpy as np


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _atomic_text(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_name(path.name + ".partial")
    with partial.open("w", encoding="utf-8", newline="\n") as stream:
        stream.write(content)
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(partial, path)


def _jsonl(path: Path) -> list[dict]:
    return [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines()
            if line]


def _manifest_output(result_dir: Path, record: Mapping) -> Path:
    configured = Path(str(record["path"]))
    path = configured if configured.is_file() else result_dir / configured.name
    if not path.is_file():
        raise FileNotFoundError(path)
    if _sha256(path) != record["sha256"]:
        raise ValueError(f"aggregation output checksum mismatch: {path}")
    return path


def _rank_correlation(left: Mapping[str, int], right: Mapping[str, int],
                      shared: set[str]) -> float | None:
    if len(shared) < 2:
        return None
    ordered = sorted(shared)
    x = np.asarray([left[key] for key in ordered], dtype=np.float64)
    y = np.asarray([right[key] for key in ordered], dtype=np.float64)
    if np.std(x) <= 1e-12 or np.std(y) <= 1e-12:
        return None
    return float(np.corrcoef(x, y)[0, 1])


def _markdown(report: Mapping) -> str:
    lines = [
        "# Multi-cocrystal pre-docking analysis",
        "",
        f"Accepted: `{str(report['accepted']).lower()}`",
        "",
        "## Molecule coverage",
        "",
        f"- Union: {report['molecules']['union']:,}",
        f"- Shared by all queries: {report['molecules']['shared_all']:,}",
        f"- Shared/union: {report['molecules']['shared_all_fraction']:.4%}",
        "",
        "| Query | Molecules | Exclusive | Shared-all coverage |",
        "|---|---:|---:|---:|",
    ]
    for query_id, row in report["queries"].items():
        lines.append(
            f"| {query_id} | {row['molecules']:,} | {row['exclusive_molecules']:,} | "
            f"{row['shared_all_fraction']:.4%} |")
    lines += ["", "## Pairwise Top-K overlap", "",
              "| Query pair | K | Intersection | Jaccard | Overlap coefficient |",
              "|---|---:|---:|---:|---:|"]
    for pair in report["pairwise"]:
        label = f"{pair['query_a']} vs {pair['query_b']}"
        for row in pair["top_k"]:
            lines.append(
                f"| {label} | {row['k']:,} | {row['intersection']:,} | "
                f"{row['jaccard']:.4f} | {row['overlap_coefficient']:.4f} |")
        correlation = pair["shared_rank_spearman"]
        correlation_text = f"{correlation:.4f}" if correlation is not None else "NA"
        lines.append(
            f"| {label} shared-rank Spearman | - | {pair['shared_molecules']:,} | "
            f"{correlation_text} | - |")
    queue = report["docking_queue"]
    lines += ["", "## Docking queue", "",
              f"- Admitted molecules: {queue['admitted_molecules']:,}",
              f"- Receptor-specific tasks: {queue['docking_tasks']:,}",
              f"- Multi-query admitted molecules: {queue['multi_query_admitted']:,}",
              "", "## Interpretation boundary", "",
              "This report measures retrieval overlap and queue construction. It does not "
              "establish binding, pose correctness, docking accuracy, scaffold diversity, "
              "or biological activity.", ""]
    return "\n".join(lines)


def analyze_multi_cocrystal_result(result_dir: Path, output_dir: Path, *,
                                    top_k: Iterable[int] = (100, 500, 1000, 5000)) -> dict:
    """Analyze query complementarity without comparing cross-query raw scores."""
    result_dir = result_dir.resolve()
    output_dir = output_dir.resolve()
    manifest_path = result_dir / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("format") != "aidd-multi-cocrystal-molecule-aggregation":
        raise ValueError("not a multi-cocrystal aggregation directory")
    if manifest.get("status") != "complete":
        raise ValueError("aggregation is not complete")
    if int(manifest.get("counts", {}).get("sites", 0)) != 1:
        raise ValueError("this overlap report requires exactly one biological site")
    requested = tuple(sorted(set(int(value) for value in top_k)))
    if not requested or requested[0] <= 0:
        raise ValueError("Top-K values must be positive")

    paths = {name: _manifest_output(result_dir, record)
             for name, record in manifest["outputs"].items()}
    evidence = _jsonl(paths["query_molecule_evidence"])
    summaries = _jsonl(paths["molecule_summary"])
    admissions = _jsonl(paths["docking_admission"])
    tasks = _jsonl(paths["docking_tasks"])

    by_query: dict[str, list[dict]] = {}
    for row in evidence:
        by_query.setdefault(row["query_id"], []).append(row)
    if len(by_query) < 2:
        raise ValueError("overlap analysis requires at least two queries")
    for rows in by_query.values():
        rows.sort(key=lambda row: (row["query_molecule_rank"], row["molecule_id"]))
    molecule_sets = {query: {row["molecule_id"] for row in rows}
                     for query, rows in by_query.items()}
    union = set().union(*molecule_sets.values())
    shared_all = set.intersection(*molecule_sets.values())

    query_report = {}
    for query_id, rows in sorted(by_query.items()):
        ids = molecule_sets[query_id]
        exclusive = {molecule for molecule in ids
                     if sum(molecule in values for values in molecule_sets.values()) == 1}
        objectives = sorted(rows[0]["objective_best"]) if rows else []
        agreement = {}
        for left, right in combinations(objectives, 2):
            same = sum(row["objective_best"][left]["global_id"] ==
                       row["objective_best"][right]["global_id"] for row in rows)
            agreement[f"{left}__{right}"] = {
                "same_conformer": same,
                "fraction": same / len(rows) if rows else 0.0,
            }
        all_same = sum(len({row["objective_best"][name]["global_id"]
                            for name in objectives}) == 1 for row in rows)
        query_report[query_id] = {
            "molecules": len(ids), "exclusive_molecules": len(exclusive),
            "shared_all_molecules": len(ids & shared_all),
            "shared_all_fraction": len(ids & shared_all) / len(ids) if ids else 0.0,
            "objective_pose_agreement": agreement,
            "all_objectives_same_conformer": {
                "count": all_same, "fraction": all_same / len(rows) if rows else 0.0},
        }

    pairwise = []
    for query_a, query_b in combinations(sorted(by_query), 2):
        ranks_a = {row["molecule_id"]: row["query_molecule_rank"]
                   for row in by_query[query_a]}
        ranks_b = {row["molecule_id"]: row["query_molecule_rank"]
                   for row in by_query[query_b]}
        shared = set(ranks_a) & set(ranks_b)
        top_rows = []
        previous = -1
        for k in requested:
            a = {molecule for molecule, rank in ranks_a.items() if rank <= k}
            b = {molecule for molecule, rank in ranks_b.items() if rank <= k}
            intersection = len(a & b)
            if intersection < previous:
                raise AssertionError("Top-K intersection must be nondecreasing")
            previous = intersection
            top_rows.append({
                "k": k, "query_a_size": len(a), "query_b_size": len(b),
                "intersection": intersection, "union": len(a | b),
                "jaccard": intersection / len(a | b) if a or b else 0.0,
                "overlap_coefficient": intersection / min(len(a), len(b))
                if a and b else 0.0,
            })
        pairwise.append({
            "query_a": query_a, "query_b": query_b,
            "shared_molecules": len(shared),
            "shared_rank_spearman": _rank_correlation(ranks_a, ranks_b, shared),
            "top_k": top_rows,
        })

    support = {(row["site_id"], row["molecule_id"]): row["query_support_count"]
               for row in summaries}
    admitted_keys = {(row["site_id"], row["molecule_id"]) for row in admissions}
    task_support: dict[tuple[str, str], set[tuple[str, str]]] = {}
    for row in tasks:
        key = (row["site_id"], row["molecule_id"])
        task_support.setdefault(key, set()).add((row["query_id"], row["receptor_id"]))
    if any(len(task_support.get(key, set())) != support[key] for key in admitted_keys):
        raise ValueError("admitted molecule task multiplicity differs from query support")
    reason_counts = Counter(reason for row in admissions for reason in row["admission_reasons"])
    task_query_counts = Counter(row["query_id"] for row in tasks)
    task_receptor_counts = Counter(row["receptor_id"] for row in tasks)
    support_distribution = Counter(support[key] for key in admitted_keys)

    report = {
        "format": "aidd-multi-cocrystal-overlap-analysis", "version": 1,
        "accepted": True, "aggregation_manifest": str(manifest_path),
        "aggregation_manifest_sha256": _sha256(manifest_path),
        "raw_cross_query_scores_compared": False,
        "molecules": {
            "union": len(union), "shared_all": len(shared_all),
            "shared_all_fraction": len(shared_all) / len(union) if union else 0.0,
        },
        "queries": query_report, "pairwise": pairwise,
        "docking_queue": {
            "admitted_molecules": len(admissions), "docking_tasks": len(tasks),
            "multi_query_admitted": sum(support[key] > 1 for key in admitted_keys),
            "support_count_distribution": dict(sorted(support_distribution.items())),
            "admission_reason_counts": dict(sorted(reason_counts.items())),
            "task_counts_by_query": dict(sorted(task_query_counts.items())),
            "task_counts_by_receptor": dict(sorted(task_receptor_counts.items())),
        },
        "input_hashes_verified": True,
    }
    if len(union) != manifest["counts"]["unique_site_molecules"]:
        raise ValueError("recomputed union does not match aggregation manifest")
    if len(admissions) != manifest["counts"]["admitted_site_molecules"]:
        raise ValueError("recomputed admission count does not match manifest")
    if len(tasks) != manifest["counts"]["docking_tasks"]:
        raise ValueError("recomputed task count does not match manifest")

    output_dir.mkdir(parents=True, exist_ok=True)
    _atomic_text(output_dir / "analysis.json",
                 json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n")
    _atomic_text(output_dir / "analysis.md", _markdown(report))
    return report
