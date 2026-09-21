"""Immutable evidence review, same-pose selection preview and docking handoff."""
from __future__ import annotations

import csv
import json
import os
from pathlib import Path
import re
import tempfile

import numpy as np

from . import library_acceptance as ev
from .expanded_wee1 import fingerprint
from .gaussian_batch import _atomic_json, _load_query, prepare_gaussian_query
from .interaction_review import OBJECTIVE
from .project_context import ensure_within


def archive(path):
    with np.load(path, allow_pickle=False) as data:
        arrays = {key: data[key] for key in data.files}
    # Retrieval stores S16 bytes; Gaussian stores U16 text. str(bytes) creates
    # a Python representation such as "b'm1'", not the actual identifier.
    for key in ("molecule_ids", "conformer_ids"):
        if key not in arrays: continue
        values = arrays[key]
        if values.ndim != 1 or values.dtype.kind not in {"S", "U"}:
            raise ValueError("Identifier arrays must contain one-dimensional bytes or text")
        try:
            decoded = [v.decode("utf-8") if isinstance(v, (bytes, np.bytes_)) else str(v) for v in values]
        except UnicodeDecodeError as exc:
            raise ValueError("Identifier bytes must be valid UTF-8") from exc
        if any(not v or "\0" in v for v in decoded):
            raise ValueError("Identifiers must be nonempty and contain no embedded NUL")
        arrays[key] = np.asarray(decoded, dtype=str)
    return arrays


def check_hashes(hashes):
    if not hashes or fingerprint(hashes) != hashes:
        raise ValueError("Upstream evidence changed; create a fresh review")


def upstream(execution, run_id, action):
    if not isinstance(run_id, str) or not re.fullmatch(r"PROMPT-[a-f0-9]{16}", run_id):
        raise ValueError("Invalid source run ID")
    runs = execution.parent.parent
    source = ensure_within(runs / run_id, runs)
    plan_path = source / "plan.json"
    plan, seal = ev.read(plan_path), ev.read(source / "plan-seal.json")
    current = ev.read(execution.parent / "plan.json")
    if any(plan[k] != current[k] for k in ("user", "project")):
        raise ValueError("Source Project ownership mismatch")
    if seal != dict(plan_sha256=ev.sha(plan_path), user=plan["user"], project=plan["project"]):
        raise ValueError("Source plan seal mismatch")
    report_path = source / "execution/report.json"
    report, marker = ev.read(report_path), ev.read(report_path.parent / "RUN_STATUS.json")
    if report.get("status") != "complete" or marker != dict(status="complete", report_sha256=ev.sha(report_path)):
        raise ValueError("Source task must have a completed report receipt")
    allowed = {action} if isinstance(action,str) else set(action)
    steps = [(sid, s) for sid, s in report["steps"].items() if s["action"] in allowed and s["status"] == "complete"]
    if len(steps) != 1:
        raise ValueError("Source task must contain exactly one matching completed step")
    sid, step = steps[0]
    receipt = ev.read(report_path.parent / (sid + ".stage.json"))
    for p in receipt["outputs"]:
        ensure_within(Path(p), source)
    check_hashes(receipt["outputs"])
    result = step["result"]
    if result != receipt["result"]:
        raise ValueError("Source step receipt mismatch")
    child = ensure_within(Path(result["report"]), source)
    return child, [plan_path, source / "plan-seal.json", report_path,
                   report_path.parent / "RUN_STATUS.json", report_path.parent / (sid + ".stage.json"),
                   *map(Path, receipt["outputs"])]


def anchor_evidence(query_path):
    manifest, query = _load_query(query_path)
    source = manifest["source"]
    paths = [query_path, query_path.with_suffix(".manifest.json")]
    for name in ("mmcif", "ccd", "query_manifest"):
        p = Path(source[name])
        if ev.sha(p) != source[name + "_sha256"]:
            raise ValueError("Crystal anchor provenance mismatch")
        paths.append(p)
    document = ev.read(source["query_manifest"])
    if document["query_id"] != source["query_id"]:
        raise ValueError("Anchor query identity mismatch")
    mapping = source.get("anchor_mapping")
    if mapping is None:
        # Legacy E034 queries lack the mapping. Rebuild only the tiny crystal query,
        # requiring exact equality to every stored array before using new metadata.
        with tempfile.TemporaryDirectory() as tmp:
            rebuilt = Path(tmp) / "query.npz"
            reconstructed = prepare_gaussian_query(Path(source["mmcif"]), Path(source["ccd"]),
                                                   Path(source["query_manifest"]), rebuilt)
            arrays = archive(rebuilt)
            if set(arrays) != set(query) or any(not np.array_equal(query[k], arrays[k]) for k in query):
                raise ValueError("Legacy query reconstruction differs; anchor mapping cannot be accepted")
            mapping = reconstructed["source"]["anchor_mapping"]
    indices = list(map(int, query["anchor_feature_indices"]))
    if set(mapping) != {a["anchor_id"] for a in document["anchors"]} or set(mapping.values()) != set(indices):
        raise ValueError("Incomplete anchor-to-feature mapping")
    rows = []
    for anchor in document["anchors"]:
        feature = mapping[anchor["anchor_id"]]
        rows.append(dict(anchor_id=source["query_id"] + "/" + anchor["anchor_id"],
                         feature_index=feature, score_column=indices.index(feature),
                         feature_class=anchor["feature_type"], ligand_atom_indices=anchor["ligand_atom_indices"],
                         atom_center=anchor["atom_center"], weight=anchor["weight"],
                         evidence=anchor["evidence"], query_id=source["query_id"],
                         evidence_class="crystal_geometry_with_reported_angle" if anchor["evidence"].get("angle_degrees") is not None else "crystal_geometry_without_reported_angle",
                         interpretation="Crystal-derived feature hypothesis; candidate match is not a validated hydrogen bond"))
    return rows, paths, query


def validate_arrays(rigid, side, anchor_indices):
    for key in ("global_ids", "molecule_ids", "conformer_ids"):
        if not np.array_equal(rigid[key], side[key]):
            raise ValueError("Candidate identities differ between pose and annotation")
    n = len(rigid["global_ids"])
    scores, assignments = side[OBJECTIVE + "__anchor_scores"], side[OBJECTIVE + "__anchor_assignments"]
    if (scores.shape != (n, len(anchor_indices)) or assignments.shape != scores.shape
            or not np.array_equal(side["query_anchor_feature_indices"], anchor_indices)
            or not np.isfinite(scores).all() or np.any((scores < 0) | (scores > 1))
            or len(np.unique(rigid["global_ids"])) != n
            or not np.isfinite(rigid[OBJECTIVE + "__objective"]).all()):
        raise ValueError("Invalid candidate/anchor arrays")


def review(search_report, output):
    search_report, output = Path(search_report), Path(output)
    report = ev.read(search_report)
    marker = ev.read(search_report.parent / "RUN_STATUS.json")
    if report["status"] != "complete" or marker != dict(status="complete", report_sha256=ev.sha(search_report)):
        raise ValueError("Search completion receipt mismatch")
    output.mkdir(parents=True, exist_ok=True)
    paths = [search_report, search_report.parent / "RUN_STATUS.json", search_report.parent / "protocol.json"]
    protocol = ev.read(paths[-1])
    queries, union, retrieval_union, retrieved_ids = [], set(), set(), set()
    for row in report["queries"]:
        exhaustive = row.get("retrieval", {}).get("mode") == "exhaustive"
        if exhaustive:
            retrieval = row["retrieval"]
            if not protocol.get("exhaustive") or retrieval.get("coverage_fraction") != 1.0 or retrieval.get("compared_conformers") != retrieval.get("total_conformers"):
                raise ValueError("Incomplete exhaustive coverage")
        elif not row["gaussian_equivalence"]["passed"] or not row["annotation_equivalence"]["passed"]:
            raise ValueError("Search equivalence did not pass")
        cfg = row["gaussian"]["config"]
        rigid_path = ensure_within(Path(row["gaussian"]["final_result"]), search_report.parent)
        label = row["query_id"].split(":")[0].lower()
        side_path = search_report.parent / label / "interaction-matches.npz"
        side_manifest_path = side_path.with_suffix(".manifest.json")
        side_manifest = ev.read(side_manifest_path)
        if side_manifest["output_sha256"] != ev.sha(side_path) or row["gaussian"]["final_result_sha256"] != ev.sha(rigid_path):
            raise ValueError("Search output checksum mismatch")
        for item in side_manifest["inputs"].values():
            if ev.sha(item["path"]) != item["sha256"]:
                raise ValueError("Annotation input changed")
            paths.append(Path(item["path"]))
        query_path = Path(cfg["query"])
        query_manifest_path = query_path.with_suffix(".manifest.json")
        if protocol["sources"].get(str(query_manifest_path.resolve())) != ev.sha(query_manifest_path):
            raise ValueError("Query metadata differs from locked search provenance")
        if side_manifest["inputs"]["query"]["path"] != str(query_path) or side_manifest["inputs"]["rigid_result"]["path"] != str(rigid_path):
            raise ValueError("Annotation lineage mismatch")
        anchors, query_paths, query = anchor_evidence(query_path)
        if any(a["query_id"] != row["query_id"] for a in anchors):
            raise ValueError("Wrong crystal query")
        r, s = archive(rigid_path), archive(side_path)
        validate_arrays(r, s, query["anchor_feature_indices"])
        schedule_path = Path(cfg["candidate_schedule"])
        if ev.sha(schedule_path) != cfg["candidate_schedule_sha256"]:
            raise ValueError("Candidate schedule changed")
        schedule = archive(schedule_path)
        if not set(map(int, r["global_ids"])).issubset(set(map(int, schedule["global_ids"]))):
            raise ValueError("Refined candidates not in retrieval")
        retrieved_molecules = dict(zip(map(int, schedule["global_ids"]), map(str, schedule["molecule_ids"])))
        if any(retrieved_molecules[int(g)] != str(m) for g,m in zip(r["global_ids"], r["molecule_ids"])):
            raise ValueError("Retrieval/refinement molecule identity mismatch")
        paths.extend([rigid_path, side_path, side_manifest_path, schedule_path, *query_paths])
        union.update(map(str, r["molecule_ids"]))
        retrieval_union.update(map(str, schedule["molecule_ids"]))
        retrieved_ids.update(map(int, schedule["global_ids"]))
        with (output / (label + "-candidate-anchors.csv")).open("w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(["global_id", "molecule_id", "conformer_id", "gaussian_score", *[a["anchor_id"] for a in anchors]])
            for i in range(len(r["global_ids"])):
                writer.writerow([int(r["global_ids"][i]), str(r["molecule_ids"][i]), str(r["conformer_ids"][i]),
                    float(r[OBJECTIVE + "__objective"][i]), *[float(s[OBJECTIVE + "__anchor_scores"][i, a["score_column"]]) for a in anchors]])
        for a in anchors:
            values = s[OBJECTIVE + "__anchor_scores"][:, a["score_column"]]
            assigned = s[OBJECTIVE + "__anchor_assignments"][:, a["score_column"]] >= 0
            a["diagnostic_counts"] = [{"minimum_score": t, "conformers": int(np.sum(assigned & (values >= t))),
                "molecules": len(set(map(str, r["molecule_ids"][assigned & (values >= t)])))} for t in (.25, .5, .75)]
        queries.append(dict(query_id=row["query_id"], anchors=anchors, retrieval=row.get("retrieval", {}),query_npz=str(query_path),
            counts=dict(retrieved_conformers=len(schedule["global_ids"]), retrieved_molecules=len(set(map(str, schedule["molecule_ids"]))),
                        refined_conformers=len(r["global_ids"]), refined_molecules=len(set(map(str, r["molecule_ids"])))),
            rigid=str(rigid_path), sidecar=str(side_path),
            artifact_catalog=side_manifest["inputs"]["artifact_catalog"]["path"],
            chemical_companion=side_manifest["inputs"]["chemical_companion"]["path"]))
    # E036 records the E034 report as a hashed source, including accepted library counts.
    library = None
    for name, checksum in protocol["sources"].items():
        p = Path(name)
        if p.name == "report.json":
            if ev.sha(p) != checksum: raise ValueError("Search source report changed")
            source_report = ev.read(p)
            if "acceptance" in source_report:
                library = {k: source_report["acceptance"][k] for k in ("library_conformers", "library_molecules")}
                paths.append(p)
    result = dict(status="complete", kind="screening_evidence", library=library, queries=queries,
        refined_molecule_union=len(union), retrieved_molecule_union=len(retrieval_union),
        retrieved_conformer_union=len(retrieved_ids), objective=OBJECTIVE, sources=fingerprint(paths),
        unsupported_classes=["hydrophobic_contacts", "pi_stacking", "salt_bridges", "water_bridges", "metal_coordination"],
        e031_changes_ranking=False, biological_quality="not_evaluated",
        policy="Diagnostic threshold counts only; no candidates selected. Anchor columns shared by mapped anchors are not independent evidence.")
    for q in queries:
        if q['retrieval'].get('mode')=='exhaustive' and (not library or q['retrieval']['compared_conformers']!=library['library_conformers']):
            raise ValueError('Exhaustive compared count differs from accepted library')
    result["classes"] = {kind: [a["anchor_id"] for q in queries for a in q["anchors"] if a["feature_class"] == kind]
                         for kind in sorted({a["feature_class"] for q in queries for a in q["anchors"]})}
    _atomic_json(output / "report.json", result)
    with (output / "anchors.csv").open("w", newline="", encoding="utf-8") as f:
        fields = ["anchor_id", "feature_class", "ligand_atom_indices", "feature_index", "protein_partner", "interaction", "distance_angstrom", "angle_degrees", "angle_status"]
        writer = csv.DictWriter(f, fieldnames=fields); writer.writeheader()
        for q in queries:
            for a in q["anchors"]:
                writer.writerow({**{k:a[k] for k in fields[:4]},
                    **{k:a["evidence"].get(k) for k in fields[4:]}})
    return result


def validate_selection(policy):
    required = {"required_anchors", "match_mode", "minimum_score"}
    if not isinstance(policy, dict) or not required <= policy.keys() or policy.keys() - required - {"max_molecules", "coarse_constraints"}:
        raise ValueError("Selection needs explicit anchors, all/any mode and minimum score")
    anchors = policy["required_anchors"]
    if not isinstance(anchors, list) or not anchors or len(anchors) > 64 or any(not isinstance(a, str) or not 1 <= len(a) <= 160 for a in anchors) or len(set(anchors)) != len(anchors):
        raise ValueError("Choose distinct anchor IDs from the evidence report")
    if policy["match_mode"] not in {"all", "any"} or type(policy["minimum_score"]) not in (float, int) or not 0 < policy["minimum_score"] <= 1:
        raise ValueError("Use all/any and a score in (0,1]; this is not a probability")
    if "max_molecules" in policy and (type(policy["max_molecules"]) is not int or not 1 <= policy["max_molecules"] <= 100000):
        raise ValueError("Invalid molecule cap")
    if 'coarse_constraints' in policy:
        from .joint_coarse import validate_constraints
        validate_constraints(policy['coarse_constraints'])
    return policy


def select_rows(r, s, columns, policy, query_id):
    """All conditions evaluated within one stored anchored-objective pose."""
    hits = (s[OBJECTIVE + "__anchor_scores"][:, columns] >= policy["minimum_score"]) & (s[OBJECTIVE + "__anchor_assignments"][:, columns] >= 0)
    passed = np.all(hits, axis=1) if policy["match_mode"] == "all" else np.any(hits, axis=1)
    if 'coarse_constraints' in policy:
        if 'joint_eligible' not in r:
            raise ValueError('Joint coarse eligibility must be evaluated before selection')
        passed &= r['joint_eligible']
    chosen = {}
    for index in np.lexsort((r["global_ids"], -r[OBJECTIVE + "__objective"])):
        if not passed[index]: continue
        mid = str(r["molecule_ids"][index])
        chosen.setdefault(mid, dict(molecule_id=mid, global_id=int(r["global_ids"][index]),
            conformer_id=str(r["conformer_ids"][index]), row_index=int(index), query_id=query_id,
            gaussian_score=float(r[OBJECTIVE + "__objective"][index])))
    return list(chosen.values()), int(passed.sum())


def preview(review_path, output, policy):
    validate_selection(policy)
    evidence = ev.read(review_path)
    if evidence.get('kind') in ('full_library_conditions','condition_funnel'):
        from .full_library_screen import preview as full_preview
        return full_preview(review_path, output, policy)
    if evidence.get("kind") != "screening_evidence": raise ValueError("Expected screening evidence")
    check_hashes(evidence["sources"])
    available = {a["anchor_id"]: (q, a) for q in evidence["queries"] for a in q["anchors"]}
    if not set(policy["required_anchors"]) <= available.keys(): raise ValueError("Unknown anchor ID; review evidence first")
    # Explicitly scope one selection to one query. Users may make separate branches.
    query_ids = {available[a][0]["query_id"] for a in policy["required_anchors"]}
    if len(query_ids) != 1: raise ValueError("Select anchors from one crystal query per preview; do not mix query frames")
    q = available[policy["required_anchors"][0]][0]
    columns = sorted({available[a][1]["score_column"] for a in policy["required_anchors"]})
    rigid = archive(q['rigid'])
    if 'coarse_constraints' in policy:
        from .joint_coarse import eligibility
        from .gaussian_batch import ArtifactCatalogReader, _load_query
        _, query = _load_query(Path(q['query_npz']))
        rigid['joint_eligible'] = eligibility(query, ArtifactCatalogReader(Path(q['artifact_catalog'])),
                                              rigid['global_ids'], policy['coarse_constraints'])
    rows, conformers = select_rows(rigid, archive(q["sidecar"]), columns, policy, q["query_id"])
    matched = len(rows)
    if "max_molecules" in policy: rows = rows[:policy["max_molecules"]]
    output = Path(output); output.mkdir(parents=True, exist_ok=True)
    result = dict(status="complete", kind="selection_preview", policy=policy, query=q,evidence_report=str(Path(review_path).resolve()),
        counts=dict(**q["counts"], matching_conformers=conformers, matching_molecules=matched,
                    selected_molecules=len(rows), selected_representatives=len(rows)), representatives=rows,
        sources={**evidence["sources"], **fingerprint([review_path])}, e031_changes_ranking=False,
        user_selection_changes_membership=True, ranking="Original Gaussian score among qualifying poses; stable global-ID ties",
        approval="Awaiting separate user export request; no docking submitted", biological_quality="not_evaluated")
    check_hashes(result["sources"])
    _atomic_json(output / "report.json", result)
    with (output / "selected-molecules.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["molecule_id", "global_id", "conformer_id", "row_index", "query_id", "gaussian_score"])
        writer.writeheader(); writer.writerows(rows)
    return result


def export(preview_path, output):
    from .gaussian_batch import ArtifactCatalogReader
    from .chemical_companion import ChemicalCompanionReader
    from .gaussian_overlay import apply_transform
    from .predocking_qc import _chemical_sdf
    selected = ev.read(preview_path)
    if selected.get("kind") != "selection_preview": raise ValueError("Expected a selection preview")
    if selected.get('representatives_jsonl'):
        from .full_library_screen import export as full_export
        return full_export(preview_path, output)
    check_hashes(selected["sources"])
    q = selected["query"]
    # Recompute selection from sealed inputs before writing an export.
    with tempfile.TemporaryDirectory() as tmp:
        review_paths = ([Path(selected['evidence_report'])] if selected.get('evidence_report') else
            [Path(p) for p in selected["sources"] if Path(p).name == "report.json" and ev.read(p).get("kind") == "screening_evidence"])
        if any(str(p.resolve()) not in selected['sources'] for p in review_paths): raise ValueError('Unsealed selection evidence')
        if len(review_paths) != 1: raise ValueError("Ambiguous evidence lineage")
        expected = preview(review_paths[0], Path(tmp), selected["policy"])
    if expected["representatives"] != selected["representatives"] or expected["counts"] != selected["counts"]:
        raise ValueError("Selection preview changed")
    output = Path(output); output.mkdir(parents=True, exist_ok=True)
    r = archive(q["rigid"])
    rows = selected["representatives"]
    if rows and os.name == "posix":
        from .expanded_wee1 import ensure_file_descriptor_limit
        ensure_file_descriptor_limit(len(ev.read(q["artifact_catalog"])["shards"]))
    # A zero-result preview has a valid empty handoff and does not open library readers.
    reader = ArtifactCatalogReader(Path(q["artifact_catalog"])) if rows else None
    chemistry = ChemicalCompanionReader(Path(q["chemical_companion"])) if rows else None
    with (output / "selected-poses.sdf").open("w", encoding="utf-8") as f:
        for row in rows:
            gid, i = row["global_id"], row["row_index"]
            candidate, chem = reader.get(gid), chemistry.get(gid)
            if any(obj.molecule_id != row["molecule_id"] or obj.conformer_id != row["conformer_id"] for obj in (candidate, chem)):
                raise ValueError("Export molecule identity mismatch")
            points = apply_transform(candidate.shape_points, r[OBJECTIVE + "__transform"][i])
            text = _chemical_sdf(points, chem, q["query_id"])
            text = text.replace("$$$$", ">  <AIDD_MOLECULE_ID>\n" + row["molecule_id"] + "\n\n$$$$")
            f.write(text)
    _atomic_json(output / "selected-ids.json", rows)
    check_hashes(selected["sources"])
    result = dict(status="complete", kind="docking_handoff", counts=selected["counts"], policy=selected["policy"],query=q,
        sdf=str(output / "selected-poses.sdf"), ids=str(output / "selected-ids.json"),
        sources={**selected["sources"], **fingerprint([preview_path])},
        outputs=fingerprint([output / "selected-poses.sdf", output / "selected-ids.json"]),
        docking_status="not_run_requires_engine_specific_preparation", e031_changes_ranking=False,
        preparation="Heavy-atom rigid poses only; hydrogens/protonation, full stereochemical annotations and receptor/grid preparation remain required")
    _atomic_json(output / "report.json", result)
    return result
