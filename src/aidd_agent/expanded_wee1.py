"""E034: query-matched WEE1 retrieval and resumable Gaussian/E031 handoff.

No registry or library writes, docking, or biological acceptance are performed.
"""
from __future__ import annotations

import argparse
import importlib
import os
from pathlib import Path
import platform
import shutil
import sys
import time

import numpy as np

from . import library_acceptance as ev
from .gaussian_batch import (
    _atomic_json, _load_query, prepare_gaussian_query,
    run_scaled_gaussian_reranking,
)
from .interaction_matching import score_interaction_matches

QUERY_SPECS = (("8bju", "8BJU:QT9:A:601", "8BJU.cif", "QT9.cif"),
               ("1x8b", "1X8B:824:A:901", "1X8B.cif", "824.cif"))


def fingerprint(paths):
    return {str(Path(p).resolve()): ev.sha(p) for p in paths}


def runtime_versions():
    # Import names work for both pip and conda FAISS distributions.
    return {name: str(importlib.import_module(name).__version__)
            for name in ("numpy", "rdkit", "Bio", "gemmi", "faiss")}


def ensure_file_descriptor_limit(shards):
    """Readers mmap many shard payloads; change only this process's soft limit."""
    import resource
    required = max(4096, 24 * shards + 512)
    soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
    if soft != resource.RLIM_INFINITY and soft < required:
        ev.require(hard == resource.RLIM_INFINITY or hard >= required,
                   f"Need at least {required} open files for {shards} shards; hard limit is {hard}. "
                   "Ask the workstation administrator to raise RLIMIT_NOFILE before running E034.")
        resource.setrlimit(resource.RLIMIT_NOFILE, (required, hard))
        ev.log(f"Raised this process's open-file soft limit from {soft} to {required}")
    return required


def checked_stage(root, name, inputs, operation):
    """Reuse only a committed stage with unchanged inputs and output bytes."""
    receipt = root / f"{name}.stage.json"
    hashes = fingerprint(inputs)
    if receipt.exists():
        saved = ev.read(receipt)
        ev.require(saved["inputs"] == hashes, f"{name}: changed stage inputs")
        ev.require(saved["outputs"] == fingerprint(saved["outputs"]),
                   f"{name}: changed completed output; use a fresh run for repair")
        ev.log(f"Reusing completed {name}")
        return saved["result"]
    ev.log(f"Starting {name}")
    result, outputs = operation()
    ev.require(hashes == fingerprint(inputs), f"{name}: input changed during execution")
    _atomic_json(receipt, dict(inputs=hashes, outputs=fingerprint(outputs), result=result))
    return result


def verify_acceptance(batch, evaluation):
    """Bind a completed E033 full acceptance to the current immutable batch."""
    marker = ev.read(evaluation / "EVALUATION_COMPLETE.json")
    report_path = evaluation / "report.json"
    ev.require(marker["report_sha256"] == ev.sha(report_path), "E033 report checksum mismatch")
    report = ev.read(report_path)
    ev.require(marker["acceptance_status"] == report["acceptance"]["status"] == "passed",
               "E033 full acceptance required (metadata-only is insufficient)")
    ev.require(marker["calibration_status"] == report["calibration_status"] == "gate_passed_on_panel",
               "E033 calibration gate has not passed")
    ev.require(report["acceptance"]["hashes"] == "full_generated_payloads"
               and report["parameters"]["metadata_only"] is False, "E033 full byte acceptance required")
    protocol = ev.read(evaluation / "protocol.json")
    ev.require(protocol["parameters"] == report["parameters"], "E033 protocol/report parameter mismatch")
    for key, path in (
        ("artifact_catalog_sha256", batch / "artifacts/catalog.json"),
        ("faiss_manifest_sha256", batch / "faiss/manifest.json"),
        ("transform_sha256", batch / "faiss/transform.npz"),
        ("query_sha256", evaluation / "queries.npz"),
    ):
        ev.require(protocol[key] == ev.sha(path), f"E033 lineage mismatch: {key}")
    ev.require(ev.read(evaluation / "acceptance.json") == report["acceptance"],
               "E033 acceptance/report mismatch")
    # Recheck current metadata, IDs and registry without another full payload scan.
    catalog, current = ev.accept_library(batch, full=False)
    for key in ("library_conformers", "library_molecules", "shards", "registered_sources"):
        ev.require(current[key] == report["acceptance"][key], f"E033/current {key} mismatch")
    ev.require(ev.sha(batch / "faiss/index.faiss") == ev.read(batch / "faiss/manifest.json")["index_sha256"],
               "Current FAISS index checksum mismatch")
    return catalog, report["acceptance"]


def prepare_query(directory, output, spec):
    from .chemistry_prep import _ccd_molecule, _mmcif_dict, enumerate_ligand_instances, standardize_parent
    from .similarity import usrcat_descriptor

    label, expected, cif_name, ccd_name = spec
    output.mkdir(parents=True, exist_ok=True)
    sources = [directory / cif_name, directory / ccd_name, directory / "query_manifest.json"]
    for path in sources:
        shutil.copyfile(path, output / path.name)
    mmcif, ccd, manifest = [output / p.name for p in sources]
    document = ev.read(manifest)
    ev.require(document["query_id"] == expected, f"Expected locked query {expected}")
    ev.require(bool(document.get("anchors")), f"{expected}: no E031 anchors")
    pdb_id, ccd_id, chain, residue = expected.split(":")
    ev.require(_mmcif_dict(mmcif).get("_entry.id", [""])[0].upper() == pdb_id,
               "Query mmCIF entry does not match locked PDB")
    ev.require(ccd_id in _mmcif_dict(ccd).get("_chem_comp.id", []), "Wrong CCD identity")
    rows = [r for r in enumerate_ligand_instances(mmcif, [ccd_id])
            if (r["chain_id"], r["residue_number"]) == (chain, residue)]
    ev.require(len(rows) == 1, "Locked query must resolve to exactly one crystal ligand")
    parent, smiles = standardize_parent(_ccd_molecule(ccd, rows[0]["atoms"]))
    vector = np.asarray(usrcat_descriptor(parent), dtype=np.float32)
    ev.require(vector.shape == (60,) and np.isfinite(vector).all(), "Invalid crystal USRCAT")
    np.save(output / "usrcat.npy", vector)
    gaussian = output / "gaussian-query.npz"
    prepare_gaussian_query(mmcif, ccd, manifest, gaussian)
    loaded, arrays = _load_query(gaussian)
    ev.require(loaded["source"]["query_id"] == expected, "Gaussian query identity mismatch")
    ev.require(len(arrays["anchor_feature_indices"]) > 0, "Gaussian query has no anchors")
    result = dict(query_id=expected, standardized_smiles=smiles,
                  exclusion_policy="external crystal query; no verified library molecule-ID self-exclusion",
                  conventions="USRCAT standardized parent; Gaussian original CCD ligand, same crystal instance",
                  sources=fingerprint(sources), directory=str(output))
    return result, [*map(lambda p: output / p.name, sources), output / "usrcat.npy",
                    gaussian, gaussian.with_suffix(".manifest.json")]


def choose_nprobe(results, query_count, gate=.95):
    for nprobe in (128, 256):
        rows = [r for r in results if r["nprobe"] == nprobe]
        if (len(rows) == query_count and len({r["query"] for r in rows}) == query_count
                and all(r["recall"]["top1000_strict"] >= gate for r in rows)):
            return nprobe
    return None


def retrieve(batch, output, catalog, acceptance, query_dirs, threads):
    import faiss

    output.mkdir(parents=True, exist_ok=True)
    corpus = ev.Corpus(catalog)
    with np.load(batch / "faiss/transform.npz", allow_pickle=False) as tr:
        mean, std = np.asarray(tr["mean"], dtype=np.float32), np.asarray(tr["std"], dtype=np.float32)
    ev.require(mean.shape == std.shape == (60,) and np.isfinite(mean).all()
               and np.isfinite(std).all() and np.all(std > 0), "Invalid frozen transform")
    raw = np.asarray([np.load(p / "usrcat.npy", allow_pickle=False) for p in query_dirs], dtype=np.float32)
    queries = np.ascontiguousarray((raw - mean) / std, dtype=np.float32)
    faiss.omp_set_num_threads(threads)
    t = time.perf_counter()
    index = faiss.read_index(str(batch / "faiss/index.faiss"))
    load_seconds = time.perf_counter() - t
    ev.require(index.ntotal == corpus.total and index.d == 60, "FAISS shape/count mismatch")
    ev.require(index.nlist >= 256, "Index cannot support locked nprobe=256 comparison")
    truth, exclusions, exact_timing = ev.exact_truth(corpus, queries, [b""] * len(queries),
                                                   mean, std, 1000, 65536)
    rows, outputs = [], []
    for qi, (spec, q, reference) in enumerate(zip(QUERY_SPECS, queries, truth)):
        truth_path = output / f"{spec[0]}-truth.npz"
        np.savez(truth_path, global_ids=reference[1], squared_l2=reference[0])
        outputs.append(truth_path)
        for nprobe in (128, 256):
            row, arrays = ev.evaluate_setting(index, corpus, q, b"", exclusions[qi], reference,
                                              mean, std, nprobe, 10000, 5, [100, 1000], [],
                                              acceptance["library_molecules"])
            candidate = output / f"{spec[0]}-np{nprobe}-candidates.npz"
            np.savez(candidate, **arrays)
            outputs.append(candidate)
            row.update(query=spec[1], candidate_file=str(candidate),
                       truth_denominators={str(k): min(k, len(reference[1])) for k in (100, 1000)})
            rows.append(row)
    result = dict(results=rows, selected_nprobe=choose_nprobe(rows, len(queries)),
                  recall_gate=.95, candidate_budget=10000, exact_reference_timing=exact_timing,
                  index_load_seconds=load_seconds, query_panel="two target-specific calibration queries; not holdout")
    path = output / "retrieval.json"
    _atomic_json(path, result)
    return result, [*outputs, path]


def timing_gate(gaussian, sidecar_seconds):
    stage = gaussian["stages"]["refine"]
    seconds = float(stage["wall_seconds_this_invocation"])
    valid = (stage["chunks_reused"] == 0 and stage["chunks_computed"] > 0
             and np.isfinite(seconds) and seconds > 0)
    ev.require(np.isfinite(sidecar_seconds) and sidecar_seconds >= 0, "Invalid E031 call timing")
    return dict(refine_wall_seconds=seconds, e031_call_wall_seconds=sidecar_seconds,
                overhead_fraction=sidecar_seconds / seconds if valid else None,
                passed=sidecar_seconds / seconds <= .10 if valid else None,
                status="measured_same_run" if valid else "unavailable_reused_or_invalid_refine_timing",
                threshold=.10, definition="full E031 call / same-query Gaussian refine stage wall time")


def run_pose_stages(root, batch, query_dir, candidates, workers):
    artifact = batch / "artifacts/catalog.json"
    chemical = batch / "chemical/catalog.json"
    query = query_dir / "gaussian-query.npz"
    gaussian_dir = root / "gaussian"
    rigid = gaussian_dir / "refine/merged-scores.npz"

    def gaussian_operation():
        result = run_scaled_gaussian_reranking(artifact, query, candidates, gaussian_dir,
                                              workers=workers, top_n_per_objective=5000,
                                              max_pair_seeds=512, resume=True)
        ev.require(result["status"] == "complete", "Gaussian run did not complete")
        return result, sorted(p for p in gaussian_dir.rglob("*") if p.is_file() and p.suffix in (".json", ".npz"))

    gaussian = checked_stage(root, "gaussian", [artifact, query, query.with_suffix(".manifest.json"), candidates],
                             gaussian_operation)
    sidecar = root / "interaction-matches.npz"

    def interaction_operation():
        before = ev.sha(rigid)
        t = time.perf_counter()
        manifest = score_interaction_matches(artifact, chemical, query, rigid, sidecar)
        elapsed = time.perf_counter() - t
        ev.require(ev.sha(rigid) == before, "E031 changed immutable rigid results")
        with np.load(rigid, allow_pickle=False) as r, np.load(sidecar, allow_pickle=False) as s:
            for key in ("global_ids", "molecule_ids", "conformer_ids"):
                ev.require(np.array_equal(r[key], s[key]), f"E031 identity/order mismatch: {key}")
            candidate_count = gaussian["config"]["candidate_count"]
            ev.require(len(r["global_ids"]) <= candidate_count, "Refined count exceeds retrieval count")
            with np.load(candidates, allow_pickle=False) as c:
                ev.require(np.isin(r["global_ids"], c["global_ids"]).all(), "Refined ID outside candidate set")
                candidate_molecules = len(np.unique(c["molecule_ids"]))
            counts = dict(refined_conformers=len(r["global_ids"]),
                          refined_molecules=len(np.unique(r["molecule_ids"])),
                          retrieved_conformers=candidate_count, retrieved_molecules=candidate_molecules,
                          refine_conformer_retention_fraction=len(r["global_ids"]) / candidate_count,
                          refine_molecule_retention_fraction=len(np.unique(r["molecule_ids"])) / candidate_molecules,
                          e031_retention_fraction=1.0)
        return dict(manifest=manifest, counts=counts, timing=timing_gate(gaussian, elapsed)), [
            sidecar, sidecar.with_suffix(".manifest.json")]

    interaction = checked_stage(root, "e031", [artifact, chemical, query, rigid,
                                              rigid.with_suffix(".manifest.json")], interaction_operation)
    return dict(gaussian=gaussian, interaction=interaction)


def render_report(report, output):
    lines = ["# E034 WEE1 expanded-library handoff", "", f"Execution: {report['status']}", "",
             "USRCAT recall and E031 cost are engineering measurements; biological and docking quality remain unassessed.",
             "No molecule cap, torsion relaxation or E031 reranking is applied.", "",
             "| Query | nprobe | Recall Top1000 | Candidates / molecules | Search p50 (s) |", "|---|---:|---:|---:|---:|"]
    for row in report["retrieval"]["results"]:
        lines.append(f"| {row['query']} | {row['nprobe']} | {row['recall']['top1000_strict']:.4f} | "
                     f"{row['retained_conformers']} / {row['retained_molecules']} | {row['warm_search_p50_seconds']:.6f} |")
    lines += ["", f"Selected nprobe: {report['retrieval']['selected_nprobe']}", "",
              "| Query | Refined conformers / molecules | Coarse (s) | Refine (s) | E031 (s) | Overhead | <=10% |",
              "|---|---:|---:|---:|---:|---:|---|"]
    for name, row in report.get("poses", {}).items():
        interaction = row["interaction"]; timing = interaction["timing"]; counts = interaction["counts"]
        ratio = timing['overhead_fraction']
        overhead = "unavailable (resumed)" if ratio is None else f"{ratio:.2%}"
        coarse = row['gaussian']['stages']['coarse']['wall_seconds_this_invocation']
        lines.append(f"| {name} | {counts['refined_conformers']} / {counts['refined_molecules']} | "
                     f"{coarse:.3f} | {timing['refine_wall_seconds']:.3f} | {timing['e031_call_wall_seconds']:.3f} | "
                     f"{overhead} | {timing['passed']} |")
    lines += ["", "Completed stage reuse preserves original timings; partial chunk reuse cannot establish the E031 latency gate.",
              "Costs exclude preflight integrity checks and exact-reference scanning unless individually labeled.",
              "Five warm repeats / two target queries are descriptive, not a service SLA. Exact distances are descriptor distances.",
              "Budget/Top-N retention does not measure chemical rejection. Review report.json for separate gate states and evidence.", ""]
    (output / "report.md").write_text("\n".join(lines), encoding="utf-8")


def check_output_location(output, batch, evaluation, query_dirs):
    for protected in [batch, evaluation, *query_dirs]:
        ev.require(not output.is_relative_to(protected) and not protected.is_relative_to(output),
                   "Output must be separate from batch, E033 report and source queries")


def resolve_resume(output, requested):
    if not output.exists():
        if requested:
            ev.log(f"No prior E034 output at {output}; starting a fresh run")
        return False
    ev.require(requested, f"Output already exists: {output}. Use --resume for the same run, or choose a new E034_OUTPUT.")
    ev.require((output / "protocol.json").is_file(),
               f"Cannot resume {output}: protocol.json is missing. Existing contents were preserved. "
               "Set E034_OUTPUT to a new directory and run without --resume; "
               "or point E034_OUTPUT to the actual prior E034 run containing protocol.json.")
    return True


def run(args):
    # Workstation entry point: exclusive existing library lock; no library writes.
    import fcntl
    batch, evaluation, output = [p.resolve() for p in (args.batch, args.e033, args.output)]
    query_dirs = [args.qt9.resolve(), args.x8b.resolve()]
    check_output_location(output, batch, evaluation, query_dirs)
    sources = [p / filename for p, spec in zip(query_dirs, QUERY_SPECS)
               for filename in (spec[2], spec[3], "query_manifest.json")]
    for path in sources:
        ev.require(path.is_file(), f"Missing locked query input: {path}")
    ev.require(args.workers > 0 and args.threads > 0, "Workers/threads must be positive")
    resume = resolve_resume(output, args.resume)
    versions = runtime_versions()
    with (batch / "run.lock").open("r") as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        ev.log("Verifying E033 provenance and current batch before execution")
        catalog, acceptance = verify_acceptance(batch, evaluation)
        ensure_file_descriptor_limit(len(catalog["shards"]))
        # Pin all implementation modules, dependency versions, paths and small
        # library manifests; never allow a resume to mix methods or query bytes.
        module_files = sorted(Path(__file__).parent.glob("*.py"))
        library_files = [batch / name for name in ("COMPLETE.json", "batch.json", "artifacts/catalog.json",
                         "chemical/catalog.json", "pharmacophore/catalog.json", "faiss/manifest.json", "faiss/transform.npz")]
        evidence_files = [evaluation / name for name in ("report.json", "protocol.json", "acceptance.json", "EVALUATION_COMPLETE.json")]
        protocol = dict(format="aidd-e034-expanded-wee1", version=1, workers=args.workers, threads=args.threads,
                        sources=fingerprint(sources + library_files + evidence_files), code=fingerprint(module_files),
                        versions=versions, python=sys.version, host=platform.node(), cpu_count=os.cpu_count(),
                        query_ids=[s[1] for s in QUERY_SPECS], budget=10000, nprobes=[128, 256], recall_gate=.95,
                        top_n_per_objective=5000, max_pair_seeds=512, excluded_molecule_ids=["", ""])
        if resume:
            ev.require(ev.read(output / "protocol.json") == protocol, "Changed E034 inputs/code/environment; use a fresh output")
        output.mkdir(parents=True, exist_ok=resume)
        _atomic_json(output / "protocol.json", protocol)
        _atomic_json(output / "RUN_STATUS.json", dict(status="running"))
        try:
            queries = []
            for directory, spec in zip(query_dirs, QUERY_SPECS):
                destination = output / spec[0] / "query"
                inputs = [directory / name for name in (spec[2], spec[3], "query_manifest.json")]
                checked_stage(output, f"{spec[0]}-query", inputs,
                              lambda d=directory, p=destination, s=spec: prepare_query(d, p, s))
                queries.append(destination)
            retrieval = checked_stage(output, "retrieval", [output / "protocol.json", *[q / "usrcat.npy" for q in queries]],
                                      lambda: retrieve(batch, output / "retrieval", catalog, acceptance, queries, args.threads))
            report = dict(format="aidd-e034-report", status="retrieval_gate_failed", acceptance=acceptance,
                          integrity_mode="inherited E033 full payload acceptance; current metadata/IDs and index bytes rechecked",
                          retrieval=retrieval, poses={}, biological_quality="not_evaluated",
                          pose_quality="not_independently_validated", e031_changes_ranking=False)
            if retrieval["selected_nprobe"] is not None:
                for spec, query in zip(QUERY_SPECS, queries):
                    candidates = output / "retrieval" / f"{spec[0]}-np{retrieval['selected_nprobe']}-candidates.npz"
                    report["poses"][spec[1]] = run_pose_stages(output / spec[0], batch, query, candidates, args.workers)
                report["status"] = "complete"
            _atomic_json(output / "report.json", report)
            render_report(report, output)
            _atomic_json(output / "EXECUTION_COMPLETE.json", dict(status=report["status"], report_sha256=ev.sha(output / "report.json")))
            _atomic_json(output / "RUN_STATUS.json", dict(status=report["status"], report_sha256=ev.sha(output / "report.json")))
            ev.log(f"E034 {report['status']}: {output / 'report.md'}")
            return 0 if report["status"] == "complete" else 2
        except Exception as exc:
            _atomic_json(output / "RUN_STATUS.json", dict(status="failed", error=str(exc)))
            _atomic_json(output / "FAILED.json", dict(error=str(exc), type=type(exc).__name__,
                         note="Historical failed attempt; current completion is determined by report and completion marker."))
            raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--batch", type=Path, required=True)
    parser.add_argument("--e033", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--qt9", type=Path, required=True, help="Directory with 8BJU.cif, QT9.cif, query_manifest.json")
    parser.add_argument("--x8b", type=Path, required=True, help="Directory with 1X8B.cif, 824.cif, query_manifest.json")
    parser.add_argument("--workers", type=int, default=16)
    parser.add_argument("--threads", type=int, default=20)
    parser.add_argument("--resume", action="store_true")
    return run(parser.parse_args())


if __name__ == "__main__":
    raise SystemExit(main())
