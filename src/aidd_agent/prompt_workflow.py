"""Project-scoped natural-language plans compiled to deterministic scientific adapters."""
from __future__ import annotations

import argparse
from contextlib import contextmanager
import hashlib
import json
import os
import platform
from pathlib import Path
import subprocess
import sys
import uuid

from . import library_acceptance as ev
from .acquisition import acquire_rcsb_mmcif
from .ai_recommendation import create_ai_review_request, import_ai_recommendation
from .expanded_wee1 import checked_stage, fingerprint
from .gaussian_batch import _atomic_json
from .prediction import write_alphafold3_input, load_model_profile, prediction_command, inspect_alphafold3_output
from .project_context import require_project_owner, require_active_project, ensure_within
from .prompt_plan import CAPABILITIES, SYSTEM_PROMPT, chat_plan, strict_json, validate_plan
from .protein_data import fetch_protein, resolve_protein, get_json
from .rcsb import search_structures
from .registry import connect
from .workflow_skills import skills_for_plan
from .prediction import validate_af3_installation


class Blocked(RuntimeError):
    pass


@contextmanager
def file_lock(path):
    with path.open("a+b") as handle:
        if os.name == "nt":
            import msvcrt
            handle.seek(0, 2)
            if handle.tell() == 0: handle.write(b"0"); handle.flush()
            handle.seek(0)
            msvcrt.locking(handle.fileno(), msvcrt.LK_NBLCK, 1)
            try: yield
            finally:
                handle.seek(0); msvcrt.locking(handle.fileno(), msvcrt.LK_UNLCK, 1)
        else:
            import fcntl
            fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
            yield


def project_root(db, user, project):
    with connect(db) as connection:
        row = require_project_owner(connection, user, project)
        require_active_project(connection, user, project)
        return Path(row["workspace_path"]).resolve()


def create_plan(db, user, project, prompt, *, llm_profile=None, response_file=None, local_plan=None):
    root = project_root(db, user, project)
    if not isinstance(prompt, str) or not 1 <= len(prompt) <= 20000: raise ValueError("Invalid prompt length")
    if sum(v is not None for v in (llm_profile, response_file, local_plan)) != 1:
        raise ValueError("Choose one LLM profile, offline fixture or local coordinator plan")
    if local_plan is not None:
        plan = validate_plan(local_plan)
        model = dict(provider="local-coordinator", model="session-bound-explicit-workflow")
    elif response_file is not None:
        plan = validate_plan(strict_json(Path(response_file).read_text(encoding="utf-8")))
        model = dict(provider="offline-fixture", model="fixture-not-a-live-LLM", fixture_sha256=ev.sha(response_file))
    else:
        plan, model = chat_plan(prompt, strict_json(Path(llm_profile).read_text(encoding="utf-8")))
        if any("source_run" in s["params"] for s in plan["steps"]):
            raise ValueError("Use separate chat evidence/selection/export actions for prior tasks")
    with connect(db) as connection:
        request_id = create_ai_review_request(connection, user, project, task_type="workflow_planning",
            subject_type="project", subject_id=project, prompt_version="e038-v1", prompt_text=prompt,
            evidence={"capabilities": CAPABILITIES}, data_classes=["user_prompt", "project_configuration"])
        recommendation = import_ai_recommendation(connection, user, request_id, provider=model["provider"],
            model_name=model["model"], model_version=None,
            response=dict(action="structured_workflow_plan", rationale=plan["summary"], confidence=0.0,
                          confidence_note="not_calibrated", evidence_refs=["capabilities"],
                          uncertainties=plan["clarifications"], requires_human_review=bool(plan["clarifications"]), plan=plan))
    directory = ensure_within(root / "runs" / ("PROMPT-" + uuid.uuid4().hex[:16]), root)
    directory.mkdir(parents=True)
    envelope = dict(format="aidd-prompt-envelope-v1", user=user, project=project, prompt=prompt,
                    workflow_skills=skills_for_plan(plan),
                    plan=plan, model=model, request_id=request_id, recommendation_id=recommendation,
                    system_prompt_sha256=hashlib.sha256(SYSTEM_PROMPT.encode()).hexdigest())
    target = directory / "plan.json"
    _atomic_json(target, envelope)
    _atomic_json(directory / "plan-seal.json", dict(plan_sha256=ev.sha(target), user=user, project=project))
    (directory / "plan.md").write_text("# Prompt workflow\n\n" + plan["summary"] + "\n\n" +
        "\n".join(f"- {s['id']}: {s['action']} {json.dumps(s['params'], ensure_ascii=False)}" for s in plan["steps"]) +
        "\n" + "\n".join(plan["clarifications"]) + "\n", encoding="utf-8")
    return target


class Services:
    mode = "live"
    fetch_json = staticmethod(get_json)
    pdb_search = staticmethod(search_structures)
    pdb_fetch = staticmethod(acquire_rcsb_mmcif)

    @staticmethod
    def run_command(argv, log):
        with log.open("a", encoding="utf-8") as stream:
            result = subprocess.run(argv, stdout=stream, stderr=subprocess.STDOUT, shell=False)
        if result.returncode:
            raise RuntimeError(f"Scientific adapter exited {result.returncode}; inspect local log {log.name}")


def load_runtime(path):
    if path is None: return {}, []
    path = Path(path).resolve()
    cfg = strict_json(path.read_text(encoding="utf-8"))
    if not isinstance(cfg, dict) or set(cfg) - {"af3_profile", "search"}: raise ValueError("Unknown runtime settings")
    inputs = [path]
    if cfg.get("af3_profile"):
        p = Path(cfg["af3_profile"])
        if not p.is_absolute(): p = path.parent / p
        cfg["af3_profile"] = str(p.resolve())
        if p.is_file():
            inputs.append(p)
            raw = ev.read(p)
            if raw.get("execution") == "apptainer":
                image = Path(raw.get("image", ""))
                # Hash once per invocation, not once for every workflow stage.
                if image.is_file(): cfg["af3_image_sha256"] = ev.sha(image)
            else:
                runner = Path(raw.get("runner", ""))
                if runner.is_file(): inputs.append(runner)
    if "search" in cfg:
        search = cfg["search"]
        if not isinstance(search, dict) or set(search) != {"batch", "e034", "workers", "coarse_chunk", "refine_chunk", "bounded_pair_seeds"}:
            raise ValueError("Invalid trusted search profile")
        for field in ("batch", "e034"):
            p = Path(search[field])
            if not p.is_absolute(): p = path.parent / p
            search[field] = str(p.resolve())
        if any(type(search[k]) is not int or search[k] <= 0 for k in ("workers", "coarse_chunk", "refine_chunk")):
            raise ValueError("Invalid search scheduling")
        if type(search["bounded_pair_seeds"]) is not bool: raise ValueError("Invalid seed engine")
    return cfg, inputs


def execute_step(step, directory, execution, results, cfg, allow_compute, services):
    action, params = step["action"], step["params"]
    if action in {"review_screening", "classify_screening", "full_library_screen", "condition_funnel", "benchmark_funnel", "select_screening", "export_screening", "prepare_docking", "run_docking"}:
        from .screening_selection import upstream, review, preview, export
        source_action = {"review_screening": "search_3d", "classify_screening":"review_screening", "full_library_screen":"classify_screening", "condition_funnel":"select_screening", "benchmark_funnel":"select_screening", "select_screening": ("review_screening","classify_screening","full_library_screen","condition_funnel"),
                         "export_screening": "select_screening", "prepare_docking":"export_screening", "run_docking":"prepare_docking"}[action]
        source, _ = upstream(execution, params["source_run"], source_action)
        output = directory / "screening"
        if action == "review_screening": summary = review(source, output)
        elif action == "classify_screening":
            if not allow_compute: raise Blocked("Enable --allow-compute to score expanded 3D features")
            from .classified_features import classify
            summary = classify(source, output)
        elif action == "benchmark_funnel":
            if not allow_compute: raise Blocked("Enable --allow-compute for a bounded hardware pilot")
            from .funnel_benchmark import run as benchmark
            if params.get('coarse_only'):
                summary = benchmark(source, output, coarse_only=True, count=10000, include_gpu=False)
            else:
                summary = benchmark(source, output)
            if summary["status"] != "complete": raise ValueError("Hardware equivalence failed; inspect screening/report.json")
        elif action == "condition_funnel":
            if not allow_compute: raise Blocked("Enable --allow-compute for a full-library condition funnel")
            from .full_library_screen import run_funnel
            summary = run_funnel(source, output)
        elif action == "full_library_screen":
            if not allow_compute: raise Blocked("Enable --allow-compute for uncapped all-conformer pose evaluation")
            from .full_library_screen import run as run_full_library
            summary = run_full_library(source, output)
        elif action == "select_screening": summary = preview(source, output, {k:v for k,v in params.items() if k != "source_run"})
        elif action in {'prepare_docking','run_docking'}:
            from .docking_preparation import prepare, run
            if action=='run_docking' and not allow_compute:raise Blocked('Enable --allow-compute for licensed preparation and docking')
            summary=(prepare if action=='prepare_docking' else run)(source,output)
        else: summary = export(source, output)
        return dict(status=summary["kind"], report=str(output / "report.json"),
                    e031_changes_ranking=False, biological_quality="not_evaluated")
    if action == "protein_fetch":
        record = fetch_protein(params["accession"], params.get("organism_id"), services.fetch_json)
    elif action == "protein_resolve":
        record = resolve_protein(params["gene"], params["organism_id"], services.fetch_json)
    else:
        record = None
    if record is not None:
        record["evidence_mode"] = services.mode
        _atomic_json(directory / "protein.json", record)
        (directory / "protein.fasta").write_text(f">{record['accession']}\n{record['sequence']}\n", encoding="utf-8")
        return dict(accession=record["accession"], organism=record["organism"], length=record["length"],
                    sequence_sha256=record["sequence_sha256"], source_url=record["source_url"])
    if action in {"pdb_search", "af3_prepare"}:
        protein = ev.read(execution / params["protein_step"] / "protein.json")
    if action == "pdb_search":
        query, candidates = services.pdb_search(protein["accession"], max_resolution=params.get("max_resolution", 3.0))
        _atomic_json(directory / "pdb-candidates.json", dict(query=query, candidates=candidates,
                      receptor_selected=False, note="UniProt matches require construct/pocket review"))
        return dict(accession=protein["accession"], candidates=len(candidates), receptor_selected=False)
    if action == "pdb_fetch":
        return services.pdb_fetch(params["pdb_id"], directory)
    if action == "af3_prepare":
        sequence = protein["sequence"]
        start, end = params.get("start", 1), params.get("end", len(sequence))
        if end > len(sequence): raise ValueError("Construct exceeds verified protein sequence")
        path = write_alphafold3_input(params["name"], [sequence[start-1:end]], directory / "af3-input.json",
                                      directory, seeds=params.get("seeds", [1]))
        payload = ev.read(path)
        for index, code in enumerate(params.get("ligand_ccd", [])):
            component = services.fetch_json(f"https://data.rcsb.org/rest/v1/core/chemcomp/{code}")
            if component.get("chem_comp", {}).get("id") != code: raise ValueError("CCD identity mismatch")
            _atomic_json(directory / f"ccd-{index}.json", component)
            payload["sequences"].append({"ligand": {"id": chr(66+index), "ccdCodes": [code]}})
        _atomic_json(path, payload)
        evidence = dict(accession=protein["accession"], source_url=protein["source_url"],
                        evidence_mode=services.mode,
                        full_sequence_sha256=protein["sequence_sha256"], start=start, end=end,
                        construct_sequence_sha256=hashlib.sha256(sequence[start-1:end].encode()).hexdigest(),
                        input_sha256=ev.sha(path), status="prepared_not_predicted")
        _atomic_json(directory / "input-provenance.json", evidence)
        return evidence
    if action in {"af3_run", "search_3d"} and not allow_compute:
        raise Blocked("Compute not enabled for this invocation; prepared inputs retained. Use --allow-compute to run.")
    if action == "af3_run":
        if not cfg.get("af3_profile") or not Path(cfg["af3_profile"]).is_file():
            raise Blocked("Configure the installed AF3 profile in the local runtime JSON")
        profile = load_model_profile(Path(cfg["af3_profile"]))
        if profile["backend"] != "alphafold3": raise ValueError("AF3 action requires alphafold3 profile")
        try:
            validate_af3_installation(profile)
        except ValueError as exc:
            raise Blocked(str(exc)) from None
        input_path = execution / params["input_step"] / "af3-input.json"
        attempt = 1
        while (directory / f"prediction-attempt-{attempt:03d}").exists(): attempt += 1
        prediction_dir = directory / f"prediction-attempt-{attempt:03d}"
        # No model-supplied argv is ever loaded or executed.
        command = prediction_command(profile, input_path, prediction_dir)
        if profile.get("execution") == "apptainer":
            prediction_dir.mkdir(parents=True)
        _atomic_json(directory / "command.json", dict(argv=command, source="local_operator_profile"))
        services.run_command(command, directory / "execution.log")
        payload = ev.read(input_path)
        manifest = inspect_alphafold3_output(prediction_dir, construct_name=payload["name"],
                    model_version=str(profile.get("version", "local_unspecified")),
                    chain_ids=[next(iter(entity.values()))["id"] for entity in payload["sequences"]], input_path=input_path)
        _atomic_json(directory / "prediction-manifest.json", manifest)
        return dict(status="prediction_completed", model=manifest["structure_path"],
                    structure_sha256=manifest["structure_sha256"], confidence=manifest["confidence"],
                    pose_quality="not_independently_validated")
    if action == "search_3d":
        search = cfg.get("search")
        if not search: raise Blocked("Configure batch and E034 evidence in the trusted local search profile")
        query = {"wee1_qt9": "8bju", "wee1_824": "1x8b", "wee1_both": "both"}[params["query"]]
        output = directory / "search"
        command = [sys.executable, "-m", "aidd_agent.fast_3d_search", "--batch", search["batch"],
                   "--e034", search["e034"], "--output", str(output), "--query", query,
                   "--workers", str(search["workers"]), "--coarse-chunk", str(search["coarse_chunk"]),
                   "--refine-chunk", str(search["refine_chunk"]), "--fresh-retrieval"]
        if search["bounded_pair_seeds"]: command.append("--bounded-pair-seeds")
        if params.get("retrieval_mode", "exhaustive") == "exhaustive": command.append("--exhaustive")
        if output.exists(): command.append("--resume")
        _atomic_json(directory / "command.json", dict(argv=command, source="allowlisted_search_adapter"))
        services.run_command(command, directory / "execution.log")
        report, marker = ev.read(output / "report.json"), ev.read(output / "RUN_STATUS.json")
        if report["status"] != "complete" or marker["status"] != "complete" or marker["report_sha256"] != ev.sha(output / "report.json"):
            raise ValueError("Search report/completion mismatch")
        return dict(status="search_completed", query=params["query"], report=str(output / "report.json"),
                    e031_changes_ranking=False, biological_quality="not_evaluated")
    raise ValueError("Unsupported action")


def run_plan(db, user, project, plan_path, *, runtime=None, allow_compute=False, services=None):
    root = project_root(db, user, project)
    path = ensure_within(Path(plan_path), root)
    seal_path = path.parent / "plan-seal.json"
    envelope, seal = ev.read(path), ev.read(seal_path)
    if envelope.get("format") != "aidd-prompt-envelope-v1" or envelope.get("user") != user or envelope.get("project") != project:
        raise ValueError("Plan ownership mismatch")
    if seal != dict(plan_sha256=ev.sha(path), user=user, project=project): raise ValueError("Plan changed after creation")
    plan = validate_plan(envelope["plan"])
    # Evidence branches consume sealed source artifacts; unrelated AF3 image hashing
    # and runtime discovery would add avoidable startup cost to every preview.
    evidence_only = bool(plan["steps"]) and all("source_run" in s["params"] for s in plan["steps"])
    cfg, config_inputs = ({}, []) if evidence_only else load_runtime(runtime)
    services = services or Services()
    execution = ensure_within(path.parent / "execution", root)
    execution.mkdir(exist_ok=True)
    with file_lock(execution / "workflow.lock"):
        protocol = dict(plan_sha256=ev.sha(path), configuration=cfg, configuration_hashes=fingerprint(config_inputs),
                        code=fingerprint(sorted(Path(__file__).parent.glob("*.py"))), service_mode=services.mode,
                        user=user, project=project, python=sys.version, platform=platform.platform())
        protocol_path = execution / "protocol.json"
        if protocol_path.exists() and ev.read(protocol_path) != protocol:
            raise ValueError("Changed workflow inputs/code/runtime; create a fresh plan")
        _atomic_json(protocol_path, protocol)
        results = {}
        report = dict(status="running", steps=results, clarifications=plan["clarifications"],
                      service_mode=services.mode, model=envelope["model"], e031_changes_ranking=False)
        _atomic_json(execution / "report.json", report)
        if plan["clarifications"]:
            report["status"] = "blocked"
        else:
            for step in plan["steps"]:
                sid = step["id"]
                directory = ensure_within(execution / sid, execution); directory.mkdir(exist_ok=True)
                dependencies = [protocol_path, path, seal_path, *config_inputs]
                for ref in ("protein_step", "input_step"):
                    if ref in step["params"]:
                        parent = execution / step["params"][ref]
                        dependencies.extend(p for p in parent.rglob("*") if p.is_file())
                def operation():
                    result = execute_step(step, directory, execution, results, cfg, allow_compute, services)
                    _atomic_json(directory / "result.json", result)
                    return result, sorted(p for p in directory.rglob("*") if p.is_file())
                try:
                    results[sid] = dict(status="running", action=step["action"])
                    _atomic_json(execution / "report.json", report)
                    if "source_run" in step["params"]:
                        from .screening_selection import upstream
                        source_action = {"review_screening": "search_3d", "classify_screening":"review_screening", "full_library_screen":"classify_screening", "condition_funnel":"select_screening", "benchmark_funnel":"select_screening", "select_screening": ("review_screening","classify_screening","full_library_screen","condition_funnel"),
                                         "export_screening": "select_screening", "prepare_docking":"export_screening", "run_docking":"prepare_docking"}[step["action"]]
                        _, source_files = upstream(execution, step["params"]["source_run"], source_action)
                        dependencies.extend(source_files)
                        if step['action']=='prepare_docking':
                            from .docking_preparation import profile_path
                            dependencies.append(profile_path())
                        for source_file in source_files:
                            if source_file.name == "report.json":
                                dependencies.extend(map(Path, ev.read(source_file).get("sources", {})))
                    result = checked_stage(execution, sid, dependencies, operation)
                    results[sid] = dict(status="complete", action=step["action"], result=result)
                    _atomic_json(execution / "report.json", report)
                except Exception as exc:
                    status = "blocked" if isinstance(exc, Blocked) else "failed"
                    # Do not include arbitrary server response bodies or credentials in reports.
                    message = str(exc) if isinstance(exc, (ValueError, Blocked)) else type(exc).__name__
                    results[sid] = dict(status=status, action=step["action"], error=message)
                    report["status"] = status
                    break
            else:
                report["status"] = "complete"
        for step in plan["steps"]:
            results.setdefault(step["id"], dict(status="not_run", action=step["action"]))
        _atomic_json(execution / "report.json", report)
        _atomic_json(execution / "RUN_STATUS.json", dict(status=report["status"], report_sha256=ev.sha(execution / "report.json")))
        lines = ["# Prompt workflow", "", plan["summary"], "", f"Status: {report['status']}",
                 f"Adapters: {services.mode}", "", "| Step | Action | Status |", "|---|---|---|"]
        lines += [f"| {sid} | {r['action']} | {r['status']} |" for sid, r in results.items()]
        lines += ["", *plan["clarifications"], "Protein/PDB evidence is not receptor acceptance; AF3 outputs are predictions."]
        (execution / "report.md").write_text("\n".join(lines)+"\n", encoding="utf-8")
        return report, execution / "report.json"


def main():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="mode", required=True)
    init = sub.add_parser("init")
    init.add_argument("--storage-root", type=Path, required=True)
    init.add_argument("--username", default="workstation")
    init.add_argument("--project-name", default="Prompt AIDD")
    for name in ("plan", "run", "ask"):
        child = sub.add_parser(name)
        child.add_argument("--db", type=Path, required=True)
        child.add_argument("--user", required=True)
        child.add_argument("--project", required=True)
        if name in {"plan", "ask"}:
            text = child.add_mutually_exclusive_group(required=True)
            text.add_argument("--prompt"); text.add_argument("--prompt-file", type=Path)
            model = child.add_mutually_exclusive_group(required=True)
            model.add_argument("--llm-profile", type=Path); model.add_argument("--response-file", type=Path)
        if name == "run":
            child.add_argument("--plan", type=Path, required=True)
        if name in {"run", "ask"}:
            child.add_argument("--runtime", type=Path)
            child.add_argument("--allow-compute", action="store_true")
    args = p.parse_args()
    if args.mode == "init":
        print(json.dumps(initialize_context(args.storage_root, args.username, args.project_name), ensure_ascii=False))
        return 0
    if args.mode in {"plan", "ask"}:
        prompt = args.prompt if args.prompt is not None else args.prompt_file.read_text(encoding="utf-8")
        path = create_plan(args.db, args.user, args.project, prompt, llm_profile=args.llm_profile, response_file=args.response_file)
        print(json.dumps(dict(status="planned", plan=str(path)), ensure_ascii=False))
        if args.mode == "plan": return 0
        args.plan = path
    report, path = run_plan(args.db, args.user, args.project, args.plan, runtime=args.runtime, allow_compute=args.allow_compute)
    print(json.dumps(dict(status=report["status"], report=str(path)), ensure_ascii=False))
    return 0 if report["status"] == "complete" else 2


def initialize_context(storage_root, username="workstation", project_name="Prompt AIDD"):
    from .registry import initialize
    from .context import create_user, slugify
    from .project_context import create_scientific_project, activate_project
    storage_root = Path(storage_root).resolve()
    database = storage_root / "registry/aidd.sqlite3"
    database.parent.mkdir(parents=True, exist_ok=True)
    initialize(database)
    with connect(database) as connection:
        name = slugify(username)
        existing = connection.execute("SELECT id FROM app_user WHERE username=?", (name,)).fetchone()
        user = existing["id"] if existing else create_user(connection, name)
        matches = connection.execute("SELECT id FROM project WHERE owner_user_id=? AND name=?", (user, project_name)).fetchall()
        if len(matches) > 1: raise ValueError("Ambiguous existing Project name")
        project = matches[0]["id"] if matches else create_scientific_project(connection, user, project_name,
                    "Prompt-driven protein, AF3 and calibrated 3D workflows", storage_root)
        context = activate_project(connection, user, project)
    context["db"] = str(database)
    target = storage_root / "users" / name / "prompt-context.json"
    target.parent.mkdir(parents=True, exist_ok=True)
    _atomic_json(target, context)
    return dict(**context, context_file=str(target))


if __name__ == "__main__":
    raise SystemExit(main())
