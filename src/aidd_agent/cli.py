from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
import time

from .importer import import_library
from .registry import connect, initialize, library_summary
from .campaign import (
    campaign_status, create_campaign, register_structure_candidates,
    select_structure, set_target,
)
from .rcsb import search_structures
from .acquisition import acquire_rcsb_mmcif
from .doctor import environment_report
from .prediction import load_model_profile, prediction_command, write_alphafold3_input
from .context import create_user
from .project_context import (
    activate_project, active_project, create_scientific_project,
    deactivate_project, initialize_storage_root, list_owned_projects, project_summary,
)
from .provenance import (
    create_workflow_run, import_result_manifest, record_project_decision, run_summary,
)
from .cluster import (
    export_slurm_bundle, import_submission_receipt, register_cluster,
    set_user_cluster_identity, validate_bundle,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="aidd-agent")
    subparsers = parser.add_subparsers(dest="command", required=True)

    init_parser = subparsers.add_parser("init", help="Initialize the local registry")
    init_parser.add_argument("--db", type=Path, required=True)

    import_parser = subparsers.add_parser("import-mol2", help="Import a batch MOL2 library")
    import_parser.add_argument("--db", type=Path, required=True)
    import_parser.add_argument("--library", required=True)
    import_parser.add_argument("--input", type=Path, required=True)
    import_parser.add_argument("--no-recursive", action="store_true")
    import_parser.add_argument("--dry-run", action="store_true")
    import_parser.add_argument(
        "--progress-every", type=int, default=10_000,
        help="Print progress after this many records; use 0 to disable",
    )

    summary_parser = subparsers.add_parser("summary", help="Summarize a registered library")
    summary_parser.add_argument("--db", type=Path, required=True)
    summary_parser.add_argument("--library", required=True)

    create = subparsers.add_parser("create-campaign", help="Create an AIDD campaign")
    create.add_argument("--db", type=Path, required=True)
    create.add_argument("--user", required=True)
    create.add_argument("--project", required=True)
    create.add_argument("--name", required=True)
    create.add_argument("--objective", required=True)

    target = subparsers.add_parser("set-target", help="Bind a biological target")
    target.add_argument("--db", type=Path, required=True)
    target.add_argument("--user", required=True)
    target.add_argument("--campaign", required=True)
    target.add_argument("--name", required=True)
    target.add_argument("--organism", required=True)
    target.add_argument("--uniprot")
    target.add_argument("--gene")
    target.add_argument("--notes", default="")

    search = subparsers.add_parser("search-pdb", help="Search and register RCSB candidates")
    search.add_argument("--db", type=Path, required=True)
    search.add_argument("--user", required=True)
    search.add_argument("--campaign", required=True)
    search.add_argument("--uniprot", required=True)
    search.add_argument("--max-resolution", type=float, default=3.0)
    search.add_argument("--method", default="X-RAY DIFFRACTION")

    choose = subparsers.add_parser("select-pdb", help="Freeze one PDB candidate")
    choose.add_argument("--db", type=Path, required=True)
    choose.add_argument("--user", required=True)
    choose.add_argument("--campaign", required=True)
    choose.add_argument("--pdb", required=True)
    choose.add_argument("--rationale", required=True)

    status = subparsers.add_parser("campaign-status", help="Show campaign state")
    status.add_argument("--db", type=Path, required=True)
    status.add_argument("--user", required=True)
    status.add_argument("--campaign", required=True)

    user = subparsers.add_parser("create-user", help="Create a platform user")
    user.add_argument("--db", type=Path, required=True)
    user.add_argument("--username", required=True)
    user.add_argument("--display-name")

    project = subparsers.add_parser("create-project", help="Create a user-owned project")
    project.add_argument("--db", type=Path, required=True)
    project.add_argument("--user", required=True)
    project.add_argument("--name", required=True)
    project.add_argument("--objective", required=True)
    project.add_argument("--description", default="")
    project.add_argument("--storage-root", type=Path, required=True)

    init_storage = subparsers.add_parser("init-storage", help="Initialize the fixed AIDD storage root")
    init_storage.add_argument("--storage-root", type=Path, required=True)

    projects = subparsers.add_parser("list-projects", help="List Projects owned by a user")
    projects.add_argument("--db", type=Path, required=True)
    projects.add_argument("--user", required=True)

    activate_project_parser = subparsers.add_parser("activate-project", help="Set a user's active Project")
    activate_project_parser.add_argument("--db", type=Path, required=True)
    activate_project_parser.add_argument("--user", required=True)
    activate_project_parser.add_argument("--project", required=True)

    active_project_parser = subparsers.add_parser("active-project", help="Show a user's active Project")
    active_project_parser.add_argument("--db", type=Path, required=True)
    active_project_parser.add_argument("--user", required=True)

    deactivate_project_parser = subparsers.add_parser("deactivate-project", help="Clear a user's active Project")
    deactivate_project_parser.add_argument("--db", type=Path, required=True)
    deactivate_project_parser.add_argument("--user", required=True)

    project_status = subparsers.add_parser("project-status", help="Show Project methods and Runs")
    project_status.add_argument("--db", type=Path, required=True)
    project_status.add_argument("--user", required=True)
    project_status.add_argument("--project", required=True)

    create_run_parser = subparsers.add_parser("create-run", help="Register an executed workflow method")
    create_run_parser.add_argument("--db", type=Path, required=True)
    create_run_parser.add_argument("--user", required=True)
    create_run_parser.add_argument("--project", required=True)
    create_run_parser.add_argument("--campaign")
    create_run_parser.add_argument("--run-type", required=True)
    create_run_parser.add_argument("--tool", required=True)
    create_run_parser.add_argument("--tool-version", required=True)
    create_run_parser.add_argument("--backend", choices=("local", "slurm"), required=True)
    create_run_parser.add_argument("--parameters-json", type=Path, required=True)
    create_run_parser.add_argument("--input-artifact", action="append", default=[])

    import_manifest_parser = subparsers.add_parser("import-result-manifest", help="Verify and register Run outputs")
    import_manifest_parser.add_argument("--db", type=Path, required=True)
    import_manifest_parser.add_argument("--user", required=True)
    import_manifest_parser.add_argument("--manifest", type=Path, required=True)

    run_status_parser = subparsers.add_parser("run-status", help="Show Run provenance and results")
    run_status_parser.add_argument("--db", type=Path, required=True)
    run_status_parser.add_argument("--user", required=True)
    run_status_parser.add_argument("--run", required=True)

    decision_parser = subparsers.add_parser("record-decision", help="Record a meaningful Project decision")
    decision_parser.add_argument("--db", type=Path, required=True)
    decision_parser.add_argument("--user", required=True)
    decision_parser.add_argument("--project", required=True)
    decision_parser.add_argument("--campaign")
    decision_parser.add_argument("--decision-type", required=True)
    decision_parser.add_argument("--subject", required=True)
    decision_parser.add_argument("--rationale", required=True)

    cluster_parser = subparsers.add_parser("register-cluster", help="Register a credential-free Cluster Profile")
    cluster_parser.add_argument("--db", type=Path, required=True)
    cluster_parser.add_argument("--profile-json", type=Path, required=True)

    identity_parser = subparsers.add_parser("link-cluster-user", help="Link a user to a cluster username")
    identity_parser.add_argument("--db", type=Path, required=True)
    identity_parser.add_argument("--user", required=True)
    identity_parser.add_argument("--cluster", required=True)
    identity_parser.add_argument("--cluster-username", required=True)
    identity_parser.add_argument("--transport", choices=("auto", "openssh", "putty"), default="auto")
    identity_parser.add_argument("--putty-session")

    bundle_parser = subparsers.add_parser("export-slurm-bundle", help="Generate a validated Slurm Bundle")
    bundle_parser.add_argument("--db", type=Path, required=True)
    bundle_parser.add_argument("--user", required=True)
    bundle_parser.add_argument("--run", required=True)
    bundle_parser.add_argument("--cluster", required=True)
    bundle_parser.add_argument("--slurm-spec-json", type=Path, required=True)
    bundle_parser.add_argument("--output", type=Path, required=True)

    validate_parser = subparsers.add_parser("validate-bundle", help="Verify a Slurm Bundle without submission")
    validate_parser.add_argument("--bundle", type=Path, required=True)

    receipt_parser = subparsers.add_parser("import-submission-receipt", help="Register an authenticated Slurm submission")
    receipt_parser.add_argument("--db", type=Path, required=True)
    receipt_parser.add_argument("--user", required=True)
    receipt_parser.add_argument("--receipt", type=Path, required=True)

    fetch_parser = subparsers.add_parser("fetch-pdb", help="Download an RCSB mmCIF into a Project")
    fetch_parser.add_argument("--pdb", required=True)
    fetch_parser.add_argument("--project-root", type=Path, required=True)

    af_parser = subparsers.add_parser("prepare-alphafold3", help="Create a validated AlphaFold 3 job input")
    af_parser.add_argument("--name", required=True)
    af_parser.add_argument("--sequence", action="append", required=True)
    af_parser.add_argument("--seed", action="append", type=int, default=[])
    af_parser.add_argument("--project-root", type=Path, required=True)
    af_parser.add_argument("--output", type=Path, required=True)
    af_parser.add_argument("--profile", type=Path)

    subparsers.add_parser("doctor", help="Report workstation dependencies")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "init":
        initialize(args.db)
        print(json.dumps({"database": str(args.db), "status": "initialized"}, indent=2))
        return 0
    if args.command == "import-mol2":
        started = time.monotonic()
        print(
            f"Starting {'validation' if args.dry_run else 'import'}: {args.input}",
            file=sys.stderr,
            flush=True,
        )

        def show_progress(report, path: Path) -> None:
            elapsed = time.monotonic() - started
            rate = report.records_seen / elapsed if elapsed else 0.0
            print(
                f"Processed {report.records_seen:,} records "
                f"({rate:,.0f} records/s), current file: {path.name}",
                file=sys.stderr,
                flush=True,
            )

        report = import_library(
            args.db,
            args.library,
            args.input,
            recursive=not args.no_recursive,
            dry_run=args.dry_run,
            progress_every=max(0, args.progress_every),
            progress=show_progress,
        )
        elapsed = time.monotonic() - started
        print(f"Completed in {elapsed:,.1f} seconds", file=sys.stderr, flush=True)
        print(json.dumps(report.to_dict(), indent=2, ensure_ascii=False))
        return 1 if report.invalid else 0
    if args.command == "summary":
        with connect(args.db) as connection:
            result = library_summary(connection, args.library)
        print(json.dumps(result, indent=2, ensure_ascii=False))
        return 0
    if args.command == "create-campaign":
        initialize(args.db)
        with connect(args.db) as connection:
            campaign_id = create_campaign(connection, args.user, args.project,
                                          args.name, args.objective)
        print(json.dumps({"campaign_id": campaign_id, "state": "draft"}, indent=2))
        return 0
    if args.command == "set-target":
        initialize(args.db)
        with connect(args.db) as connection:
            target_id = set_target(connection, args.user, args.campaign, name=args.name,
                                   organism=args.organism, uniprot_id=args.uniprot,
                                   gene_name=args.gene, notes=args.notes)
        print(json.dumps({"target_id": target_id, "state": "target_set"}, indent=2))
        return 0
    if args.command == "search-pdb":
        query, candidates = search_structures(args.uniprot,
                                              max_resolution=args.max_resolution,
                                              method=args.method)
        with connect(args.db) as connection:
            inserted = register_structure_candidates(connection, args.user, args.campaign,
                                                      candidates, query)
        print(json.dumps({"inserted": inserted, "candidates": candidates},
                         indent=2, ensure_ascii=False))
        return 0
    if args.command == "select-pdb":
        with connect(args.db) as connection:
            select_structure(connection, args.user, args.campaign, args.pdb, args.rationale)
        print(json.dumps({"pdb_id": args.pdb.upper(), "state": "structure_selected"}, indent=2))
        return 0
    if args.command == "campaign-status":
        with connect(args.db) as connection:
            result = campaign_status(connection, args.user, args.campaign)
        print(json.dumps(result, indent=2, ensure_ascii=False))
        return 0
    if args.command == "create-user":
        initialize(args.db)
        with connect(args.db) as connection:
            user_id = create_user(connection, args.username, args.display_name)
        print(json.dumps({"user_id": user_id, "username": args.username}, indent=2))
        return 0
    if args.command == "create-project":
        initialize(args.db)
        with connect(args.db) as connection:
            project_id = create_scientific_project(
                connection, args.user, args.name, args.objective,
                args.storage_root, args.description,
            )
        print(json.dumps({"project_id": project_id, "owner_user_id": args.user}, indent=2))
        return 0
    if args.command == "init-storage":
        print(json.dumps(initialize_storage_root(args.storage_root), indent=2))
        return 0
    if args.command == "list-projects":
        with connect(args.db) as connection:
            result = list_owned_projects(connection, args.user)
        print(json.dumps({"projects": result}, indent=2, ensure_ascii=False))
        return 0
    if args.command == "activate-project":
        with connect(args.db) as connection:
            result = activate_project(connection, args.user, args.project)
        print(json.dumps(result, indent=2, ensure_ascii=False))
        return 0
    if args.command == "active-project":
        with connect(args.db) as connection:
            result = active_project(connection, args.user)
        print(json.dumps({"active_project": result}, indent=2, ensure_ascii=False))
        return 0
    if args.command == "deactivate-project":
        with connect(args.db) as connection:
            deactivate_project(connection, args.user)
        print(json.dumps({"user_id": args.user, "active_project": None}, indent=2))
        return 0
    if args.command == "project-status":
        with connect(args.db) as connection:
            result = project_summary(connection, args.user, args.project)
        print(json.dumps(result, indent=2, ensure_ascii=False))
        return 0
    if args.command == "create-run":
        parameters = json.loads(args.parameters_json.read_text(encoding="utf-8"))
        with connect(args.db) as connection:
            run_id = create_workflow_run(
                connection, args.user, args.project, run_type=args.run_type,
                tool_name=args.tool, tool_version=args.tool_version,
                parameters=parameters, execution_backend=args.backend,
                campaign_id=args.campaign, input_artifact_ids=args.input_artifact,
            )
        print(json.dumps({"run_id": run_id, "status": "created"}, indent=2))
        return 0
    if args.command == "import-result-manifest":
        with connect(args.db) as connection:
            result = import_result_manifest(connection, args.user, args.manifest)
        print(json.dumps(result, indent=2, ensure_ascii=False))
        return 0
    if args.command == "run-status":
        with connect(args.db) as connection:
            result = run_summary(connection, args.user, args.run)
        print(json.dumps(result, indent=2, ensure_ascii=False))
        return 0
    if args.command == "record-decision":
        with connect(args.db) as connection:
            decision_id = record_project_decision(
                connection, args.user, args.project,
                decision_type=args.decision_type, subject=args.subject,
                rationale=args.rationale, campaign_id=args.campaign,
            )
        print(json.dumps({"decision_id": decision_id}, indent=2))
        return 0
    if args.command == "register-cluster":
        initialize(args.db)
        profile = json.loads(args.profile_json.read_text(encoding="utf-8"))
        with connect(args.db) as connection:
            cluster_id = register_cluster(connection, **profile)
        print(json.dumps({"cluster_id": cluster_id, "name": profile["name"]}, indent=2))
        return 0
    if args.command == "link-cluster-user":
        with connect(args.db) as connection:
            set_user_cluster_identity(
                connection, args.user, args.cluster, args.cluster_username,
                args.transport, args.putty_session,
            )
        print(json.dumps({"user_id": args.user, "cluster_id": args.cluster,
                          "cluster_username": args.cluster_username}, indent=2))
        return 0
    if args.command == "export-slurm-bundle":
        spec = json.loads(args.slurm_spec_json.read_text(encoding="utf-8"))
        with connect(args.db) as connection:
            result = export_slurm_bundle(
                connection, args.user, args.run, args.cluster, args.output, spec,
            )
        print(json.dumps(result, indent=2, ensure_ascii=False))
        return 0
    if args.command == "validate-bundle":
        print(json.dumps(validate_bundle(args.bundle), indent=2, ensure_ascii=False))
        return 0
    if args.command == "import-submission-receipt":
        with connect(args.db) as connection:
            result = import_submission_receipt(connection, args.user, args.receipt)
        print(json.dumps(result, indent=2, ensure_ascii=False))
        return 0
    if args.command == "fetch-pdb":
        print(json.dumps(acquire_rcsb_mmcif(args.pdb, args.project_root), indent=2))
        return 0
    if args.command == "prepare-alphafold3":
        path = write_alphafold3_input(
            args.name, args.sequence, args.output, args.project_root,
            seeds=args.seed or (1,),
        )
        result = {"input": str(path)}
        if args.profile:
            profile = load_model_profile(args.profile)
            result["command"] = prediction_command(profile, path, path.parent / "output")
        print(json.dumps(result, indent=2))
        return 0
    if args.command == "doctor":
        report = environment_report()
        print(json.dumps(report, indent=2))
        return 0 if report["ready"] else 1
    raise AssertionError("unreachable")


if __name__ == "__main__":
    raise SystemExit(main())
