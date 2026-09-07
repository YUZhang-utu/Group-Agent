from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
import time

from .importer import import_library
from .registry import connect, initialize, library_summary
from .campaign import (
    campaign_status, create_campaign, register_predicted_structure,
    register_structure_candidates, select_receptor_candidate, select_structure, set_target,
)
from .rcsb import search_structures
from .acquisition import acquire_alphafold_db_mmcif, acquire_rcsb_mmcif
from .structure_compare import (
    comparison_set_status, create_campaign_comparison_set,
    create_receptor_ensemble, generate_comparison_pymol_review,
    generate_receptor_ensemble_pymol_review, receptor_ensemble_status,
)
from .ai_recommendation import (
    ai_review_status, create_ai_review_request,
    create_ligand_query_review_request, create_receptor_eligibility_review_request,
    import_ai_recommendation,
)
from .doctor import environment_report
from .prediction import (
    inspect_alphafold3_output, load_model_profile, prediction_command, write_alphafold3_input,
    write_prediction_input,
)
from .ligand_workflow import (
    campaign_ligand_status, compare_campaign_ligands, register_campaign_ligands,
    run_hierarchical_library_search, select_query_ligand,
)
from .chemistry_prep import build_library_indices, prepare_campaign_ligands
from .anchor_extraction import enrich_manifest_with_reduce
from .pharmacophore_index import (
    build_pharmacophore_index, compile_pharmacophore_query,
    load_external_l1_ids, search_pharmacophore_index,
)
from .pharmacophore_validation import run_pharmacophore_validation
from .gaussian_batch import (
    prepare_gaussian_query, run_staged_gaussian_reranking,
    score_gaussian_candidates,
)
from .reduce_pocket import (
    build_reduce_het_dictionary, export_query_pocket_pdb, run_reduce, validate_reduce_run,
)
from .context import (
    activate_task, active_task, create_task, create_user, deactivate_task, list_tasks,
)
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
    create.add_argument("--task", required=True)
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

    import_prediction = subparsers.add_parser(
        "import-alphafold3-result",
        help="Inspect and register a completed AlphaFold 3 structure candidate")
    import_prediction.add_argument("--db", type=Path, required=True)
    import_prediction.add_argument("--user", required=True)
    import_prediction.add_argument("--campaign", required=True)
    import_prediction.add_argument("--output-dir", type=Path, required=True)
    import_prediction.add_argument("--input", type=Path)
    import_prediction.add_argument("--construct", required=True)
    import_prediction.add_argument("--model-version", required=True)
    import_prediction.add_argument("--chain", action="append", default=[])
    import_prediction.add_argument("--runtime-json", type=Path)
    import_prediction.add_argument("--manifest-output", type=Path)

    select_receptor = subparsers.add_parser(
        "select-receptor", help="Freeze an experimental or predicted receptor candidate")
    select_receptor.add_argument("--db", type=Path, required=True)
    select_receptor.add_argument("--user", required=True)
    select_receptor.add_argument("--campaign", required=True)
    select_receptor.add_argument("--candidate", required=True)
    select_receptor.add_argument("--rationale", required=True)

    status = subparsers.add_parser("campaign-status", help="Show campaign state")
    status.add_argument("--db", type=Path, required=True)
    status.add_argument("--user", required=True)
    status.add_argument("--campaign", required=True)

    comparison = subparsers.add_parser(
        "create-structure-comparison", help="Create a Campaign PDB comparison set")
    comparison.add_argument("--db", type=Path, required=True)
    comparison.add_argument("--user", required=True)
    comparison.add_argument("--campaign", required=True)
    comparison.add_argument("--name", required=True)
    comparison.add_argument("--pdb", action="append", required=True)
    comparison.add_argument("--reference-pdb", required=True)
    comparison.add_argument("--pocket-residues", required=True,
                            help="Comma-separated residue numbers")
    comparison.add_argument("--chain", action="append", default=[],
                            help="PDB=chain; repeat for structure-specific chains")
    comparison.add_argument("--rationale", required=True)

    ensemble = subparsers.add_parser(
        "create-receptor-ensemble",
        help="Create an all-candidate experimental/predicted receptor ensemble")
    ensemble.add_argument("--db", type=Path, required=True)
    ensemble.add_argument("--user", required=True)
    ensemble.add_argument("--campaign", required=True)
    ensemble.add_argument("--name", required=True)
    ensemble.add_argument("--reference-pdb", required=True)
    ensemble.add_argument("--pocket-residues", default="")
    ensemble.add_argument("--chain", action="append", default=[])
    ensemble.add_argument(
        "--exclude-pdb", action="append", default=[],
        help="PDB=reason; retain provenance but omit an ineligible receptor")
    ensemble.add_argument("--rationale", required=True)

    ensemble_status = subparsers.add_parser(
        "receptor-ensemble-status", help="Show a unified receptor ensemble")
    ensemble_status.add_argument("--db", type=Path, required=True)
    ensemble_status.add_argument("--user", required=True)
    ensemble_status.add_argument("--ensemble", required=True)

    ensemble_pymol = subparsers.add_parser(
        "prepare-receptor-ensemble-pymol",
        help="Generate all-structure and ligand-pocket PyMOL review")
    ensemble_pymol.add_argument("--db", type=Path, required=True)
    ensemble_pymol.add_argument("--user", required=True)
    ensemble_pymol.add_argument("--ensemble", required=True)
    ensemble_pymol.add_argument("--project-root", type=Path, required=True)
    ensemble_pymol.add_argument("--output", type=Path, required=True)

    comparison_status_parser = subparsers.add_parser(
        "structure-comparison-status", help="Show a structure comparison set")
    comparison_status_parser.add_argument("--db", type=Path, required=True)
    comparison_status_parser.add_argument("--user", required=True)
    comparison_status_parser.add_argument("--comparison", required=True)

    review = subparsers.add_parser(
        "prepare-pymol-review", help="Generate a PyMOL script for a comparison set")
    review.add_argument("--db", type=Path, required=True)
    review.add_argument("--user", required=True)
    review.add_argument("--comparison", required=True)
    review.add_argument("--project-root", type=Path, required=True)
    review.add_argument("--output", type=Path)

    ai_request = subparsers.add_parser(
        "create-ai-review", help="Create an audited model-facing recommendation request")
    ai_request.add_argument("--db", type=Path, required=True)
    ai_request.add_argument("--user", required=True)
    ai_request.add_argument("--project", required=True)
    ai_request.add_argument("--campaign")
    ai_request.add_argument("--task-type", required=True)
    ai_request.add_argument("--subject-type", required=True)
    ai_request.add_argument("--subject-id", required=True)
    ai_request.add_argument("--prompt-version", required=True)
    ai_request.add_argument("--prompt-file", type=Path, required=True)
    ai_request.add_argument("--evidence-json", type=Path, required=True)
    ai_request.add_argument("--data-class", action="append", required=True)

    receptor_ai = subparsers.add_parser(
        "create-receptor-eligibility-review",
        help="Create a provider-neutral LLM review of all Campaign receptor candidates")
    receptor_ai.add_argument("--db", type=Path, required=True)
    receptor_ai.add_argument("--user", required=True)
    receptor_ai.add_argument("--campaign", required=True)
    receptor_ai.add_argument("--requirements-json", type=Path, required=True)
    receptor_ai.add_argument("--computed-evidence-json", type=Path)
    receptor_ai.add_argument("--prompt-version", default="receptor-eligibility-v1")

    ligand_ai = subparsers.add_parser(
        "create-ligand-query-review",
        help="Create a provider-neutral LLM review of registered query ligands")
    ligand_ai.add_argument("--db", type=Path, required=True)
    ligand_ai.add_argument("--user", required=True)
    ligand_ai.add_argument("--campaign", required=True)
    ligand_ai.add_argument("--requirements-json", type=Path, required=True)
    ligand_ai.add_argument("--comparison-json", type=Path)
    ligand_ai.add_argument("--prompt-version", default="ligand-query-selection-v1")

    ai_import = subparsers.add_parser(
        "import-ai-recommendation", help="Import a structured model recommendation")
    ai_import.add_argument("--db", type=Path, required=True)
    ai_import.add_argument("--user", required=True)
    ai_import.add_argument("--request", required=True)
    ai_import.add_argument("--provider", required=True)
    ai_import.add_argument("--model", required=True)
    ai_import.add_argument("--model-version")
    ai_import.add_argument("--response-json", type=Path, required=True)

    ai_status_parser = subparsers.add_parser(
        "ai-review-status", help="Show an AI review request and recommendation")
    ai_status_parser.add_argument("--db", type=Path, required=True)
    ai_status_parser.add_argument("--user", required=True)
    ai_status_parser.add_argument("--request", required=True)

    ligand_import = subparsers.add_parser(
        "register-campaign-ligands", help="Register filtered co-crystal ligand instances")
    ligand_import.add_argument("--db", type=Path, required=True)
    ligand_import.add_argument("--user", required=True)
    ligand_import.add_argument("--campaign", required=True)
    ligand_import.add_argument("--manifest", type=Path, required=True)

    ligand_prepare = subparsers.add_parser(
        "prepare-campaign-ligands", help="Extract retained mmCIF ligands using authoritative RCSB CCD bonds")
    ligand_prepare.add_argument("--db", type=Path, required=True)
    ligand_prepare.add_argument("--user", required=True)
    ligand_prepare.add_argument("--campaign", required=True)
    ligand_prepare.add_argument("--structures-dir", type=Path, required=True)
    ligand_prepare.add_argument("--output-dir", type=Path, required=True)
    ligand_prepare.add_argument("--ccd-cache", type=Path, required=True)

    index_build = subparsers.add_parser(
        "build-library-indices", help="Build production Morgan and conformer-level USRCAT indices")
    index_build.add_argument("--db", type=Path, required=True)
    index_build.add_argument("--library", required=True)
    index_build.add_argument("--output-dir", type=Path, required=True)

    reduce_export = subparsers.add_parser(
        "export-reduce-pocket", help="Export a query ligand and complete 5A pocket residues as PDB")
    reduce_export.add_argument("--mmcif", type=Path, required=True)
    reduce_export.add_argument("--output", type=Path, required=True)
    reduce_export.add_argument("--ccd-id", required=True)
    reduce_export.add_argument("--chain", required=True)
    reduce_export.add_argument("--residue", required=True)
    reduce_export.add_argument("--radius", type=float, default=5.0)

    reduce_run = subparsers.add_parser(
        "run-reduce-pocket", help="Hydrogenate and optimize a query pocket using Linux Reduce")
    reduce_run.add_argument("--profile", type=Path, required=True)
    reduce_run.add_argument("--pocket-pdb", type=Path, required=True)
    reduce_run.add_argument("--output-dir", type=Path, required=True)

    reduce_validate = subparsers.add_parser(
        "validate-reduce-pocket", help="Verify ligand/pocket hydrogenation and Reduce warnings")
    reduce_validate.add_argument("--output-dir", type=Path, required=True)
    reduce_validate.add_argument("--ccd-id", required=True)
    reduce_validate.add_argument("--chain", required=True)
    reduce_validate.add_argument("--residue", required=True)

    reduce_dictionary = subparsers.add_parser(
        "build-reduce-het-dictionary", help="Build a query-specific Reduce HET dictionary from CCD")
    reduce_dictionary.add_argument("--ccd", type=Path, required=True)
    reduce_dictionary.add_argument("--ccd-id", required=True)
    reduce_dictionary.add_argument("--output", type=Path, required=True)

    reduce_enrich = subparsers.add_parser(
        "enrich-query-anchors", help="Add Reduce angles and projection points to a query manifest")
    reduce_enrich.add_argument("--query-manifest", type=Path, required=True)
    reduce_enrich.add_argument("--ccd", type=Path, required=True)
    reduce_enrich.add_argument("--hydrogenated-pdb", type=Path, required=True)
    reduce_enrich.add_argument("--output", type=Path, required=True)

    pharmacophore_build = subparsers.add_parser(
        "build-pharmacophore-index",
        help="Build reusable shard-local 3D pharmacophore pair postings from artifacts")
    pharmacophore_build.add_argument("--artifact-catalog", type=Path, required=True)
    pharmacophore_build.add_argument("--output-dir", type=Path, required=True)
    pharmacophore_build.add_argument("--bin-width", type=float, default=0.5)
    pharmacophore_build.add_argument("--max-distance", type=float, default=20.0)

    pharmacophore_compile = subparsers.add_parser(
        "compile-pharmacophore-query",
        help="Compile one co-crystal anchor manifest without reading the library")
    pharmacophore_compile.add_argument("--query-manifest", type=Path, required=True)
    pharmacophore_compile.add_argument("--output", type=Path, required=True)

    pharmacophore_search = subparsers.add_parser(
        "search-pharmacophore-index",
        help="Run loose/balanced/strict partial 3D motif retrieval and union external L1 IDs")
    pharmacophore_search.add_argument("--index-catalog", type=Path, required=True)
    pharmacophore_search.add_argument("--query-plan", type=Path, required=True)
    pharmacophore_search.add_argument("--external-l1-ids", type=Path)
    pharmacophore_search.add_argument("--output", type=Path, required=True)

    pharmacophore_validate = subparsers.add_parser(
        "validate-pharmacophore-retrieval",
        help="Build/reuse, query, search, and validate the complete pre-docking admission stage")
    pharmacophore_validate.add_argument("--artifact-catalog", type=Path, required=True)
    pharmacophore_validate.add_argument("--pharmacophore-index-dir", type=Path, required=True)
    pharmacophore_validate.add_argument("--query-manifest", type=Path, required=True)
    pharmacophore_validate.add_argument("--output-dir", type=Path, required=True)
    pharmacophore_validate.add_argument("--external-l1-ids", type=Path)
    pharmacophore_validate.add_argument("--faiss-index-dir", type=Path)
    pharmacophore_validate.add_argument("--mmcif", type=Path)
    pharmacophore_validate.add_argument("--ccd", type=Path)
    pharmacophore_validate.add_argument("--ccd-id")
    pharmacophore_validate.add_argument("--faiss-k", type=int, default=100_000)
    pharmacophore_validate.add_argument("--nprobe", type=int, default=256)
    pharmacophore_validate.add_argument("--bin-width", type=float, default=0.5)
    pharmacophore_validate.add_argument("--max-distance", type=float, default=20.0)

    gaussian_query = subparsers.add_parser(
        "prepare-gaussian-query",
        help="Package one locked co-crystal ligand for artifact Gaussian reranking")
    gaussian_query.add_argument("--mmcif", type=Path, required=True)
    gaussian_query.add_argument("--ccd", type=Path, required=True)
    gaussian_query.add_argument("--query-manifest", type=Path, required=True)
    gaussian_query.add_argument("--output", type=Path, required=True)

    gaussian_score = subparsers.add_parser(
        "score-gaussian-candidates",
        help="Rerank stable artifact IDs without changing L1 admission")
    gaussian_score.add_argument("--artifact-catalog", type=Path, required=True)
    gaussian_score.add_argument("--query", type=Path, required=True)
    gaussian_score.add_argument("--candidate-ids", type=Path, required=True)
    gaussian_score.add_argument("--output", type=Path, required=True)
    gaussian_score.add_argument("--sigma", type=float, default=1.0)
    gaussian_score.add_argument("--cutoff", type=float, default=4.5)
    gaussian_score.add_argument("--pair-tolerance", type=float, default=2.0)
    gaussian_score.add_argument("--axial-samples", type=int, default=6)
    gaussian_score.add_argument("--max-pair-seeds", type=int, default=512)

    gaussian_staged = subparsers.add_parser(
        "run-staged-gaussian-reranking",
        help="Run resumable parallel PCA coarse scoring and Top-N pair refinement")
    gaussian_staged.add_argument("--artifact-catalog", type=Path, required=True)
    gaussian_staged.add_argument("--query", type=Path, required=True)
    gaussian_staged.add_argument("--candidate-ids", type=Path, required=True)
    gaussian_staged.add_argument("--output-dir", type=Path, required=True)
    gaussian_staged.add_argument("--stage", choices=("coarse", "refine", "all"),
                                 default="all")
    gaussian_staged.add_argument("--workers", type=int, default=16)
    gaussian_staged.add_argument("--chunk-size", type=int, default=1000)
    gaussian_staged.add_argument("--top-n-per-objective", type=int, default=5000)
    gaussian_staged.add_argument("--sigma", type=float, default=1.0)
    gaussian_staged.add_argument("--cutoff", type=float, default=4.5)
    gaussian_staged.add_argument("--pair-tolerance", type=float, default=2.0)
    gaussian_staged.add_argument("--axial-samples", type=int, default=6)
    gaussian_staged.add_argument("--max-pair-seeds", type=int, default=512)
    gaussian_staged.add_argument("--progress-every", type=int, default=1)
    gaussian_staged.add_argument("--no-resume", action="store_true")

    ligand_status = subparsers.add_parser(
        "campaign-ligands", help="List Campaign ligands and the immutable query lock")
    ligand_status.add_argument("--db", type=Path, required=True)
    ligand_status.add_argument("--user", required=True)
    ligand_status.add_argument("--campaign", required=True)

    ligand_compare = subparsers.add_parser(
        "compare-campaign-ligands", help="Compare Campaign ligands in 2D and available 3D")
    ligand_compare.add_argument("--db", type=Path, required=True)
    ligand_compare.add_argument("--user", required=True)
    ligand_compare.add_argument("--campaign", required=True)

    ligand_select = subparsers.add_parser(
        "select-query-ligand", help="Irreversibly lock the Campaign similarity query ligand")
    ligand_select.add_argument("--db", type=Path, required=True)
    ligand_select.add_argument("--user", required=True)
    ligand_select.add_argument("--campaign", required=True)
    ligand_select.add_argument("--ligand", required=True)
    ligand_select.add_argument("--rationale", required=True)

    ligand_search = subparsers.add_parser(
        "search-similar-ligands", help="Run hierarchical Morgan/USRCAT library search")
    ligand_search.add_argument("--db", type=Path, required=True)
    ligand_search.add_argument("--user", required=True)
    ligand_search.add_argument("--campaign", required=True)
    ligand_search.add_argument("--library", required=True)
    ligand_search.add_argument("--morgan-index", type=Path, required=True)
    ligand_search.add_argument("--usrcat-index", type=Path)
    ligand_search.add_argument("--conformer-map-json", type=Path)
    ligand_search.add_argument("--two-d-pool", type=int, default=1000)
    ligand_search.add_argument("--limit", type=int, default=100)

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

    create_task_parser = subparsers.add_parser("create-task", help="Create a Project-scoped task")
    create_task_parser.add_argument("--db", type=Path, required=True)
    create_task_parser.add_argument("--user", required=True)
    create_task_parser.add_argument("--project", required=True)
    create_task_parser.add_argument("--name", required=True)
    create_task_parser.add_argument("--objective", required=True)

    list_tasks_parser = subparsers.add_parser("list-tasks", help="List tasks accessible to a user")
    list_tasks_parser.add_argument("--db", type=Path, required=True)
    list_tasks_parser.add_argument("--user", required=True)

    activate_task_parser = subparsers.add_parser("activate-task", help="Set a user's active task")
    activate_task_parser.add_argument("--db", type=Path, required=True)
    activate_task_parser.add_argument("--user", required=True)
    activate_task_parser.add_argument("--task", required=True)

    active_task_parser = subparsers.add_parser("active-task", help="Show a user's active task")
    active_task_parser.add_argument("--db", type=Path, required=True)
    active_task_parser.add_argument("--user", required=True)

    deactivate_task_parser = subparsers.add_parser("deactivate-task", help="Clear a user's active task")
    deactivate_task_parser.add_argument("--db", type=Path, required=True)
    deactivate_task_parser.add_argument("--user", required=True)

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

    afdb_parser = subparsers.add_parser(
        "fetch-alphafold-db", help="Download an AlphaFold DB model into a Project")
    afdb_parser.add_argument("--uniprot", required=True)
    afdb_parser.add_argument("--project-root", type=Path, required=True)

    af_parser = subparsers.add_parser("prepare-alphafold3", help="Create a validated AlphaFold 3 job input")
    af_parser.add_argument("--name", required=True)
    af_parser.add_argument("--sequence", action="append", required=True)
    af_parser.add_argument("--seed", action="append", type=int, default=[])
    af_parser.add_argument("--project-root", type=Path, required=True)
    af_parser.add_argument("--output", type=Path, required=True)
    af_parser.add_argument("--profile", type=Path)

    prediction_parser = subparsers.add_parser(
        "prepare-structure-prediction", help="Create a backend-native protein prediction input")
    prediction_parser.add_argument(
        "--backend", choices=("alphafold2", "alphafold3", "boltz2", "chai1"), required=True)
    prediction_parser.add_argument("--name", required=True)
    prediction_parser.add_argument("--sequence", action="append", required=True)
    prediction_parser.add_argument("--seed", action="append", type=int, default=[])
    prediction_parser.add_argument("--project-root", type=Path, required=True)
    prediction_parser.add_argument("--output", type=Path, required=True)
    prediction_parser.add_argument("--profile", type=Path, required=True)

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
            campaign_id = create_campaign(connection, args.user, args.project, args.task,
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
    if args.command == "import-alphafold3-result":
        initialize(args.db)
        runtime = (json.loads(args.runtime_json.read_text(encoding="utf-8"))
                   if args.runtime_json else {})
        manifest = inspect_alphafold3_output(
            args.output_dir, construct_name=args.construct,
            model_version=args.model_version, chain_ids=args.chain or ["A"],
            input_path=args.input, runtime=runtime,
        )
        if args.manifest_output:
            args.manifest_output.parent.mkdir(parents=True, exist_ok=True)
            args.manifest_output.write_text(
                json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
        with connect(args.db) as connection:
            candidate_id = register_predicted_structure(
                connection, args.user, args.campaign, manifest)
        print(json.dumps({"candidate_id": candidate_id, "manifest": manifest},
                         indent=2, ensure_ascii=False))
        return 0
    if args.command == "select-receptor":
        initialize(args.db)
        with connect(args.db) as connection:
            select_receptor_candidate(
                connection, args.user, args.campaign, args.candidate, args.rationale)
        print(json.dumps({"candidate_id": args.candidate,
                          "state": "structure_selected"}, indent=2))
        return 0
    if args.command == "campaign-status":
        with connect(args.db) as connection:
            result = campaign_status(connection, args.user, args.campaign)
        print(json.dumps(result, indent=2, ensure_ascii=False))
        return 0
    if args.command == "create-structure-comparison":
        try:
            pocket_residues = [int(value.strip()) for value in args.pocket_residues.split(",")
                               if value.strip()]
        except ValueError as exc:
            raise ValueError("Pocket residues must be comma-separated integers") from exc
        chains = {}
        for specification in args.chain:
            if specification.count("=") != 1:
                raise ValueError("Each chain must use PDB=chain syntax")
            pdb_id, chain_id = specification.split("=", 1)
            chains[pdb_id.upper()] = chain_id
        with connect(args.db) as connection:
            comparison_id = create_campaign_comparison_set(
                connection, args.user, args.campaign, args.name, args.pdb,
                args.reference_pdb, pocket_residues, chains, args.rationale,
            )
        print(json.dumps({"comparison_set_id": comparison_id,
                          "campaign_state": "structures_review"}, indent=2))
        return 0
    if args.command == "create-receptor-ensemble":
        try:
            pocket_residues = [int(value.strip()) for value in args.pocket_residues.split(",")
                               if value.strip()]
        except ValueError as exc:
            raise ValueError("Pocket residues must be comma-separated integers") from exc
        chains = {}
        for specification in args.chain:
            if specification.count("=") != 1:
                raise ValueError("Each chain must use candidate=chain syntax")
            key, chain = specification.split("=", 1)
            chains[key.upper() if len(key) == 4 else key] = chain
        exclusions = {}
        for specification in args.exclude_pdb:
            if specification.count("=") != 1:
                raise ValueError("Each exclusion must use PDB=reason syntax")
            pdb_id, reason = specification.split("=", 1)
            exclusions[pdb_id.upper()] = reason
        with connect(args.db) as connection:
            ensemble_id = create_receptor_ensemble(
                connection, args.user, args.campaign, args.name,
                args.reference_pdb, pocket_residues, chains, args.rationale,
                exclusions=exclusions)
        print(json.dumps({"ensemble_id": ensemble_id,
                          "campaign_state": "structures_review"}, indent=2))
        return 0
    if args.command == "receptor-ensemble-status":
        with connect(args.db) as connection:
            result = receptor_ensemble_status(connection, args.user, args.ensemble)
        print(json.dumps(result, indent=2, ensure_ascii=False))
        return 0
    if args.command == "prepare-receptor-ensemble-pymol":
        with connect(args.db) as connection:
            result = generate_receptor_ensemble_pymol_review(
                connection, args.user, args.ensemble, args.project_root, args.output)
        print(json.dumps(result, indent=2, ensure_ascii=False))
        return 0
    if args.command == "structure-comparison-status":
        with connect(args.db) as connection:
            result = comparison_set_status(connection, args.user, args.comparison)
        print(json.dumps(result, indent=2, ensure_ascii=False))
        return 0
    if args.command == "prepare-pymol-review":
        with connect(args.db) as connection:
            result = generate_comparison_pymol_review(
                connection, args.user, args.comparison, args.project_root, args.output)
        print(json.dumps(result, indent=2, ensure_ascii=False))
        return 0
    if args.command == "create-ai-review":
        evidence = json.loads(args.evidence_json.read_text(encoding="utf-8"))
        prompt = args.prompt_file.read_text(encoding="utf-8")
        with connect(args.db) as connection:
            request_id = create_ai_review_request(
                connection, args.user, args.project, campaign_id=args.campaign,
                task_type=args.task_type, subject_type=args.subject_type,
                subject_id=args.subject_id, prompt_version=args.prompt_version,
                prompt_text=prompt, evidence=evidence, data_classes=args.data_class,
            )
        print(json.dumps({"request_id": request_id, "status": "pending"}, indent=2))
        return 0
    if args.command == "create-receptor-eligibility-review":
        requirements = json.loads(args.requirements_json.read_text(encoding="utf-8"))
        computed = None
        if args.computed_evidence_json:
            computed = json.loads(
                args.computed_evidence_json.read_text(encoding="utf-8"))
        with connect(args.db) as connection:
            request_id = create_receptor_eligibility_review_request(
                connection, args.user, args.campaign, requirements=requirements,
                computed_evidence=computed, prompt_version=args.prompt_version)
            result = ai_review_status(connection, args.user, request_id)
        print(json.dumps(result, indent=2, ensure_ascii=False))
        return 0
    if args.command == "create-ligand-query-review":
        requirements = json.loads(args.requirements_json.read_text(encoding="utf-8"))
        comparison = None
        if args.comparison_json:
            comparison = json.loads(args.comparison_json.read_text(encoding="utf-8"))
        with connect(args.db) as connection:
            request_id = create_ligand_query_review_request(
                connection, args.user, args.campaign,
                selection_requirements=requirements,
                comparison_summary=comparison, prompt_version=args.prompt_version)
            result = ai_review_status(connection, args.user, request_id)
        print(json.dumps(result, indent=2, ensure_ascii=False))
        return 0
    if args.command == "import-ai-recommendation":
        response = json.loads(args.response_json.read_text(encoding="utf-8"))
        with connect(args.db) as connection:
            recommendation_id = import_ai_recommendation(
                connection, args.user, args.request, provider=args.provider,
                model_name=args.model, model_version=args.model_version,
                response=response,
            )
        print(json.dumps({"recommendation_id": recommendation_id,
                          "request_id": args.request}, indent=2))
        return 0
    if args.command == "ai-review-status":
        with connect(args.db) as connection:
            result = ai_review_status(connection, args.user, args.request)
        print(json.dumps(result, indent=2, ensure_ascii=False))
        return 0
    if args.command == "register-campaign-ligands":
        records = json.loads(args.manifest.read_text(encoding="utf-8"))
        if isinstance(records, dict):
            records = records.get("ligands")
        with connect(args.db) as connection:
            ligand_ids = register_campaign_ligands(
                connection, args.user, args.campaign, records)
        print(json.dumps({"inserted": len(ligand_ids), "ligand_ids": ligand_ids}, indent=2))
        return 0
    if args.command == "prepare-campaign-ligands":
        with connect(args.db) as connection:
            result = prepare_campaign_ligands(
                connection, args.user, args.campaign, args.structures_dir,
                args.output_dir, args.ccd_cache)
        print(json.dumps(result, indent=2, ensure_ascii=False))
        return 0
    if args.command == "build-library-indices":
        with connect(args.db) as connection:
            result = build_library_indices(connection, args.library, args.output_dir)
        print(json.dumps(result, indent=2, ensure_ascii=False))
        return 0
    if args.command == "export-reduce-pocket":
        result = export_query_pocket_pdb(args.mmcif, args.output, args.ccd_id,
                                         args.chain, args.residue, args.radius)
        print(json.dumps(result, indent=2, ensure_ascii=False))
        return 0
    if args.command == "run-reduce-pocket":
        result = run_reduce(args.profile, args.pocket_pdb, args.output_dir)
        print(json.dumps(result, indent=2, ensure_ascii=False))
        return 0
    if args.command == "validate-reduce-pocket":
        result = validate_reduce_run(args.output_dir, args.ccd_id, args.chain, args.residue)
        print(json.dumps(result, indent=2, ensure_ascii=False))
        return 0 if result["accepted"] else 2
    if args.command == "build-reduce-het-dictionary":
        result = build_reduce_het_dictionary(args.ccd, args.ccd_id, args.output)
        print(json.dumps(result, indent=2, ensure_ascii=False))
        return 0
    if args.command == "enrich-query-anchors":
        result = enrich_manifest_with_reduce(args.query_manifest, args.ccd,
                                             args.hydrogenated_pdb, args.output)
        print(json.dumps({"output": str(args.output), "anchors": len(result["anchors"])}, indent=2))
        return 0
    if args.command == "build-pharmacophore-index":
        result = build_pharmacophore_index(
            args.artifact_catalog, args.output_dir,
            bin_width=args.bin_width, max_distance=args.max_distance)
        print(json.dumps(result, indent=2, ensure_ascii=False))
        return 0
    if args.command == "compile-pharmacophore-query":
        result = compile_pharmacophore_query(args.query_manifest, args.output)
        print(json.dumps(result, indent=2, ensure_ascii=False))
        return 0
    if args.command == "search-pharmacophore-index":
        external = load_external_l1_ids(args.external_l1_ids)
        result = search_pharmacophore_index(
            args.index_catalog, args.query_plan, args.output, external_l1_ids=external)
        print(json.dumps(result, indent=2, ensure_ascii=False))
        return 0
    if args.command == "validate-pharmacophore-retrieval":
        result = run_pharmacophore_validation(
            args.artifact_catalog, args.pharmacophore_index_dir,
            args.query_manifest, args.output_dir,
            external_l1_path=args.external_l1_ids,
            faiss_index_dir=args.faiss_index_dir, mmcif=args.mmcif,
            ccd=args.ccd, ccd_id=args.ccd_id, faiss_k=args.faiss_k,
            nprobe=args.nprobe, bin_width=args.bin_width,
            max_distance=args.max_distance)
        print(json.dumps(result, indent=2, ensure_ascii=False))
        return 0 if result["accepted"] else 2
    if args.command == "prepare-gaussian-query":
        result = prepare_gaussian_query(
            args.mmcif, args.ccd, args.query_manifest, args.output)
        print(json.dumps(result, indent=2, ensure_ascii=False))
        return 0
    if args.command == "score-gaussian-candidates":
        result = score_gaussian_candidates(
            args.artifact_catalog, args.query, args.candidate_ids, args.output,
            sigma=args.sigma, cutoff=args.cutoff,
            pair_tolerance=args.pair_tolerance, axial_samples=args.axial_samples,
            max_pair_seeds=args.max_pair_seeds)
        print(json.dumps(result, indent=2, ensure_ascii=False))
        return 0
    if args.command == "run-staged-gaussian-reranking":
        result = run_staged_gaussian_reranking(
            args.artifact_catalog, args.query, args.candidate_ids,
            args.output_dir, stage=args.stage, workers=args.workers,
            chunk_size=args.chunk_size,
            top_n_per_objective=args.top_n_per_objective,
            sigma=args.sigma, cutoff=args.cutoff,
            pair_tolerance=args.pair_tolerance,
            axial_samples=args.axial_samples,
            max_pair_seeds=args.max_pair_seeds,
            resume=not args.no_resume,
            progress_every=args.progress_every)
        print(json.dumps(result, indent=2, ensure_ascii=False))
        return 0
    if args.command == "campaign-ligands":
        with connect(args.db) as connection:
            result = campaign_ligand_status(connection, args.user, args.campaign)
        print(json.dumps(result, indent=2, ensure_ascii=False))
        return 0
    if args.command == "compare-campaign-ligands":
        with connect(args.db) as connection:
            result = compare_campaign_ligands(connection, args.user, args.campaign)
        print(json.dumps(result, indent=2, ensure_ascii=False))
        return 0
    if args.command == "select-query-ligand":
        with connect(args.db) as connection:
            select_query_ligand(
                connection, args.user, args.campaign, args.ligand, args.rationale)
            result = campaign_ligand_status(connection, args.user, args.campaign)
        print(json.dumps(result["query_lock"], indent=2, ensure_ascii=False))
        return 0
    if args.command == "search-similar-ligands":
        conformer_map = None
        if args.conformer_map_json:
            conformer_map = json.loads(args.conformer_map_json.read_text(encoding="utf-8"))
        with connect(args.db) as connection:
            result = run_hierarchical_library_search(
                connection, args.user, args.campaign, args.library,
                args.morgan_index, usrcat_index_path=args.usrcat_index,
                conformer_to_molecule=conformer_map, two_d_pool=args.two_d_pool,
                limit=args.limit)
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
    if args.command == "create-task":
        initialize(args.db)
        with connect(args.db) as connection:
            task_id = create_task(
                connection, args.user, args.project, args.name, args.objective,
            )
        print(json.dumps({"task_id": task_id, "project_id": args.project}, indent=2))
        return 0
    if args.command == "list-tasks":
        with connect(args.db) as connection:
            result = list_tasks(connection, args.user)
        print(json.dumps({"tasks": result}, indent=2, ensure_ascii=False))
        return 0
    if args.command == "activate-task":
        with connect(args.db) as connection:
            result = activate_task(connection, args.user, args.task)
        print(json.dumps(result, indent=2, ensure_ascii=False))
        return 0
    if args.command == "active-task":
        with connect(args.db) as connection:
            result = active_task(connection, args.user)
        print(json.dumps({"active_task": result}, indent=2, ensure_ascii=False))
        return 0
    if args.command == "deactivate-task":
        with connect(args.db) as connection:
            deactivate_task(connection, args.user)
        print(json.dumps({"user_id": args.user, "active_task": None}, indent=2))
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
    if args.command == "fetch-alphafold-db":
        print(json.dumps(acquire_alphafold_db_mmcif(
            args.uniprot, args.project_root), indent=2))
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
    if args.command == "prepare-structure-prediction":
        profile = load_model_profile(args.profile)
        if profile["backend"] != args.backend:
            raise ValueError("Profile backend does not match --backend")
        path = write_prediction_input(
            args.backend, args.name, args.sequence, args.output,
            args.project_root, seeds=args.seed or (1,),
        )
        print(json.dumps({"backend": args.backend, "input": str(path),
                          "command": prediction_command(profile, path, path.parent / "output")},
                         indent=2))
        return 0
    if args.command == "doctor":
        report = environment_report()
        print(json.dumps(report, indent=2))
        return 0 if report["ready"] else 1
    raise AssertionError("unreachable")


if __name__ == "__main__":
    raise SystemExit(main())
