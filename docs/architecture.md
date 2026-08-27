# Architecture

## Multi-user Project boundary

Shared models, docking engines, MD engines, compound libraries, execution
environments, and content-addressed compute caches belong to the platform.
Scientific state belongs to exactly one Project:

```text
user -> project -> optional campaign -> executed runs and decisions
```

Every scientific command requires an authenticated user and an active Project.
Campaign commands additionally verify that the Campaign belongs to that active
Project. Access to another user's private Project is rejected before any scientific
state is read or written.

Each Project receives an independent English-only workspace under a configured
storage root:

```text
D:/AIDD/users/{username}/projects/{project-id}-{project-slug}/
  project.yaml
  research-state.yaml
  research_log.md
  findings.md
  experiments.md
  benchmarks.md
  target/
  campaigns/
  runs/
  results/
  data/
  reports/
```

Shared artifacts may be referenced by multiple Projects, but Project decisions,
parameters, outputs, approvals, and reports are never implicitly merged.

Methods are optional. A model, docking, rescoring, MD, or custom-analysis Run is
created only when that method actually executes. Absence of a Run means the
method was not performed. A skip decision is recorded only when its rationale
is scientifically meaningful.

## Run and Artifact provenance

```text
Project
  -> Workflow Run
      -> input Artifacts
      -> tool/version/parameters/compute key
      -> local or Slurm execution metadata
      -> output Artifacts
      -> metrics
```

Large files remain on the configured filesystem. SQLite stores ownership,
paths, SHA-256 hashes, sizes, formats, lineage, metrics, and execution status.
Completed local and Slurm jobs produce a versioned `result_manifest.json`,
which is verified before metadata registration.

## Control and data boundaries

The future LangGraph controller stores only campaign IDs, artifact references,
QC summaries, decisions, and pending approvals. MOL2 structures, docking poses,
trajectories, and experimental measurements remain in local stores.

## Screening graph

```text
locate_user
  -> list_or_create_project
  -> activate_project
  -> optionally_create_campaign
  -> resolve_target
  -> search_structures
  -> human_select_structure
  -> prepare_receptor
  -> define_pocket
  -> human_confirm_pocket
  -> select_library
  -> choose_optional_methods
      -> optional_docking
      -> optional_model_prediction
      -> optional_rescoring
      -> optional_MD
  -> aggregate_by_molecule
  -> rank_unique_candidates
  -> human_select_md
  -> prepare_and_validate_md
  -> human_approve_submission
  -> monitor_md
  -> analyze_md
  -> human_select_synthesis
```

The candidate count is applied after conformer and pose aggregation. Twenty
requested candidates therefore means twenty unique molecule IDs.

## Planned adapters

1. RCSB structure search and download.
2. Receptor preparation and pocket definition.
3. PLANTS docking using the existing local installation and Slurm scripts.
4. Versioned local property/model adapters.
5. MD preparation, validation, submission, monitoring, and analysis.

Every expensive Run uses a compute key derived from input artifact hashes,
tool/model version, parameters, and execution environment so successful results
can be reused safely.

## Cluster execution

The desktop control plane exports credential-free Slurm Bundles. A client on
Windows, macOS, or Linux completes the cluster's interactive authentication,
uploads the Bundle, invokes the cluster runner, and returns a submission
receipt. See `docs/cluster-bridge.md`.
