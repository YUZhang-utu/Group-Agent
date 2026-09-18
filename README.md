# Group-Agent: AIDD Macrocycle Agent

An automated, local-first, resumable AIDD workflow for collaborative drug design.
The first milestone
imports batch MOL2 libraries into a versioned SQLite registry while keeping all
structures and metadata local.

The platform is multi-user and Project-first. Models, docking, rescoring, MD,
libraries, and compute caches are shared; scientific state is isolated by user
and Project. A Project must be activated before scientific work can begin.
Campaigns are optional branches inside a Project.

## Current capabilities

For molecule-level pose review and E031 equivalence/performance validation on the
completed E034 results, use [E035 workstation instructions](to_human/E035_WORKSTATION_RUN.md).
Default scale is 100,000 distinct library conformers (600,000 scoring comparisons);
the optional million-conformer tier remains an engineering stress test, not docking validation.

For the expanded-library WEE1 retrieval → Gaussian refinement → E031 workstation
run, see [E034 execution and resume instructions](to_human/E034_WORKSTATION_RUN.md).
The E033 calibration panel has passed per the supplied workstation report;
E034 target-specific results are still pending.

- Streams one or many `@<TRIPOS>MOLECULE` records from each MOL2 file.
- Groups terminal names such as `compound_conf0`, `compound_conf1`, and
  `compound_conf2` under one molecule.
- Treats a molecule without a terminal conformer suffix as conformer `0`.
- Validates declared atom/bond counts and required MOL2 sections.
- Assigns stable IDs and records SHA-256 provenance.
- Is idempotent: re-importing the same record does not duplicate it.
- Produces a JSON import report and a queryable SQLite database.
- Creates audited Campaigns, binds a biological target, searches RCSB/PDB,
  registers candidates, and freezes a human-reviewed receptor choice.
- Records only methods that actually run, using generic Run, Artifact,
  input/output lineage, metric, and decision metadata.
- Verifies result manifests and rejects outputs outside the active Project.
- Exports credential-free Slurm Bundles for OpenSSH on macOS/Linux/Windows and
  PuTTY Plink/PSCP on Windows, then imports content-bound submission receipts.
- Downloads validated RCSB mmCIF structures with checksum metadata.
- Generates validated AlphaFold 3 inputs and credential-free execution commands.
- Builds persistent Morgan 2D and USRCAT 3D ligand-search indices.

## Quick start

```powershell
cd D:\agent\projects\aidd_macrocycle_agent
python -m aidd_agent.cli init-storage --storage-root D:\AIDD
python -m aidd_agent.cli init --db D:\AIDD\registry\aidd.sqlite3
python -m aidd_agent.cli import-mol2 --db data\aidd.sqlite3 --library macrocycles-v1 --input D:\path\to\mol2
python -m aidd_agent.cli summary --db data\aidd.sqlite3 --library macrocycles-v1
python -m aidd_agent.cli create-user --db D:\AIDD\registry\aidd.sqlite3 --username alice --display-name "Alice"
python -m aidd_agent.cli create-project --db D:\AIDD\registry\aidd.sqlite3 --user USR-... --name "WEE1 Macrocycle Screening" --objective "Discover macrocycle binders" --storage-root D:\AIDD
python -m aidd_agent.cli activate-project --db D:\AIDD\registry\aidd.sqlite3 --user USR-... --project PRJ-...
python -m aidd_agent.cli create-task --db D:\AIDD\registry\aidd.sqlite3 --user USR-... --project PRJ-... --name "Structure Selection" --objective "Select target structures"
python -m aidd_agent.cli activate-task --db D:\AIDD\registry\aidd.sqlite3 --user USR-... --task TSK-...
python -m aidd_agent.cli create-campaign --db D:\AIDD\registry\aidd.sqlite3 --user USR-... --project PRJ-... --task TSK-... --name campaign-001 --objective "Screen the ATP site"
python -m aidd_agent.cli set-target --db D:\AIDD\registry\aidd.sqlite3 --user USR-... --campaign CAM-... --name "Target name" --organism "Homo sapiens" --uniprot P12345
python -m aidd_agent.cli search-pdb --db D:\AIDD\registry\aidd.sqlite3 --user USR-... --campaign CAM-... --uniprot P12345 --max-resolution 3.0
python -m aidd_agent.cli create-structure-comparison --db D:\AIDD\registry\aidd.sqlite3 --user USR-... --campaign CAM-... --name ATP-pocket-review --pdb 1ABC --pdb 2XYZ --reference-pdb 1ABC --pocket-residues 41,54,80,83,104 --chain 1ABC=A --chain 2XYZ=A --rationale "Compare pocket completeness, conformation, and bound ligands"
python -m aidd_agent.cli prepare-pymol-review --db D:\AIDD\registry\aidd.sqlite3 --user USR-... --comparison CMP-... --project-root D:\AIDD\users\alice\projects\PRJ-...
python -m aidd_agent.cli create-ai-review --db D:\AIDD\registry\aidd.sqlite3 --user USR-... --project PRJ-... --campaign CAM-... --task-type structure_selection --subject-type campaign --subject-id CAM-... --prompt-version structure-review-v1 --prompt-file prompt.md --evidence-json evidence.json --data-class public_structure_metadata --data-class computed_summary
python -m aidd_agent.cli import-ai-recommendation --db D:\AIDD\registry\aidd.sqlite3 --user USR-... --request AIR-... --provider openai --model MODEL-NAME --response-json recommendation.json
python -m aidd_agent.cli select-pdb --db D:\AIDD\registry\aidd.sqlite3 --user USR-... --campaign CAM-... --pdb 1ABC --rationale "Relevant state and ligand; complete pocket"
python -m aidd_agent.cli create-run --db D:\AIDD\registry\aidd.sqlite3 --user USR-... --project PRJ-... --run-type docking --tool PLANTS --tool-version 1.2 --backend slurm --parameters-json parameters.json
python -m aidd_agent.cli import-result-manifest --db D:\AIDD\registry\aidd.sqlite3 --user USR-... --manifest D:\AIDD\users\alice\projects\PRJ-...\runs\RUN-...\result_manifest.json
```

Use `--dry-run` to validate a library without writing to the database. By
default files are referenced in place and are never modified. Large imports
print progress every 10,000 records; change this with `--progress-every`.

## Architecture boundary

Codex/LangGraph will control decisions and resumable workflow state. Local
adapters will run structure preparation, docking, prediction, and MD. Heavy
artifacts and experimental measurements do not enter LLM state.

See `docs/architecture.md` for the planned screening graph.
See `docs/end-to-end-roadmap.md` for the computational screening, docking,
prediction, laboratory automation, and closed-loop learning roadmap.
See `docs/local-pharmacophore-retrieval.md` for the once-per-library local 3D
feature-pair index and per-co-crystal loose/balanced/strict query workflow.
See `docs/gaussian-artifact-reranking.md` for artifact-backed Gaussian reference
reranking and its current atom-centered feature boundary.
See `docs/billion-scale-tiered-search.md` for the fixed-budget shard-native path
that retains slim coarse and detailed outputs.
See `docs/multi-cocrystal-aggregation.md` for conformer-to-molecule collapse,
site-scoped multi-query fusion, protected query quotas, and docking task export.
See `docs/cluster-bridge.md` for cross-platform authenticated Slurm submission.
See `docs/desktop-deployment.md` for workstation, RDKit, PyMOL, and AlphaFold setup.

## E036 fast 3D execution
See [workstation instructions](to_human/E036_FAST_3D_RUN.md) for bounded parallel
validation chunks and whole-index WEE1 search with finer Gaussian scheduling.
The completed E035 million-conformer evidence is summarized in
[to_human/E035_WORKSTATION_1M_RESULT.json](to_human/E035_WORKSTATION_1M_RESULT.json).
