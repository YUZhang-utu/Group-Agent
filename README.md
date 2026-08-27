# AIDD Macrocycle Agent

Local-first foundation for a resumable AIDD screening agent. The first milestone
imports batch MOL2 libraries into a versioned SQLite registry while keeping all
structures and metadata local.

The platform is multi-user and Project-first. Models, docking, rescoring, MD,
libraries, and compute caches are shared; scientific state is isolated by user
and Project. A Project must be activated before scientific work can begin.
Campaigns are optional branches inside a Project.

## Current capabilities

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
python -m aidd_agent.cli create-campaign --db D:\AIDD\registry\aidd.sqlite3 --user USR-... --project PRJ-... --name campaign-001 --objective "Screen the ATP site"
python -m aidd_agent.cli set-target --db D:\AIDD\registry\aidd.sqlite3 --user USR-... --campaign CAM-... --name "Target name" --organism "Homo sapiens" --uniprot P12345
python -m aidd_agent.cli search-pdb --db D:\AIDD\registry\aidd.sqlite3 --user USR-... --campaign CAM-... --uniprot P12345 --max-resolution 3.0
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
See `docs/cluster-bridge.md` for cross-platform authenticated Slurm submission.
See `docs/desktop-deployment.md` for workstation, RDKit, PyMOL, and AlphaFold setup.
