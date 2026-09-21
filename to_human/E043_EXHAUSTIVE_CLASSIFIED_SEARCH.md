# Exhaustive 3D screening and classified crystal query features

This extends E042 without requiring PLIP. Candidate 3D feature matching is not
candidate-complex interaction analysis. Docking preparation is a separate later
step, never an implicit dependency of search or classification.

## Update and start

```bash
git pull --ff-only origin main
bash scripts/setup_prompt_profiles.sh
export AIDD_LLM_CONFIG_DIR=/mnt/local/hand/yuzhang/aidd/config/llm
export AIDD_RUNTIME_PROFILE=/mnt/local/hand/yuzhang/aidd/config/prompt-runtime-apptainer.local.json
export AIDD_DOCKING_PROFILE=/mnt/local/hand/yuzhang/aidd/config/docking.local.json
bash scripts/run_chat_agent.sh --allow-compute
```

Stop an old server first. Retain API keys in existing environment variables.
Restarting with new code requires new tasks, not resuming old code-bound execution.
Previously completed searches can still be attached and reviewed.

## 1. Compare every conformer before refinement

> Run exhaustive 3D search for wee1_both. Compare every library conformer descriptor before selecting candidates. Keep E031 annotation-only and do not dock.

New chat searches default to `retrieval_mode=exhaustive`. The low-level E036 CLI
requires `--exhaustive`; its historical default remains unchanged. The exhaustive
pass reads each descriptor shard in chunks of 65,536 rows, compares every stored
conformer and merges exact Top-10,000 per query with stable global-ID ties.
It does not visit only selected FAISS partitions. Memory is bounded by chunk and
candidate budgets. Both queries currently make separate streaming passes.

Gaussian refinement then uses the existing per-objective Top-N allocation.
Exhaustive descriptor retrieval is not exhaustive docking, activity recall or
proof that discarded candidates cannot bind. Changed candidates cannot be
compared for output equivalence against the old approximate candidate set.
Coverage must equal the accepted library conformer count; actual descriptor
scanning time is included in query execution. Do not promise previous ANN latency.

## 2. Review and classify crystal-derived evidence

> Review search task SEARCH_TASK_ID and show the full-library compared count and retained counts.

Or `/evidence SEARCH_TASK_ID`. Then:

> Classify the crystal-derived 3D query features from review task REVIEW_TASK_ID. Show interaction categories and candidate feature-match counts. Do not build candidate complexes or run docking.

Or `/classify REVIEW_TASK_ID`. You can classify the existing ANN search now to
inspect the table, but that does not retroactively make that search exhaustive.

The classification task produces `interactions.html`, `interaction-classes.csv`,
`report.json` and classified match NPZs. Open the HTML locally in a browser. Each
category contains the query, exact selection ID, ligand atom indices, protein
partner, crystal geometry, evidence level and illustrative match counts.

| Category | Query evidence and matching scope |
|---|---|
| Hydrogen bond | Existing HBD/HBA crystal anchors; preserve original angle limitations |
| Hydrophobic | Indexed ligand hydrophobe near a protein carbon from a conservative residue template |
| Pi stacking | Indexed ligand aromatic plane and complete protein aromatic-ring template; distance, normal angle and offset rules |
| Cation-pi | Indexed aromatic/positive-ionizable feature near a complementary protein template; distance hypothesis |
| Salt bridge | Indexed ionizable group near complementary ASP/GLU or LYS/ARG group; protonation assumption is explicit |
| Water bridge | Ligand polar feature near a crystallographic water and a nearby protein polar atom; two-distance hypothesis only |
| Metal coordination | Indexed acceptor near a crystallographic metal; proximity hypothesis, not validated coordination chemistry |
| Halogen bond | Explicitly unsupported: this library schema has no halogen-specific feature channel |

These are transparent native rules, not a PLIP-equivalent implementation. Defaults
are recorded in the report. Charge templates do not determine protonation states;
water orientation, metal coordination geometry and alternative locations require
review. Missing rings/waters/metals or no detected hypothesis do not establish
chemical absence. For pi stacking, both near-parallel and near-perpendicular plane
angles are allowed under the recorded offset rule. The report is exploratory.

The six existing indexed feature families are donor, acceptor, positive/negative
ionizable, hydrophobe and aromatic. No re-precompute of the library is required.
Crystal templates propose relevant query features; the scorer compares their
positions/types/directions against stored candidate poses. It does not claim
that each candidate makes the named protein contact. Multiple residue hypotheses
using the same ligand feature share a score column and are not independent proof.
Expanded assignment competition can change diagnostic feature scores; original
E031 files and Gaussian ranks remain immutable. This is not an E031 equivalence run.

## 3. Select on the classified table, then export

> From classification task CLASSIFICATION_TASK_ID, preview molecules matching all of EXACT_ANCHOR_ID_1 and EXACT_ANCHOR_ID_2, with minimum score 0.5 and no molecule cap. Do not export yet.

Replace IDs with the actual table entries. The threshold is an illustrative
user-selected setting, not a validated cutoff. `all` conditions must hold in one
stored pose. Select anchors from one crystal query per preview; cross-query
conditions are not merged. Counts separate conformers, molecules and cap effects.

> Confirm selection task SELECTION_TASK_ID and export its selected poses and IDs for docking preparation. Do not run docking.

Or `/export SELECTION_TASK_ID`. Results remain accessible through `/results`.
Only compact summaries and up to 64 anchor IDs per query are included in LLM
context; the HTML/CSV/JSON contain all evidence. Supply exact IDs for longer lists.

## 4. Derive the pocket and prepare a docking workflow

> Prepare docking inputs from export task EXPORT_TASK_ID using the crystal ligand pocket. Show its center and dimensions. Do not execute docking yet.

Or `/prepare_docking EXPORT_TASK_ID`. No pre-existing grid is required.
The workflow removes the selected reference ligand from the first crystal model,
retains other residues/waters/cofactors for review, exports the native ligand and
derives the site from its bounding box midpoint with 5 A padding. It writes a
`pocket.json`, receptor PDB, native/selected SDFs, Glide input files and command plan.
The provisional inner box is 10 A and outer dimensions at least 20 A; larger
selected ligands and nearby alternate sites require review. No command runs yet.

This is preparation input generation, not a claim that hydrogens, missing residues,
protonation or stereochemistry have already been resolved. Inspect retained
cofactors, default installed preparation chemistry, and the exported heavy-atom
ligand limitations. Preparation cannot infer every biological choice from a pocket.

## 5. Explicitly execute preparation and gated docking

> Run the prepared Glide workflow from task PREPARATION_TASK_ID. Execute protein and ligand preparation, build the grid, validate reference redocking, and only then dock the selected candidates if that gate passes.

Or `/run_docking PREPARATION_TASK_ID`. The trusted local profile supplies executable
paths; the LLM cannot supply shell commands. Execution uses PrepWizard, LigPrep,
Glide and structconvert. Logs and output hashes are retained at each stage.
Reference validation uses best-of-returned-pose symmetry-aware heavy-atom RMSD
in the receptor frame, without realignment. The provisional 2 A gate is an
engineering setting, not a claim of biological enrichment or top-ranked-pose quality.
Chemically incompatible reference poses fail validation. Candidate docking is not
submitted when reference validation fails. Candidate pose quality remains unvalidated.

The installed 2025-1 CLI/options, licenses and real engine outputs must be checked
on the workstation; local tests use mocked processes. PLANTS execution is not
implemented in this adapter. A cancelled wrapper may leave vendor job-control
jobs running; inspect the vendor job manager before retrying. Failed docking tasks
require a fresh task, because partially written engine outputs are not reused.

## Evidence and validation boundaries

Local tests compare exhaustive retrieval to dense Top-K, test category geometry,
escape report content, retain selection/export provenance and block candidate
docking when a mock reference gate fails. Real full-library timing, RDKit/Biopython
crystal preparation, live model routing and licensed docking remain workstation
acceptance. Halogen-specific indexing and independent interaction/biological
validation remain explicit gaps; the five stages are not all scientifically validated.

Engineering references: [PLIP definitions](https://github.com/pharmai/plip/blob/master/DOCUMENTATION.md)
informed the distinction between feature hypotheses and geometric contacts, but
PLIP is not called. [Glide grid documentation](https://learn.schrodinger.com/public/python_api/2025-3/glide.html)
describes receptor grids. [RDKit CalcRMS](https://www.rdkit.org/docs/cppapi/namespaceRDKit_1_1MolAlign.html)
documents in-place symmetry-aware RMSD. These do not validate this implementation.
