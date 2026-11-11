# Scientific P0 delivery and workstation testing

## Delivery scope

This release adds three project-scoped model tools: `pose_qc`,
`screening_benchmark`, and `af3_batch`. They use the existing sealed task queue,
project ownership, cancellation, task reports and dependency workflows. It extends
`molecule_prepare` with explicit chemical-state policies. Restart Chat after updating.

| P0 | Implemented | Acceptance boundary |
|---|---|---|
| Structure quality | Partner-chain fix in legacy QC; assembly-aware ligand checks; saved candidate-pose QC; AF3 QC reports | Available geometry checks, not experimental validity |
| Full-library evaluation | Saved-ranking benchmark, positive attrition, exact-reference recovery, conditional enrichment, measured search timing breakdown | Production speed/recall and assay-label quality require workstation data |
| Molecular preparation | Explicit fragment/tautomer/user-state policies, parent/state lineage, adaptive conformer budgets, macrocycle settings | Finite sampling; no automatic pH model or validated macrocycle coverage |
| AF3 batches | Chat and project workflow entry, independent compounds, shared compatible MSA/templates, sequential resource lock, failed-child resume, sample confidence and pose consistency | Local simulated runtime checks do not establish AF3/GPU success |
| Real acceptance | Versioned golden-prompt runner, environment checks, exact result/input assertions, task IDs and restart receipts | User will run live provider/AF3/production library acceptance |

## Update and start

In the existing workstation repository, with the Agent Python environment active:

```bash
git pull --ff-only origin feature/structure-guided-chat
python -m pip install -e .
```

Stop and restart the existing Chat server with the same runtime configuration,
provider configuration and storage root. Preserve all project data. For example,
replace the configuration paths with your existing paths:

```bash
bash scripts/run_chat_agent.sh \
  --storage-root /mnt/local/hand/yuzhang/aidd/prompt-workspace \
  --runtime /absolute/path/to/your/runtime.json \
  --llm-config-dir /absolute/path/to/your/llm-config-directory \
  --allow-compute
```

Keep any existing `--input-root` / `--result-import-root` arguments needed for local
receptor/reference files. Do not replace an installed AF3 profile with an example.
No AF3 weights, databases or GPU software are installed by this release.

## Manual checks through Chat

Use one target Project. Replace all placeholder IDs with IDs actually returned by
the application. A queued receipt confirms submission only. Inspect the task panel
or `/status TASK_ID` and `/results TASK_ID` before proceeding.

### 1. Molecular preparation

```text
Prepare CC(=O)C with 3 conformers and seed 47. Preserve its chemical state.
Report the parent molecule, generated state identifiers, actual conformer count,
optimization failures and saved structure paths. Use molecule_prepare.
```

Then explicitly test a state decision:

```text
Prepare CC[NH3+].[Cl-]. Explicitly select the largest fragment and also include
the user-supplied state CCN. Keep the original input and record every chemical
transformation. Use max_states 2 and conformer_policy adaptive.
```

Expected: one input record, two chemical states, original salt provenance and
removed fragment retained. No automatic protonation-state prediction is performed.
Tautomer enumeration is opt-in (`tautomer_policy enumerate`); the parent state is
retained and truncation is recorded. Adaptive budgets are exploratory heuristics,
not proof of conformational coverage. Original 3D coordinates are preserved.
Check `preparation_status`, `failed_molecules`, `molecules.json` and `failures.csv`.

### 2. Full local library and saved ranking

Attach a single 3D SDF or MOL2 ligand, then send:

```text
Use the attached 3D ligand to search every stored conformer of our configured
local library with ligand_library. Export the top 10 unique molecules as original
and posed MOL2. Keep the complete ranking locally and report coverage, timing,
export count and any shortfall. Do not substitute an ANN budget search.
```

Expected: `full_coverage: true`, actual conformer/molecule counts, `ranking.sqlite`,
original/posed MOL2, and timing fields. `fresh_latency_eligible` is false for a
resumed invocation. Timings include their stated scope; no speedup is asserted.

```text
Import the complete ranking from task SEARCH_TASK_ID as a candidate set named
Full reference. Then create a separate set containing its first 100 molecules.
Keep the complete reference set and reuse saved results without another scan.
```

Supply independent assay labels as CSV `molecule_id,label`, with labels `0` or `1`.
The IDs must match the local library; unknown/duplicate labels are rejected.

```text
Use screening_benchmark with reference_set_id FULL_SET_ID and stage_set_ids
[TOP_SET_ID], top_k 100, and the attached label CSV. Report positive losses,
positive recall, exact-reference recovery and enrichment where justified.
```

Without labels, exact score-ranking recovery is still available. Enrichment is
null unless the entire reference universe has explicit binary labels. An unlabelled
molecule is not treated as inactive. A filtered reference is only a candidate-panel
benchmark, not a whole-library validation. Review assay conditions and held-out
splits before making predictive claims. Catalog identity does not establish that
two rankings used comparable target/reference protocols.

### 3. Assembly-aware pose QC

Attach the intended receptor assembly, including both chains for an interface
pocket. Use the saved candidate set:

```text
Use pose_qc on candidate set TOP_SET_ID, count 10, with the attached receptor
assembly. Evaluate all supplied receptor chains. Report chemical identity,
specified tetrahedral geometry, bond-length anomalies, internal clashes and
receptor clashes separately. Do not interpret QC as affinity.
```

If the receptor still contains its native ligand, explicitly name it for exclusion
using `exclude_residues: [["A", "901", "LIG"]]` with the actual chain/residue/name.
Otherwise it is included in collision checks. Optional `required_contacts` uses
explicit chain/residue/distance rules in the same candidate pose. Contacts are
distance checks, not hydrogen-bond or metal coordination assignments.

Alternatively attach one posed SDF/MOL2 ligand and a receptor and ask for pose QC.
All supplied atoms, including other chains/cofactors/waters, are considered unless
explicitly excluded. The software does not infer missing biological assemblies.
`passed_available_checks` is distinct from `experimental_status: not_assessed`.
Missing receptor evidence remains null; it is never a receptor-QC pass.

### 4. AF3 monomer, dimer and independent ligand batch

First reuse your previously successful exact-sequence single AF3 prompt. Check the
input JSON, CIF, `sample-summary.json`, `pose-quality.json` and execution log.
Missing/ambiguous structure or atom correspondence is reported, never fabricated.

For a batch, reuse a completed project task with exactly one AF3 prepared input:

```text
Run af3_batch named CIP2A_candidates using the exact protein assembly from
protein_task_id PRIOR_AF3_TASK_ID. Predict separately with ligands
c1cc4OCOc4(cc1CN2CCN(CC2)C(Nc3ccc(cc3)Cl)=S) and CCO.
Use exactly model seeds: [47,48,49,50,51]. Enable shared_msa.
These are independent complexes, not two co-ligands in the same structure.
Keep all results locally. This is a runtime check, not evidence that either binds.
```

To test a homodimer, provide the complete sequence directly and explicitly request
`protein_copies 2`. For a heteromer, supply each sequence in `protein_sequences`;
repeat entries for intended stoichiometry. Reusing an existing task preserves its
assembly and does not multiply copies again. Batch inputs support at most 25
protein chains plus one ligand, and 1..100 seeds. Do not invent a dimer merely
because the adapter supports it.

To use screened compounds:

```text
Run af3_batch named screened_candidates using protein_task_id PRIOR_AF3_TASK_ID,
set_id TOP_SET_ID and count 5. Use shared_msa and exactly model seeds: [47,48].
Run each candidate as an independent complex and retain all confidence reports.
```

Expected: one prepared input per ligand; compatible protein features computed once;
separate child attempt directories; exact seeds/stoichiometry; completed CIF hashes;
all recorded sample confidence and the recorded best-ranking row. The sample report
also reports seed coverage: mismatched recorded seeds fail import; absent ranking
evidence is explicitly not verified and is never inferred from the input alone.
Pose consistency uses exploratory 2-Angstrom greedy clusters after whole-protein C-alpha alignment,
with atom-name correspondence and no symmetry correction. It is not validated
binding-mode discovery. AF3 confidence and sample ranking are not affinity rankings.

The native and Apptainer profiles use the same trusted command builder. Missing
configuration preserves prepared inputs and blocks computation. Changing sealed
code/runtime requires a new plan. Do not resume an old pre-upgrade task with changed
code or overwrite its artifacts.

For a known child failure, `/resume BATCH_TASK_ID` reuses verified completed children
and retries failed children. For interrupted GPU dispatch, first verify the old
process is stopped on the workstation, then reconcile that child locally:

```bash
PYTHONPATH=src python -m aidd_agent.af3_batch \
  --state /absolute/path/to/execution/work/batch-state.json \
  --job ligand-00002 --confirm-process-stopped
```

Use `--job shared-data` for an uncertain shared-feature dispatch. Then resume the
same task in Chat. This operator-only command is not available to the cloud model.
Changed model/input/confidence files invalidate reuse. Cancellation retains local
artifacts. The per-user host lock serializes AF3 tasks; occupied resources report
blocked and can be resumed after the active task exits.

## Automated golden-prompt runner

Start with environment checks only:

```bash
bash scripts/run_p0_acceptance.sh \
  --suite configs/p0-cpu-acceptance.json \
  --output /absolute/path/to/p0-check \
  --storage /absolute/path/to/p0-isolated-workspace \
  --runtime /absolute/path/to/your/runtime.json \
  --config-dir /absolute/path/to/your/llm-config-directory
```

Add `--run` to execute the live-provider CPU preparation case. The default provider
is deepseek; use `--provider gpt` for the configured GPT profile. API keys remain in
environment variables and are not written to acceptance reports.

Copy `configs/p0-live-acceptance.example.json` to a local, untracked suite and replace
all `REPLACE_` values with your selected target, exact sequences and reference files.
Run it with the same wrapper and `--run`. It includes uploaded full-library search,
explicit-target search, AF3 monomer and independent-ligand homodimer cases. Use a
dedicated acceptance workspace, not the running production Chat workspace. Keep the
target/reference files inside that workspace's project or its explicitly approved
input roots when using the manual Chat workflow.
For automated explicit-target cases, pass `--input-root /directory/containing/inputs`
to the acceptance wrapper; this option may be repeated for distinct input roots.

Each case checks actual task status, actions, exact result fields, input seeds and
stoichiometry where specified, and available payload hashes. `acceptance.json`
retains task IDs and local report paths. Re-running the same unchanged suite resumes
those tasks rather than resubmitting prompts. An uncertain prompt submission blocks
for reconciliation. Timeout interrupts this runner's worker and retains task IDs;
uncertain GPU children still need the operator reconciliation above. The suite does
not claim browser visual acceptance or universal natural-language reliability.

## Feedback checklist

E137 adds original-output review for an existing completed `af3_run` task:

```text
Use af3_review on task TASK_ID. Read the original AF3 ranking table and sample
files without launching inference or rewriting existing reports. Report observed
seeds, sample/model counts, the best recorded seed and sample, its CIF path, and
missing evidence. Distinguish metric dimensions from chain-label correspondence.
```

Both bare and job-prefixed sample files are supported. Sample count is observed,
not inferred from seed count. Follow pagination for long sample/ranking tables.

For each manual case, record the exact prompt, task ID, terminal status, expected
versus actual behavior and the relevant report. For a failure, include the task
error and the last approximately 40 lines of its local execution log, without keys.
For scientific review, include input identity, coverage/retention, actual output
count, QC failures, AF3 sample consistency and measured timing scope. Keep private
documents out of cloud chat; the existing Knowledge libraries boundary is unchanged.

The local development host has CPU chemistry/structure dependencies and an RTX 5070
Laptop GPU (8151 MiB), but no configured AF3 runtime, production catalog or provider
credential. Live AF3, cloud-provider and production-scale quality/performance remain
workstation acceptance tasks, as requested by the user.
