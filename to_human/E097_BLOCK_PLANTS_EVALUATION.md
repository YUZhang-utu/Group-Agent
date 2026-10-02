# E097: reproducible block sampling and PLANTS evaluation

Protocol recorded before implementation/fixture evaluation, 2026-10-02.

Sample 100 conformers uniformly without replacement from every regular block of
E094, E095 and E096, using a recorded seed and independent per-scheme/block random
streams. Expected slots: 38,400 + 118,800 + 118,900 = 276,100. Deduplicate by source
conformer ID, never by molecule. Preserve every scheme/block association. The
132,885-member special pool and deferred unsupported chemistry are outside this panel.

Hypothesis (exploratory): joint spatial partitions may produce more useful
within-block docking score distributions than backbone-only or broad-property
partitions. Block count alone cannot establish that improvement. Compare identical
receptors, engine settings and score definitions; a reused conformer/receptor result
must have identical provenance. Report missing scores and per-block sample fractions.
Population summaries weight block sample means by block population. Equal per-block
sampling is not uniform library sampling. Cross-scheme samples and same-molecule
conformers are dependent; do not report naive significance or validated search recall.

Default PLANTS mode is rigid_ligand=1 with ChemPLP: preserve input heavy-atom
conformations while docking their placement. The PLANTS 1.2 manual section 1.9
permits donor-group flexibility. This is not a guarantee that all hydrogen positions
remain fixed. Flexible docking is an explicit alternative and changes interpretation.
Prepared receptor MOL2 and binding-site center/radius can now be generated from
explicitly adopted experimental pocket representatives. Source ligand states are
preserved, with the review status recorded separately. No receptor, site or preparation
quality is inferred from the target name. Real PLANTS execution and receptor validation remain workstation
acceptance steps; local tests use synthetic fixtures.

The LLM selects typed actions; trusted runtime configuration supplies filesystem
paths. Follow-up tasks consume owned sealed task reports. Sampling, preparation,
execution and analysis are separate resumable stages. No unassessed block is rejected.

## Automatic receptor workflow

The user supplied `~/pdb/1/SPORES_64bit` and the working invocation
`SPORES_64bit --mode complete input.pdb output.mol2`. The example PDB 6Y4Q is
syntax evidence only, not an adopted MDM2 receptor. The installed binary hash and
command are recorded. Older SPORES documentation uses different mode names; the
configured user-tested `complete` mode takes precedence. The
[SPORES publication](https://pubmed.ncbi.nlm.nih.gov/20882397/) describes its protein
and ligand preparation role; this adapter does not claim pKa enumeration.

1. `protein_from_pdb` resolves RCSB polymer/UniProt identity. Multiple distinct targets
   pause for entity selection rather than guessing.
2. `receptor_assess` reuses the existing experimental-structure census and pocket
   analysis. The short recommendation distinguishes insufficient evidence, unresolved
   differences, discrete representatives and coverage sets. Literature interpretation
   remains available through the domain agent's existing literature tool.
3. Human acceptance invokes `adopt_receptors`: record the selected pocket IDs, then
   prepare those representatives with SPORES. No adoption occurs while merely asking
   for a recommendation.
4. Keep the chosen protein chain; remove other chains, the bound ligand, waters and
   other components as requested. Record all removed residue types and flag nearby
   non-water components. Missing heavy atoms are not reconstructed by this adapter.
5. Compute the native reference ligand's heavy-atom geometric centroid and an enclosing
   radius with 5 A padding. This is not the whole-protein centroid or a mass-weighted
   center. Macrocycle access/box suitability still needs target-specific evaluation.
6. Confirm the SPORES output retains the heavy-atom inventory and native coordinates
   within 0.05 A. Write protein MOL2, site evidence and `plants-profile.json` automatically.

There is no need to manually prepare receptor files or enter XYZ coordinates for this
route. Ambiguous native ligands or non-exportable chain identifiers produce specific
diagnostics. Distinct prepared receptors get separate docking jobs/results; reuse is
only across schemes for the same conformer and receptor under identical settings.

## Workstation setup

Keep any running old computation in its original checkout. For the completed E096
checkout, update the existing `feature/structure-guided-chat` branch with a fast-forward
pull. Do not reset local modifications.

```bash
cd /mnt/medchem_taltio/wrk/yu_agent/Group-Agent
git pull --ff-only origin feature/structure-guided-chat
export PYTHONPATH="$PWD/src${PYTHONPATH:+:$PYTHONPATH}"
BASE=/mnt/local/hand/yuzhang/aidd/mc_runs/e081-001/audit-stage/blocks-20000
CFG="$BASE/e097-config"
OUT="$BASE/e097-block-evaluation"
python scripts/configure_e097_blocks.py --base "$BASE" --output "$CFG"
```

The configurator discovers PLANTS from the existing environment/PATH/home installation.
If discovery fails, sampling is still available. Supply the actual executable using
`--plants /absolute/path/to/PLANTS1.2_64bit` and a fresh config output directory. SPORES
defaults to the user's path above. Use `--runtime /path/to/existing/runtime.json` when
generating the new configuration to preserve existing search, AF3 and viewer settings.
No existing runtime file is overwritten.

To start the existing chat with these capabilities, use the same storage root and LLM
configuration directory as the existing session, and pass the new runtime:

```bash
bash scripts/run_chat_agent.sh \
  --storage-root /path/to/your/existing/chat-workspace \
  --runtime "$CFG/runtime.json" \
  --llm-config-dir /path/to/your/existing/llm-config \
  --allow-compute
```

Example conversational requests (other languages also go through the LLM):

- "Use PDB ID XXXX to identify the target and assess whether one or multiple receptors are appropriate. Show me the recommendation first."
- "I accept the proposed receptor representatives. Prepare them using SPORES."
- "Sample 100 conformers from every block in E094, E095 and E096. Reuse duplicate conformers and export them."
- "Prepare PLANTS jobs for that sample using the receptor task we just completed."
- "Run the prepared docking jobs and analyze the results by block and scheme."

`XXXX` is a placeholder for a real user-selected ID. Recommendations alone never run
SPORES or docking. An explicit request for the entire evaluation can queue sample,
prepare, dock and analyze as one typed plan using an existing prepared receptor task.
The existing selected LLM provider handles intent; no fixed phrase matching is required.
The domain agent receives the same workflow contract and can inspect owned result reports.

## Independent batch entry points

For users who prefer to start sampling immediately without Chat:

```bash
nohup bash scripts/run_e097_block_samples.sh "$CFG/sampling.json" "$OUT" \
  > "$BASE/e097-launch.log" 2>&1 &
tail -f "$OUT/sampling.log"
```

Outputs: `samples/report.json`, `samples/block_samples.csv`, `samples/samples.sqlite`,
`samples/exports.json` and `samples/ligands/*.mol2`. These are globally deduplicated
chunks; the CSV links every conformer to every sampled scheme/block. Chunks of at
most 100 are scheduling units, not new chemical blocks. The last chunk of a source
file may be smaller. All 276,100 slots remain represented in the manifest.

Standalone prepared receptor and job execution (paths must reference actual reports):

```bash
python -m aidd_agent.plants_receptors \
  --adoption /path/to/adopted-pocket/report.json \
  --profile "$CFG/receptor-tools.json" --output "$OUT/receptors"
python -m aidd_agent.block_plants prepare \
  --source "$OUT/samples/report.json" \
  --profile "$OUT/receptors/plants-profile.json" --output "$OUT/prepared"
# First verify one actual engine job and its pose identity mapping.
python -m aidd_agent.block_plants run \
  --source "$OUT/prepared/report.json" --output "$OUT/docking" --max-jobs 1
# Same output resumes and reuses successful jobs; no repeated docking for duplicates.
python -m aidd_agent.block_plants run \
  --source "$OUT/prepared/report.json" --output "$OUT/docking"
python -m aidd_agent.block_plants analyze \
  --source "$OUT/docking/report.json" --output "$OUT/analysis"
```

Chat-owned stages stay in their Project run directories, not this standalone OUT.
Standalone CLI artifacts are not automatically imported as owned chat tasks; use the
chat sampling route when conversational continuation is desired. This avoids making
unowned filesystem paths selectable by the model.

Analysis outputs include `block_summary.csv`, `block_scores.csv`, `review.md` and a
machine-readable report. Each row retains CID, receptor, pose file and record index
for later N-E rescoring. Failed/missing poses remain visible and never receive zero
scores. Population-weighted scheme means are withheld until all required samples
are scored. No N-E or binding-affinity model is executed in this version.

## Validation and practical limits

Local synthetic tests cover deterministic per-block sampling, deduplication, original
coordinate preservation, source mutation rejection, failed-job retry, sealed result
reuse, score-to-pose mapping, weighted summaries, owned task continuation, PDB target
ambiguity, receptor cleaning, centroid generation and mocked SPORES/PLANTS execution.
Real binary formats, full-library runtime, MDM2 receptor suitability and live-provider
intent quality are separate workstation acceptance checks. An optimizer RNG seed is
not claimed: the seed recorded here governs panel sampling only.
