# E101: fair block comparison and reserved rescoring interface

Protocol dated 2026-10-05. Status: design and local contract validation only.
No new docking, model inference, measured enrichment or affinity result exists.
The user confirmed N-E is their own model; code and an I/O example are pending.
This protocol refines E082 using the actual E097 conformer-level experiment.

## Question and evidence

Can a partition route a fixed computation budget to useful molecules more
effectively than uniform sampling, and does independent rescoring improve this?
User-reported E097 completion: 8445 attempted jobs, zero failed jobs, 821619
scored conformer/receptor pairs, first-job gate passed. These are a pasted receipt,
not independent verification of remote files. There are 273873 distinct sampled
conformer CIDs across three receptors. They are not 273873 unique molecules.
Use E098 attachment and sealed analysis before interpreting remote output.

E094/E095/E096 contain 384/1188/1189 regular blocks and used 100 conformers per
block. Consequently raw scheme sampling counts and costs differ. The special
pool of 132885 conformers was excluded. No full-library recall claim follows.
Successful export does not resolve the pending ligand-chemistry review.

## Frozen comparison design

1. Preserve source reports, hashes, CIDs, molecule IDs, receptor identities,
   microstates, pose indices and block membership for every scheme. Validate
   bond orders, stereochemistry, ring closure, clashes and atom mapping before
   conversion for another scorer. Keep rejected cases in the coverage report.
2. Report existing ChemPLP distributions, missingness and sampling fractions.
   Population-weighted means are descriptive; alternative partitions of one
   population cannot establish superiority simply by changing its mean score.
3. For paired comparisons map the same evaluated candidates to all three full
   membership tables. The union of block-stratified samples is not a uniform
   sample. Report panel-conditional estimates unless inclusion probabilities
   support justified population weighting. Do not silently extrapolate.
4. Split by molecule, never by conformer alone, before learning block priorities.
   Freeze discovery/held-out split and tie-breaking before inspection. Prefer an
   independent uniform held-out panel for a population claim; use retrospective
   disjoint molecules for an explicitly exploratory first comparison.
5. Compare routing policies at equal evaluated conformer/receptor counts AND
   report actual wall time/GPU time. Include unblocked uniform sampling and a
   size-stratified baseline. Keep a fixed receptor panel and conformer budget per
   molecule. Include a prespecified random exploration allocation.
6. Primary routing outcome: held-out top-score-tail recall versus computation
   budget, with unique-molecule coverage and chemical diversity. A conformer-level
   diagnostic is reported separately. Fix the reference scorer and reference
   population; using ChemPLP defines a ChemPLP proxy, not active-compound recall.
7. Use paired molecule-cluster bootstrap intervals; keep all conformers and all
   scheme memberships of a molecule together. No superiority claims until the
   budget, split, sample support and uncertainty are explicit.

No good hit among 100 sampled conformers proves no absence of rare hits: under
an independent 1% hit-rate illustration, the probability of zero hits is
0.99**100, approximately 36.6%. This is an illustration, not a fitted model.

## Scoring ladder and pilot

- ChemPLP: retain the completed baseline, lower is better. Current run retained
  one pose per conformer/receptor, so rescoring cannot recover discarded poses.
- N-E: reserved user-owned model interface; no meaning, units or accuracy inferred.
- GNINA: first public baseline candidate, using existing coordinates in score-only
  mode. Keep CNNscore (pose quality) and CNNaffinity (predicted affinity) separate.
  Pin binary/version, weights and configuration. Demonstrate coordinate/atom
  preservation before scaling; input compatibility does not prove macrocycle
  accuracy. Neither output is experimental binding affinity.
- PIGNet2: optional second public structure-based model after the first baseline,
  subject to noncanonical chemistry and macrocycle validation. Not a claim of
  superiority to GNINA or N-E.
- Boltz-2: optional eligible subset, not the default whole-library scorer. Official
  affinity guidance limits support to small-molecule ligand chains, permits up to
  128 atoms under its RDKit RemoveHs convention, and discourages sizes significantly
  above 56. Measure our size distribution first. Cofolded complexes require new
  pose/receptor provenance and cannot be treated as unchanged PLANTS poses.
- Peptide-specific research models remain candidates, not approved defaults:
  DeepPpIScore needs direct review of noncanonical residue and cyclization support.
  Recent cyclic-complex rescoring work evaluates pose correctness, not measured Kd.

Pilot before scale: freeze a molecule-stratified representative panel spanning
size, ring chemistry, N-methylation, charge and ChemPLP score range, plus available
positive/negative controls. Include random candidates, not only docking leaders.
First test a small technical panel, then choose and record scientific pilot size,
seed, model hashes and budget BEFORE inference. Evaluate the same panel for each
model; retain unsupported/failure status and both shared-support performance and
whole-panel coverage. No numerical scientific pilot size is fixed in this draft.

## Enrichment and affinity endpoints

With experimentally supported binary labels on a frozen unique-molecule panel:
`n_top = ceil(0.01 * N)` and `EF1 = (A_top / n_top) / (A_all / N)`.
Freeze molecular deduplication, aggregation, label thresholds and deterministic
tie handling. Unknown activity is not inactivity; without active controls EF is
undefined. Declare label coverage and distinguish decoys from tested inactives.
Report EF1%, EF5%, top-k hit counts/recall and BEDROC (primary alpha=20, fixed
before evaluation), with PR-AUC and ROC-AUC secondary. Report actual top-set sizes
and intervals. Small labeled panels may be too imprecise for an EF1 comparison.

Without labels, report score-tail recovery under an explicitly named reference
score. Do not name it experimental EF1 or replace labels with a model prediction.
Keep receptor-specific tables first; freeze ensemble aggregation on validation
controls. Never average ChemPLP, CNN outputs and another model's raw units.

For continuous experimental affinity, use matched-target/construct, assay-aware
Ki or Kd data, fixed molar/log conversion, Spearman correlation and MAE/RMSE.
Keep Ki, Kd, IC50 and censored observations distinct. Audit train/test ligand and
target overlap; report a macrocycle/noncanonical subset separately. Model-only
agreement establishes neither affinity accuracy nor biological enrichment.

## Reserved implementation boundary

`src/aidd_agent/rescoring_contract.py` defines PoseInput, ModelSpec, PoseScore and
PoseRescorer plus batch validation. One model signal per invocation; zero-based
pose indices must be explicitly converted from source indexing. Model configuration
and weights carry hashes. Output covers every input, with null failed values.
This is an internal contract only: no N-E/GNINA executor, Chat action, subprocess,
file verifier or automatic ranking change is registered. Existing sealed PLANTS
runs are unchanged. Coordinate-changing methods require a separate future contract.

N-E handoff needed: callable/script, environment, weights, real input/output pair,
units and score direction, supported chemistry, and coordinate-change behavior.
The user has not yet established availability of experimental affinity labels.

## Primary references checked 2026-10-05

- GNINA code: https://github.com/gnina/gnina
- GNINA 1.3 paper: https://pmc.ncbi.nlm.nih.gov/articles/PMC11874439/
- PIGNet2 implementation: https://github.com/mseok/PIGNet2
- Boltz affinity applicability: https://github.com/jwohlwend/boltz/blob/main/docs/prediction.md
- Early recognition / BEDROC: https://pubmed.ncbi.nlm.nih.gov/17288412/
- DeepPpIScore: https://www.nature.com/articles/s41401-025-01659-8
- Cyclic pose rescoring preprint, not affinity validation:
  https://www.biorxiv.org/content/10.64898/2026.08.20.746104v2
