# E105: macrocycle block-prioritized search across targets

Date: 2026-10-06. Status: proposed protocol; no benchmark launched.
Current MDM2 results are exploratory development evidence. Freeze a version of
this protocol before confirmatory execution; record amendments separately.

## Question and hypotheses

Can target-independent macrocycle backbone/side-chain/spatial partitions,
prioritized by a small target-specific docking/rescoring panel, improve the
cost-versus-retrieval frontier over whole-library 3D search, random allocation,
property tranches and molecular active learning?

- H1: E095/E096 prioritize unseen promising molecules better than random blocks
  and property-only partitions at matched total online cost.
- H2: spatial side-chain information in E096 adds value beyond E095 chemistry
  descriptors and E094 backbone partitioning. A negative result is informative.
- H3: Top-10 distinct-molecule means predict held-out block yield better than
  whole-block means; Top-5, quantiles and uncertainty-aware allocation are
  secondary ablations, not changes selected after inspecting test targets.
- H4: gains persist across unrelated targets and do not result solely from
  molecular size, lipophilicity, block population or conformer multiplicity.

The contribution is a testable macrocycle-specific allocation method. Partitioning
and representative prescreening are established ideas, especially AdaptiveFlow.
Chat automation is reproducibility infrastructure, not proof of scientific novelty.

## Evidence boundary and counting

User receipts show 821619 successful EquiScore conformer/receptor pairs and 8283
block/receptor ranking rows. They do not establish experimental activity or
unseen-library recall. The three current MDM2 structures are receptor states of
one target, not three independent targets. The latest search receipt is a blocked
configuration retry; no fresh running/completed search receipt was supplied here.

Current production delivery remains at most 100000 conformers, multiple per
molecule allowed. Always report unique molecules, conformers and receptor pairs.
For evaluation, molecule-level recall is primary; conformer-level recovery is
secondary. Keep a fixed conformer ensemble across methods. Do not confuse the
Top-10 molecule block statistic with the conformer delivery budget.

## Staged design

1. Finish and inspect the existing MDM2 continuation, without rerunning its old
   PLANTS or EquiScore panels. This remains development, not held-out evidence.
2. Construct a common tractable reference library, initially about 100000 distinct
   molecules with their frozen conformer ensembles if feasible. Sample from the
   whole library before block prioritization, not from current leading blocks.
   Inventory actual molecules/conformers and estimate cost before fixing size.
3. Exhaustively evaluate this reference set with one fixed 3D search engine and
   a matched docking protocol. Score every reference molecule if claiming docking
   Top-K recall. A whole-library search followed by docking only its own Top-K is
   not an exhaustive docking oracle; label that comparison as workflow overlap.
4. Develop on MDM2, then freeze settings. Select 4-6 additional unrelated targets
   using an explicit inventory of resolved ligand-bound structures, query quality,
   macrocycle relevance, activity labels and assay feasibility. Include different
   pocket classes; exact targets remain unselected until this audit. Treat additional
   structures of a target as nested replicates. If a target class lacks known
   macrocycle activity, distinguish technical stress testing from activity validation.
5. Scale the frozen method to the real large library. Full unpartitioned 3D search
   is a computational baseline where feasible; do not assert exhaustive docking
   recall on a library whose full docking scores are unavailable.
6. Plan prospective assays on 1-2 feasible targets after computational evaluation,
   preserving equal unique-compound budgets and selection records across arms.

## Comparators and allocation

| Arm | Purpose |
| --- | --- |
| Whole-library 3D search, same queries and final docking budget | Direct workflow comparator |
| Exhaustive docking on the tractable reference set | Computational Top-K reference, not activity truth |
| Random conformers and random blocks, matched spend | Allocation baseline; repeat seeds |
| Property-only tranches with representative prescreening | AdaptiveFlow-style close baseline; distinguish faithful reproduction from adaptation |
| E094 / E095 / E096 | Isolate partition representation contributions |
| Molecular active learning, e.g. MolPAL | Competitive allocation baseline under the same docking engine |

Within each partition compare ChemPLP, EquiScore and union selection. Union arms
must pay for their extra selected population; equal numbers of blocks are not
equal budgets. Report population-budget curves in addition to Top-5/Top-10 block
examples. Allocate blocks under a frozen cost rule. Count deduplicated physical
evaluations once and attribute shared work transparently, including baseline setup.

Report both total online wall/resource cost and docking-call counts. Include pilot
docking, EquiScore inference, model training, file mapping/I/O, search, export and
final docking. The completed 821619-pair panel is not free in an end-to-end cost
comparison. Report offline library preparation/partitioning separately and amortize
it explicitly over a stated number of targets. Match hardware, concurrency and
input state. Fix chemistry, conformer generation, reference templates, alignment,
search settings and docking receptor/site across comparator arms.

## Leakage, sampling and uncertainty

- Split pilot/evaluation by molecule, keeping all conformers of a molecule in one
  split even if they occupy different blocks. Exclude pilot molecules from primary
  held-out yield and report them separately in total production yield.
- Freeze target-independent partition definitions before target-specific scoring.
  Use scaffold/sequence-held-out analyses appropriate to cyclic peptides; specify
  identity rules rather than assuming standard small-molecule scaffolds suffice.
- Exclude query ligands and predefined close analogues from novel-hit evaluation.
  Audit overlap with learned scoring training data where available; disclose unknowns.
- Equal 100-conformer panels contain unequal distinct-molecule counts. For the
  confirmatory Top-10 experiment use a fixed distinct-molecule sampling budget
  where possible and report undersized blocks separately. Do not silently reuse
  current conformer panels as if they satisfied this condition.
- Repeat pilot sampling with five fixed seeds as the initial design; budget this
  before execution. Bootstrap molecules, not individual conformers. Report paired
  differences per target, confidence intervals and failure cases; targets are the
  units for cross-target generalization. Do not pool raw scores across receptors.
- Whole-block means and Top-N tail means answer different questions. Tail means
  are potentially noisy under small samples; measure stability and winner bias.

## Metrics and decision rule

Primary computational endpoint: held-out molecule recall of the reference docking
Top-1% versus total cost. Also show Top-0.1%, Top-K recovery, time to a fixed recall,
novel scaffold/sequence coverage, block rank stability and molecule size biases.
Search-score recall and docking-score recall must be separately named. No failed
score is zero-imputed; report coverage and exclude incomplete reference panels
from claims of exact recall.

With actual activity labels, report molecule-level PR-AUC and EF1%:
EF1% = (actives in top 1% / molecules in top 1%) / (all actives / all molecules).
Docking Top-1% recovery is not experimental EF1%. Define rounding/ties upfront.

Advance if frozen partitions improve the cost-recall frontier over random and
property-tranche comparators on held-out targets without collapsing diversity.
Report superiority, equivalence or losses against active learning honestly. Final
effect-size and precision targets require a feasibility/variance pilot and must
be frozen before the untouched-target experiment. No effect size is claimed now.

## Prospective validation and manuscript scope

Balance selected and control arms by number of distinct molecules tested and
availability/synthesis filters. Preserve selection before assay results. Obtain
identity/purity confirmation, concentration-response measurements, an orthogonal
binding or functional assay as appropriate, and counterscreens for interference
and selectivity. For a strong lead, pursue cellular target engagement and structural
or mutational support where feasible. Do not equate a docking score with Kd/Ki.

A strong manuscript would establish an interpretable macrocycle allocation
advantage across targets plus prospective novel hits. Journal suitability depends
on effect sizes, novelty relative to AdaptiveFlow and others, and biological depth;
experimental hits alone do not guarantee a high-impact publication.

## Sources

- Cecchini et al., *AI-enhanced adaptive virtual screening of large libraries for
  ligand discovery*, Nature Biotechnology, published 1 September 2026.
  https://doi.org/10.1038/s41587-026-03217-x
  Closest novelty challenge: property partitions and target-guided prescreening.
- Graff, Shakhnovich and Coley, *Accelerating high-throughput virtual screening
  through molecular pool-based active learning*, Chemical Science (2021).
  https://doi.org/10.1039/D0SC06805E
  Competitive budget-aware molecular selection baseline.
- *Synthon-based ligand discovery in virtual libraries of over 11 billion
  compounds*, Nature (2022). https://doi.org/10.1038/s41586-021-04220-9
  Hierarchical chemical-space exploration; its synthon library assumptions differ
  from the fixed macrocycle conformer library here.
