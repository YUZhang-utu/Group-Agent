# E053: general-target reference, pocket and consensus contract

Status: requirements and acceptance contract; not an implemented pipeline.
MDM2 is the first validation project, not the source of generic defaults.
This document supersedes treating an E052 contact-overlap proposal as a
validated same-pocket cohort. Preserve E050 jobs and all sealed E052 outputs.

## Objective and scope

Provide a chat-driven workflow that verifies target and site, derives independent
interaction consensus and shape-template sets, validates retrieval quality, and
screens an entire configured library. The user reviews evidence and decisions;
routine work must not require scripts, JSON editing, or manually entered paths.
Unknown identity/site/quality is an explicit state, never silent inclusion.

## 1. Immutable project reference

The target specification contains organism/taxonomy, canonical accession,
isoform/variant if applicable, sequence checksum, intended domain and exclusions.
Never substitute a homolog, paralog, other organism or isoform automatically.

The user can supply a PDB ID or a local co-crystal structure through the future
project attachment adapter. Resolve it to a specific model, biological assembly,
target entity and chain instance, ligand entity/instance, residue identifiers,
alternate-state policy and coordinate checksum. A PDB ID alone is not a unique
pocket. If several plausible complexes remain, display them and ask the user
which site is intended; do not silently choose chain A or the first ligand.

An unambiguous reference may be resolved automatically and displayed for review.
If no reference is supplied, propose verified candidates for human site selection.
An uploaded structure does not bypass target or pocket checks. Aligned/truncated
exports preserve local provenance and original atom mappings; missing native
metadata must be recovered from a verified source or reported as unresolved.

Freeze a reference frame, target-residue map, pocket-lining residue set and actual
ligand envelope. Changing the reference/site creates a new design version and
invalidates dependent membership decisions rather than overwriting old results.

## 2. Candidate discovery is not admission

PDB search results are candidate entries only. Admission is per ligand-target
complex instance, not per PDB title or entry. A matching protein can occur in an
entry while its ligand binds a different chain or an entirely different protein.

Required checks, with stored evidence:

1. Target metadata: accession, organism, entity and construct associations.
2. Sequence correspondence: verify the actual target chain/construct against the
   project sequence; map label/auth numbering, insertion codes, substitutions,
   tags, missing segments and chain breaks. Evaluate intended-domain and pocket
   coverage separately from full-sequence coverage. Never require a domain-only
   construct to cover the full-length protein or conceal a pocket mutation.
3. Local coordinates: demonstrate that the ligand contacts the verified target
   chain instance. Do not combine contacts from unrelated chains to pass a gate.
4. Assembly context: distinguish the intended biological/interface site from
   contacts to crystal symmetry mates or packing neighbors. Multi-chain target
   sites require an explicit assembly definition; do not always discard interfaces.
5. Identity classification: correct target, explicit alternate cohort, rejected
   mismatch, or unresolved. Contradictory annotations are not settled by LLM vote.

Homolog/antitarget structures can be useful in a separate comparison cohort but
cannot silently enter target consensus denominators. An explicit scope change
may create such a cohort; human approval must not relabel a mismatch as a match.

## 3. Same-site gate before interaction consensus

First establish sequence-based residue correspondences, then perform a protein-
based fit of the relevant domain/pocket into the project reference frame. Store
the atom/residue pairs, transform, alignment coverage, residuals and outliers.
Do not align ligands to each other to prove that they occupy the same pocket.

Assess site agreement with all of:

- Mapped pocket-residue correspondence and comparable observed coverage.
- Local protein alignment quality, with core-fit and pocket displacements shown
  separately when a genuine conformational change is present.
- Ligand location in the aligned receptor frame: pocket/envelope overlap and
  separation, accounting for ligand size and subpocket occupancy.
- Chain/assembly context and whether the ligand belongs to the reference site,
  an adjacent site, an allosteric site or an unresolved interface.

Three shared contact residues alone are insufficient to admit a complex. No
single fixed RMSD/distance threshold is declared universal in this contract.
Gate policies must be explicit, versioned and validated against same-site and
wrong-site controls. Candidates near a boundary remain pending review. A genuine
alternative receptor conformation can form a separate site-state branch, not an
averaged pocket that no experimental structure exhibits.

## 4. Quality gate and ligand chemistry

Report experimental method/resolution, ligand occupancy, alternate states,
heavy-atom identity/completeness, available density-fit metrics, pocket omissions,
mutations, protonation/tautomer uncertainty and local model-quality concerns.
Missing density metrics mean unassessed, not passed. A resolution cutoff alone
cannot establish reference quality. Record unsupported checks explicitly.

Small molecules, peptides and macrocycles are discovered and counted. Polymer
ligands require complete validated chemical graphs, including nonstandard
residues, stereochemistry, termini and covalent/cyclization links. Sequence
similarity is diagnostic only and cannot substitute for chemical/3D similarity
or certify a prepared search query. Partial evidence may remain visible but is
ineligible for a hard constraint dependent on the missing information.

## 5. Independent consensus model and shape templates

Only identity/site/quality-admitted complex instances can enter consensus.
Cluster compatible receptor states and ligand binding modes before pooling.
Show stratified evidence across ligand classes; class membership alone neither
proves compatibility nor requires incompatible modes to be combined.

The pocket model uses stable pocket-anchor IDs independent of any ligand's
feature-array index. Each anchor records interaction class, mapped receptor
residues, spatial distribution, direction distribution, uncertainty, supporting
PDB/ligand/contact IDs, eligible denominator and independent chemotype support.
Coordinates are distributions/tolerances derived after validated alignment,
not simple averages of unrelated contacts or conflicting modes.

For every anchor, distinguish observation, observed absence and unobservable
data. Missing residues, ambiguous waters and unavailable atom/angle evidence
must not become negative observations. Count distinct structures and also report
chemotype-weighted support, preventing replicated chains or one chemical series
from dominating. Freeze the weighting and denominator definitions in the design.

LLM proposals must cite actual admissible evidence IDs. Deterministic validation
enforces a configured support floor and minimum independent support for automatic
mandatory proposals. A 1/12 observation in the relevant comparable cohort cannot
be upgraded to mandatory by LLM plausibility. High frequency still does not prove
energetic necessity. Evidence-poor cases abstain or remain optional; any manually
specified exceptional requirement is separately labeled with its justification.

Select shape templates independently for quality, chemical/shape diversity and
coverage of compatible modes. An interaction-typical ligand need not be a shape
template; an unusual shape need not define the pocket's mandatory anchors.
Maintain template-to-site-state compatibility and never form an impossible union
of constraints across mutually exclusive modes.

## 6. Constraints, scoring and adaptive thresholds

Support mandatory anchors, alternative groups, weighted optional anchors, hard
exclusion regions and soft exclusion penalties. Distinguish physical collisions
from evidence-based chemical undesirability. Absence of a contact does not imply
an exclusion region. Antitarget selectivity requires comparative evidence and
its own cohort; it cannot be inferred merely by placing a forbidden sphere.

All requirements/rewards/penalties are evaluated on the same physical pose.
Do not combine one pose's Gaussian score with another pose's anchor rewards.
Hard constraints are eligibility rules and cannot be offset by optional rewards.
Normalize scoring terms and prevent duplicate contacts sharing one feature from
inflating the score. Retain the term breakdown and model/template provenance.

Derive coarse eligibility baselines from validated reference/known-ligand
distributions by compatible ligand/mode branch, with explicit uncertainty and
margins. Do not combine small-molecule and peptide size ranges indiscriminately.
A user-facing permissiveness control may choose a validated parameter schedule;
it must never relax target/site/provenance checks or modify mandatory chemistry
silently. Throughput tuning is constrained by independent positive retention.
If recall and capacity cannot both be met, report the conflict and use batching
or a revised design rather than quietly discarding known positives.

## 7. Quality panel and full-library reporting

Separate structural preparation controls from independent activity controls.
PDB binding poses are not automatically potency labels. Activity records require
target/organism/assay/units/relation and data-quality checks. Keep holdout
chemotypes separate from reference selection and threshold fitting where possible.
Prepare controls through the same chemistry/conformer pipeline as library inputs;
crystal-coordinate self-recovery alone is insufficient.

Controls absent from the library are tagged validation inputs, not counted as
library candidates or inserted silently into the catalog. Report unique-molecule
retention at every gate, cumulative retention, denominators, preparation failures,
lost-positive IDs/reasons and uncertainty. Do not count an unavailable control as
a successful pass. Record diversity and throughput alongside this minimum metric.

Full-library reporting distinguishes conformers, generated poses, Gaussian-scored
poses and molecules. Multi-template screening merges molecular membership across
valid branches while retaining the actual pose, template, site-state and anchor
evidence. Union across references never means simultaneous anchor satisfaction.

## 8. Chat review checkpoints

1. Show resolved target and reference complex/site; clarify genuine ambiguity.
2. Show discovery/admission counts and rejected/pending reasons. Provide reference
   overlays, residue correspondence and alignment/contact evidence for inspection.
3. Show mode-specific consensus frequency tables and separate shape-template
   coverage; accept natural-language edits or explicit adoption.
4. Show validation-panel retention, frozen rules and resource expectations.
5. Execute the requested full library when resources are available, without a
   mandatory small-library pilot. Quality-control panels are separate checks,
   not a substitute for the user's requested full-library run.
6. Show funnel metrics, lost controls and candidate combinations for human review.

Do not ask the user to hand-correct routine identifiers or write code. LLM text
explains measured evidence; deterministic code establishes membership and counts.

## 9. Acceptance cases and implementation order

Identity/site admission tests must cover: a title-only false match; a matching
target chain with ligand bound to another protein; paralog/antitarget; wrong
organism; pocket mutant; truncated but valid domain; inserted/relabelled residue;
same target but different site; neighboring pocket with three shared contacts;
packing-only contact; multi-chain biological site; alternate states; duplicate
crystal copies; and a genuine same-pocket conformational change.

Consensus tests cover missing-as-unobservable denominators, repeated-chemotype
bias, incompatible modes, 1/12 mandatory rejection, source citation mismatch and
same-feature double counting. Reference/template edits must remain independent.
Scoring tests cover same-pose terms, hard-gate precedence and distinct exclusions.
Validation tests cover held-out positives, unavailable chemistry and control/library
count separation. Generic fixtures use no MDM2-specific IDs, residues or cutoffs.

Implementation order: reference registration -> chain-level identity gate ->
protein-aligned site gate -> quality/cohort table -> independent pocket-anchor
schema -> consensus evidence checks -> template selection -> scored constraints ->
control panel/adaptive schedule -> multi-reference full-library execution and chat.
MDM2 supplies real positive/ambiguous cases; generic negative fixtures and another
structurally different target are needed before claiming cross-target validation.

## Current implementation gap

E051/E052 provide metadata discovery, partial sequence mapping, exploratory
contact overlap, small-molecule diversity and single-query screening. They do
not yet satisfy this contract. In particular, E052's 68 organic ligands and 46
polymer sequence/link variants are exploratory proposals, not E053-accepted
same-pocket consensus members. Preserve the old reports and issue a new cohort
report after the stronger gates are implemented.
