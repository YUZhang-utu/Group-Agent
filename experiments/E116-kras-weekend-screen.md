# E116 KRAS Switch-II native-state weekend screen

2026-10-23, prospective exploratory protocol. User requests KRAS G12D and G12V
noncovalent Switch-II full-library searches on a 24-core, 64-GB workstation.
Do not claim mutant selectivity, biological binding, or a Monday completion date.

Prepare experimentally observed, GDP-bound native ligand/receptor pairs separately
for each allele. Primary references under review: 7RPZ/6IC for G12D and
12FF/A1DB0 for G12V. Audit sequence against human KRAS P01116, including engineered
background substitutions not listed in deposited mutation descriptions. Do not
silently mutate coordinates or transfer the old MDM2 Ile61 rule. Keep GDP/Mg
as receptor-context atoms; exclude crystallization additives and waters.

Use existing AIDD feature generation, directional one-to-one assignment and
bounded seed generation (512 pair seeds, six axial samples, 2 A tolerance).
Every catalog conformer is eligible for evaluation; chunking is execution, not
top-block selection. Exact finite seed-set coverage is not exhaustive pose space.
All tests must be satisfied on the same seed and paired native receptor.

Exploratory KRAS baseline: at least three distinct native-contact residues and
at least one matched native donor/acceptor feature, directional match >=0.5;
whole-ligand Gaussian shape >=0.50 and color >=0.30; heavy-atom minimum distance
>=2 A; contact-envelope coverage >=0.65 and each supported spatial region >=0.50.
These numerical floors are inherited diagnostic settings, not validated KRAS
activity cutoffs. Whole color remains baseline; E115 grouped coverage must not be
silently adopted without its returned evidence. Preserve counts before/after color
and geometric near misses so zero hits can be diagnosed. No lowering to force hits.

Before library search, recover each reference from a fixed rigid perturbation with
the actual seed generator and all gates. Abort on a failed positive self-control.
This checks geometry and search recovery, not enrichment or sensitivity to other
chemotypes. Audit reference provenance, atom completeness and engineered sequence.

Run a deterministic spread sample for timing before full search. Record wall time,
workers, sample size, library size, survival counts and linear extrapolation with
explicit uncertainty. Default 20 workers shared by both tasks, single-thread BLAS.
Execute resumable bounded chunks, interleaving targets to avoid starving one target.
Checkpoint each completed chunk, seal input/code identities, preserve all passing
conformers and a bounded diagnostic near-miss export. Resume only matching inputs.
Budget exhaustion produces partial status and explicit coverage, never full status.

Verification: synthetic shape/contact/identity tests and resume/seal failure tests;
native reference recovery from downloaded structures locally where feasible.
Actual full-library timings/results can only be established on the workstation.

Pre-library amendment after the user's explicit request for effective hard
screening: require three distinct matched reference feature indices in addition
to three residues, a native polar contact to Glu62 or Asp69, and a separate native
contact to His95 or Tyr96. This prevents arbitrary three-contact matches from
meeting the rule. Derived from the native reference inventory, before viewing any
KRAS library scores; still exploratory and chemotype-biased. Compare sequence to
KRAS4B P01116-2, not the alternative canonical KRAS4A C-terminal isoform.
