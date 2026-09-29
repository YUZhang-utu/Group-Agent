# E092: evidence-gated pocket states and coverage representation

Design discussion, not an implemented change. User redirected attention while
workstation boundary routing runs; do not interrupt or rerun that process.

Local reinspection of E085 MDM2: 2,211 pairs; global 1-IoU dominates the current
max in 1,334 pairs, half-local in 877, no exact ties. Correlation 0.756996.
Global 10/50/90 percentiles: 0.13465/0.18988/0.23891. Half-local:
0.12282/0.18418/0.24609. The max is deterministic, continuous but nonsmooth,
and suppresses the smaller component, not a random coin flip. Marginal overlap
alone is not evidence of unstable selection; perturbation tests must assess it.

Same-entry different-chain comparisons: 9 pairs; max 0.232905, median 0.127902,
p95 0.217212. User reports max 0.214 for different-entry same-ligand pairs; the
exact same-ligand matching definition and those pairs have not been verified in
this turn. These controls include genuine biological/crystallographic variation
as well as representation uncertainty. Maxima are not established hard noise
floors. Any replacement distance requires recalibrating these distributions.

## Proposed changes

1. Replace winner-takes-larger fusion with a versioned nonnegative weighted sum
   of zero-preserving normalized global, local and chemical dissimilarities.
   Prefer robust frozen calibration scales over per-batch ranks/min-max. Do not
   force equal marginal distributions: that can amplify low-signal channels.
   Correlated global/local channels require redundancy and weight sensitivity
   checks. Local material occlusion also receives an explicit diagnostic flag
   so averaging does not hide a functionally important small subpocket change.

2. Separate numerical reproducibility controls (grid origin, spacing, alignment,
   atom/alternate-location treatment) from near-repeat experimental controls
   (same entry chains, same ligand across entries). Match target sequence,
   construct/site, ligand identity and structural quality. Report distributions,
   effective independent samples and uncertainty, rather than one maximum.
   Do not enforce a cutoff lower bound from nine correlated chain pairs.

3. Treat a robust nontrivial plateau as a prerequisite for publishing discrete
   states, not a theorem deciding physical continuity. Examine exact dendrogram
   merge intervals, not only a coarse cutoff grid. Exclude trivial k=1 and k=n
   plateaus. Require stable partitions under PDB-entry grouped resampling,
   representation perturbations and reasonable descriptor weights. Count
   stability alone is insufficient. Support should exceed near-repeat variation
   and be examined against a suitable continuous/null baseline if one can be
   justified. Plateau widths/criteria must be fixed before judging MDM2, not
   chosen to recover a desired state count.

4. Reporting modes: discrete_states_supported; no_robust_partition_resolved;
   insufficient_evidence. The latter two use coverage representatives. Do not
   equate no plateau with proof that no metastable/discrete states exist.
   Finite PDB observations are biased samples, not equilibrium populations.

5. Coverage mode reports representative count, weighted coverage within an
   explicit tolerance, median/p95/max residual distance, and uncovered structures.
   Start from a quality-qualified weighted medoid. Select additions by reduction
   in uncovered experimental support / weighted residual error; compare with
   farthest-first rather than letting a lone extreme artifact dominate selection.
   Weight each PDB entry once across chains, with further redundancy assessment.
   Preserve low-support but credible unusual structures in explicit review;
   frequency alone must not delete possible cryptic or ligand-specific pockets.
   Stop at a declared resolution/coverage or computation budget, not a desired k
   or assumed elbow. Representative count does not mean number of states.

6. Downstream consensus differs by reporting mode. Discrete mode permits
   within-state consensus after adoption. Coverage mode computes contacts per
   representative/local neighborhood, preserving alternatives; do not call
   nearest-representative bins independent physical states or intersect all
   representative contacts into universal mandatory interactions.

## Current MDM2 implication

E085 smooth defaults yield 28/11/3/1 groups at .20/.25/.30/.35 and k=1 thereafter.
The k=1 interval is trivial. Existing sparse sweeps do not establish a stable
nontrivial partition; nor can they prove that no such interval exists. The old
implementation is weighted average linkage cut by distance, not fixed-k, but its
reporting lacked this evidence gate. Provisionally report no robust partition
resolved and a coverage representation, pending metric/control recalibration.

## Primary methodological sources

- U. von Luxburg, Clustering Stability: An Overview, arXiv:1007.1075:
  https://arxiv.org/abs/1007.1075
- U. von Luxburg et al., Clustering: Science or Art?, PMLR 27 (2012):
  https://proceedings.mlr.press/v27/luxburg12a.html

These sources support caution about interpreting stability as an automatic
physical cluster-count criterion. The proposed molecular evidence gate and
coverage workflow are project design choices, not claims directly validated
by those papers. No new pocket metric or clustering was implemented this turn.
