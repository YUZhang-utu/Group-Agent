# E084: PDB pocket states before consensus

Hypothesis: aligned local cavity overlap with a modest chemical-field contribution
can consolidate incidental coordinate changes while preserving materially different
pocket states. Cluster counts are observations, not fixed expectations.

Protocol: implement shared-grid cavity masks and residue-template chemical fields,
average-linkage clustering with a permissive exploratory cutoff and sensitivity
counts, quality-aware weighted medoids and explicit adoption receipts. Source
hashes and target/site identity must remain verifiable. PDB multiplicity is support,
never equilibrium occupancy. Do not force a maximum cluster count or drop small
clusters. Do not admit missing sidechains as apparent cavity expansion.

Tests before rollout: rigid transform invariance; unchanged/small perturbation
merging; material occlusion separation; changed chemical annotation detection;
singleton preservation; duplicate support control; changed-source adoption refusal;
unadopted consensus refusal; state-restricted contact input. Use existing local
MDM2 experimental structures for exploratory acceptance if available. Report that
run separately from synthetic engineering checks. No docking recall/affinity claim.

Exploratory result: the final weighted-linkage implementation admitted 68 chains
from 60 PDB entries in the existing E052 collection, yielding 6 clusters at 0.40,
12 at 0.30 and 3 at 0.50. 180 candidate chains/entries remain held for review,
including unmapped target entries and unsupported assembly contexts. This is not
a count of 180 unique rejected PDB structures. Initial 0.30/default full-residue
quality checks overpartitioned the sample and held remote missing atoms. The
permissive cutoff and pocket-local atom scope were selected during development.
Same-PDB chain weighting changed the provisional unweighted 5-cluster result to
6; retain this distinction. No user adoption or downstream docking has run.
