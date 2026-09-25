# Deferred macrocycle block design

Status: design proposal. Production clustering and incremental assignment are not
implemented. Chat/AF3/PyMOL acceptance and the seed audit have priority.
Existing offline block replay utilities evaluate supplied block labels only.

## Identity and records

Keep molecule_id, the exact molecule name, conformer_id (conf1/conf2/conf3), source
MOL2 record and content hash. User-defined cis variants retain separate molecule
IDs. Conformers can occupy different blocks; final screening budgets deduplicate
by molecule. Never silently rename, merge stereoisomers or reconstruct missing
hydrogen/stereochemical metadata as if verified.

## Hierarchy and descriptor versions

1. Extract the designated macrocycle backbone from building-block connectivity,
   with an auditable atom map. Store ring size and connectivity class. Fused,
   bridged and multiple-ring cases need an explicit definition or an unassigned
   queue, not an arbitrary smallest-ring choice.
2. Preserve the ordered omega cis/trans pattern within each mapped backbone as a
   hard grouping constraint. Near-boundary/atypical omega gets an explicit state.
   Use an invariant cyclic correspondence without merging enantiomers or changing
   chemically inequivalent building-block order.
3. Record a backbone flexibility descriptor: eligible torsion count, amide count,
   N-methylation, constrained/proline-like units and crosslinks. A rotatable-bond
   count is not the true number of independent degrees of freedom in a closed ring.
   Test its value as a descriptor before adding another hard partition.
4. Compare backbone torsions using sin/cos and structural alignment. Benchmark phi
   sign and residue-conditioned phi/psi basin labels as alternate partitions.
   A geometric boundary is not evidence of an energy barrier; three conformers
   cannot establish transition kinetics or thermal populations. Include boundary
   handling and report fragmentation caused by extra hard partitions.
5. Add side-chain descriptors using building-block identity and position, typed
   pharmacophore points, attachment-to-feature vectors, distance from backbone,
   radial/tangential/normal projections, charge and donor/acceptor/aromatic/
   hydrophobic categories. Chemical type is part of the representation, so a donor
   and acceptor at the same coordinate remain different features. Preserve atom
   mapping and chemistry provenance for every derived point.

## Coordinate frames

A centroid plus principal axes is a useful initial frame, but eigenvector signs
are ambiguous and near-equal eigenvalues can swap axes. Store the construction
rule, handedness, eigenvalue gaps and quality flag. Use chemically mapped backbone
landmarks or reference alignment to resolve ambiguity, with invariant pairwise
feature distances as a fallback. Never permit mirror reflection during alignment.
For puckered rings, record best-fit plane residuals and local normals; do not treat
a global plane normal as a reliable universal orientation. A local backbone frame
at each attachment is an additional benchmark variant.

## Capacity and incremental insertion

Target 10k, 20k or 30k conformer records per block, not exactly equal populations at
the cost of chemical coherence. Within hard groups use capacity-constrained splits;
allow small residual blocks and an outlier queue. Record actual sizes and unique
molecule counts. Never fill a short block with another ring size/cis state.

Freeze descriptor schema, atom mapping, cluster prototypes, metric, hard groups,
assignment radius, capacity rule and version. New molecules compute descriptors
and join eligible nearby blocks; full blocks split into versioned children or an
overflow queue. Out-of-domain rows remain unassigned for review. Existing membership
and search receipts do not change silently. Adding molecules invalidates any claim
that an old ranking is the complete ranking of the enlarged library.

## Adaptive exploration experiment

For every block, sample 100 or 500 unique molecules and score their relevant stored
conformers with an identical search protocol. Use a fixed upper-tail fraction, or
compare fixed top-k only at the same sample size. Report top-k mean, maximum,
effective k, unique sample size and uncertainty. Five observations out of 100 and
five out of 500 estimate different tails and must not be compared as identical.

Prioritize promising blocks, but reserve a prespecified random exploration budget
for low-ranked/unvisited blocks. A missed rare high scorer remains possible; a low
sample tail is not proof that a block contains no useful molecules. Track costs at
the molecule/conformer/template evaluation level and cache repeat molecules.

Use the completed full-library ranking as a retrospective scoring reference, not
biological ground truth. Freeze a tuning subset and evaluate on held-out molecules,
chemistry strata and repeated random seeds. Compare uniform exploration, current
retrieval, backbone-only, backbone plus chemistry and full side-chain-geometry
variants at equal scoring budgets. Report molecule recall@100000, retained
contact-score distribution, unique backbone diversity, uncovered rare families,
scored fraction, wall time, indexing cost and incremental assignment stability.
If the reference ranking came from an ANN pool, restrict recall claims to that pool.

First decision: does chemistry-aware blocking improve held-out top-molecule recall
at the same compute budget? Only then add live adaptive scheduling.
