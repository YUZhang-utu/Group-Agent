# E077 full source audit

Protocol before execution: inspect every CSV row and every MOL2 record in the two
supplied shards. Use a disk-backed name index with duplicate tracking, strict
SMILES/MOL2 identity comparisons, conformer suffix accounting and full input hashes.
Recognize only the observed residue alphabet A,W,L,F,V,P,G with optional d and
nme/NMe. Validate labels against graph-derived peptide units using stereochemical
amino-acid reference fragments, not string position or a universal D=R rule.
Recover a main ring only when units form one unambiguous closed peptide cycle.
Keep unknown names, multiple cycles, identity mismatches and missing records
explicit. Check all conformers, not only the first conformer of each molecule.

This is a full audit of the supplied sources, not of unseen workstation shards.
No production registry or screening data will be rewritten. Report count closure,
all issue categories and examples before recommending partition strata. The
experiment is exploratory; it does not validate energy barriers or recall.

## Implementation correction during execution

The first attempt was interrupted after the MOL2 residue comparison exposed a
representation mismatch: cleavage retained MOL2 no-implicit-H flags and produced
uncapped radical fragments. Whole-molecule identity had not failed. Explicitly
cap the analysis-only peptide endpoints, then compare to reference fragments.
An explicit-H/no-implicit-H regression and 200 real MOL2 records, including
proline-containing cases, passed. Preserve the first attempt as superseded;
only `data/e077-full-audit-v2` may support the final source conclusions.

## Final results

Full supplied-source scan completed in 420.83 seconds, four workers. All 100000
CSV molecules and 299999 MOL2 conformers pass the implemented identity and
residue-sequence checks. Main-ring size is 18 throughout. Duplicate CSV names,
duplicate conformer names and CSV names without MOL2 are all zero. Independent
canonical-SMILES grouping also finds zero duplicate chemical-structure groups.
99999 molecules have three conformers; one has conf1/conf2 only:
`c--A-Wnme-PdFnme-dW-dLnme-c`.

27658 names have multiple observed omega patterns; 16653 have a definite cis/trans
change at a matched bond after excluding ambiguous cyclic alignments. Boundary
angles occur in 16171 conformers; repeated-sequence cyclic alignment is ambiguous
for 18 conformers. These are descriptor/coverage considerations, not chemical
identity failures. Do not infer amide state from the residue D-prefix.

28 focused tests passed. Exploratory conclusion: source-wide peptide mapping is
feasible here and eliminates the old generic-ring ambiguity, but future blocking
must preserve per-conformer omega and cyclic symmetry. Whole workstation library
coverage and ranking enrichment have not been measured.
