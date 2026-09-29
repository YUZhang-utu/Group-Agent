# E086: capacity and boundary fragmentation inspection

Protocol (2026-09-29): user requests inspection before regrouping. Preserve
cis/trans distinctions. A 20,000-conformer limit is not required for a logical
class; operational batching can be separate. Do not modify the frozen model.

Hypotheses: small pre-existing hard groups dominate fragmentation; boundary
patterns may create many low-population classes. Count original tree roots,
leaves and capacity-induced extra leaves exactly. Attribute root populations
to boundary-containing versus definite cis/trans patterns using one descriptor
per hard group. This is conditional on the already passed whole-model validator.
Do not rerun MOL2 chemistry or scan every descriptor. Do not infer a causal
boundary-only contribution from association: chirality and size also stratify.

Important code finding: boundary is the entire interval 30 < abs(omega) < 150
degrees, not a narrow tolerance at a cutoff. Never silently force all such
records into cis/trans or discard them. Angular and compatibility assessment
will follow the fragmentation counts before selecting a new assignment rule.

Validation: synthetic population with known extra capacity leaf and boundary
root; byte-identical input databases; local small existing model smoke run.
Whole-library counts remain pending workstation execution. The supplied user
receipt reports 24,663,736 described conformers in 16,531 leaves, capacity 20,000;
it does not report the number of original hard groups.

## Local check

Synthetic integrity/nonmutation test passed. Existing typed fixture gave
24 conformers, 10 leaves, 9 roots, 3 boundary classes containing 3 conformers.
These confirm diagnostic execution only, not the whole-library hypothesis.
