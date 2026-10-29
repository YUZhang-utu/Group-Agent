# E125: General prompt-driven AF3 and receptor screening

Date: 2026-10-29. Protocol recorded before tests.

Implement direct literal sequences, heteromeric/homomeric protein chains and
explicit seeds or deterministic seed counts in the existing AF3 adapter. Preserve
ligand identities and existing monomer plans. Add an explicit receptor/reference
workflow that builds coordinate-backed queries for the existing full-library
engine without target-name dispatch. Bound local file access to configured roots,
seal inputs and use the existing Chat queue, compute gate and progress receipts.

Test literal sequence provenance, chain IDs, seeds, invalid dependencies and
runtime invocation offline. Test arbitrary receptor labels, contact-rule mapping,
native-reference acceptance, mismatched coordinates, preserved cofactors, file
ownership, compute gating, full scan dispatch and result paths. Use synthetic
fixtures and cached structures; no production full-library scan or GPU inference.

Success requires executable generic adapters, not only external report import.
No claim of universal target recall or biological validation follows from these
software tests. Publish only task-related changes to the existing GitHub branch.

Results: confirmatory for the specified software contracts. Initial generic AF3,
query and prompt/structure suite: 66 passed. The next regression run exposed an
empty pre-download evidence-map guard on PDB-only input and a prompt-size guard
inconsistent with the larger default; both were corrected. The Chat request test
was updated to assert preservation of original user evidence, which prevents
router-generated sequences from being treated as user-supplied sequences.
Affected follow-up suite: 65 passed. Final release validation is recorded in
research_log.md. Dispatch fixtures simulate external inference/search; they do
not establish school GPU success or production-library throughput.

Release acceptance: 155 tests passed in 40.85 seconds in the isolated delivery
checkout, including both E124/E125 and existing AF3/Chat/full-library/container
regressions. English guard passed for 719 maintained files; JavaScript syntax,
scoped diff whitespace and Chat CLI help checks passed. No production compute.
