# E124: Core workflow integration

Date: 2026-10-29
Status: implemented; offline integration checks passed

Question: Can Chat expose existing KRAS export and AF3 batch evidence without
rerunning scientific calculations, while supporting explicitly requested AF3
homomers and keeping block routing out of the default development path?

Protocol:
1. Add session-owned immutable report attachments with explicit import roots.
2. Validate report schemas, completion counters and reported artifact hashes;
   never infer that referenced remote payloads have been verified locally.
3. Exercise deterministic commands, HTTP access, domain-tool contracts and UI
   integration using offline fixtures. Test foreign paths/sessions, changed
   sources, inconsistent exports and AF3 host checks versus predictions.
4. Test AF3 protein-copy input generation, ligand chain IDs, invalid copy counts
   and preservation of existing monomer behavior with offline services.
5. Run relevant regression tests and the repository English-content check.

Acceptance: no compute jobs on report import/read; no cross-session result
disclosure; no host check described as inference; explicit same-complex ligand
semantics; old workflows remain usable. No GPU inference, library search or
binding validation is part of this experiment. Results are integration evidence,
not scientific validation or a performance benchmark.

Results: confirmatory for the bounded integration requirements. Actual E123
returned-report import matched 254364/983854 exports with zero jobs launched.
AF3 copy counts 1/2/8 preserved constructs, seeds, literal SMILES and unique IDs.
Initial homomer fixtures omitted the SMILES from the user request and were correctly
rejected; fixtures were corrected without relaxing the literal-SMILES guard.
The initial test invocation needed explicit PYTHONPATH=src in this local environment.

Related audit fix: AF3 batch launch/artifact-validation failures now produce failed
receipts, total job counts and separate process/validation return codes. This was
covered with missing model, invalid confidence JSON and launch-error fixtures.

Validation records: 31 plan/AF3 tests; 56 result/domain/full-library/structure/budget
tests; 30 batch/result tests before adding the pending-job case. Suites overlap.
Final checks are recorded in research_log.md. Frontend syntax and HTTP API were
checked; live browser interaction, provider planning and school GPU execution were
not validated. No new scientific screening or inference was performed.
