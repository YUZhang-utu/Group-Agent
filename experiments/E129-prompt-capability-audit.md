# E129: prompt capability audit

Date: 2026-11-11. Question: can Chat dispatch full-library molecular searches and
AF3 protein predictions including homodimers and heterodimers?

Protocol: trace the supplied Chat tool schema, workflow validation, input
generation, queue worker and engine calls. Exercise homodimer and heterodimer
requests with and without ligand through the real Chat worker using a scripted
planner and a synthetic external runner. Verify exact chain sequences, distinct
IDs, seeds, one task, retained input and completed execution receipts. Exercise
natural-language tool dispatch for generic target search with a scripted provider;
run existing query/full-catalog and provenance suites. No real LLM/GPU/production
library is accessed. These tests verify plumbing, not scientific predictions or
live-provider interpretation accuracy.

No production behavior change is planned unless a reproducible defect is found.

## Results

Confirmed that target_search continued planning after a successful queue receipt;
the one-response scripted provider reproduced a second unwanted model call.
Added deterministic submission acknowledgement and same-turn pending/mutation
guards. The first dispatch test also exposed a test helper naming collision,
corrected before reproducing the product behavior.

110 tests passed across the audit, generic workflows, prompt execution, container,
full-library, preselection, guided search, domain-agent and task-reliability suites.
One additional pending-operation guard test passed separately (111 distinct
passing cases total). The English guard passed 739 files. Four AF3
composition tests used real worker subprocesses with synthetic planner/runtime.
No live provider or GPU/production-library run was performed. Scientific inference
generality and screening speed/recall remain workstation validation questions.
