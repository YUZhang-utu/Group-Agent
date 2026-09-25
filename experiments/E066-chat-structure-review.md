# E066: Chat structure review acceptance protocol

Question: can existing scientific outputs be reviewed through Chat with preserved
structure identity/alignment, and can a literal ligand SMILES reach the configured
AF3 input without invented chemistry?

Before execution: retain original jobs and sealed artifacts. New predictions use
fresh tasks. Use temporary local fixture projects for adapter tests. No remote
screening, AF3 or GUI execution is part of local unit validation.

Local checks: strict SMILES/chirality/charge preservation; rejection of invented or
disconnected chemistry; confidence explanations use actual fields; session/project
ownership; recorded rigid transforms; source hash mismatch; typed contact ledger
filtering; allowlisted PyMOL commands; downloadable receipt artifact paths; existing
Chat regressions; English-only maintained content.

Workstation exploratory acceptance: follow E066_CHAT_STRUCTURE_REVIEW.md. Reproduce
one previously successful protein/SMILES AF3 input and run it; inspect summary/model
identity. Open a completed consensus, inspect saved alignment, export contact CSV,
PNG and PSE, and compare selected distances manually. Test one live provider with
the documented English prompts. Record versions, task IDs and failures.

Success claims are separate: passing local tests establishes implementation
contracts only. Real AF3 compatibility, desktop control, scientific contact
validity and screening performance are separate outcomes. Macrocycle clustering
is deferred and is not an E066 acceptance requirement.

Local confirmatory result (2026-09-24): full regression 515 passed, 4 skipped;
English-content guard passed; JavaScript syntax check passed. RDKit 2026.3.6 was
loaded from a workspace-only test dependency directory. GUI, live provider and
real AF3 workstation tests were not executed. No measured speed or recall claim.
