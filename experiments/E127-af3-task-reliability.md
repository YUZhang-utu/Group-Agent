# E127: AF3 task planning and failure visibility

## Protocol, recorded before implementation

User supplied a real CIP2A literal-sequence/SMILES request, seeds 47..51, and a
trace showing task 0a928af2a73a479b queued, two subsequent malformed conversational
responses, then a generic background failure notification. The background error
is not present in the supplied trace; request it without resubmitting the job.

Hypotheses: a deterministic queue acknowledgement avoids additional model failures
after a successful submission; bounded structured-plan repair can recover transient
format errors before any workflow is launched; persistent phase-specific failure
receipts make pre-execution errors visible without a missing execution report.

Keep original sequence/SMILES/seed provenance and validated scientific adapters.
Do not silently replace chemistry, truncate the construct, retry inference, or
submit another workstation job. Test the supplied input through preparation with
the existing offline/local services and verify exact seed, chain and ligand values.
Test malformed-plan recovery/exhaustion, no repeated mutations, background failure
notifications and failure evidence accessible from the existing status/report tools.
Treat local engineering fixtures as confirmatory only for those contracts. Live
runtime success and the user's actual failure cause require returned evidence.

## Returned evidence and results

The user subsequently returned `Invalid plan fields`, confirming failure at the
top-level workflow-plan validator before any scientific report. The exact invalid
response remains unavailable. Added an explicit schema, actionable field feedback,
two bounded repair retries, deterministic queue receipts and phase-aware failure
logs/notifications. Added discoverable input paths, compact step cards and a test
mode banner. Preserved labelled full sequence/construct and explicit seeds.

Confirmatory local checks: 119 core/planner/AF3/Chat cases passed (33.68 seconds),
67 additional scientific-adapter/chain/result cases passed (6.46 seconds), then
27 reliability/HTTP cases passed (8.24 seconds), including one new recovery test.
The latter run overlaps the first; distinct total 187. Four Markdown checks and
JavaScript syntax passed. The actual queued worker prepared the supplied CIP2A
JSON after schema repair, with exact chemistry, 905 residues and five seeds, then
stopped at the compute gate. A failed planning task resumed under its original ID
and retained prepared input artifacts. No GPU or live model used; no scientific
accuracy/throughput conclusion. Guide: to_human/E127_TASK_RELIABILITY.md.
