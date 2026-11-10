# E126: local-first conversational workbench

## Protocol (recorded before implementation and checks)

Complete streaming answers, safe Markdown display and local file attachments.
Keep scientific payloads on the Agent host; provide the model with owned file
references, not automatic file bodies. Preserve the existing task queue and
input/code-bound execution receipts. Audit independent-ligand AF3 reuse and
full-library screening/structure handoff without launching production calculations.

Hypotheses: streaming can expose answer deltas and tool progress without executing
partial model output; uploaded structures can enter existing typed workflows
without moving the library to a provider; restart-safe AF3 receipts can avoid
repeating completed jobs without accepting changed inputs or missing outputs.

Verification: provider SSE fixtures (including fragmented JSON and invalid output),
authenticated HTTP streaming/upload tests, safe Markdown rendering checks,
session isolation and byte/hash preservation, existing Chat/domain/planner tests,
AF3 repeated-run and changed-input checks, and generic full-library regression.
Use small synthetic fixtures and existing local dependencies only. Record failures
and fixes. No GPU throughput, production recall or biological accuracy claim follows
from these checks. Inspect the actual browser when a supported surface is available.

Scientific audit: distinguish catalog coverage from finite pose sampling, exact
filters from heuristic gates, conformers from unique molecules, and original source
coordinates from posed structures. Report export readiness and missing adapters.

Results are confirmatory for the listed engineering contracts; scientific workflow
findings are exploratory code audits until validated on the configured workstation.

## Results, 2026-11-10

Implemented provider SSE/NDJSON answer streaming, typed tool progress, safe Markdown
subset and local bounded uploads. Added AF3 completed-job/shared-data recovery,
atomic reports and a bundle lock. Removed eager scientific imports from basic Chat
snapshots and avoided unchanged polling rerenders.

- Chat/domain/uploads/AF3: 57 tests passed in 12.53 seconds in the final run.
- Generic query/full-library/guided/results/campaign/planner/structure/budget:
  138 passed and two environment-dependent queue failures in the first run.
  Both failed tests passed in 13.71 seconds with the existing scientific venv;
  their child-worker fixtures overwrite PYTHONPATH, so root Python lacked RDKit.
  Total: 197 passing Python cases across these targeted suites, not one full-suite run.
- Markdown: four Node tests passed; JavaScript syntax and diff whitespace checked.
- Added a subprocess regression proving ordinary Chat snapshots do not load RDKit.
- Early new-test failures were fixture assertions: JSON-escaped Windows paths and
  an existing appended Agent trace. Assertions were corrected to preserve the
  production behavior. Slow native imports ultimately completed.
- Visual browser acceptance remains pending; browser access was unavailable and
  the isolated preview service failed to start. No live model or GPU job was run.

Confirmatory engineering result: actual streaming events precede completion,
invalid answer evidence clears previews, uploads preserve bytes and ownership,
and resume skips only verified completed jobs. Exploratory scientific audit:
generic search retains local candidates but lacks a direct MOL2 export adapter;
source offset indexing and target-specific rejection/recall calibration are the
next priorities. See the guide for the concrete handoff and validation proposals.

## Release follow-up

Added a portable temporary UI server plus a full HTTP upload/prompt/real-read-only-
tool/stream/persistence test. Final affected suite: 58 passed in 14.54 seconds;
four Markdown tests passed. Previous broader scientific tests remain applicable.
The scientific-venv rerun stalled in an existing case after the integration cases;
the validated root-Python/RDKit overlay completed. The preview server now starts,
but browser access/minimized-window recovery failed; no visual pass is claimed.
GitHub delivery was explicitly authorized. Installation and acceptance steps are
recorded in to_human/E126_UPDATE_AND_TEST.md.
