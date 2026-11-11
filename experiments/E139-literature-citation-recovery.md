# E139: resolve literature citation contracts

User workstation feedback confirms public document listing and three retrieved
passages on pages 6 and 9, followed by synthesis citation validation failure.
The actual model answer is unavailable. Missing citations, grouped citations and
unknown IDs cannot be distinguished from the old combined error. Do not assume
which one occurred or that scientific summarization succeeded.

Protocol: expose an explicit evidence-to-document/page reference table, accept
unambiguous grouped evidence citations and exact retrieved passage/page citations,
reject unknown IDs and unread pages, and allow one bounded citation correction
using the same public evidence with no tools. Never silently append sources to an
uncited answer. Separate missing and unknown citation errors in safe diagnostics.
Test actual serialized provider requests, repeated failure, transport failures,
time budget, metadata-only sources, unknown citations and privacy isolation.

## Confirmatory local results (2026-11-11)

165 tests passed in 73.31 seconds across literature recovery, knowledge privacy,
domain-agent behavior, streaming/uploads, prompt workflows, project guidance,
scientific P0 and chat-agent behavior. The captured native recovery transport also
passes the two-response case: citation-free draft, then a cited answer; both calls
exclude the private-library sentinel and neither dispatches tools or calculations.
Grouped/whitespace evidence references and exact passage/page aliases resolve;
unknown IDs, unread pages and metadata-only evidence are rejected. Repeated invalid
citations stop after one correction; transport failure is not retried; exhausted
time budget prevents another request. English-content and whitespace checks pass.

These are synthetic provider tests, not live-model or scientific-accuracy evidence.
The user's PDF and cloud provider require a new workstation acceptance run. The
previous workstation error cannot reveal which unsupported citation form was used.
