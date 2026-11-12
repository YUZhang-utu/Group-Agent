# E142: saved-review retrieval without repeated JSON planning

User reports a malformed model response after public knowledge access, with zero
documents/passages in the progress message. The remote trace and exact tool calls
are unavailable. Code inspection confirms that saved-review list/index receipts
do not count as text reading and still require another model-generated JSON step
before any text can reach the synthesis fallback. The old progress wording conflates
these stages. Do not claim the remote import or review index is empty.

Protocol: for explicit read-only saved-review summary requests, coordinate bounded
local list/index/read operations without model planning. Limit selected records,
characters, old-trace indexing and existing turn time/steps; recheck public sources
on each read. Synthesize once as text from actual retrieved summaries, with the
existing citation correction allowance. Keep computation, general requests and
mixed-action prompts on their existing route. Distinguish listing, indexing,
summary reading and final-answer failures in receipts. Test native serialized
transport, empty/stale sources, legacy indexing, pagination, caps and no dispatch.

## Confirmatory local validation (2026-11-12)

212 integrated Python tests passed in 85.63 seconds across saved-review flow,
literature reading/recovery, review memory, privacy, domain dispatch, guidance,
chat/streaming, prompt workflows and scientific P0. Five safe Markdown rendering
tests passed. The English-content and whitespace guards passed.

Native-path tests make any JSON planner invocation fail immediately, then inspect
the actual serialized text-synthesis request through a synthetic transport. The
saved-review request succeeds with one text call, actual read evidence citations,
no tools/JSON response constraint and no private sentinel. Old-review adoption,
empty/archived sources, indexing opt-out, direct review IDs, failure-stage messages,
partial coverage and hard read/index caps pass. No calculation is dispatched.

The remote failure cause and exact operations remain unverified without its trace.
Live cloud-provider reliability, actual-PDF scientific accuracy and workstation
browser acceptance are not established by these tests.
