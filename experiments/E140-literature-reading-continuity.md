# E140: literature follow-up context and page completion

Workstation feedback reports complete indexed page-1 retrieval, partial Methods
search results, and a response claiming the prior summary was unavailable. This is
user-returned acceptance evidence, not independently inspected remote output.

Code inspection confirms that the planner truncates each history message to 1500
characters and the recovery synthesizer receives no previous answer. Search also
accepts document_id/page arguments without applying them as filters.

Protocol: retain a bounded complete previous same-document answer from this
conversation as comparison context, never fresh scientific evidence; carry it into
normal planning and recovery. Apply document/page search filters. Automatically
follow pagination only for a selected explicit page within existing turn budgets.
Separate complete contiguous page reads from search fragments and partial reads;
prompt the planner once to expand relevant Methods hits before an early answer.
Expose the actual coverage even if the model stops early. Do not claim full paper
or section comprehension, reuse old evidence IDs, access internal knowledge, or
start scientific calculations. Test scope, paging holes/failure, comparison
context, budget exhaustion, privacy and unchanged non-literature dispatch.

## Confirmatory results (2026-11-12)

179 related Python tests passed in 71.16 seconds, including 14 new reading cases.
The synthetic extraction fixture exercises the real local index and domain loop:
an explicitly selected multi-part page completes with only two planner calls;
pagination receipts preserve offsets and no scientific task is created. Document
and page search filters, empty/out-of-scope results, invalid filters, a bounded
Methods nudge, gaps, stalled/failed cursors and turn-budget partial coverage pass.
The actual serialized recovery request contains the full previous summary tail
past the old 1500-character cutoff, excludes the private-library sentinel, and
does not promote previous-turn E-number citations into current evidence. Other
conversations are not automatically incorporated into comparison context.

English-content and whitespace checks pass. No live cloud-provider invocation,
remote trace access, full-PDF scientific verification, or production calculation
was performed. Explicit-page completeness is an indexed-text receipt; full-section
selection/comprehension remains model behavior requiring workstation acceptance.
