# E133: operational project direction and evidence-backed learning

Protocol (2026-11-11): implement versioned project missions, a rendered PROJECT.md,
bounded task context, evidence-linked feedback and persistent task cards. Validate
typed mission restrictions against actual plans before execution. Detect equivalent
scientific plans independently of prompt wording and retain reuse lineage; never
silently re-execute failed/uncertain work or reuse changed/missing artifacts.
Keep user changes separate from model proposals, observations from interpretations,
and requirement acceptance separate from successful task execution. Provide direct
project controls plus natural-language tools, and derive project experience from
actual task transitions and user-sourced feedback. Test aging context, conflict,
restart, revision, ownership, drift, deduplication, explicit replication and artifact
invalidation with real small local workflows. No live-provider/GPU quality claim.

Confirmatory software results: 192 integrated Python cases passed, including
14 new mission/learning cases. Four Markdown tests, JavaScript syntax and the
English-content guard (759 maintained files) passed. The regression set covers
project sharing/isolation, real worker execution, AF3 plan/input contracts,
streaming/uploads, target search, reusable candidates, dependency continuation,
report import and confidence parsing. Scripted provider tests exercise exact
feedback capture and delivery of pinned corrections to a later conversation.
HTTP tests exercise direct mission changes, stale revision rejection, project
isolation, revision-bound user acceptance and restart persistence.

Real local checks prepared ethanol, reused its unchanged artifact inventory across
conversations/restarts, rejected a changed report or missing SDF, and preserved
explicit independent repeats. A two-conformer library was searched and ranked,
then a changed export count blocked a repeat scan; importing that ranking and
exporting the larger candidate set completed through worker subprocesses. This
does not measure production-library speed, live-provider decisions or GPU science.
Task links cannot be rewritten after execution starts; changed missions block an
old task's resumption. Live AF3 computation is not implicitly cached by this guard.

An existing literal-SMILES path incorrectly rejected an empty external-input hash
map even though the literal request was sealed. Fixed that path and its prepared
query inspection. Added legacy-note rendering so the document and pinned context
both retain existing project constraints. This final rendering-only addition is
covered by a focused guidance regression run recorded with the delivery.

Final focused guidance regression: 15 passed in 22.58 seconds after the
legacy-note rendering addition. Across the integrated and focused runs, 193
distinct Python cases passed. No browser visual pass, live-provider or production
science acceptance is claimed.
