# E135: usable, bounded document retrieval and privacy regressions

Protocol: retain the E134 document-level boundary and same-project public cloud
use. Reject internal chat commands before storing either query or response. Add
stable list/search/read pagination, source and passage integrity checks, reversible
local-only archive controls, idempotent same-library imports and protection against
publishing a known identical internal document. Improve Chinese literal retrieval
with local bigram terms without claiming semantic embeddings. Upgrade existing
indexes in place, retaining IDs and originals. Provide readable local result cards,
page controls and import feedback, without writing internal UI queries/results into
chat. Verify HTTP authentication, malicious/invalid requests, migration/restart,
long-document coverage, tampering and actual serialized model request contents
using a local fake transport. Run affected regression tests; no live-provider,
production retrieval-quality, GPU or browser visual-acceptance claim.

Confirmatory results: all 222 integrated Python cases passed, including 29
knowledge/privacy cases. Four Markdown tests, JavaScript syntax, whitespace and
the English-content guard (765 maintained files) passed. The fake transport
captures actual serialized chat-completion HTTP request bodies: public evidence
is present; rejected internal command queries, internal titles, IDs and content
are absent. Same-project public work continues after internal imports.

Pagination covers a 23-document catalog and a long multi-passage document without
missing/duplicated citation IDs. Tests verify stable identity and archive state
on repeated upload, archive exclusion even for direct model reads, restore,
exact-content internal/public restriction, original/passage/metadata mutation
detection, verified source re-upload repair, old-index migration/restart and
Chinese terms within longer sentences. Invalid paging is rejected before opening
even an empty library. Legacy metadata/passages receive local seals during the
trusted migration; this is not retroactive proof against prior host compromise.

These are confirmatory software checks. No live provider, production corpus
retrieval-quality measurement, embedding model, GPU experiment or browser visual
acceptance was performed.
