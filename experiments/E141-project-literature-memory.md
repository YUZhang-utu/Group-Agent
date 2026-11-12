# E141: paper-level citations and reusable project literature reviews

User requests simpler paper-level references instead of tool receipt IDs and asks
to finish project knowledge reuse across conversations. Protocol: render deduplicated
paper labels and links, keep page/reading details collapsed and internal receipts in
the trace, support locally edited public bibliographic labels without inventing
authors or dates. Preserve strict internal evidence validation.

Save bounded, versioned model-generated reviews only from public literature tools,
with source hashes, coverage and provenance. List/search/read them within the owned
project across chats; introduce only a bounded index into planning, then retrieve
summaries as needed. Never recursively treat summaries as fresh paper readings.
Validate every source before exposing a saved review, including archive, internal
reclassification, source mutation and project isolation. Review storage must not
change the project mission or imply experimental verification. Test serialization,
restart, multiple papers, safe rendering, privacy and synthesis recovery.

## Confirmatory validation (2026-11-12)

195 integrated Python tests passed in 74.79 seconds across the literature stack,
knowledge privacy, domain loop, project guidance, chat/streaming, prompt workflows
and scientific P0. Five Markdown tests pass, including collapsed source details,
safe external links and inert HTML/script URLs. Chat JavaScript syntax, the English
content guard and whitespace checks pass. Additional eligibility checks cover
uncited answers and mixed-tool conversations not becoming public review memory.
The final review-specific suite passed all 18 cases in 4.65 seconds; it overlaps
the integrated suite and is not an additive count.

Captured synthetic provider requests verify saved-review synthesis and exclusion
of private source text/IDs. Real local index tests verify cross-chat/restart reuse,
source archive/reclassification/mutation exclusion, project isolation, immutable
records, exact bounded summary paging, idempotent adoption and paper deduplication.
No live cloud-provider quality, scientific extraction accuracy or workstation
browser visual acceptance is claimed. Guide:
to_human/E141_PAPER_CITATIONS_AND_PROJECT_REVIEWS.md.
