# E138: public literature pagination and answer recovery

User reports an exhausted agent turn with malformed model JSON after importing a
public document. Remote trace contents are unavailable; no exact provider failure
cause or successful import is assumed. Code inspection confirms that model-visible
knowledge pages can exceed the evidence envelope and lose their structured cursor,
and invalid final JSON consumes the same bounded loop as evidence reading.

Protocol: bound model knowledge results by whole passages, preserve next_offset;
retain an explicit read ledger; repair only unambiguous whole JSON fences; offer a
tool-free plain-text synthesis recovery using already retrieved public evidence;
on failure return exact read progress and specific failure stage. Do not resend
mutating tools, weaken private-library gates, invent unread findings or claim full
text coverage from snippets. Verify native transport payloads, privacy, exhaustion,
successful recovery, invalid evidence citations, and unchanged source documents.

## Validation (2026-11-11)

Confirmatory local regression: 149 tests passed across literature recovery,
knowledge privacy, domain-agent dispatch, streaming/uploads, prompt workflows,
project guidance, scientific P0 and chat-agent behavior (73.23 seconds). After
adding explicit step-budget and metadata-only coverage cases, all 12 literature
recovery tests passed (2.64 seconds). These runs overlap; counts are not additive.
The English-content guard passed and git diff --check reported no whitespace errors.

The native text-recovery HTTP payload was inspected through a captured synthetic
transport. It contains no tools or JSON response-format constraint. Tests confirm
that only public-library passages enter recovery, private library contents remain
excluded, no calculation is submitted, and existing imports remain intact.
No live provider request, remote trace inspection or user-PDF accuracy evaluation
was performed. Workstation acceptance remains pending; the user-facing protocol
is to_human/E138_PUBLIC_LITERATURE_RECOVERY.md.
