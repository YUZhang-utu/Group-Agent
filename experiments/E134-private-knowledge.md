# E134: public and internal project knowledge

Protocol: separate public and internal knowledge storage and retrieval at the
document level. User clarified that importing internal material MUST NOT block
cloud use for the rest of that project. Model-facing retrieval opens only the
public store; internal raw text, titles, summaries and search results stay in a
separate local UI path without chat/memory insertion. Public classification
requires direct user confirmation and cannot be set by model tools. Validate
project isolation, private metadata absence, model tool denial, continued public
cloud context, restart persistence and authenticated local retrieval. Keep this
first retrieval implementation local and lexical; no embedding or cloud OCR.
No claim about private text manually pasted into ordinary cloud chat, preexisting
reports, arbitrary plugins, OS access controls or full host/network isolation.

Confirmatory software results: 208 integrated Python tests passed, including 15
new knowledge/privacy cases. The same project continues scripted cloud-assisted
public retrieval after an internal import; captured model context contains public
passages and no private text/title/ID. Internal search/list/read are rejected
before index access, guessed internal IDs do not resolve in the public store,
private file registration outside the project root is denied, and document-borne
instructions do not unlock private retrieval. Local authenticated HTTP retrieval
does not insert results into chat or project memory. Cross-chat persistence,
restart and project isolation passed. Local text-PDF extraction preserves its
actual page number; a scanned blank PDF reports local OCR as missing without
any external fallback. Four Markdown cases, JavaScript syntax, English763files
and whitespace checks passed. Providers were scripted; no live-provider,
production-corpus retrieval-quality, network penetration or browser visual claim.

The focused regression exposed an existing cancellation publication race: terminal
job status became visible before the process handle was cleared. Reordered the
cleanup before status publication; the integrated cancellation test passed.
