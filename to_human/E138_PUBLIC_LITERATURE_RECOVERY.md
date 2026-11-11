# Public literature reading and answer recovery

## What changed

Public model retrieval now returns complete passages within the agent evidence
budget. Paging remains structured: use next_offset, not offset plus the requested
limit. The local Knowledge libraries UI keeps its normal paging behavior.

The agent records which document pages/passages it actually read. Two consecutive
invalid model responses stop the formatting retry loop. After successful public
reading, a separate tool-free Markdown synthesis can recover an answer using only
the already retrieved evidence. It cannot dispatch scientific work. Public reading
reserves time for synthesis; the turn remains bounded. If synthesis also fails,
the answer reports reading progress, continuation cursors and the failure phase.
It does not claim that a successful import failed or invent unread results.

Entire fenced JSON objects are accepted after strict parsing and schema validation.
Fragments, guessed values and arbitrary JSON repairs are not accepted. Provider
diagnostics record finish reason and output length without logging raw responses
or credentials. Internal-library gating is unchanged; neither private data nor
private lookup metadata is included in the recovery context.

## Update

```bash
git pull --ff-only origin feature/structure-guided-chat
```

Restart the existing Chat server with the same configuration and storage directory,
then refresh the browser. Existing imported PDFs and indexes are preserved. No
additional parser, model weights or new optional dependency is required by this fix.

## First verify import independently of the model

Inside the intended target Project, enter:

```text
/knowledge {"operation":"list","collection":"public","limit":5}
```

This local command should return the imported document ID, title/source, hash and
indexed-page information. If the document is absent, check the selected project,
public/internal selector, archive status and the original import result. The model's
formatting error alone does not prove a PDF parsing failure.

## Then test one paper

Use the actual returned document ID:

```text
Review public document DOCUMENT_ID in this project. Read its abstract and the
results relevant to CIP2A, following project_knowledge pagination as needed.
Write a concise summary of at most 400 words covering the target construct,
experimental methods, main findings and limitations relevant to this project.
Cite the document ID and PDF pages. Separate reported evidence from interpretation.
Explicitly state which sections were not read or could not be parsed.
Do not launch calculations or change the project goals.
```

An indexed page is not proof that all its content, tables or figures was extracted.
The importer remains local text extraction and lexical retrieval; scanned pages,
complex tables and figures may require a separate local parsing/OCR workflow.

If the turn reports incomplete reading, continue the same document using the
returned cursor and ask for the still-unread relevant sections. Do not re-import
the PDF. A short answer with explicit partial coverage is not a full-paper review.
This change improves runtime recovery; it does not establish scientific extraction
accuracy. Compare the returned claims, units, methods and page references against
the original PDF during workstation acceptance.

## Verification boundary

Tests cover whole-passage paging, cursor completion, malformed model JSON, the
native fallback path with a captured provider HTTP request, known evidence IDs,
no tool dispatch during synthesis, internal-data exclusion and preservation of
existing public documents. A synthetic transport verifies software behavior only;
live provider and the user's actual PDF still require workstation testing.
