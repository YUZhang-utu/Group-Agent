# Paper citations and reusable project reviews

## Reader-facing sources

Public-literature answers display paper names linked to their recorded source.
Multiple tool receipts for the same paper collapse to one paper in the source list.
The Sources and reading details disclosure contains page coverage; tool evidence
IDs, document UUIDs and full execution receipts remain in the Agent trace panel.
This presentation does not remove internal evidence validation. Failure messages
can still expose exact cursors needed for diagnosis.

In Knowledge libraries, select Public literature, list documents and expand
Citation details on a document. Supply a short title and, optionally, an author
label and publication year. Save citation details. These labels are user-edited,
not independently verified bibliographic metadata. Without them, the existing
document filename is used; no author or publication date is fabricated. Existing
chat messages retain their historical display; new answers use current labels.

## Project review records

Eligible completed paper reviews are saved automatically in the owning project's
local public knowledge index. Each record retains its interpretation, original
document hashes, page coverage, originating conversation and trace. Records are
immutable versions; they do not rewrite the project mission or become accepted
scientific facts. Mixed scientific-tool conversations and answers with unavailable
source dependencies are not automatically stored as public literature reviews.

The planner receives a bounded index of relevant/recent review records, not all
paper text or full conversation histories. It can list/search records and read
their summaries progressively. This is local retrieval of project experience, not
fine-tuning. Reusing a review does not count as rereading the original PDF. New
claims and verification should return to original passages. A newer review does
not automatically settle conflicting interpretations.

Every retrieval checks the original sources. If a source is archived, reclassified
as internal, modified or unavailable, dependent reviews are withheld from model
retrieval. Previously displayed chat text or provider disclosures are not erased.
No internal documents or internal bibliographic labels are indexed into these
public review records. Ordinary chat text is outside document-library isolation.

## Existing work and local browsing

Knowledge libraries / Public literature has Saved public reviews and Index previous
public reviews controls. Indexing inspects bounded batches of previous completed
review traces from this owned project, verifies their public source dependencies,
and adopts eligible records without a PDF reread, new model request or calculation.
Repeat the index action if another batch remains. Repeated indexing is idempotent.
Saved summaries can be inspected locally without putting their text into chat.

The same operation is available through an ordinary prompt:

```text
Index our previous public-literature reviews in this project, without rereading
PDFs or launching calculations. Then list the available paper reviews.
```

After that, open a new conversation in the SAME project:

```text
Use this project's saved public-literature reviews to summarize what we have
learned about CIP2A, the limitations and unresolved questions. Cite the papers
by name. Distinguish saved interpretations from freshly checked source evidence.
Do not launch calculations or reread original papers unless needed to resolve
missing evidence; explain any such gap first.
```

## Update and acceptance

Pull feature/structure-guided-chat, restart the existing Agent server with the same
storage directory, and refresh the browser to load the new JavaScript. Do not
re-import documents. Check that paper links replace E-number citations in final
answers, source details start collapsed, old eligible reviews are indexed, and a
new same-project conversation reuses them with the correct limitations.

Automated tests exercise citation presentation, safe Markdown, actual captured
synthesis payloads, local index persistence/restart, project isolation, source
revocation/mutation, idempotent adoption and bounded summary pagination. Synthetic
model transports do not establish live-provider quality or scientific correctness;
browser visual acceptance and actual-paper review remain workstation checks.
