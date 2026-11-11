# Knowledge library update

E135 completes everyday library controls while retaining document-level privacy:
public documents can support the cloud model and internal documents remain in the
local library interface. The rest of the same project remains cloud-enabled.

## Update the workstation

From the existing checkout on `feature/structure-guided-chat`:

```bash
git pull --ff-only origin feature/structure-guided-chat
python -m pip install -e ".[knowledge]"
```

Use the Python environment that runs the Agent. The optional knowledge extra adds
local PDF extraction; UTF-8 text retrieval does not require it. Restart the Chat
server and reload the browser. Keep the existing workspace database and document
directories. Existing indexes upgrade locally in place on first use; document IDs,
original files and source references are preserved. Do not delete the workspace.

## Using the local panel

Open **Knowledge libraries** inside the intended target project. Imports default
to **Internal / group-only**. Public literature requires a public source URL and
explicit classification confirmation. The model cannot import, reclassify, archive
or restore documents through its knowledge tool.

The panel now shows readable result cards with title, document/passage citation,
source, page and text. Use **Read document**, **Read this page**, a page number,
and **Previous / Next** to inspect more than the initial results. A blank page
selector reads across the document. List and search results also paginate.
Switching libraries clears the displayed results and local query. Requests remain
bound to the project that opened the dialog.

**Archive** reversibly removes a document from normal retrieval. **Include archived
documents** is a local view control, and **Restore** makes an eligible document
retrievable again. Model retrieval always excludes archived records, including
direct reads using an old document ID. Sources remain on disk. Archiving does not
erase passages already sent in earlier chats or provider requests.

Re-uploading the same bytes into the same library preserves the existing document
ID, original citation metadata and archive state. It rebuilds derived passages
from the matching source instead of creating another document. This can repair a
changed/missing source copy or passage index when the re-upload matches the sealed
source hash. Changed citation metadata requires explicit local investigation;
untrusted replacement metadata is not silently accepted.

An identical byte-for-byte document already recorded as internal in this project
cannot be uploaded or restored into its public library. Importing an internal copy
of an existing identical public document archives the public copy and returns an
explicit notice. Earlier disclosures cannot be undone. This is an exact-content
guard, not detection of all rewritten, summarized or cross-project copies.

## Public retrieval through chat

```text
Search the public knowledge library for permeability assays. Read the relevant
pages, follow pagination if necessary, and cite the document and page for each
conclusion. Keep hypotheses separate from the reported observations.
```

The model tool returns `has_more` and `next_offset`; `limit` is at most 20.
`read` accepts an optional page and an offset, so long text documents and pages
are no longer limited to their first few passages. Source bytes, saved metadata
and selected passage hashes are checked before returning evidence. Missing or
changed material gives an explicit error rather than an unchanged-looking quote.

Search is local FTS5/BM25 with English words and Chinese character-pair terms.
It supports literal Chinese terms inside longer sentences without an external
tokenizer. This is lexical retrieval, not trained semantic embeddings, and no
retrieval-quality improvement is claimed without a corpus benchmark. PDF parsing
is local; scanned PDFs still require separately supplied local OCR text.

## Privacy boundary

Internal search/list/read remains local-UI-only. Neither private queries nor their
results are automatically inserted into chat history, notes, project memory or
cloud requests. In particular, an accidental `/knowledge` internal query is rejected
before its text is saved as a chat message. This closes the previous command-history
path while keeping public commands usable.

The protection applies to documents managed through Knowledge libraries. Ordinary
chat text, general-purpose attachments and preexisting scientific reports are not
automatically classified. Do not paste internal material into ordinary cloud chat.
There is no general semantic DLP detector or private-document LLM analysis in this
release. Local seals detect changes relative to trusted records; they do not protect
a compromised host or erase existing external copies.

## Verification

See [E135 protocol/results](../experiments/E135-knowledge-hardening.md). Coverage
includes long-document/catalog pagination, archive/restore, source/passage/metadata
changes, duplicate identity, original-index migration and restart, Chinese literal
matching, and actual serialized provider request bodies captured by a fake local
transport. This is software validation, not a live-provider security audit,
production retrieval-quality benchmark or browser visual acceptance.
