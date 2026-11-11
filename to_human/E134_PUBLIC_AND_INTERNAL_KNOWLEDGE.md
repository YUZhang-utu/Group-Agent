# Public and internal knowledge in one project

E134 applies confidentiality to documents and retrieval routes. Importing internal
group material does **not** disable cloud use for the project. Public literature
and normal authorized scientific work remain available in the same project.

## Import and use

Pull `feature/structure-guided-chat`, restart Chat and reload the browser. In the
target project's sidebar, open **Knowledge libraries**. This dedicated interface
is separate from ordinary chat attachments.

The import selector defaults to **Internal / group-only**. Choose a PDF or UTF-8
TXT/MD/CSV/TSV/JSON file, up to 20 MiB. Internal files can be listed and searched
in this local interface. Results are not inserted into conversations, project
notes, the project dossier or model-visible memory. Importing or searching does
not make an LLM request.

For public literature, select **Public literature**, supply its public source URL,
and explicitly confirm that the document, title and source contain only public
information. This is a user classification, not automatic verification of the
URL or document's confidentiality. No URL is fetched during import. The existing
literature-search tool can still discover online titles and abstracts; its results
are not automatically imported into this persistent library.

For public library retrieval through chat:

```text
Search the public project knowledge library for permeability assay methods.
Compare the retrieved evidence and cite each document ID, source and page.
```

The direct chat command is also public-only:

```text
/knowledge {"operation":"search","collection":"public","query":"permeability assay","limit":5}
```

Private queries belong in the Knowledge libraries panel, whose search selector
chooses **Internal / group-only**. Private query text and result text remain out
of chat history. The API's local user controls are authenticated by the existing
Chat bearer token and same-origin checks. Model tools cannot import documents,
change their classifications or invoke the internal retrieval route.

## Storage and model boundary

| Store | Location | Model access |
| --- | --- | --- |
| Public | `<project_root>/knowledge-public/` | Public search, list and page reads |
| Internal | `<chat_root>/internal-knowledge/<project_id>/` | Denied; local user interface only |

Each store has its own SQLite full-text index and original source copies under
`documents/<document_id>/`. Records include a document ID, title, source reference,
source SHA-256, page and passage number. Internal storage is outside the project
artifact root that model-facing report/file tools can access. Guessing an internal
document ID in a public read finds no record. Requests for the internal collection
are rejected before its index is opened, so no titles, counts, snippets or summaries
are returned through that model tool. A private retrieval request in `/knowledge`
is also denied because normal chat responses become future model context.

Both libraries are project-scoped and persist across conversations and restart.
They are not recursively added to PROJECT.md or project memory. Public passages
retrieved by the model are supplied as context and can therefore reach the
configured cloud provider. Internal passages are displayed only through the
separate local UI/API; there is no automatic model-generated summary of them.

## Retrieval and operational limits

Retrieval currently uses local SQLite FTS5/BM25 keyword matching, bounded passages
and source/page citations. It is not embedding-based semantic search. Search
quality for paraphrases, multilingual terminology and long Chinese phrases has
not been benchmarked. No external embedding, reranking or OCR service is used.

Text PDFs are parsed locally with the optional `pypdf` package. To enable PDF
imports in the Agent's Python environment:

```bash
python -m pip install -e ".[knowledge]"
```

UTF-8 text import needs no PDF package. Scanned/image-only or unparseable PDFs
produce an explicit local-extraction error; no cloud OCR fallback occurs. The
importer limits PDFs to 500 pages and documents to 5,000 passages.

This boundary applies to documents imported through Knowledge libraries. It does
not automatically classify existing scientific reports, ordinary attachments,
chat messages or manually written project notes. Pasting internal information
into ordinary cloud chat, or explicitly importing it as public, can expose it.
There is no universal content-based DLP detector or retroactive removal from
earlier provider requests. Likewise, this change does not isolate arbitrary
external plugins, OS processes, third-party runtimes or a compromised local host.
Keep workspace backups and filesystem access controlled; SQLite and original
files are not newly encrypted by this feature.

Internal model reasoning and automatic scientific workflows over private document
contents are not part of this first retrieval boundary. Such additions need a
separate local processing path and explicit control of any derived outputs.

## Validation

See [E134 protocol and results](../experiments/E134-private-knowledge.md). Tests
exercise continued public chat in a project that contains internal documents,
denial before index access, guessed IDs, project isolation, restart, local HTTP
authentication, absence from history/memory, public prompt injection attempting
internal retrieval, public source/page provenance and local PDF extraction.
Scripted providers inspect the constructed context; no live-provider or network
penetration-test claim is made. UI JavaScript is checked without a browser visual
acceptance claim.
