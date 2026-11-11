# Target or uploaded ligand to local-library candidates and MOL2

Date: 2026-11-11. This corrects the overly narrow receptor/reference-only framing
in the E129 audit. A standalone 3D ligand does not require a receptor for ligand
similarity search.

## The two entries

1. Target name/accession: resolve target identity, retrieve experimental complexes,
   analyze ligand diversity and pocket states, and prepare actual cocrystal ligand
   queries using the existing structure workflow. Inspect those queries and select
   the intended reference for whole-library 3D similarity and MOL2 export. Multiple
   distinct binding sites or ambiguous target identities must be resolved using
   evidence; do not silently invent a pocket or pick a ligand name.
2. Upload: attach one nondegenerate 3D conformer in a single-molecule SDF or MOL2.
   The new ligand_library tool resolves the session-owned attachment, prepares
   chemical features and directly searches the configured local library. A receptor
   is optional for this ligand-similarity purpose. A 2D file or SMILES table still
   needs a separate explicit conformer-preparation step; uploading alone does not
   make every format directly searchable.

The existing target pipeline is staged and retains review points. This change
does not establish unattended target-name-to-export completion across multiple
asynchronous tasks or independently validate live-provider routing. It does
connect prepared queries and uploaded ligands to real search/export execution.

## Copyable prompts

Attach a 3D SDF/MOL2 and write:

```text
Use the attached 3D molecule as the reference for our local-library 3D search.
Use ligand_library to score every stored conformer in the configured library.
Rank unique molecules by their best 3D shape/chemical-feature similarity and
export the top 1000 unique candidates as original and aligned MOL2 files.
Keep the complete ranking, molecule names, scores and provenance locally.
Report the task ID, processed conformer count, unique molecule count and output paths.
Do not require a receptor and do not run docking.
```

For target-driven discovery:

```text
For human <TARGET_GENE_OR_UNIPROT>, find experimental protein-ligand complexes
and analyze their bound ligands, structural quality, binding sites and diversity.
Use the existing structure workflow to prepare evidence-backed cocrystal ligand
queries. Report distinct site/reference choices if they are ambiguous.
The objective is to search the entire configured local molecular library and
export the top 1000 unique 3D-similarity candidates as original and aligned MOL2.
Keep results locally and report actual task IDs and output paths for each stage.
Do not ask me to manually provide a receptor that you can retrieve from the PDB.
Do not run docking.
```

Once the intended prepared query is available:

```text
Inspect the prepared ligand queries in task <TASK_ID> using ligand_library.
Use query <EXACT_RETURNED_QUERY_ID> to search every stored conformer in our
local library and export the top 1000 unique candidates as original and aligned MOL2.
```

## Execution and scientific meaning

The trusted tool takes either attachment_id or task_id/query_id. Its execution
plan is local-only; the LLM cannot invent unrestricted file paths. Input hashes,
library identity/lineage, code and settings are recorded. File ownership and
compute gating remain enforced. `search.batch`, `workers`, `refine_chunk` and
`--allow-compute` use the existing host runtime.

The new path reuses the Gaussian scorer and chunk scheduler. All configured
catalog conformers are scored with the existing finite rigid-pose seeds. There
is no candidate retrieval Top-K or new internal-conformer generation. The rank
is best unweighted Gaussian shape/color score per unique molecule, with stable
ties. The complete ranking remains in SQLite. The export count limits delivery,
not library coverage; shortfall is reported. Source catalog and registry molecule
coverage must agree before full coverage is reported.

The original-record exporter verifies record hashes and chemical identity, keeps
original names/bonds/hydrogens/charges, and writes original plus transformed MOL2.
Output includes `molecules.csv`, `molecules.jsonl`, `original/`, `posed/`, and a
per-file integrity manifest. Chat results and task cards expose these paths.
The single-template RRF compatibility field is monotonic in Gaussian rank; it
does not indicate multi-template consensus evidence.

This lane ranks ligand similarity, without receptor contacts, clash screening,
binding affinity or activity validation. For contact-first multi-reference
screening, the existing consensus/budget workflow remains available and already
exports MOL2. It queries the full-library ANN index, then scores retrieved
molecules; it must not be labelled all-conformer scoring. The existing
target_search contact-funnel path remains separate. Whole-catalog scoring can
be more expensive than indexed retrieval; production speed/recall were not
measured in this change. Existing results are retained; no production search is
automatically launched by updating the code.

## Verification

Small-library tests actually build conformer and chemical artifacts, score 3D
poses, check complete catalog coverage, rank molecules and export original/posed
MOL2. The upload path also runs through the real Chat worker subprocess. Prepared
cocrystal query discovery is represented by a local fixture; no live PDB/provider
or production-library run was performed. Tests cover deduplication, stable ties,
source tampering, attachment ownership, 2D rejection, compute gating, preserved
names, chunk reuse and immediate queue receipts.

Pull feature/structure-guided-chat and restart using the same storage/runtime/
provider configuration. Do not change code during an active sealed scientific run.
