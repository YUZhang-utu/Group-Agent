---
name: aidd-protein-preparation
description: Retrieve verified UniProt protein sequences and RCSB structure evidence for this AIDD project, prepare explicit constructs and preserve target identity and provenance.
---

Use an owned, active Project. Read the
[protein and prompt guide](../../../to_human/E038_PROMPT_PROTEIN_AF3.md) for exact
CLI/profile setup. Existing adapters are in `src/aidd_agent/protein_data.py`,
`rcsb.py`, `acquisition.py` and `prompt_workflow.py`.

## Resolve identity before preparing structures

- Prefer an explicit canonical UniProt accession with the intended organism.
- Otherwise use `protein_resolve` with a gene and taxonomy ID. It requires exactly
  one reviewed match. If ambiguous, present the candidates and request the missing
  accession; do not choose the first result or invent a sequence.
- `protein_fetch` checks returned accession, optional organism, sequence length
  and hash. Preserve `protein.json`, `protein.fasta` and source URLs.
- A requested construct uses a1-based inclusive start/end interval. Preserve the
  full-sequence hash, interval and construct hash. Do not infer a domain boundary
  or replace nonstandard residues without an explicit scientific decision.

## Structure evidence

`pdb_search` consumes a preceding verified protein step; its default is X-ray
structures at resolution<=3Angstrom. `pdb_fetch` downloads an explicit PDB ID into
the Project with provenance. The prompt adapter does not select a receptor.

A UniProt association alone does not establish the desired chain/domain, complete
pocket, ligand state or suitable conformation. Review the actual structure and
existing comparison/campaign records before selection. Preserve missing residues,
alternate states, construct differences and source confidence as unresolved facts.
Do not infer protein sequence from a ligand search result.

Sequence and raw structure artifacts stay local; the planner receives capability
descriptions and the user prompt, not automatic uploads of these artifacts or
experimental measurements. Record failures and unresolved identity separately from
successful evidence retrieval.
