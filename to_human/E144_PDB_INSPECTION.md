# Explicit PDB inspection

The new `pdb_inspect` scientific action takes `pdb_id` only. It does not require
an isoform decision, a binding-site reference or a receptor-adoption decision.
The active project stores the downloaded mmCIF, acquisition metadata and
`inspection-report.json` under the new task's execution directory.

## Workstation acceptance

Pull the delivery branch and restart Chat with the existing storage and runtime
configuration. In the target project, create a new conversation and send:

```text
Create a new structure-inspection task using pdb_inspect with pdb_id 4TQ9.
Preserve earlier blocked tasks; do not resume or modify their plans.
Inspect the deposited sequence, chain identifiers, coordinate coverage,
assembly annotations, nonwater ligands and ions, database mappings and sequence
differences. Save the structure and inspection report in this project.
Report missing metadata and unresolved isoform assignment without guessing.
Do not adopt a receptor, prepare molecules, run AF3, dock or search the library.
```

After completion, ask:

```text
Read the completed inspection report for task TASK_ID. Summarize its deposited
evidence and limitations. Distinguish asymmetric-unit chains from assembly
annotations and nucleotide/ion components from screening references.
Report the saved structure and report paths. Record a concise sourced project
note with this task ID. Do not run another inspection or calculation.
```

The task report persists automatically. A narrative project note is a separate
outer-agent operation and must have its own successful receipt. In a new chat
within the same project, ask to retrieve this task and its saved report.

## Interpretation limits

- Sequences are deposited canonical polymer sequences, not fetched isoform references.
- Coordinate coverage uses per-model label residue indices against the deposited
  construct. It is not full-length UniProt coverage. Missing indexing yields null.
- Database alignment boundaries and sequence differences are deposited annotations.
  They can support a reported mutation but are not an independent alignment check.
- Assembly rows, generator expressions and operations are preserved. No new
  assembly coordinates or biological validation are produced.
- Nonwater nonpolymer instances preserve model, author/label chain, residue and
  insertion identifiers. Polymer-associated modified residues remain polymers.
- No ligand role or receptor is adopted automatically. Unreported fields are not
  evidence that the corresponding feature is absent.

This release adds the previously missing inspection capability. It does not fix
the separate general blocked-plan clarification/replanning lifecycle or implement
independent KRAS isoform mapping. Old sealed blocked plans remain unchanged.
