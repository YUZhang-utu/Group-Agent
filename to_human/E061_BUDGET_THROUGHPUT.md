# Budget throughput diagnosis

User observation: about 2500/660000 chunks in 600 seconds. Same-rate extrapolation
is 158400 seconds (44 hours), excluding later merge and MOL2 export. Template and
molecule costs differ, so this is provisional. Larger chunks reduce file/task
overhead but do not reduce the number of conformers, templates or poses evaluated.

The new implementation uses completion-order bounded scheduling rather than
waiting for the first submitted task. Contact-first consensus ranking evaluates
Gaussian only on the exact maximum-contact ties within each conformer. All
original seeds, assignments, physical clash checks, conformers and templates are
retained. All-zero contacts still require Gaussian across all tied poses.
The per-conformer global Gaussian maximum is now null with an explicit scope
annotation; representative contact, Gaussian and composite terms remain available.

Local equivalence checks cover unique, tied, zero-contact and masked seeds plus
chunk receipt counts. They do not establish a workstation speedup. Seed generation
and contact assignment remain expensive. If Gaussian is only a small fraction,
lazy Gaussian cannot yield an order-of-magnitude end-to-end acceleration.

Read timing receipts while the current run continues, without editing its checkout:

```bash
python -m aidd_agent.budget_profile --run /absolute/path/to/execution/guided/guided/search
```

This command requires the new module in a separate updated checkout if the old
checkout is still running. Alternatively, inspect one existing
`search/chunks/00/*.receipt.json` and `search/protocol.json` directly. Supply
`worker_seconds`, `counts`, workers, chunk_conformers and assignment_backend.
The first available receipts may be biased toward an early template.

Do not pull changes into the active checkout or resume an old sealed run under
new code. Keep the active run and artifacts unchanged until its timing breakdown
is assessed. A new optimized run needs a fresh Chat task/output directory; old
receipts cannot be silently adopted under new code fingerprints.

For larger work reductions, changing templates-per-molecule, conformer budgets,
seed budgets or retrieval size changes the search and requires a separately
reported quality/recall comparison. GPU availability alone does not accelerate
the current CPU seed generation and one-to-one contact assignment.
