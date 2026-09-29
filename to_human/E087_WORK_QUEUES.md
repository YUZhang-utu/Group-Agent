# Uncapped classes and pooled evaluation planning

This change addresses scheduling overhead without inventing chemical equivalence.
The old databases and descriptor coordinates remain untouched. A new sidecar
maps old leaf IDs to their original hard classes and to evaluation queues.
The 20,000-member capacity does not apply to logical classes in this plan.

Rules:

1. Definite cis/trans classes with >=100 conformers get their own queue.
2. Smaller definite classes share a queue only at identical peptide length and
   ordered cis/trans pattern. Original chirality-bearing class IDs remain separate.
3. Boundary-containing classes share review queues by peptide length. They are
   not assigned a definite cis/trans state or merged into a physical family.

The 100-member threshold is a workload policy, not a geometry or activity cutoff.
Queues can contain different chemical identities, explicitly retained in
classes.csv. Every old block appears exactly once in old_block_mapping.csv.
The source conformer population and each class population must balance exactly.

Each queue receives 100/200 planned conformer slots. Allocation is balanced
across classes, deterministic and nested between budgets. If there are more
classes than slots, the remaining classes are marked unassessed_do_not_reject.
No unsampled class is discarded. This is an initial pilot budget, not exhaustive
class validation. Boundary queues with many classes will need additional rounds
or increased budgets before broad conclusions are possible.

These are **conformer-slot plans**, not the old molecule-based sample manifests.
No ligand coordinates, PLANTS input, unique-molecule sample, affinity score or
recall result is produced. Existing search/docking dispatch is not automatically
changed. The plan is a prerequisite for a subsequent coordinate export and
scoring run. Do not accidentally run old per-leaf sampling to use the new plan.

## Workstation command

Update the feature/structure-guided-chat branch first. From the repository root:

```bash
export PYTHONPATH="$PWD/src${PYTHONPATH:+:$PYTHONPATH}"
BUILD=/mnt/local/hand/yuzhang/aidd/mc_runs/e081-001/audit-stage/blocks-20000/backbone
DIAG="${BUILD}-fragmentation-review"
python -m aidd_agent.conformer_work_queues \
  --diagnostic "$DIAG" \
  --output "${BUILD}-work-queues-v1"
```

DIAG must be the directory holding the E086 diagnostic report.json and
hard_groups.csv. If a different --output was used for E086, adjust DIAG.
Use a fresh output directory. The diagnostic and model receipt hashes are
checked, and the current tree is checked against recorded class populations.
Large SQLite hashes are not repeated; the prior successful full validator
remains required. No MOL2 scan or full descriptor scan is needed.

Outputs:

- report.json: exact queue counts, class/population conservation and budgets.
- classes.csv: original identity and new queue for every class.
- queues.csv: population and number of classes per queue.
- old_block_mapping.csv: complete old leaf-to-class-to-queue mapping.
- sampling_plan.csv: nested class quotas and explicit unassessed flags.

## Evidence

Workstation E086: 24,663,736 conformers; 16,531 leaves; 15,236 hard classes.
12,143 leaves hold 143,112 conformers (~0.58%). Boundary-containing classes
hold 1,384,339 conformers (5.61%) across 11,301 groups. These two populations
are not interchangeable, and their intersection was not supplied.

Local tests: synthetic rare, large, cis/trans and boundary classes; nested budgets;
source nonmutation; tampered CSV, model receipt and tree population rejection.
Five tests passed including E086. Existing 24-conformer typed fixture changed
from 10 old leaves to 9 retained classes and 2 evaluation queues. This is only
a smoke test, not the full-library queue count. Workstation execution pending.

Scientific follow-up: recover actual omega angles for boundary records, check
compatibility without losing chirality/alignment, and evaluate held-out geometry
and screening recall before any physical regrouping or rejection. Boundary is
30 < abs(omega) < 150 degrees, so automatic nearest-state assignment is avoided.
