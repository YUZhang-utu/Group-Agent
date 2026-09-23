# Protein-contact-only screening policy

The user excludes water-mediated interactions from the current non-explicit-water
workflow. New recommendations omit water_bridge and unmodeled metal_coordination
anchors from model context; local validation rejects their selection. Adoption and
consensus funnel entry reject such roles, even at zero weight. Source evidence is
preserved and historical proposals remain readable by diagnostic tools. This policy
does not assert that all crystallographic waters are erroneous, nor validate direct
contact hypotheses as energies. Ligand polar features remain in shape/color scoring.

An explicit migration command creates a NEW proposal without mediators. It removes
affected mandatory roles and alternative/optional group members, drops empty groups,
and records the removed IDs. Review any removed hard conditions. It preserves an
existing fixed budget, or freezes the old weight sum for a legacy weighted mean.
Legacy dimension coefficients are normalized to sum to one without changing their
relative contribution. All other weights and templates remain unchanged. An empty
remaining anchor set is rejected. It does not adopt or run the proposal.

For the current 20-mode proposal, migration removes six water modes, retaining 14
selected anchors and budget 9.75. A fresh LLM recommendation starts a new budget;
use migration for controlled before/after comparisons. The original explicit pose
threshold remains unchanged; a null threshold is recalculated upon adoption.

After pulling the feature branch, activate the existing environment and set
PYTHONPATH to src. Run:

```bash
python -m aidd_agent.contact_policy \
  --recommendation /absolute/path/to/original/recommendation/report.json \
  --output /absolute/path/to/new/protein-only-proposal
python -m aidd_agent.evidence_audit \
  --recommendation /absolute/path/to/new/protein-only-proposal/report.json \
  --output /absolute/path/to/new/protein-only-audit
```

Historical water relations may still appear in the full evidence diagnostics, but
no water anchor is selected for scoring or eligibility in the revised proposal.
Do not resume an older adopted water-containing funnel under new code. In Chat,
restart the same workspace and request a new recommendation, or explicitly edit
the existing recommendation to remove water IDs from all roles while preserving
its budget. An unedited /adopt of the old water-containing proposal is rejected.

Validation: 82 focused tests passed, including history readability, mediator
rejection, partial-family preservation, fixed-budget preservation and empty-design
rejection. The saved local proposal migration removed the expected six IDs.
