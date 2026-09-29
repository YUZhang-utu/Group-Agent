# Inspect block fragmentation before changing assignments

User direction: retain cis/trans distinctions; a logical class does not require
a 20,000-member cap. First inspect existing outputs, then choose adjustments.
This diagnostic never rebuilds or edits the model or descriptors.

From the repository root on the workstation (after updating the feature branch):

```bash
export PYTHONPATH="$PWD/src${PYTHONPATH:+:$PYTHONPATH}"
BUILD=/absolute/path/to/the/completed/build
python -m aidd_agent.conformer_block_diagnostics \
  --model "$BUILD/capacity-20000" \
  --descriptors "$BUILD/descriptors.sqlite" \
  --output "${BUILD}-fragmentation-review"
```

Use a fresh output directory. Return report.json; hard_groups.csv provides each
class population and omega pattern; block_size_histogram.csv provides exact
leaf sizes. Report bins are disjoint successive intervals, not cumulative.

Key values:

- blocks_without_capacity_limit: exact existing root count if capacity splitting
  is removed while keeping current hard classes unchanged.
- capacity_extra_blocks: leaves minus roots, the extra blocks due to the cap.
- boundary_attribution: populations, root and leaf counts with/without boundary.
- maximum_uncapped_group_size: largest class if current capacity splits disappear.

The tool reads the small tree table and one indexed descriptor per hard group.
It does not scan all source MOL2 records or recompute chemistry. It relies on
the completed model's prior successful full integrity validation and records
receipt hashes; it does not rehash the large SQLite files.

Current hard groups encode variant, peptide unit count, ordered chirality and
ordered omega classes. Boundary means 30 < abs(omega) < 150 degrees. Its broad
range includes strongly nonplanar states; it is not merely cutoff jitter.
Boundary-associated fragmentation is not proof of boundary-only causation.
No automatic conversion to cis/trans is performed.

Next decision: eliminate capacity splitting at the logical-class level if
appropriate; retain independent batch limits for memory and docking execution.
For boundary records, inspect angular distributions and compatibility with
definite patterns before choosing assignment, overlapping retrieval, or a
traceable uncertain pool. Preserve handedness and unambiguous cis/trans states.

Local validation: synthetic known-count test passed; input bytes unchanged.
Small existing typed model: 24 conformers, 10 leaves, 9 original classes;
3 boundary-containing classes hold 3 conformers. These are local test data,
not the workstation's 24,663,736-conformer result. Workstation counts pending.
