# E058 pilot remediation protocol

Exploratory, before implementation measurements. User workstation baseline:
256 molecules, 808 conformers, 11 templates, 226 hits, 65.3243 s total,
36.5514 s scan. Scaffold reconstruction: 87 successes and 139 failures.

Hypotheses: original MOL2 aromatic fallback and/or lost hydrogen metadata explain
reconstruction failures; optional ANY eligibility and broad cohort coarse bounds
explain permissive union; deterministic assignment compilation reduces CPU work.

1. Read existing poses without rescanning. Count template hits, unique additions,
contact masks and score distributions; failures remain identifiable.
2. Recover a bounded set of original registry records with content hashes. Compare
heavy atom numbers, charges, aromatic flags and bonds; inspect original loader mode
and hydrogen counts. Never repair chemistry by guessing hydrogen or bond orders.
3. Compare p53 and ligand-defined regions on identical aligned coordinates. Record
agreement, boundary, unoccupied and unknown; missing definitions are not negatives.
No hard gate or new weights until reference outcomes and selections are reviewed.
4. Compile the existing assignment loop without fastmath, keeping tie ordering.
Require exact assignment/score equivalence on ties and randomized matrices, then
repeat the exact saved sample on workstation with all templates and unchanged rules.
Measure compile/startup separately where possible. Never extrapolate partial runs.

Local tests are implementation checks. Actual failing molecules and throughput
require workstation execution; no biological recall claims.
