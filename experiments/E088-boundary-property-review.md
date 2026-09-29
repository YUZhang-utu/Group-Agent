# E088: boundary candidate coverage and continuous steric properties

Protocol: 2026-09-29. User wants >90-95% useful ordinary-block coverage, with a
small exhaustively searched special set. Do not tune angular cutoffs to force
coverage. Distinguish workload pooling from scientifically supported membership.

Recover omega from stored sin/cos, inspect all boundary descriptors without
MOL2 recomputation. Recover the omitted explicit chirality signature by exact
verification of the original hard-group hash over R/S/achiral sequences with a
bounded enumeration budget. Unresolved signatures remain special. Use only
directed cyclic rotations, keep known cis/trans positions unchanged, and target
only existing definite classes. Exact 90-degree states remain unresolved.
Report candidate upper-bound coverage at tolerances 35,45,60,75 degrees from
the nearest cis/trans centre, no production assignment. Every boundary CID
must be accounted for exactly once. No bridge/union of definite states.

Continuous property review uses side-chain size, charge, cyclic N and N-methyl
annotations plus donor/acceptor/aromatic/hydrophobic feature counts. Add actual
coordinate VDW-expanded local-axis extents, proximal branching, spread and reach
as steric proxies. Keep property and steric differences separate until calibrated.
These are not union volumes or binding energies. Do not require exact sequences.

Property data availability on the workstation is unknown. Check sibling
chemistry/typed receipts. Existing typed feature moments are partial geometry,
not a complete side-chain occupancy model. New steric extraction is for bounded
review samples, not an unannounced full-library re-extraction.

Tests: angle roundtrip/state validation, handedness and definite-state rejection,
cyclic rotation, midpoint and unresolved-signature handling. Small existing
backbone fixture smoke run. Steric profile rigid transform invariance and
response to changed side-chain coordinates. Workstation evidence pending.

## Local results

25 focused tests passed, including offline validator exit handling. Existing
24-conformer backbone fixture: all 9 chirality signatures recovered, 3 boundary
records accounted for; 2 candidates at 35 degrees and 3 at 45. This is exploratory
fixture evidence only. Workstation run and actual special-set size are pending.
