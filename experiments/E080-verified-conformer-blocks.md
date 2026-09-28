# E080: verified conformer blocks

Protocol, 2026-09-28, before implementation validation.

Implement a source-backed, audit-gated block builder with directed peptide units,
cyclic correspondence, phi/psi/omega descriptors and explicit chemical/typed
side-chain variants. Keep original molecule and conformer identity. Support an
optional exact provenance join to a specified production library; otherwise mark
source-only IDs clearly. Do not change production screening or registry data.

Use disk-backed capacity-bounded median trees; freeze split rules, centroids,
radii, descriptor version and membership hashes. Compare capacities 10000,
20000 and 30000; small blocks are valid. Incremental proposals must preserve the
frozen model, enforce capacity and route unsupported/outlying records to review
or overflow rather than silently changing previous memberships.

Generate deterministic per-block unique-molecule samples at 100 and 500, exporting
all relevant conformers in that block. Keep cross-block molecule repetitions
explicit and cache at conformer identity. Retrospective tail/recall evaluation
requires a complete user-supplied conformer score reference, not invented scores.

Confirmatory checks: proper rigid motion and atom-order invariance; proline and
cyclic symmetry; chemical feature channel separation; distinct conformer blocks;
capacity, determinism, frozen routing and incremental overflow; source/hash/registry
failure handling; unique-molecule sampling across blocks and fair tail fractions.
Exploratory validation: reuse E078 mixed real inputs and their existing audit;
do not rerun E077. Report scientific validation and workstation execution as pending.

Results: 46 focused confirmatory tests passed. Exploratory real mixed-source replay
passed for backbone, chemistry and typed variants: 24/24 conformers, 8 molecules,
no review records; 10 blocks at diagnostic capacity 4 and 9 at each requested large
capacity. Evidence: data/e080-final-{backbone,chemistry,typed}. This small panel
does not measure large-library runtime, clustering quality or screening recall.
