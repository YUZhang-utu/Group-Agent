# E087: separate uncapped classes from evaluation queues

Protocol 2026-09-29, before execution. The E086 workstation result has 16,531
leaves and 15,236 hard groups. Capacity adds only 1,295 leaves. 12,143 leaves
hold 143,112 conformers. Boundary-containing classes hold 1,384,339 conformers
in 11,301 groups. These are user-supplied results, not locally recomputed.

Hypothesis: unnecessary per-class sampling overhead can be removed without
unvalidated chemistry merging. Preserve every hard-group ID and all source
memberships. Remove capacity subdivisions from logical classes via a sidecar.
Pool boundary classes by peptide length for uncertain-work scheduling only.
Pool small definite classes (<100 conformers) by length and exact omega pattern.
Other definite classes keep separate evaluation queues. Chirality remains in
the original class ID, even when execution is pooled. No physical state claim.

Output deterministic nested conformer-slot allocations at budgets 100 and 200
per queue, balanced across constituent classes, recording zero-slot classes as
unassessed. No unassessed class may be rejected. This is a planning artifact,
not actual ligand export or unique-molecule sampling. Report full conservation,
mapping, reduced queue counts and sampling workload. Existing adapters remain
unchanged until a reviewed export/dispatch stage is connected.

Validate synthetic rare/large/cis/trans/boundary classes, capacity leaf mappings,
nonmutation and corrupt inputs. Run small existing fixture. Workstation run
needed for actual queue counts; do not extrapolate them from E086 aggregates.
