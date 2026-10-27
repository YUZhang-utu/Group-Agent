# KRAS potential-candidate retention policy

Date: 2026-10-27. Decision record; no production filtering changes implemented here.

The objective is conservative reduction of the existing direct-contact pool, not
quota filling. Preserve required contact definitions. Failed bounded rigid search
is unresolved evidence, not a molecule-level rejection.

## Evidence categories

- Pose-supported: E120 finds one pose meeting all original proximity rules and
  the 2 A heavy-atom floor. Source-coordinate validation is still required at export.
  This is not validated affinity or a complete interaction-quality assessment.
- Unresolved potential: a molecule belongs to the original direct-contact pool,
  but lacks an E120 accepted pose, including unprocessed molecules. Retain its ID
  for flexible docking; do not label the original overlapping pose acceptable.
- Necessary-condition failure: reject only when a documented necessary condition
  fails across the explicitly permitted representations. No such additional
  molecule-level rejection is established by current optimizer failures.

The existing counts.sqlite remains the authoritative potential pool: G12D
2,632,299 IDs and G12V 4,080,670 IDs in prior completed evidence. E120 ranking.sqlite
contains a subset and does not delete failures from that source pool. Do not
overwrite either file or edit code used by an active run.

The 600-second user-provided progress contains 6 G12D and 501 G12V eligible IDs;
these are progress counts, not final export counts. Optimizer failure, a clash
in one pose, low whole-ligand color, or low heuristic rank must not by themselves
remove an ID from the potential pool. A top-N budget is prioritization, not proof
that excluded molecules cannot bind.

G12V residue 12 currently means steric compatibility, not a required favorable
Val12 contact. G12D requires the configured Asp12 donor proximity. Shared rules
are typed contacts to His95, Asp69, Glu62, Met72 and Arg68, not arbitrary residue
proximity.

## Remaining implementation and validation

Provide a separate unresolved-pool export/manifest rather than repurposing the
strict E120 export as high-recall screening. Any further reduction needs a stated
necessary-condition argument or a measured false-negative audit with expanded
pose/conformer searches. Existing pool membership itself arose from finite seeds;
retention within this pool does not guarantee recall over the full library.
