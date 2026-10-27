# E120: reuse the direct pool with contact-constrained rigid refinement

User retains every original mandatory contact, asks for useful candidates today,
and wants up to 1000000 molecule IDs per target in 100000-record MOL2 batches.
Do not rerun full direct counting or replace the contact policy with random sampling.

Existing inputs: sealed E116 references/library; completed direct-count SQLite;
E119 capped boundary poses including global/conformer IDs and transforms. First
refine these saved poses directly. Further work uses a deterministic spread
subset/order of existing library conformers filtered by the direct molecule pool.
Excluded molecule IDs receive no seed generation. The plan records coverage,
including an explicit candidate/conformer and wall-time budget.

Refinement changes only six rigid-body degrees of freedom. It does not change
internal bond lengths, angles or torsions, protonation, topology or receptor atoms.
Optimize atom-distance violations and all mandatory contact-distance violations
together, with analytic rigid-body derivatives. Use bounded incremental Euler
rotations (+/- 0.65 radians each) and translations (+/- 4 A each). Target a 2.4 A
soft receptor distance margin and 0.05 A inside each contact upper bound. Original
2.0 A hard exclusion and direct-contact masks remain the final acceptance tests.
These distances remain gross-clash/proximity proxies, not hydrogen-bond energy or
angle validation. No claim of a validated docking pose or binding affinity.

Use up to three lowest-overlap direct-contact seed poses per conformer; retain
already feasible poses. Score only independently rechecked, same-pose contact
and geometry survivors. Do not add features from separate poses. Save start/end
metrics, transform, conformer ID, source lineage, optimizer evaluations and
feasibility independently of optimizer success. Failed optimization is not evidence
that the molecule cannot bind. Reuse E118 soft ranking/diversity and source export.

Native translated/rotated perturbation controls and synthetic gradient/contact
loss tests precede release. A bounded gate checks saved E119 examples; continue
only if at least one real-library pose is feasible, and disclose per-target counts.
This gate is selected-example evidence, not a representative enrichment benchmark.

Batch default: four-hour cooperative compute budget, drain in-flight chunks,
checkpoint all committed results, then export actual eligible candidates even
from a partial run under explicitly labelled budgeted scope. Never silently claim
full coverage or fill a shortfall with failed poses. Export can take additional
I/O time. Budget and total target numbers are not completion guarantees.
