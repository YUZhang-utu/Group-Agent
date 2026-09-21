# E048: explicit joint coarse eligibility

Hypothesis: requested anchor feasibility alone admits too many conformers.
Adding explicit whole-ligand size, shape-extent and typed-feature coverage
predicates before seed generation can reduce expensive evaluations.

These are new user-selected eligibility predicates, not proven necessary bounds
on Gaussian similarity. No defaults, percentile cutoffs, Top-K or forced rejection
quota are permitted. Existing anchor-only policies retain their meaning.

Validate rigid-motion invariance, threshold boundaries, invalid inputs, chemical
coverage multiplicity, policy propagation, and selection/funnel consistency.
Compare unfiltered reference poses with optimized survivors under the SAME new
policy. Include current positives; zero-hit panels cannot validate hit retention.
Measure pre-seed rejection and wall latency separately. A 90% rejection target is
an unmeasured performance objective, never an asserted outcome.

Pocket-derived classified features can be explicitly required in the anchor list.
Receptor clashes need placed poses; these predicates do not test clashes, pocket
occupancy, docking, activity or all possible conformational orientations.
