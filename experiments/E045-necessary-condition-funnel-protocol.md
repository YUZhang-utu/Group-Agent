# E045: Conservative condition-driven full-library funnel

Trigger: user E044 log shows 18 additional 2,048-row chunks in 429 seconds,
36,864 / 429 = about 85.9 conformers/second, for QT9. Straight local-rate
extrapolation gives about 82.6 hours for remaining QT9 IDs; this is not a full
runtime prediction. The unfiltered all-pose design does not meet fast screening.

Hypothesis: an explicit same-pose feature rule admits a conservative rigid-invariant
prefilter. Type/direction-kind compatibility, distinct feature assignment and
pair-distance bounds can reject impossible conformers before Gaussian poses.
No descriptor Top-K, Gaussian Top-N or molecule budget may decide membership.

For feature score spatial * angular >= t, angular <= 1 implies positional error
<= sigma * sqrt(-2 log(t)), also limited by the scoring cutoff. Triangle inequality
then bounds the discrepancy between corresponding feature-pair distances by the
sum of positional radii. Use a conservative numerical margin. Arc consistency and
bipartite matching are necessary conditions only, not proofs of a feasible pose.
ANY rules and single-feature diagnostics cannot use conjunctive pair constraints.

Reuse an explicit selection preview as the rule definition. Never silently turn
all listed contact categories into required anchors, or invent an acceptance
threshold. E044 unconstrained per-feature diagnostics cannot be replaced by a
conjunctive rule without changing the question. Subsequent counts are conditional
on the selected rule, not unconditional counts for every feature.

Validation: analytic impossible/possible fixtures, alias deduplication, matching
collisions, threshold boundaries, random rigid transforms with independently
computed passing assignments, and exact selected-ID equivalence against an
unfiltered small library. Skipped rows have explicit no-pose markers; full coverage
counts include rejected rows, but matching counts exclude them. Frozen legacy E044
chunks may be reused only with matching inputs, scoring code, parameters and
feature columns. Log any contradiction between the necessary-condition bound and
a previously passing pose as a hard failure. Record measured speedup as pending
until real candidate survival and timing are available.
