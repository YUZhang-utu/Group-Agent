# E090: matched target analysis of workstation E089 panel

Protocol recorded before analysis. Verify user-supplied CSV SHA256
5167b31a3896c342b6ec1dceadc2e668474e5812e2327f636b747ed16cfb297d.
Reproduce pooled counts and medians. Collapse each query/target/rotation over
its reference molecules using the median, then subtract that target's median
within-class control. Aggregate matched deltas equally per query, preserving
all rotations and targets. Report missing controls explicitly.

Explore targets represented in all three angular bands as a common-target
comparison. Report reference counts and same-molecule pairs. A diagnostic
envelope can count how many query-target rotations are within all four observed
control maxima, but this is not an acceptance rule or recall estimate: there
are at most three reference-reference pairs per class, selected nonrandomly.
Do not infer whole-library coverage from this class-balanced panel.

Output raw matched tables and summary, including query counts for envelope
existence; contrast existence of a compatible observed candidate with all
candidates being compatible. No automatic minimum-score assignment or merging.

## Result

Completed hash-verified workstation table analysis. See E090_PAIRED_PANEL_FINDINGS.md.
Synthetic matched-control/equal-query arithmetic test passed. Results exploratory,
with no inferred population acceptance, recall or final class membership.
