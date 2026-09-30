# E091 workstation receipt review - 2026-09-30

The user reports successful completion. Aggregate population conservation passes:
23,279,397 unchanged definite + 1,251,454 provisionally assigned boundary +
132,885 special = 24,663,736 admitted conformers. Special fraction is 0.538787%,
regular provisional coverage 99.461213%, meeting the <=1% engineering target.
The local routing implementation SHA256 matches the reported implementation.
Remote SQLite files and their reported hashes have not been independently checked.

121,630 of the special records (91.53%) are outside the 45-degree policy or
unresolved. Only 11,255 remain for other reported reasons. Do not relax the
angular or chemistry policy to force further reduction. Sparse/failed model
counts are model counts, not conformer counts or final block counts.

Wall time is 4.49 hours; source feature extraction accounts for 94.05% of it.
Routing itself took 205 seconds. This is build time, not search time.

Next: validate the CID sidecar join, one-to-one coverage, target identities,
allowed directed rotations, provenance and assigned/special disjointness; emit
final logical-class and execution-queue population summaries; benchmark exhaustive
special-set search and compare held-out queries against complete reference search.
Do not reuse old block radii as bounds for expanded memberships. No rebuilding
of source features is indicated solely by this successful receipt.

The 1,150,072 source chemistry-review records remain outside this admitted-set
denominator. They have not been resolved by E091. PLANTS/MDM2 and N-E evaluation
remain later steps, and the affinity model is still undecided.

Earlier candidate-band counts exclude some no-target cases; this report bands
all boundary rows by deviation. The totals differ by 90, matching the reported
no_compatible_definite_class population. Exact band allocation requires checking
assignment rows; no corruption conclusion follows from the aggregate difference.
