# E057 bounded consensus funnel pilot

Before execution: sample unique registry molecule IDs uniformly without replacement
using a recorded RNG seed. Include ALL catalog conformers of each selected molecule
and ALL selected templates. Start with 2000 molecules, not a promised 100000 in
30 minutes. Reuse production preselection compute, matching, thresholds and
scaffold grouping. Persist sampled identities and stage receipts.

Measure sample construction, worker initialization, scan, persistence, merge and
grouping separately; distinguish summed worker time from wall time. Enforce a
supervised wall budget and terminate owned workers on timeout. Partial samples are
informative for profiling only: do not extrapolate survivor rates from fast-finishing
tasks. Complete samples get exploratory population estimates with uncertainty.

No automatic threshold fitting to output capacity. No biological recall claim.
This run excludes live active-control retrieval, full-library integrity byte scans,
docking and unimplemented dual-region scoring. Source/hash and metadata checks
remain; document sampled-integrity scope. Production full-run costs may differ.
