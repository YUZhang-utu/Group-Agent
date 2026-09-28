# E079: frozen source inventory and block admission preparation

Protocol written before execution, 2026-09-28.

Hypothesis: an explicit source inventory and a separate audit admission receipt
can prevent partial-library and ambiguous-identity records from silently becoming
production block assignments.

Implement a read-only source inventory with SHA-256 hashes, optional same-stem
CSV coverage and orphan CSV reporting. Keep outputs outside input directories.
For completed audit artifacts, verify recorded output hashes and stream a
per-conformer readiness manifest. Retain original source identities. Quarantine
failed records, duplicate conformer names and same-name chemistry conflicts.
Represent graph-only alignment, equivalent rotations and boundary omega as
separate flags, without discarding accepted geometry solely for these flags.

Compare capacity lower bounds at 10000, 20000 and 30000. These are planning
counts, not fitted clusters. Registry joining and descriptor extraction remain
required before production assignment.

Confirmatory tests: missing CSV, orphan CSV, output placement, stale audit hash,
incomplete audit, cross-record conflicts and exact row accounting. Exploratory
replay: reuse the completed E078 mixed fixture, without rerunning the two-shard
audit. No screening, docking, external execution or production registry changes.

Results: 34 focused confirmatory tests passed. Exploratory E078 receipt reuse
retained all 24 geometry candidates and three boundary-omega flags. The source
inventory correctly reports two MOL2 files, one companion CSV and one missing
companion. Evidence: data/e079-mixed-readiness and data/e079-mixed-inventory.
No production assignments or whole-workstation coverage claims are made.
