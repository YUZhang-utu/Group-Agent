# E078 optional CSV audit

Protocol: discover all MOL2 files recursively and same-stem CSVs when present.
Freeze the input list. Keep the explicit-list strict mode for compatibility, but
enable MOL2-only extraction in directory discovery mode. Retain cross-reference
conflicts; never turn duplicate/invalid CSV entries into a no-CSV pass.

Verify molecular graph, stereochemistry and geometry from MOL2; use graph-only
peptide mapping if a name cannot be verified, recording that limitation. Report
CSV-backed versus MOL2-only coverage and same-name chemical inconsistencies.
Test discovery, CSV conflict behavior and MOL2-only geometry. Run a bounded real
mixed-source integration test (12 conformers from each real source, one CSV
present and one absent), not another full unchanged-source audit.

Results: 30 focused tests passed. Mixed real fixture completed with 24 accepted
conformers: 12 CSV-backed and 12 MOL2-only. No name/chemistry conflict was hidden;
all MOL2-only conformers retained peptide ring and geometric omega information.
Evidence data/e078-mixed-audit/report.json. Exploratory adapter validation only;
full workstation source discovery and execution remain pending.
