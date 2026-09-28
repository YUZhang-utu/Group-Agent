# E076: name-to-SMILES source audit

Question: do the supplied CSVs resolve the chemical identities of the previously
examined MOL2 conformers by exact molecule name, independently of file order?

Protocol: stream every CSV row, retain entries only for names in the first 1000
records of each MOL2 file, and record duplicate-name ambiguity. Parse matched
SMILES strictly, compare canonical heavy-atom graphs, hydrogen/charge information
and specified stereochemistry against strictly sanitized MOL2 records. Persist
source hashes, row numbers, original names and comparison statuses. Do not rewrite
structures or transfer atom indices by name or row position. Unspecified stereo
is not evidence of an exact stereochemical match.

This is an exploratory source-prefix audit, not full-library chemical validation,
building-block token assignment or an energy-basin benchmark.

## Results

Both CSV files have 50000 rows. All 668 names in the 2000-conformer prefix have
one matching row across the supplied CSVs; all 2000 structures return identity_match.
No unspecified potential stereo was reported in those matched structures.
The 250 E075 ring ambiguities therefore remain extraction issues after this
identity check, not evidence of invalid source chemistry.
Focused regression: 11 tests passed, covering stereo/charge/isotope mismatches,
aromatic hydrogen retention, atom reordering, duplicate names and prior blocking.
Exploratory outcome only; no whole-library identity claim. Output:
data/e076-csv-audit/report.json and identities.jsonl.
