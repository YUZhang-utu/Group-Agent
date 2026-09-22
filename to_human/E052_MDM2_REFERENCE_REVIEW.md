# MDM2 project: crystal census and diverse references

The local project inputs are D:/agent/MDM2. Original files are preserved.
The generated review is D:/agent/MDM2/analysis/e052/review.html, with raw API
snapshots, downloaded original coordinates and CSV/JSON alongside it.

Use the new branch Chat with this natural-language request:

> For human MDM2, use PDB 5C5A as the reference pocket. Census all experimental
> structures by resolution, check ligand-bound coordinates, compare ligand
> similarities and propose chemically diverse references. Include peptide and
> macrocyclic peptide complexes as a separate track. Do not start a library scan.

This invokes protein verification followed by `structure_diversity`. It does not
require user-written code, file paths in the model-generated plan, or manual API
queries. The standalone HTML is listed in the task outputs; `/results` shows the
census and bounded reference summaries. Local folder upload/registration is not
yet a generic web-chat feature: this project's eight supplied files were inventoried
by the workspace agent, and public original coordinates were fetched by PDB ID.

## Actual exploratory result

RCSB search for verified human MDM2 Q00987 returned 147 experimental entries:
138 X-ray, 8 solution NMR, 1 electron microscopy. X-ray resolution bins contain
19, 73, 32, 9 and 5 entries at <=1.5, (1.5,2.0], (2.0,2.5], (2.5,3.0] and >3.0 A.
These are target-associated entry counts, not all validated drug complexes.
All 133 X-ray entries <=3 A were downloaded. Five remain excluded from contact
analysis because of unresolved strict target mapping; failures are explicit.

The organic nonpolymer track has 68 unique ligands after the stated occupancy,
heavy-atom count, resolution and shared-pocket contact checks. The initial eight
farthest-first references cover 11/68 at Morgan Tanimoto >=0.6. The greedy full
cover uses 44 references; it is neither a proof of minimum panel size nor a
screening-sensitivity estimate. The proposal is for human inspection.

The peptide track contains 46 unique deposited sequence/link variants at the
same pocket among structures <=2.5 A and partner chains <=50 residues. Its
similarity matrix is CCD-residue-token edit similarity, not chemical graph or
3D similarity. The user's 7NUS, 8GCG and 9CDZ are polymer ligands and must not
be dropped by a nonpolymer-only query. Current polymer query preparation remains
unvalidated; these examples are not marked search-ready.

## Next boundary

This stage completes the census and reference-review step. The existing guided
search consumes one adopted query. Multi-reference library execution, union of
candidate membership with retained reference provenance, and validated peptide
chemical/query preparation still need implementation. Do not describe an eight
reference proposal as an already executed multi-reference full-library search.
The E050 workstation task was not modified or restarted.
