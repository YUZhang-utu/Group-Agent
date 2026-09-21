# E042: human-selected screening evidence and docking handoff

Question: can an existing completed search be reviewed, classified, explicitly
selected and exported through persistent chat without repeating search or silently
changing the E031 ranking policy?

Implement three separate stages: evidence review; selection preview with explicit
anchor IDs, all/any policy and score threshold; human-requested export of that
preview. Resolve prior tasks within the same session and owned Project. Record
input hashes and reject modified upstream artifacts. Counts distinguish conformers,
source-grouped molecules, per-query results and cross-query union.

Anchors are crystal-ligand features with recorded receptor-partner evidence. Map
manifest anchors to exact Gaussian feature indices, including legacy reconstruction
checked against stored query arrays. Do not label feature matches as verified
candidate hydrogen bonds. Unsupported interaction classes remain unassessed.
All required anchors must match in the same conformer/objective pose. Deduplicate
molecules only after selecting qualifying poses; retain the best Gaussian-scoring
qualifying pose, with stable ID ties. No cross-pose or cross-query anchor mixing.

Confirmatory fixture tests: exact anchor mapping, same-pose AND versus OR, zero
matches, threshold validation, stable molecule deduplication, source tampering,
session isolation, separate export intent and real artifact-derived responses.
Run regression, English and syntax checks. Live workstation acceptance must reuse
or run one calibrated WEE1 search, review evidence, preview explicit conditions,
then request export; inspect SDF/IDs/counts. Docking execution is outside this
handoff until prepared receptor/grid and ligand preparation are validated.

No live library, model API or docking engine is available in the local fixture
environment. Record those stages as pending, not passed. Repository checkpoints
provide continuity; no scheduler tool is available.
