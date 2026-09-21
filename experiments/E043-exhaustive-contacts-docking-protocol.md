# E043 exhaustive screening, interaction classification and docking validation

Implement an exhaustive standardized-USRCAT coarse pass with audited row coverage,
then reuse Gaussian refinement and annotation without claiming equivalence to
different approximate candidates. Distinguish retrieval budgets from chemistry.

Use PLIP XML as the typed rule-based contact evidence contract. Preserve eight
classes, exact ligand instance, residues, atom indices and geometry. Crystal,
candidate and docking evidence must be separate. Missing water/metal evidence and
unprocessed candidates must not be treated as experimentally absent interactions.
Create an English browsable HTML/CSV report for human review and explicit selection.

Docking uses trusted operator profiles and prepared inputs, never model-generated
shell commands or invented grids. Preserve original IDs, execution logs and output
hashes. Execution completion is separate from reference-pose validation. Missing
engines or prepared inputs must block honestly. The user has been asked for local
PLIP and receptor/grid preparation details; independent development continues.

Tests: streaming exhaustive Top-K against dense reference including ties and all
row coverage; typed contact XML parsing and ligand isolation; malformed/missing
evidence, report escaping, source checksums and explicit selections; docking argv
and failure handling using mocked executables. Live whole-library, PLIP and licensed
docking checks remain workstation work. No service SLA or biological claim.

Continuity uses repository checkpoints; this environment has no scheduler tool.

## User correction before implementation of contact classification

The user clarified that this is ligand-based 3D query matching, not per-candidate
complex/docking interaction analysis. PLIP is not installed and must not become a
default dependency. Implement native crystal-derived contact hypotheses for the
existing indexed feature families, then match these features in stored candidate
poses. Keep unsupported classes explicit rather than calling all features verified
contacts. Derive pocket center/extent from the selected crystal ligand; generate
reviewable protein/ligand preparation jobs before docking. No prepared inputs exist.
