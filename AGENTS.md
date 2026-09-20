# AIDD repository instructions

Maintain all authored code comments, prompts, examples, skills, reports and
documentation in English. Preserve exact scientific identifiers, sequences,
measurements and provenance. Do not rewrite ignored raw datasets, external
environments or unrelated workspace projects as part of a language cleanup.

Repository skills are discoverable under `.agents/skills/`:

- `aidd-3d-search`: calibrated library retrieval, Gaussian/E031 and performance evidence.
- `aidd-protein-preparation`: target identity, verified sequences and structure evidence.
- `aidd-alphafold3`: verified AF3 inputs, configured execution and prediction import.
- `aidd-prompt-workflows`: provider setup, structured plans and Project execution.

Read only the relevant skill and references. These skills guide existing adapters;
they do not expand available scientific methods or external execution permissions.
Use the current `research-state.yaml` and dated results to distinguish implementation,
local fixture validation and real workstation evidence. Preserve append-only logs.

Check maintained content with `python scripts/check_english.py`. Use the existing
test suite for executable changes. Do not rerun heavy scientific validation solely
for documentation changes; code-hashed runs must use matching code or fresh outputs.
