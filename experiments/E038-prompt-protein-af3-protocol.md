# E038 prompt-driven protein, AF3 and calibrated 3D workflows

2026-09-19. Protocol before tests.

Implement a provider-configurable JSON planner, deterministic allowlisted executor,
active-Project ownership, local credential environment lookup, public UniProt/RCSB
evidence and AF3 input/run adapters. Model never supplies executable commands,
filesystem paths, protein sequences, library budgets or credentials. References
must resolve to preceding typed steps. Ambiguous target resolution must block,
not silently choose a sequence. Keep E031 annotation-only and calibrated WEE1
search boundaries explicit; AF3 outputs are predictions, not pose/biology gates.

Modes: offline fixture planner for reproducible smoke tests, live compatible chat
endpoint, plan-only and run-plan. Operator selects endpoint/model and trusted
local AF3/search profiles. Default prompt execution prepares artifacts; explicit
compute opt-in permits expensive AF3/search actions. Preserve hashes, reports,
stage receipts, and reject changed inputs/code during resume. API errors cannot
leak credentials. Cloud context contains user prompt and tool schema only, never
automatically collected sequences, library data or experimental measurements.

Tests: mocked HTTP and model responses, strict schema and dependencies, unsupported
commands/paths, ownership, ambiguity/species/sequence identity, AF3 protein/CCD
input, shell-free adapter execution, failed/blocked status, output tamper/resume.
No live LLM/AF3/GPU speed or scientific quality claim from offline tests.
Workstation smoke script must run without credentials/network/AF3; optional live
tests have explicit prerequisites and show blocked/not-run instead of success.

References checked: OpenAI Chat API, UniProt REST API, RCSB Data API, official
google-deepmind/alphafold3 docs/input.md. AF3 basic dialect v1 protein/CCD inputs
are used for broad installed-version compatibility; MSA/template generation is
left to configured AF3 pipeline, not fabricated by LLM.
