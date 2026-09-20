# Repository workflow skills

Open the repository as the agent workspace to discover `.agents/skills/`. The
skills travel with the Git checkout; no machine-specific global installation is
required. Invoke a skill by its name or let the agent select it for a matching task.

| Skill | Use it for | Runtime actions |
|---|---|---|
| `aidd-3d-search` | Calibrated WEE1 retrieval, Gaussian/E031 and performance comparison | `search_3d` |
| `aidd-protein-preparation` | Verified target identity, sequence and structure evidence | `protein_fetch`, `protein_resolve`, `pdb_search`, `pdb_fetch` |
| `aidd-alphafold3` | AF3 input preparation, installed execution and prediction import | `af3_prepare`, `af3_run` |
| `aidd-prompt-workflows` | Model configuration, structured plans and Project recovery | Orchestrates the other groups |

Examples for a repository-aware agent:

```text
Use $aidd-3d-search to review the E037 report and run the calibrated QT9 search.
Use $aidd-protein-preparation to retrieve human WEE1 evidence without selecting a receptor.
Use $aidd-alphafold3 to prepare full-length P30291 with ATP using my installed profile.
Use $aidd-prompt-workflows to test my compatible model endpoint and inspect its plan.
```

## How the two interfaces connect

Codex reads the relevant SKILL.md for procedures and evidence boundaries. The
application's compatible LLM receives compact skill/action groups from
`src/aidd_agent/workflow_skills.py` alongside its existing capability schema.
Validated plans record the applicable skill names for audit. Neither interface
executes Markdown as code or accepts model-supplied command lines.

Skills organize existing functionality; they do not turn WEE1 calibration into
arbitrary-target acceptance or AF3 confidence into biological validation. Runtime
ownership, typed dependencies, profile checks, receipts and compute controls remain
enforced by the existing adapters. Existing user authorization should be honored
without adding redundant approval steps.

## English-only maintained content

All authored prompts, comments, examples, skills and guides use English. Planner
summaries and clarifications are instructed to use English, with a Han-text rejection
guard. User-provided inputs and external scientific data retain their identity and
provenance; no automatic translation of biological sequences or raw evidence occurs.

```bash
python scripts/check_english.py
```

The check includes tracked and nonignored untracked UTF-8 files, including these
skills. It excludes ignored raw data, external environments and Git history. The
guard detects Han text, not every possible non-English language; authoring review
remains necessary. Historical reports are translated with their original dates,
measurements and evidence limits preserved. This does not require rerunning a
scientific benchmark. Code-hashed old runs still require matching code or a fresh
output directory after an update.

For workstation commands, use the [E038 guide](../to_human/E038_PROMPT_PROTEIN_AF3.md)
and the experiment-specific guides linked from each skill.
