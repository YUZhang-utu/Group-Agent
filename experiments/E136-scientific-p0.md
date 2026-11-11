# E136: scientific workflow P0 completion

## Authorized scope and protocol

Complete five engineering workstreams: assembly-aware pose QC; reproducible
screening evaluation; explicit chemical-state/conformer preparation; project-owned
AF3 batch execution with shared protein features; executable live acceptance.
Keep existing project isolation, immutable requests and cancellation semantics.
Do not infer affinity, activity or experimental acceptance from execution success.

## Confirmatory software checks

1. A clash on a second receptor chain must be detected. Invalid chemistry and
   absent evidence must not receive a scientific pass.
2. Ranking evaluation must reject duplicate/unknown IDs and report positive
   attrition, enrichment and exact-reference recovery separately from timing.
3. Chemical transformations require explicit policy, retain parent provenance,
   and account for failed records and actual conformers.
4. AF3 batch inputs preserve literal sequences, stoichiometry, SMILES and seeds.
   Shared features must match every protein chain. Restart reuses only sealed
   complete children; interrupted dispatch cannot silently launch twice.
5. Acceptance runs use explicit cases, collect local task/input/output evidence,
   and distinguish fixture, preflight and live scientific execution.

## Workstation acceptance

Requires an operator-configured provider, AF3 installation/weights/databases/GPU,
production catalog, target/reference inputs and independently supplied labels.
Inspect available configuration without printing credentials. If unavailable,
record blocked cases and preserve a runnable suite; do not claim live acceptance.
No arbitrary library truncation threshold or predictive superiority is assumed.

## Results

Implemented three model-facing tools, molecular policy extensions and golden-prompt
acceptance. CPU checks use actual RDKit geometry, tiny catalog full scoring, real
worker subprocesses, candidate export and saved-ranking evaluation. AF3 subprocess
responses are explicit fixtures; they do not establish GPU inference success.

Development preflight found chemistry/structure dependencies and an RTX 5070 Laptop
GPU with 8151 MiB; provider credential was unset, no AF3 profile or production catalog
was configured. The user explicitly chose to run live acceptance on the workstation
and return feedback after GitHub delivery. No live-provider or production-scale
scientific validation is claimed.

Final confirmatory regression: 199 Python cases passed in 98.52 seconds across
scientific P0, molecule/candidate workflows, planner/queue, project guidance,
AF3 adapters/shared-MSA/confidence, assembly QC, uploads and knowledge privacy.
Earlier failures identified explicit-state hydrogen normalization, a missing export
directory and legacy single-chain expectations; these were corrected before the
final integrated run. English-content and whitespace checks passed.
