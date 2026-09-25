# E065: retrospective block exploration

Protocol recorded before execution, 2026-09-24.

Hypothesis: an equal initial molecular sample per block, prioritized by the mean
of its top-k contact scores, recovers more high-ranking molecules per evaluated
molecule than uniform sampling. A random exploration budget protects against
rare high-scoring members of initially low-scoring blocks.

Replay a fully evaluated panel. Blocks must be defined without target scores,
using documented scaffold/side-chain groups or, separately, storage shards.
Include all evaluated molecules, including those without surviving poses. Never
treat absent calculations as failures. Scores must share a scoring/search budget.
Use a frozen full-run molecular rank for reference retrieval when available;
otherwise the reference is the maximum molecular contact score, not RRF.

Defaults: 32 initial molecules per block, top-5 mean, 32 additional molecules per
round, 25% total molecular evaluation budget, 20% random exploration probability,
five deterministic seeds. These are exploratory settings, not production gates.
Compare adaptive exploration against uniform sampling at exactly equal molecular
cost. Report full-reference top-N recall, missed IDs and per-block allocation.
Unvisited members never influence adaptive priority. Recompute the top-k mean
after each sampled batch. No automatic block rejection or parameter adoption.

This retrospective experiment excludes pose recomputation. Molecular cost is not
CPU/wall time. It measures only the completed panel, not ANN-excluded molecules
or biological activity. Tune on one panel and validate on separate blocks/seeds
before implementing live scheduling. Small blocks use min(k, sampled count).
