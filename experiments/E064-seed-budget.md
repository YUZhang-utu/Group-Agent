# E064: bounded seed kernels and search-width audit

Protocol written before execution, 2026-09-23.

Hypotheses: bounded batched numerical work reduces seed-generation overhead;
larger generated prefixes may improve contact-first molecular ranking; a survivor
target can allocate additional search to clashing conformers but needs a hard cap.

Preserve the running workstation job. Defaults remain the reference kernel and
512 generated pair seeds. Experimental survivor targets count PCA and pair poses
passing the existing physical mask, not independent pose basins. No centroid gate.
Use a fresh output, a fixed random molecule panel with all its stored selected
conformers, all templates, and the same ranking policy. Compare generated caps
128/256/512/1024/2048 and survivor targets 100/200 with cap 4096. These are
exploratory budgets, not adopted scientific thresholds.

Validate ordered prefix identity (including duplicates, antiparallel and degenerate
pairs), zero-survivor termination, bounded batches, and unchanged default scoring.
Record per-conformer tested/surviving/scored counts and stop reason, worker stage
times, rank overlap against the largest generated reference and molecular score
changes. High-budget rankings are search references, not biological ground truth.
ANN-excluded molecules and rare high-score blocks need a separate sampling audit.

Block-tail analysis reports unique-molecule max, top-k mean, sample size and an
explicit optional high-score fraction. It never discards blocks or asserts recall.
Benchmark fixtures are local evidence; workstation speed and ranking stability
remain unmeasured until the supplied command is run.
