# E046: batched pose scoring and hardware validation

Hypothesis: the unchanged E045 seed search spends substantial time in repeated
Python validation, duplicate color-distance calculations and per-seed dictionaries.
Batch seeds, reuse color kernels, and select winners before building records.
Compare NumPy and optional CuPy float64 execution without changing seed budgets,
objectives, requested anchors, thresholds or the all-library population.

User evidence: 55,296 checked rows, all displayed chunks retained 2,048/2,048
by necessary bounds. The user stopped the run. This is a failed efficiency test,
not evidence that all conformers meet the final condition.

Add conservative relative-direction bounds where query directions exist. Validate
passing rigid poses survive; never substitute a stricter user rule. Diagnose each
rejection class and timing. Use a bounded, spread-out real-library pilot before
another all-library execution. Report exact winning seeds/transforms, numeric
errors, E031 assignments, per-anchor threshold decisions and throughput. Reject
accelerator acceptance on any discrete disagreement. GPU measurements include
synchronization and host transfers. GPU availability and performance require the
workstation; local NumPy fixtures cannot certify CUDA.

Keep the reference backend. New outputs pin backend; old output protocols cannot
silently resume changed code. GPU uses one owner process, not 24 competing CUDA
contexts. CPU workers respect affinity and cap default concurrency at 24. Pilot
success is sample evidence only, not biological validation or a full-run SLA.
