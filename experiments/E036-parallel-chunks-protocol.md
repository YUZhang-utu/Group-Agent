# E036: bounded parallel chunks and Gaussian scheduling

Protocol fixed before implementation tests, 2026-09-18.

Hypothesis: independent E035 stress chunks can run in separate processes with
unchanged samples, scores, assignments and global ranks. Gaussian coarse has only
five 2000-row tasks for 16 workers in E034; smaller tasks can improve utilization,
while smaller refine tasks can reduce tail imbalance. Neither claim is a measured
speedup until workstation wall times are available.

Implement bounded in-flight work (2 x workers), per-process readers, deterministic
assembly, original hash-bound receipts, serial fallback, progress and ETA. Keep
E031 annotation-only, budgets, Top-N, seed caps and scoring equations unchanged.
Preserve old completed evidence; new code uses new output directories.

Tests: real spawn workers on synthetic mmap fixtures, serial/parallel array and
global-rank equality, resume and tamper rejection, bounded scheduling and failures.
Record cumulative worker compute separately from elapsed wall time. No online
latency claims from stress tests, reused chunks or integrity/reference scans.

Workstation Gaussian comparison reuses E034 candidates and queries. Compare all
Gaussian output arrays, identities and transforms with original outputs, then
batched E031 sidecars. Initial candidate configuration: coarse chunks 500,
refine chunks 64, 16 workers, same scoring parameters as E034. Treat timing
against historical E034 as exploratory; preserve full-call timing and reuse flags.
This is the next measured step toward full-library indexed retrieval plus 3D
refinement; it does not establish activity enrichment or new pose quality.
