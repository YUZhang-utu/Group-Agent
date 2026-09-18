# E037 bounded pair-seed generation and consolidated workstation validation

2026-09-18, before implementation tests.

Hypothesis: stopping deterministic pair-seed generation after the requested number
of unique seeds removes unused work while preserving the exact prefix of the old
unbounded generator. No seed-budget reduction or score change is permitted.

Keep default unbounded public behavior; expose opt-in bounded generation in scaled
Gaussian. Compare seed objects/order/matrices exactly against old unbounded slicing
on seeded random, degenerate, symmetric and cap-boundary cases. Compare complete
Gaussian output archives for both engines. Measure synthetic generator timings
separately from full query time; no extrapolation to workstation speedup.

Workstation suite runs paired reference/bounded full-index pipelines using both old
and fine chunk layouts: old/reference, fine/reference, fine/bounded, old/bounded.
Reverse order in repeat 2 to reduce order confounding. All keep 10000 candidates,
5000 Top-N,512 seeds and E031 annotations. Compare every output against E034; failures
block successful suite status. Store individual protocols/reports and a resumable
summary, separating fresh runs from partial reuse and original completed timing.
Do not rerun million stress validation by default. Worker profiles are an optional
separate diagnostic run; profiled timings never count as latency measurements.
