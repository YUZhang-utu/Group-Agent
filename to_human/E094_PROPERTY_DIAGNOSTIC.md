# E094 property split diagnostic

Question: why did 384 regular parent blocks yield 384 property blocks?
Confirmed implementation limitations: 69-channel mean/std/max summary loses directed residue positions; fixed RMS threshold 0.35 can dilute sparse channel differences; PCA proposes only one median-near cut; rejected proposal reasons are not persisted by v1. These do not prove the workstation failure mechanism.

Use scripts/diagnose_e094_properties.py standalone with NumPy. Copy it to the workstation export-tools directory; no git pull or change to frozen implementation is required. Run with --root pointing to final-blocks-e094, --output pointing to a fresh sibling diagnostic directory, and --top 10 (default). --top 0 replays all parents. Source SQLite files are opened read-only. The script reads cached property vectors without parsing MOL2 or recalculating chemistry. RAM scales with the largest inspected parent because exact PCA replay uses full matrices. It writes metadata statistics before replay, then updates report.json after every block. Do not run multiple writers on one diagnostic directory.

All property-block size statistics come from the full small CSV. Root proposals use every cached vector for the selected parents in CID order, matching v1 float32 arithmetic. Diagnostics include population-floor failures, no eligible projection cut, actual contrast below threshold, active channels, separate chemistry/steric RMS and largest channel shift. Large databases are not rehashed; existing successful E094 validation is required. Do not change membership or tune a threshold until reviewing the evidence. Root pass despite no original split merits environment/policy investigation.

Three local regression tests pass: size-floor versus constant features, single-channel RMS dilution, and agreement with the existing partition decision. Initial pytest invocation lacked PYTHONPATH; rerun with src on PYTHONPATH passed. English guard passed.

Next design question: distinguish order-invariant broad composition from position-sensitive side-chain chemistry on the ring. Cached summaries can support broad property analysis but cannot reconstruct directed positional patterns. A positional descriptor would require original per-unit features, with cyclic correspondence and preserved backbone/class identities; do not silently equate mixed execution pools with a homogeneous backbone state.

## Published checkout commands

Run from the existing Group-Agent checkout and the same scientific environment used for E094:

```bash
git pull --ff-only origin feature/structure-guided-chat
BASE=/mnt/local/hand/yuzhang/aidd/mc_runs/e081-001/audit-stage/blocks-20000
OUT="$BASE/e094-property-diagnostic-$(date +%Y%m%d-%H%M%S)"
OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 python scripts/diagnose_e094_properties.py \
  --root "$BASE/final-blocks-e094" \
  --output "$OUT" \
  --top 10
printf 'Diagnostic report: %s/report.json\n' "$OUT"
```

Return the diagnostic report.json. No old build command needs to be rerun.
