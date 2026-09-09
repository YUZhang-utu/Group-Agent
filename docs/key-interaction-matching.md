# Key-interaction matching sidecar

E031 measures how well every already rigidly aligned candidate reproduces the
co-crystal ligand's key interaction anchors. It performs no new pose search,
torsion sampling, receptor clash filtering, or mutation of accepted Gaussian
scores and ranks.

For each named rigid objective pose, the scorer builds a typed spatial and
directional anchor-by-candidate-feature matrix. A deterministic maximum-weight
one-to-one assignment gives every query anchor a score in `[0,1]`. The reported
`interaction_match_score` is the query-anchor-weighted mean. The NPZ sidecar
retains the assignment and contribution of every anchor.

Run both accepted WEE1 queries on the Linux workstation:

```bash
cd /mnt/medchem_taltio/wrk/yu_agent/Group-Agent
git pull --ff-only origin main
export PYTHONPATH="$PWD/src${PYTHONPATH:+:$PYTHONPATH}"
bash scripts/run_e031_key_interaction_matching.sh "$PWD"
```

The default output is
`/mnt/local/hand/yuzhang/aidd/e031-key-interaction-matches-v1`. Each query has
an NPZ sidecar, JSON manifest, and `/usr/bin/time -v` report. The terminal
summary reports absolute wall time, score distribution, matched-anchor count,
Spearman association with the accepted anchored Gaussian objective,
Top-100/500/1000 overlap, and the ten strongest interaction matches with their
original Gaussian ranks.

Input locations can be overridden without editing the script using
`E031_QT9_QUERY`, `E031_QT9_RIGID`, `E031_1X8B_QUERY`, `E031_1X8B_RIGID`,
`ARTIFACT_CATALOG`, `CHEMICAL_COMPANION`, and `E031_OUTPUT_DIR`.

This is an evidence lane, not an activated replacement ranking. Its usefulness
requires multi-target held-out enrichment or redocking calibration. The first
schema covers the currently materialized direct hydrogen-bond anchors; other
interaction classes need their own target-independent extraction validation.
