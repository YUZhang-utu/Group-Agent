# Resume checkpoint - E031 implemented, real workstation timing pending

Date: 2026-09-09

## Resume point

Resume from GitHub `main` after code commit `f6cd02d` (`feat(e031): score key
interaction coverage`). The next scientific action is the real Linux
workstation run of E031 for both accepted co-crystal queries. Do not redesign
the search or add pose relaxation before collecting that output.

The active question is:

> Can a target-independent, explainable key-interaction coverage sidecar add
> useful evidence to the existing fast rigid 3D retrieval result while costing
> no more than 10% of rigid-refinement wall time?

## Locked architecture boundary

The workflow is deliberately layered:

1. Broad reusable FAISS/pharmacophore retrieval protects recall.
2. Existing Gaussian coarse and fixed-Top-N rigid refinement perform the 3D
   pose search.
3. E031 evaluates the already stored rigid poses once and reports how well each
   candidate reproduces the co-crystal ligand's key interactions.
4. Receptor vdW annotation is restricted to a small docking handoff.
5. Torsion beams, 13 rigid micro-seeds, minimization, and docking remain
   downstream pose-preparation work.

E031 must not add new pose seeds, change admission, overwrite Gaussian scores,
or introduce WEE1-specific rules. Shape and ordinary color stay as separate
earlier evidence lanes. The new interaction score is not yet an activated
replacement ranking.

## Accepted real workstation evidence

### E029 chemical companion

- Catalog:
  `/mnt/local/hand/yuzhang/aidd/chemical-companion-v1/catalog.json`
- `split_0001`: 149,999 conformers, 226.737 seconds.
- `split_0002`: 150,000 conformers, 244.290 seconds.
- Total: 299,999 conformers in 471.027 summed shard-seconds, approximately
  636.9 conformers/second.
- Exact source/output hashes, contiguous global IDs, boundary reads, and no
  `.partial` directories passed.
- Identical-command reuse completed in 1.49 seconds at 43,008 KB maximum RSS;
  catalog and shard manifest hashes remained unchanged.

### E028 chemistry-aware pre-docking QC

- Output:
  `/mnt/local/hand/yuzhang/aidd/e028-predocking-qc-chemical-v1`
- 200 selected tasks and 200 valid chemical SDFs, 100 per query.
- Retrieval membership and admission order were unchanged.
- All 200 poses had close receptor points; 195 had severe center-distance
  points. These are pose annotations, not retrieval rejection evidence.
- 1X8B/824 vdW penalty min/median/max:
  13.942 / 50.197 / 83.300; median clashing-atom fraction 0.582; median severe
  atoms 25.5.
- 8BJU/QT9 vdW penalty min/median/max:
  6.900 / 17.337 / 49.025; median clashing-atom fraction 0.389; median severe
  atoms 15.0.

Human visual review established that both query sets are in the correct pocket
and coordinate frame. For 8BJU, the LOW pose looked usable and conflicts were
mainly terminal; no ring penetration was observed. For 1X8B, LOW was mainly a
terminal conflict, while MEDIAN/HIGH retained close macrocycle placement.
Therefore the problem is not a global frame failure, but the rigid retrieval
poses are still not docking-ready.

## E030 pivot retained as a constraint

The tested 13 deterministic micro-seeds remain a reusable primitive, but no
production batch caller uses them. Combining 13 seeds with a two-torsion beam
can reach roughly 403 evaluations per pose, over three million states for a
7,704-member refinement set. This was intentionally stopped before integration
because the user requires a general, fast retrieval stage.

## E031 implementation now on GitHub

Important commits, oldest to newest:

- `d7f7b9e`: pivoted from WEE1 pose repair to general fast reranking.
- `d45725a`: locked the final key-interaction match definition.
- `f6cd02d`: implemented and tested the E031 sidecar and workstation runner.

Implemented command:

```text
python -m aidd_agent.cli score-key-interaction-matches
```

For every named stored rigid objective pose, it computes:

```text
interaction_match_score = sum_i(weight_i * assigned_match_i)
                          / sum_i(weight_i)
```

Each pair match is bounded to `[0,1]` and combines exact feature type, Gaussian
spatial agreement, and signed or axial direction agreement. A deterministic
maximum-weight one-to-one assignment prevents one candidate feature from
satisfying several query anchors. Unmatched anchors score zero. Extra candidate
features do not penalize query coverage.

The independent NPZ sidecar retains:

- global, molecule, and conformer IDs in original order;
- one interaction score per original objective pose;
- assigned candidate feature index for every query anchor;
- every per-anchor contribution.

Its manifest records input hashes, parameters, wall time, score distribution,
median matched-anchor count, Spearman association with the original Gaussian
objective, and Top-100/500/1000 overlap. The runner also hashes each rigid result
before and after execution and fails if it changed.

Offline verification passed 123/123 tests. A local Windows in-memory kernel
measurement for 7,704 candidates x three objective poses (23,112 assignments,
three anchors and twelve candidate features) took 3.72 seconds, or 6,217
assignments/second. This is exploratory only: it excludes real artifact and
companion mmap reads and cannot accept the Linux performance gate.

Current interaction-class boundary: the materialized query anchors cover direct
hydrogen bonds. Salt bridges, aromatic stacking, cation-pi, hydrophobic contacts,
and metal coordination require separate target-independent extractors and
validation before being claimed as supported.

## Immediate workstation continuation

Activate the existing workstation environment, then run:

```bash
cd /mnt/medchem_taltio/wrk/yu_agent/Group-Agent
git pull --ff-only origin main
git rev-parse --short HEAD
export PYTHONPATH="$PWD/src${PYTHONPATH:+:$PYTHONPATH}"
bash scripts/run_e031_key_interaction_matching.sh "$PWD"
```

The revision must be `f6cd02d` or a later descendant containing it.

Default inputs are:

- artifact catalog:
  `$PWD/data/e019_artifacts/catalog.json`
- chemical companion:
  `/mnt/local/hand/yuzhang/aidd/chemical-companion-v1/catalog.json`
- QT9 query/result:
  `$PWD/data/e019_query_8bju/gaussian-query-v1.npz` and
  `$PWD/data/e020_qt9_validation/gaussian-staged-v1/refine/merged-scores.npz`
- 1X8B query/result:
  `$PWD/data/e026_query_1x8b/gaussian-query-v1.npz` and
  `$PWD/data/e026_1x8b_validation/gaussian-staged-v1/refine/merged-scores.npz`

If the runner reports `Required E031 input not found`, do not move or rebuild
data. Record the exact missing path; each default has an `E031_*` environment
override documented in `docs/key-interaction-matching.md`.

Default outputs are written under:

```text
/mnt/local/hand/yuzhang/aidd/e031-key-interaction-matches-v1
```

Resource reports can be collected with:

```bash
grep -E 'Elapsed|Maximum resident|File system inputs|File system outputs' \
  /mnt/local/hand/yuzhang/aidd/e031-key-interaction-matches-v1/*.time.txt
```

Return both query summaries and these resource lines. The QT9 runner compares
its E031 wall time against the accepted 77.99-second rigid-refine baseline and
prints `latency_gate_le_10_percent`. The locked threshold is 7.799 seconds.

## Decision after the real run

Accept the E031 latency hypothesis only if the same-workstation QT9 sidecar wall
time is at most 7.799 seconds. Regardless of speed, do not promote the score to
a production ranking based on the two WEE1 structures. First inspect:

- interaction score min/median/max and matched-anchor distribution;
- Spearman rho versus anchored Gaussian score;
- Top-100/500/1000 overlap;
- top interaction candidates' original Gaussian ranks and per-anchor scores;
- cold/warm page faults and maximum RSS.

If the score is cheap and non-redundant, the next research step is multi-target
held-out enrichment or redocking calibration. If it is slow, optimize companion
access/batching without changing the scientific definition. If it is redundant
or unstable, retain it as annotation or disable it. Receptor-specific pose
repair still stays downstream in every case.

## Files to read when resuming

1. `AGENTS.md`
2. `research-state.yaml`
3. this checkpoint
4. `experiments/E031-general-fast-chemical-reranking-protocol.md`
5. `docs/key-interaction-matching.md`
6. `findings.md`
7. the tail of `research_log.md`, `experiments.md`, and `benchmarks.md`

The authoritative large data remain on the Linux workstation under
`/mnt/local/hand/yuzhang/aidd`; they are not expected in GitHub.
