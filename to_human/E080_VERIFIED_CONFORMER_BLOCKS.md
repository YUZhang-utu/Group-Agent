# Verified conformer blocks: implementation and workstation handoff

E081 update: use `E081_GITHUB_FULL_LIBRARY_RUN.md` for the GitHub deployment and
mandatory validation workflow. The E080 shell wrapper now validates every variant
and stops on failed checks or unresolved review records. The previous local E080
ZIP is a historical snapshot; use the published source branch for current code.

## What is implemented

E080 connects the completed source audit to actual descriptor extraction, fitted
blocks, molecule-aware samples, incremental proposals and retrospective replay.
It replaces the E075 generic-ring path for this workflow. The old E075 pilot is
retained for reproducibility. No production search defaults or registry data change.

1. Validate audit output hashes and all supplied MOL2 source hashes. Every record
   must match the audited path, index, name and content hash. Every audited record
   must be visited exactly once. Audit conflicts stay in `review.jsonl`.
2. Extract directed N-CA-C(=O) peptide units, including proline side rings. Preserve
   source atom IDs and verify the main-ring atom set against the audit. Derive
   phi, psi and omega from exact source coordinates; do not infer artifact order.
3. Canonicalize only directed cyclic rotations. Chemical residue signatures choose
   the starting residue; equivalent sequences use omega state and geometric
   tie-breaks. Reverse traversal and mirror reflection are not allowed. Preserve
   symmetry and unverified-name flags. Near-symmetric tie-break stability still
   needs benchmarking.
4. Fit independent capacity-constrained median trees for 10000, 20000 and 30000
   conformers. SQLite stores descriptors and memberships; statistics use batches
   of 1024 rows and sorting can spill to disk. Hard strata preserve residue count,
   ordered alpha-CIP pattern and omega cis/trans/boundary pattern. These conservative
   strata can fragment the library; measure this before relaxing them.
5. Freeze descriptor schema, feature definition hash, input fingerprint, split
   rules, centroids, maximum fitted radii, capacity and exact memberships. The
   fitted radius is an engineering admission boundary, not an energy threshold.

The implementation is disk-backed, not a billion-molecule performance claim.
Serial RDKit extraction, repeated source hashing, SQLite sorting and multiple
descriptor copies consume real time and disk. Interrupted builds require a new
output. No automatic checkpoint resume or parallel descriptor extraction yet.

## Three descriptor variants

| Variant | Representation |
| --- | --- |
| backbone | Per-residue sin/cos of phi, psi and omega |
| chemistry | Backbone plus side-chain atom count, charge, N-methyl/proline indicator and typed feature counts |
| typed | Chemistry plus per-type mean and second moment of side-chain feature positions in local right-handed N-CA-C frames |

Feature channels: donor, acceptor, aromatic, hydrophobe, lumped hydrophobe,
positive-ionizable and negative-ionizable. Families remain separate, including
overlapping feature hypotheses. Side-chain features must lie fully within the
residue after removing backbone N/CA/C/O. Geometry is normalized by 5 angstrom,
feature count by 4 and side atom count by 10. These fixed versioned weights are
benchmark choices, not fitted scientific optima. Moments are a lossy summary,
not a guarantee of exact side-chain spatial correspondence or interaction validity.

## Identity

Without `--registry` and `--library-id`, IDs are explicitly `SOURCE-*` identifiers.
They are not production IDs. The optional registry join requires exactly one
match on library, molecule name, conformer name/index, source path/record index,
content/topology hashes and atom/bond counts. A mismatch goes to review; there is
no name-only fallback. Moved files require an explicit provenance migration.
The registry is opened read-only. Gaussian catalogs are not rewritten or assumed
to have the source coordinate precision or atom order.

## Workstation commands

Confirmed full source root: `/mnt/local/hand/yuzhang/aidd/mc_data`.
First deliver the E075-E080 source changes and install the project into the
scientific Python environment. A pull alone is insufficient while changes remain
uncommitted/unpushed. Use `AIDD_PYTHON=/path/to/python` when needed.

Run the complete audit preparation once, after source inventory review:

```bash
bash scripts/run_e079_macrocycle_prepare.sh /mnt/local/hand/yuzhang/aidd/mc_audits/e079-run-001
```

Then fit all three variants and all three capacities:

```bash
bash scripts/run_e080_macrocycle_blocks.sh \
  /mnt/local/hand/yuzhang/aidd/mc_audits/e079-run-001/audit \
  /mnt/local/hand/yuzhang/aidd/mc_blocks/e080-run-001
```

Append the actual registry SQLite path and library ID to join production identities.
Those values have not been provided, so the command above intentionally uses
source-only identities. Every output root must be new and outside source trees.
The wrapper currently rereads sources once per variant. It does not run docking
or online search. Review `counts`, every review record, block-size histograms and
source coverage before consuming memberships. Completion may include review rows.

## Sampling and evaluation

Each fit exports samples at 100 and 500 unique molecules per block, including all
their conformers in that block. Small blocks export all molecules. A molecule may
appear in several blocks; receipts distinguish block/molecule pairs, unique
molecules and conformer evaluations. Samples are deterministic and nested for a
fixed seed. Reuse scores by conformer ID across budgets.

`conformer_block_replay` accepts a complete conformer score CSV with columns
`conformer_id,score,status`. Larger scores must mean better results. Blank scores
require `no_surviving_pose`. Missing, extra or duplicate identities are errors.
Do not join external scores by name alone or mix scoring protocols.

```bash
python -m aidd_agent.conformer_block_replay \
  --model /path/to/typed/capacity-10000 --scores /path/to/complete_scores.csv \
  --output /path/to/replay-100.json --initial 100 --fraction 0.25 \
  --tail-fraction 0.05 --exploration 0.2 --top 100000
```

Repeat with 500 and multiple seeds on a fixed held-out panel, using equal actual
conformer evaluation budgets. Initial sampling must fit within the budget. Rank
blocks by a fixed upper-tail fraction of in-block molecule maxima; report effective
k. The comparator samples uniform conformers at the same actual scoring cost.
Report both molecule discovery recall and deduplicated molecule ranking recall.
The reference is the supplied complete panel, not biological ground truth or an
unseen library. Replay has an explicit one-million-conformer in-memory limit by
default; use bounded panels or deliberately adjust it after memory sizing. Tail
uncertainty estimation and live scheduling are not implemented.

## Incremental proposals

After separately auditing and extracting new conformers with the same schema:

```bash
python -m aidd_agent.conformer_block_store \
  --model /path/to/frozen/capacity-10000 --output /path/to/new-proposals \
  --descriptor-db /path/to/new-build/descriptors.sqlite \
  --schema-report /path/to/new-build/report.json
```

Existing conformer IDs, unknown strata, descriptor mismatches and out-of-radius
points are not assigned. Full leaves produce `capacity_overflow`; accepted entries
are proposals only. The original tree and memberships remain byte-identical.
Do not merge independent proposal batches without rechecking aggregate occupancy.
Splitting full blocks into a new model version is a later explicit rebuild.

## Local evidence and remaining acceptance

46 focused tests passed, including atom permutation and proper rigid-motion
invariance, proline, equivalent rotations, separate chemical channels, disk
statistics across chunk boundaries, deterministic capacity splits, exact registry
joins, incremental overflow, overlapping-molecule sampling and cost-matched replay.

All three variants extracted all 24 real mixed-fixture conformers (8 molecules),
with no review records. Capacity 4 produced 10 blocks; 10000/20000/30000 each
produced 9 small hard-stratum blocks. Evidence: `data/e080-final-backbone/`,
`data/e080-final-chemistry/`, `data/e080-final-typed/`. These tiny fixtures do not
compare large-capacity quality or throughput. Earlier `data/e080-typed-mixed/`
was an intermediate implementation run; final evidence is under `e080-final-*`.

Still pending: delivery and execution on the complete workstation library, real
registry joining, full-library resource measurement, completed score reference,
held-out recall/uncertainty comparisons and any live search integration. Do not
reject blocks during production screening based solely on this implementation.
