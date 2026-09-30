# Final work blocks without tiny regular groups

E091 reported 24,663,736 admitted conformers: 23,279,397 unchanged definite,
1,251,454 provisional boundary assignments and 132,885 special (0.538787%).
The separate 1,150,072 chemistry-review records are not part of this stage.

## Work-block policy

The default minimum regular work-block size is **5,000 conformers**. There is
no upper-capacity split. Large target classes stay intact; small classes pool
by directed cis/trans pattern and length first. Remaining small pools attach
to an existing block, preferring matching pattern, then matching length. A
mixed execution pool retains every original class/chirality identity. It is
not a declaration that different cis/trans states are chemically equivalent.
Every regular work block meets the floor unless the entire regular library
is smaller than the floor, in which case it remains one explicit exception.
Special members remain in one exhaustive pool, never silently redistributed.

The output is an enumerable SQLite overlay, not just an aggregate plan. It
contains every boundary override and a class-to-block map for all unchanged
ordinary records. Indexed joins recover CID, MID, old leaf, original class,
provisional target, rotation and final work block. Original databases are
read-only. Old centroid/radius bounds must not reject these enlarged blocks.

## Run the backbone completion first

From the existing Group-Agent checkout and scientific environment:

```bash
git pull --ff-only origin feature/structure-guided-chat
export PYTHONPATH="$PWD/src${PYTHONPATH:+:$PYTHONPATH}"
BUILD=/mnt/local/hand/yuzhang/aidd/mc_runs/e081-001/audit-stage/blocks-20000/backbone
OUT=/mnt/local/hand/yuzhang/aidd/mc_runs/e081-001/audit-stage/blocks-20000/final-blocks-e094
bash scripts/run_e094_work_blocks.sh base "$BUILD" "$OUT"
```

The runner assumes the sibling paths used in E086/E088/E091:
`backbone-fragmentation-review`, `backbone-boundary-candidates-v1`, and
`backbone-calibrated-routing-v1`. If the actual routing directory differs,
use the module CLI with `--routing` set to the actual completed directory.
No E091 rerun is required. Use a fresh output root. Hashing large input files
and enumerating the full membership take time; progress prints every 4 GiB
hashed, 50,000 boundary rows and 1,000,000 validation records.

Inspect:

- `backbone/report.json`: number of work blocks, minimum/maximum size,
  mixed-pool counts, special fraction and population conservation.
- `backbone/blocks.csv`: work blocks, population, original class count and
  preserved omega patterns. `classes.csv` maps every definite target class.
- `backbone-validation.json`: require `structural_gate: passed`.
- `backbone/work_blocks.sqlite`: exact indexed boundary overlay and class map.

This metadata/identity stage does not recompute raw coordinate descriptors.
It verifies model/candidate/routing hashes and every boundary CID and target.
The independent validator enumerates all work memberships.

## Side-chain properties with the same size floor

After base validation passes:

```bash
export AIDD_PROPERTY_WORKERS=8
bash scripts/run_e094_work_blocks.sh properties "$BUILD" "$OUT"
```

Or use mode `all` with a fresh output to run both phases sequentially. Keep
the workstation RDKit version and feature definitions used by the source
build. The script caps nested BLAS/OpenMP threads to one per process. Increase
workers only according to available CPUs/RAM; runtime is not yet benchmarked
for the full library.

E091 features are reused for covered regular records. Missing ordinary records
need source-verified profiles; the existing 1.38M cached profiles cannot stand
in for the entire 24.5M regular population. This is the expensive stage.
Completed profile batches are committed to `profiles/profiles.sqlite`. A rerun
of the properties command resumes after an extraction interruption with the
same code and input signatures. It never claims an interrupted cache is complete.
Changing code/inputs requires a new profile output. Refinement itself requires
a fresh `properties` directory; an interrupted refinement is not resumable.

The descriptor uses continuous per-side-chain chemical properties and actual
coordinate-based local bulk/branching/extent proxies. Mean, standard deviation
and maximum summarize the directed positions without demanding exact residue
identity. It is broad composition/steric similarity, not sequence equivalence
or a protein clash calculation. The full original identities remain traceable.

Refinement is optional and adaptive: it only splits a parent if both children
have >=5,000 members and normalized centroid contrast >=0.35. At most eight
children are allowed per parent; this is a ceiling, not a requested cluster
count. Identical/weakly differing profiles and small tails remain together.
The contrast threshold is an exploratory engineering default, not a validated
activity or physical-state boundary. Full profile coverage is mandatory.

Inspect `properties/report.json`, `properties/blocks.csv` and
`properties-validation.json`. Require `property_gate: passed`, conserved
regular/special populations and zero below-floor regular blocks at full scale.

## Export complete block identities

Choose an actual block ID from the corresponding CSV:

```bash
python -m aidd_agent.work_block_delivery export \
  --blocks "$OUT/backbone" --block-id ACTUAL_BLOCK_ID \
  --output "$OUT/selected-block-members.jsonl"
```

Add `--properties "$OUT/properties"` when selecting a refined property block.
`--block-id special-exhaustive` exports the preserved special pool. This is
an identity manifest, not a ligand coordinate export, docking run or production
registry join. No unassessed class is rejected on the basis of pooled samples.

## Validation and practical boundary

Local tests exercise tail pooling, preserved labels, invalid CIDs/rotations/
scores/targets, hash changes, full enumeration, property coverage, size floors,
ties, and property exports. A real 24-conformer replay preserves all 21 regular
and 3 special records, reuses three E091 profiles, computes 18 missing profiles,
and retains 11 committed profiles across an intentional interruption. It uses
a test minimum of five because this fixture is small; production default is
5,000. See `E094_LOCAL_VALIDATION.json`.

Full workstation completion is established only by its new reports and passed
validation receipts. Search recall, special-pool timing, PLANTS/MDM2, N-E scores
and affinity models are subsequent evaluations. The former 20,000 cap is gone;
runtime batch sizes remain independent of logical work-block identity.
