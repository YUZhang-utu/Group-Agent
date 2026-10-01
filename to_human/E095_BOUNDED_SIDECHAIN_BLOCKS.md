# E095: bounded backbone and side-chain work partitions

E094 retained 384 regular blocks. Root replay of the largest ten showed real
feature variation but raw RMS contrasts about 0.10, below the uncalibrated 0.35
gate. E095 uses the existing 69-channel broad side-chain cache without reparsing
MOL2 or changing frozen E094 code, chemical identities or special membership.

## Fixed engineering policy

- Keep every child inside its original backbone work parent. Mixed execution
  pools still retain their underlying distinct chemical/backbone identities.
- Fit one global feature transform to all regular cached profiles. Drop channels
  with global standard deviation below 0.02. Standardize with a 0.10 floor.
  Give chemistry and steric feature groups equal total squared-distance weight.
  Correlated mean/std/max channels are not whitened; this remains exploratory.
- Fit a PCA projection and median threshold on at most 8,192 deterministic fit
  conformers. Keep about 20% of CIDs as an internal check using CRC32 routing.
- Require at least 5% dispersion reduction in fit and check sets, check gain at
  least 75% of fit gain, and weighted standardized centroid distance at least 0.35.
  This distance is not the old 69D raw RMS and the numerical thresholds are not
  validated biological boundaries.
- Every new child has at least 5,000 members AND at least 20% of its original
  parent. Keep at most four children per parent. With 384 parents there can be
  no more than 1,536 regular blocks. This is a ceiling, not a target or capacity.
- Reject unhelpful proposals and retain the parent. Do not discard rare members,
  move them to the special pool, or create tiny independent tails.

Checks split **conformers**, not unique molecules. Shared molecules may occur
on both sides; fitting global scales and reusing adaptive checks also prevent
claiming an unbiased final evaluation. Independent molecule-held-out retrieval,
docking and N-E evaluation are still required. These broad summaries do not
recover side-chain positions on the ring and cannot establish pocket fit.

## Workstation execution

Both previous jobs have completed. From the Group-Agent checkout, in the same
scientific Python environment:

```bash
git pull --ff-only origin feature/structure-guided-chat
BASE=/mnt/local/hand/yuzhang/aidd/mc_runs/e081-001/audit-stage/blocks-20000
bash scripts/run_e095_property_blocks.sh \
  "$BASE/final-blocks-e094" "$BASE/final-blocks-e095-properties"
```

The first pass hashes the existing profile database, verifies the parent receipt
and computes global scales. Then it processes parents one at a time and commits
each complete parent transaction. It finishes by enumerating all new memberships
against the cached source parent and checking complete population conservation.
It does not rescan the 100-GiB frozen conformer model. The supplied passed E094
receipt and its source-owned cache are prerequisites, not new coordinate checks.

If interrupted, confirm the old process has stopped, then repeat the same command.
It resumes committed parents with unchanged source hashes, code, policy and NumPy
version. Never run two writers against the same output. Full source hashes are
rechecked on resume. A changed implementation requires a fresh output directory.
RAM scales with the largest parent; the existing maximum is 373,795 profiles.
Runtime has not been measured on the full workstation library.

## Review results

- report.json: added/total blocks, minimum/maximum populations, decision reasons,
  preserved special population, membership gate, source and output hashes.
- blocks.csv and property_blocks.sqlite: enumerable child assignments compatible
  with the existing work-block identity exporter.
- transform.json: global scales and channel weights.
- split-decisions.json: accepted and rejected proposals, fit/check gains, sizes,
  projection rules and effects. Reaching the child ceiling stops further trials.

Require membership_gate=passed_against_complete_cached_profiles, unchanged regular
and special counts and blocks_below_minimum=0. A complete run may legitimately
retain some or all parents. Do not lower thresholds merely to force subdivisions.
Search dispatch/rejection and side-chain positional modeling are not delivered by
this stage. The chemistry-review population remains deferred.
