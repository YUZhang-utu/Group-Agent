# E096: backbone, chemistry and spatially coupled conformer work blocks

## Scope and inputs

Keep E094 and E095 unchanged. E095 completed with 1,188 regular blocks, 5,011 to
94,822 members, and 132,885 unchanged special conformers (user-supplied receipt).
E096 is an alternative joint description built under the original 384 E094
backbone work parents. It does not recursively split all 1,188 E095 children.

The source is 24,530,851 accepted regular conformers. The separate chemistry-review
population remains deferred and the existing special-exhaustive pool is unchanged.
Cached E094 source locators and 69 broad property/steric channels are reused.
Original MOL2 coordinates must be read again because spatial positions cannot be
recovered from the old mean/std/max summaries. Every selected record reproduces
its source hash, names, units, hard group and stored backbone descriptor before
its spatial features are accepted. Use the same RDKit/feature definitions as the
original build. Strict chemistry failures stop extraction; no silent omission.

## Joint description (317 float32 channels)

- 69 existing broad chemistry and local steric channels.
- 24 actual backbone shape channels: directed torsion summaries, alpha-carbon
  chord distances, backbone principal spreads and signed local volume.
- 224 chemically typed spatial channels. For each directed residue anchor,
  construct a proper frame from N, CA and carbonyl C. Associate feature centers
  with that residue, its next residue, next-two residue and opposite residue.
  For each association, retain occupancy and signed coordinate mean/spread in
  the anchor frame. Seven RDKit families plus the full side-heavy-atom cloud
  are included. Residue-associated donor/acceptor centers include backbone sites.

This ties chemistry to actual geometry across residues instead of merely
counting side-chain types. Proper whole-molecule rotation, translation, directed
cyclic origin and atom renumbering do not create artificial differences.
Reflection/reverse-traversal alignment is never used. Existing chemical states
and chirality identities remain available in the parent overlay.

These are finite spatial summaries, not a collision-free encoding of all atom
positions. They do not compute hydrogen-bond vectors, solvent accessibility,
energy/strain, protonation ensembles or receptor contacts. Such quantities must
not be inferred from the partition. The original source remains authoritative
for detailed atom-wise comparison. Mixed parent execution pools are not a claim
of identical backbone geometry.

## Fragmentation and internal checks

One global transform assigns distance budgets of 20% broad chemistry, 10% local
sterics, 35% backbone geometry and 35% typed spatial geometry. These are fixed
engineering defaults, not validated scoring weights. Global standard deviation
below 0.02 disables a channel; the scaling denominator has a 0.10 floor. Reports
expose active channels and mean feature counts so inactive chemistry is visible.

Every new child must contain at least 5,000 conformers and at least 20% of its
original E094 parent. Maximum four children per parent: at most 1,536 regular
blocks for the current 384 parents. No forced split, maximum capacity, discarded
tail or expansion of the special pool. Required internal dispersion improvement
is 5% on both fit and check with check gain at least 75% of fit gain. The standardized
centroid effect must reach 0.35. Accepted and rejected proposals are recorded.

Unlike E095, check assignment uses **source molecule ID**, keeping all conformers
of a molecule on the same side across parents. Source registry identity is trusted;
this is not a scaffold-disjoint test. Global scaling and adaptive reuse of checks
still make this internal evidence, not an unbiased final recall/docking evaluation.

## Overnight workstation command

Run in the existing scientific environment, from the Group-Agent checkout:

```bash
git pull --ff-only origin feature/structure-guided-chat
BASE=/mnt/local/hand/yuzhang/aidd/mc_runs/e081-001/audit-stage/blocks-20000
export AIDD_SPATIAL_WORKERS=8
nohup bash scripts/run_e096_joint_spatial_blocks.sh \
  "$BASE/final-blocks-e094" "$BASE/backbone" "$BASE/final-blocks-e096-joint" \
  > "$BASE/final-blocks-e096-joint-launch.log" 2>&1 &
echo $! > "$BASE/final-blocks-e096-joint.pid"
```

The Linux runner uses flock to prevent simultaneous writers. Do not pull new code
while this job runs. A fresh output root is required for changed code or inputs.
After confirming an interrupted process is stopped, repeating the identical
command resumes committed extraction batches and completed partition parents.
Worker count can change on resume; feature definitions and numerical versions cannot.
The old E094/E095 databases are opened read-only. Logs append; final manifests
contain source/code/output hashes. An incomplete run cannot claim full coverage.

Progress:

```bash
tail -n 20 "$BASE/final-blocks-e096-joint-launch.log"
cat "$BASE/final-blocks-e096-joint/profiles/progress.json"
```

Startup verifies hashes and plans source locators before progress.json appears.
Extraction commits every 256 conformers. The expensive stage reads raw MOL2 and
computes features; it is not the six-minute cached E095 refinement. Full-library
runtime has not been measured, and next-morning completion is not guaranteed.
The float32 vector payload alone is about 29 GiB; SQLite rows, indexes, membership,
temporary sorting and journal space add overhead. Retain the source files.
Use a local output filesystem with adequate free space (approximately 100 GiB
headroom is a conservative planning allowance, not a measured disk requirement).

## Results to return

- profiles/report.json: complete spatial coverage, width, source identities and hashes.
- blocks/report.json: membership gate, block sizes/count, feature-channel statistics.
- blocks/split-decisions.json: split/rejection reasons and internal gains.
- blocks/blocks.csv and blocks/property_blocks.sqlite: exact enumerable memberships.

Final validation enumerates every new membership against the complete verified
spatial cache and checks parent ownership and conserved populations. It does not
reuse old distance bounds or authorize skipping unsampled blocks. The next step
remains held-out 3D retrieval, PLANTS/MDM2 and N-E evaluation.
