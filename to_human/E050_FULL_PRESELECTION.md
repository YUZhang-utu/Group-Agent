# Direct full-library preselection

User-authorized scope: all 25,813,808 catalog conformers, no pilot prerequisite,
Top-K, Top-N, or hit cap. This is a new output, never an ALL-run resume.

Synchronize the existing Group-Agent checkout from GitHub before launching.
The change adds the preselection module and launcher; it does not replace the old scoring modules.
Use the same activated workstation environment (NumPy, RDKit and project deps).

```bash
cd /mnt/medchem_taltio/wrk/yu_agent/Group-Agent
git pull --ff-only origin main
BASE=/mnt/local/hand/yuzhang/aidd/e050-full-preselection
mkdir -p "$BASE"
OUTPUT=$(mktemp -d "$BASE/run-$(date +%Y%m%d-%H%M%S)-XXXXXX")
SELECTION=/mnt/local/hand/yuzhang/aidd/prompt-workspace/users/workstation/projects/prj-252fa94197d6-prompt-aidd/runs/PROMPT-3e015383d6124ab4/execution/screening/screening/report.json
nohup bash scripts/run_e050_full_preselection.sh \
  --selection "$SELECTION" --output "$OUTPUT" \
  --workers 22 --chunk-size 2048 > "$OUTPUT/run.log" 2>&1 &
printf '%s\n' "$OUTPUT" > "$BASE/LATEST_RUN.txt"
echo "$OUTPUT"
```

Only run the launch block once. A family lock prevents simultaneous E050 jobs
under the same BASE; every run also has an output lock. Do not reinstall/update
code during a run. After interruption, verify the old process has stopped and
reuse the exact OUTPUT and settings; verified completed E050 chunks resume.

```bash
OUTPUT=$(cat /mnt/local/hand/yuzhang/aidd/e050-full-preselection/LATEST_RUN.txt)
tail -n 40 "$OUTPUT/run.log"
cat "$OUTPUT/progress.json"
cat "$OUTPUT/suite-report.json"
```

Initial full integrity checking precedes progress.json creation. Missing progress
during that stage is not grounds to launch another job. Errors are in run.log.

The selection supplies anchor IDs, score 0.5 and existing coarse thresholds;
match_mode becomes ANY in the new protocol without modifying the old selection.
All original seeds are retained in Gaussian winner competition. Optimistic ANY
seed filtering precedes exact anchor assignment; exact matches preserve each
observed pose-specific bitmask, not only the Gaussian winner. Gaussian scores are
conformer ranking metadata, not scores of every stored alternative pose.

Stages report input/output conformers, output unique molecules and rejections.
Worker timings separately measure artifact reads, size, type coverage, extents,
chemical reads/anchor bounds, seed generation, seed feasibility, Gaussian,
exact assignment and persistence. Worker time sums are not parallel wall time.
Overall wall time includes integrity checking and aggregation/grouping.

Outputs:
- suite-report.json: coverage, counts, stage reductions, timings, output hashes.
- candidate-poses.jsonl: one actual pose per molecule/exact anchor combination,
  with global ID, transform, assignments and scores. No molecule export cap.
- clusters.csv: exact non-stereochemical Murcko scaffold groups and display member.
- cluster-members.csv: all members; union masks across poses are explicitly not
  simultaneous interactions. Reconstruction failures are singleton groups.
- candidates.sqlite: indexed per-molecule/per-mask pose records and group membership.

This is scaffold grouping, not approximate fingerprint/3D similarity clustering.
No new Gaussian cutoff is invented, so that scoring stage can have zero rejection.
Observed stage selectivity is not guaranteed. The existing whole-ligand coarse
criteria remain additional eligibility restrictions; their recall cost is not
measured by this full run. Completion demonstrates full execution and measured
latency, not unpruned equivalence, docking performance or biological affinity.
