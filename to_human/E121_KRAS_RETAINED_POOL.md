# Existing KRAS pool: priority and unresolved docking inputs

This handoff uses only the completed direct-contact molecule pool. It does not
rerun full-library pose generation or load all library geometries. Registry
lookups resolve original structures for pool IDs. Raw MOL2 source files may
still need sequential reading to reach selected records; this is export I/O,
not a new chemical screen.

All IDs are retained without a top-N quota. E120 successes are priority only
after exact-source contact/distance recheck; failures are downgraded to unresolved.
Other IDs are unresolved, not rejected. Invalid source parsing or missing mapping
is quarantined explicitly. No claim of affinity or comprehensive recall is made.

## Stop and reuse

Stop the identified old job group, including its wrapper and worker processes.
Do not use a broad killall/pkill python command. Keep all old output files. The
new script snapshots committed SQLite results even if report.json says running.
The old in-progress chunk may be lost; committed chunks remain available.

## Run on the workstation

```bash
cd /mnt/medchem_taltio/wrk/yu_agent/Group-Agent
conda activate /mnt/medchem_taltio/wrk/yu_agent/Group_Agent/envs/aidd_workstation
export PYTHON="$(command -v python)"
git pull --ff-only origin feature/structure-guided-chat
nohup bash scripts/run_kras_retained_pool.sh \
  "$HOME/kras-direct-count-v1" "$HOME/kras-refine-v1/batch" \
  "$HOME/kras-retained-v1" \
  > "$HOME/kras-retained-v1.log" 2>&1 < /dev/null &
```

Progress: `tail -n 10 "$HOME/kras-retained-v1.log"`.
Mapping logs every 10,000 target-specific pool members. Export logs each source
completion. Completion requires export/report.json status complete, not merely
plan/report.json. Mapping and export duration are not yet measured on the real
workstation. Millions of MOL2 records may require substantial disk space/time.

Each target has export/TARGET/priority/ and export/TARGET/unresolved/ directories,
with molecules.csv and batch_NNN.mol2 files (at most 100,000 molecules each).
These are original coordinates, not docking-prepared or receptor-aligned poses.
Do not infer ten files per target: preserving the entire pool requires more.
The unresolved representative is deterministic, not optimized or scored.

The plan database contains all original membership, selected source mappings and
mapping issues. Final accounting requires priority + unresolved + quarantined
to equal each original target pool. Tiers are disjoint per target. Cross-target
overlap is expected. This version does not further reduce the unresolved pool;
additional pruning needs separately validated necessary conditions.

Use fresh output directories. If export fails after mapping completes, the plan
can be reused with `export --plan ... --output NEW_DIRECTORY`; export itself is
not resumable. Failed exports retain a running report and must not be called
complete. Source hashes and lineage mismatches stop the job.
