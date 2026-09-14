# E032 — expanded library one-time preprocessing

Date: 2026-09-11. Status: protocol written; Linux execution pending.

User-reported workstation inventory: 365 MOL2 files, 277.08 GiB, 24 usable CPUs,
65,274,164 kB total / 60,728,176 kB available memory, 6255.46 GiB free disk.
Remote source revision reported by the user: fb26c18. This Windows checkout is
d161d1c; no remote update or deployment is implied.

Objective: prepare reusable target-independent artifacts, chemistry companions,
pharmacophore postings and a USRCAT FAISS index. E031 query-dependent scoring,
docking, pose relaxation and prediction are out of scope.

Hypothesis: per-source checkpoints and chunked vector operations permit the
expanded library to be processed without changing the accepted two-shard baseline
or loading all molecular coordinates into RAM. This is a scaling/engineering
validation, not biological validation of new rankings.

Protocol:

1. Snapshot the source list, size and mtime; require stable completed uploads.
2. Use a NEW output directory. Validate old manifests/files and preserve the
   old global/molecule/conformer IDs; old source files must match their hashes.
3. Reconstruct a private build registry from the old two sources and IDs, then
   import new files transactionally. Fail on naming/index/topology conflicts;
   do not silently rename chemistry or omit invalid records.
4. Build new v1 shards with 20 workers, sequential source scheduling. Reuse
   completed outputs after integrity checks. Retain interrupted partials in
   a run-local recovery directory and rebuild only that unfinished shard.
5. Build E029 chemical companions and pharmacophore postings per source.
6. Train a shared IVF-PQ template using a reproducible uniform sample across
   ALL old/new shards; encode in bounded chunks. Save independent FAISS shards,
   then merge once. Do not store 365 cumulative full-index snapshots.
7. Run a one-new-source pilot before the full batch. Reuse the pilot artifacts
   in the subsequent all-source run. Indexing begins only after all inputs are
   prepared; the pilot does not train a biased partial-library codebook.

Acceptance: exact old ID prefix; immutable old manifests/files; unique contiguous
global ranges; source and output hashes; matching v1/chemistry/FAISS counts;
no unfinished state in published directories; successful identical-command
reuse; recorded wall time, failures and free space. Sampled exact-vs-approximate
recall must be calibrated after build before production search-speed claims.

Runtime and success for 277 GiB are unknown until measured. Historical 299,999
conformer timings are not an end-to-end estimate for this chemically broader batch.
