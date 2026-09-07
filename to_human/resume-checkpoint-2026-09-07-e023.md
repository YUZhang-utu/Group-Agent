# Resume checkpoint — E023 staged parallel Gaussian reranking

## Verified implementation

- E022 real QT9/1,000-candidate runs succeeded on the university workstation.
- Pair mode: 219.70 seconds wall; PCA-only: 16.59 seconds wall; both used about
  one logical CPU and 55 MB RSS.
- E023 implements all-candidate PCA coarse scoring followed by pair refinement
  of the union of three objective Top-N lists.
- Query/candidate self-overlaps are cached instead of recomputed per pose.
- Worker processes open independent read-only artifact maps; chunk completion
  order does not alter merged candidate order.
- Chunk NPZ and manifests are written atomically and checked by configuration,
  ID-slice, and output hashes. Identical commands resume; corrupt or missing
  chunks alone are recomputed.
- Offline one-vs-two-worker, selection, corruption-repair, and configuration
  mismatch tests passed. Complete suite: 92 passed at `4945d4d`.

## Next workstation validation

1. Pull the latest GitHub `main` and confirm 92 tests.
2. Run the locked 1,000 IDs with 1 and 16 workers in separate output dirs and
   compare merged arrays and wall time.
3. Start a 16-worker coarse run, interrupt after several chunks, rerun the same
   command, and confirm `resume_reused` is nonzero.
4. Complete the 299,999-conformer PCA stage.
5. Run refinement with `top_n_per_objective=5000` and inspect objective overlap
   and poses before exclusion-volume or docking work.

Exact commands and output paths are in
`docs/gaussian-artifact-reranking.md`.

## Boundary

The full coarse result preserves every admitted conformer. Top-N union controls
only allocation of expensive pair-seed compute. It is not a hard biological
filter. Artifact v1 still does not support projected candidate color,
exclusion-volume scoring, or terminal-torsion refinement.
