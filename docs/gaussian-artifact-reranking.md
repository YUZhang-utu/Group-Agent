# Artifact-backed Gaussian reranking

This stage consumes stable global conformer IDs and the immutable artifact
catalog. It never reopens the source MOL2 library and never removes an admitted
L1 candidate.

## 1. Build the co-crystal query once

```bash
aidd-agent prepare-gaussian-query \
  --mmcif "$QUERY_DIR/8BJU.cif" \
  --ccd "$QUERY_DIR/QT9.cif" \
  --query-manifest "$QUERY_DIR/query_manifest.json" \
  --output "$QUERY_DIR/gaussian-query-v1.npz"
```

This writes `gaussian-query-v1.npz` and
`gaussian-query-v1.manifest.json`. The package contains crystal heavy-atom
coordinates, RDKit atom-centered features, types, anchor weights, and hashes.

## 2. Score an immutable candidate set

The candidate file may be the `pharmacophore-hits.npz` produced by the E020
validator because it contains `global_ids`:

```bash
aidd-agent score-gaussian-candidates \
  --artifact-catalog "$ARTIFACT_CATALOG" \
  --query "$QUERY_DIR/gaussian-query-v1.npz" \
  --candidate-ids "$VALIDATION_OUTPUT/pharmacophore-hits.npz" \
  --output "$VALIDATION_OUTPUT/gaussian-scores-v1.npz"
```

For an initial runtime check, make a preregistered prefix without changing the
full L1 artifact:

```bash
python -c "import numpy as np; p='$VALIDATION_OUTPUT/pharmacophore-hits.npz'; o='$VALIDATION_OUTPUT/gaussian-runtime-1000.ids.npy'; z=np.load(p); np.save(o,z['global_ids'][:1000])"
```

Then pass `gaussian-runtime-1000.ids.npy` to `--candidate-ids`. The scorer emits
one row for every input ID in the same order. It stores separate best transforms
and complete overlap primitives for shape-only, unweighted atom-centered joint,
and anchor-weighted atom-centered joint objectives.

## Current boundary

Artifact schema v1 contains heavy-atom centers and atom-centered pharmacophore
features. It does not contain atomic radii/topology or directional candidate
projection points. Therefore this implementation is a deterministic reference
reranker; it does not claim projected-color, exclusion-volume, terminal-torsion,
or production-scale native-engine performance.

## Staged parallel production path

Use one persistent output directory. Scientific parameters are locked in
`run-manifest.json`; rerunning the identical command reuses every valid chunk.
A changed input or parameter requires a new output directory.

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
STAGED_OUTPUT="$VALIDATION_OUTPUT/gaussian-staged-v1"
```

Run the complete PCA-only coarse stage with 16 worker processes:

```bash
python -m aidd_agent.cli run-staged-gaussian-reranking \
  --artifact-catalog "$ARTIFACT_CATALOG" \
  --query "$GAUSSIAN_QUERY" \
  --candidate-ids "$PHARMA_HITS" \
  --output-dir "$STAGED_OUTPUT" \
  --stage coarse --workers 16 --chunk-size 1000 \
  --top-n-per-objective 5000 --max-pair-seeds 512
```

After coarse completion, run refinement with the same locked parameters:

```bash
python -m aidd_agent.cli run-staged-gaussian-reranking \
  --artifact-catalog "$ARTIFACT_CATALOG" \
  --query "$GAUSSIAN_QUERY" \
  --candidate-ids "$PHARMA_HITS" \
  --output-dir "$STAGED_OUTPUT" \
  --stage refine --workers 16 --chunk-size 1000 \
  --top-n-per-objective 5000 --max-pair-seeds 512
```

Power loss or interruption requires no special repair: rerun the same command.
Valid chunks are skipped; missing, partial, checksum-invalid, or ID-invalid
chunks are recomputed. Main outputs are:

- `coarse/merged-scores.npz`: every admitted conformer;
- `refine/selected-candidates.npz`: per-objective Top-N union and provenance;
- `refine/merged-scores.npz`: pair-refined selected conformers;
- `run-manifest.json`: locked configuration and completion state.

`--stage all` runs both stages in one invocation. `--no-resume` deliberately
recomputes all chunks but retains atomic writes.
