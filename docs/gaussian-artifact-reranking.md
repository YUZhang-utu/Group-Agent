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
