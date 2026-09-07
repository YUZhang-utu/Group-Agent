# Reusable local 3D pharmacophore retrieval

This stage finds library conformers containing part of a co-crystal ligand's
three-dimensional interaction pattern before docking. It complements whole-
molecule USRCAT/FAISS recall; it does not replace it and cannot remove an L1
candidate.

## Cost boundary

```text
Library version (expensive, once)
  MOL2 -> immutable coords/features/meta shards -> USRCAT FAISS
       -> shard-local pharmacophore-pair postings

Co-crystal query (small, once per query)
  structure + CCD -> query_manifest -> Reduce enrichment -> query-plan.json

Repeated search (fast)
  query plan + immutable postings + optional FAISS IDs -> tiered candidate NPZ
```

The pharmacophore index reads `meta.bin` and `feats.bin`; it does not read or
recompute the source MOL2 collection. Every artifact shard has an independent
index. When a new library shard arrives, only the new shard index is built.

## Build once for one library version

```bash
aidd-agent build-pharmacophore-index \
  --artifact-catalog /path/to/artifacts/catalog.json \
  --output-dir /path/to/indices/pharmacophore-pairs-v1 \
  --bin-width 0.5 \
  --max-distance 20.0
```

Repeating the command reuses matching completed shard indices. A different
source-manifest hash or parameter set is rejected instead of silently
overwriting the index.

## Compile one co-crystal query

Use the Reduce-enriched manifest when available:

```bash
aidd-agent compile-pharmacophore-query \
  --query-manifest /path/to/query_manifest.reduce.json \
  --output /path/to/queries/8BJU-QT9/pharmacophore-query-v1.json
```

This command reads no library data. The query plan contains compatible anchor
pairs, invariant distances, three search profiles and a canonical query hash.
Translation and rotation of the crystal coordinates do not change the query.

## Search all evidence tiers at once

```bash
aidd-agent search-pharmacophore-index \
  --index-catalog /path/to/indices/pharmacophore-pairs-v1/catalog.json \
  --query-plan /path/to/queries/8BJU-QT9/pharmacophore-query-v1.json \
  --external-l1-ids /path/to/queries/8BJU-QT9/faiss-global-ids.npy \
  --output /path/to/queries/8BJU-QT9/pharmacophore-hits.npz
```

`--external-l1-ids` accepts NPY, NPZ (`global_ids` or `ids`), JSON, or a text
file and is optional. Every supplied ID is retained in the output union.

The NPZ arrays are:

- `global_ids`: sorted union candidate IDs;
- `source_flags`: bit 1 pharmacophore, bit 2 external L1;
- `highest_tier`: 0 external-only, 1 loose, 2 balanced, 3 strict;
- `matched_pair_counts`: one column per profile.

The adjacent `.manifest.json` records input hashes, profile definitions and
counts. Profiles are nested:

| Profile | Pair tolerance | Required distinct query pairs |
|---|---:|---:|
| loose | +/- 2.0 A | 1 |
| balanced | +/- 1.0 A | up to 2 |
| strict | +/- 0.5 A | up to 3 |

The required count is capped by the query's available pair count. One anchor
alone has no rotation/translation-invariant 3D geometry and therefore does not
fabricate a local-spatial hit; other broad recall channels remain available.

## Interpretation

This is an admission channel, not a binding score. A loose hit states only that
at least one compatible feature pair occurs at a similar separation somewhere
in a stored conformer. Exact feature assignment, rigid alignment, Gaussian
shape/color, query-biased Tversky, protein exclusion volume and limited torsion
refinement are later reranking stages over this unchanged union.

The exact reference primitives for the next stage are implemented in
`aidd_agent.gaussian_overlay`. They preserve raw self/cross overlaps, explicit
query/candidate Tversky direction and the transform used for each score. Batch
connection to real candidate artifacts remains a separate validation step.

## Complete QT9 workstation validation

After pulling and reinstalling the package, run the checked-in wrapper with the
seven real paths used by the existing E019 run:

```bash
bash scripts/validate_qt9_pharmacophore.sh \
  /path/to/artifacts/catalog.json \
  /path/to/incremental-faiss-index-directory \
  data/e019_query_8bju/query_manifest.reduce.json \
  data/e019_query_8bju/8BJU.cif \
  data/e019_query_8bju/QT9.cif \
  /path/to/reusable-indices/pharmacophore-pairs-v1 \
  data/e020_qt9_pharmacophore_validation \
  QT9
```

The first execution builds the reusable pair index. Later QT9 executions reuse
it; other co-crystal ligands reuse the same index and compile only a new query.
The wrapper also generates a fresh 100,000-ID FAISS baseline using `nprobe=256`,
unions it with all partial-motif hits, validates that every FAISS ID survives,
and writes:

- `faiss-global-ids.npy`;
- `pharmacophore-query-v1.json`;
- `pharmacophore-hits.npz`;
- `pharmacophore-hits.manifest.json`;
- `validation-report.json`;
- `validation.log`.

Acceptance requires `accepted: true`, no failed checks, valid stable-ID ranges,
nested tiers, manifest/array count agreement, and complete preservation of the
FAISS baseline. Candidate counts and timing are measurements, not fixed pass
thresholds in the first real run; record them before setting performance gates.
