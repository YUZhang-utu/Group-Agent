# E021 Gaussian overlay kernel and rigid-seed protocol

Status: mathematical kernel implementation completed -- offline confirmatory
validation passed; real QT9 candidate integration pending

## Hypothesis

A deterministic Gaussian overlap kernel with explicit operand direction,
feature typing, weights, and objective-specific rigid transforms can provide a
reproducible exact L2 scoring reference. Pharmacophore-pair rigid seeds should
recover a known translated/rotated local motif without changing the admitted
candidate set.

## Locked definitions

- Coordinates are row vectors in angstrom.
- A row-major homogeneous 4x4 matrix maps candidate coordinates into the query
  frame.
- Gaussian point-pair contribution is `exp(-d^2 / (2 sigma^2))`.
- Cross- and self-overlap use the identical approximation order and cutoff.
- Shape uses all supplied heavy-atom centers with unit weights.
- Color permits same-type pairs only. Alternative projected points divide one
  parent feature's weight rather than duplicating its functional weight.
- Tanimoto is `AB / (AA + BB - AB)`.
- Query-biased Tversky is `AB / (alpha * AA + beta * BB)`, where A is the
  query, B the candidate, alpha=0.95 and beta=0.05 by default.
- Atom-centered, projected, shape-only, and joint objectives retain independent
  transforms. Cross-pose score mixing is invalid.

## Protocol

1. Implement weighted Gaussian self/cross overlap with optional distance cutoff.
2. Implement same-type color overlap, Tanimoto, and directed Tversky.
3. Implement validated homogeneous transforms and weighted Kabsch alignment.
4. Generate proper-rotation pair seeds from compatible feature-type pairs with
   distance tolerance and sampled rotations about the aligned pair axis.
5. Deduplicate numerically equivalent transforms deterministically.
6. Score every seed without changing input candidate identity or admission.
7. Preserve primitive cross/self values alongside normalized scores.

## Confirmatory acceptance criteria

- An identical point cloud has ShapeTanimoto and query-Tversky equal to 1 within
  numerical tolerance.
- Translation/rotation does not change self-overlap.
- Incompatible feature types have zero color cross-overlap.
- Weighted alternative points preserve total parent functional weight.
- Kabsch and at least one compatible pair seed recover locked synthetic geometry.
- Every returned rotation is proper (`det(R)=+1`) and every transform is finite.
- Query/candidate operand order is visible in the output schema.
- Invalid dimensions, negative weights, nonpositive sigma, invalid cutoff, or
  alpha/beta values fail explicitly.
- The complete dependency-light suite remains green.

## Later workstation experiment

After offline mathematical validation, score a preregistered QT9 candidate
subset and compare shape-only, unweighted color, anchored atom-centered color,
and projected color. Runtime/enrichment claims require that real run and are not
inferred from synthetic tests.

## Offline results (2026-09-07)

- Implemented weighted Gaussian self/cross overlap with consistent cutoff and
  calculation order for shape and typed color.
- Implemented explicit query/candidate Tanimoto and directed query-biased
  Tversky primitives while retaining all raw cross/self overlaps.
- Implemented parent-feature alternative-point weight splitting, homogeneous
  transforms, proper weighted Kabsch alignment, and compatible pair seeds with
  six axial rotations and deterministic transform deduplication.
- Synthetic rigid transforms were recovered; identical shape/color scored one;
  incompatible color types scored zero; all generated rotations had determinant
  +1; invalid inputs failed explicitly.
- Nine focused tests and the complete 87-test dependency-light suite passed.
- Candidate-artifact batch orchestration, objective-specific optimization and
  real QT9 ranking remain pending and are not claimed by this result.
