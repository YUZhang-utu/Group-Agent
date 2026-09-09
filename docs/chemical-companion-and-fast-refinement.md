# One-time chemical companion and fixed-budget refinement

Artifact v1 remains the immutable broad-recall source. E029 adds a sidecar keyed
by the same global conformer IDs. Building it requires one sequential pass over
each source MOL2 shard; later co-crystal queries mmap only their retained Top-N.

The companion stores heavy-atom elements, charges, aromatic/chiral flags,
heavy-heavy bonds, feature atom membership, donor/acceptor directions, aromatic
normals, and bounded terminal torsions. It does not alter FAISS or
pharmacophore membership.

```bash
python -m aidd_agent.cli build-chemical-companion \
  --library LIB-AFA68EE6888C \
  --artifact-catalog /path/e019_artifacts/catalog.json \
  --output-root /mnt/local/hand/yuzhang/aidd/chemical-companion-v1 \
  --workers 16 \
  --max-moving-atoms 12
```

If a v1 shard records a source path from another machine, provide its relocated
source explicitly without changing registry identity:

```bash
--source shard-name=/current/path/source.mol2
```

Repeat `--source` for each relocated shard. Completed shards are reused. An
existing `.SHARD.partial` is never deleted automatically; inspect it after an
interruption before deciding whether to rebuild. Shard-internal resume is not
implemented: after preserving any diagnostics, move the partial directory out
of the output root and rerun. Missing sources, source SHA-256 mismatches, and
registry/artifact identity mismatches are validated before a partial directory
is created.

The immutable artifact-v1 catalog supplies global, conformer, and molecule
identity in source-record order. Exact whole-file source SHA-256 plus per-row
artifact shape checks make a Windows build registry unnecessary on the Linux
workstation. `--db` remains optional only as redundant validation when the
exact registry used to create artifact v1 is deliberately available; never
point it at an unrelated active workstation registry.

The implemented query-time primitives are:

- signed donor/acceptor and axial aromatic directional Gaussian overlap;
- element-aware vdW penetration and soft exclusion annotations;
- deterministic beam refinement over at most two smallest terminal torsions,
  with zero angle always retained.

Hard thresholds remain prohibited until WEE1 redocking and held-out enrichment
calibrate them.
