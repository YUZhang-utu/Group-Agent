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
  --db /mnt/medchem_taltio/wrk/yu_agent/runtime/registry/aidd.sqlite3 \
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
interruption before deciding whether to resume or rebuild.

The implemented query-time primitives are:

- signed donor/acceptor and axial aromatic directional Gaussian overlap;
- element-aware vdW penetration and soft exclusion annotations;
- deterministic beam refinement over at most two smallest terminal torsions,
  with zero angle always retained.

Hard thresholds remain prohibited until WEE1 redocking and held-out enrichment
calibrate them.
