# Pre-docking aligned-pose and pocket QC

This stage applies the immutable candidate-to-query transforms from the
multi-cocrystal aggregation to artifact v1 heavy-atom coordinates. It annotates
obvious protein overlap before docking without changing retrieval or admission.

```bash
python -m aidd_agent.cli run-predocking-pocket-qc \
  --artifact-catalog /path/artifacts/catalog.json \
  --aggregation-dir /path/e026_aggregation \
  --output-dir /path/e028_predocking_qc \
  --receptor 8BJU-prepared-v1=/path/8BJU.cif \
  --receptor 1X8B-prepared-v1=/path/1X8B.cif \
  --top-n-per-query 100
```

Outputs are `pose-qc.jsonl`, one generic-element PDB point cloud per selected
task, `review.pml`, and a hash-bearing `manifest.json`. The review script loads
all point clouds disabled and enables the first pose; use `enable pose_*` or the
PyMOL object panel to inspect additional candidates.

The exported PDB files contain generic `X` points because artifact v1 has no
atomic identity or bond topology. They are visualization/QC artifacts only and
must never be supplied to a docking engine. A versioned, once-built artifact v2
companion is required for chemically correct SDF export, projected features and
terminal-torsion refinement.

For the accepted real 8BJU/1X8B experiment, run:

```bash
bash scripts/run_e028_predocking_qc.sh \
  /mnt/medchem_taltio/wrk/yu_agent/Group-Agent
```

The default receptor chain is the auth chain encoded in each query ID. Other
polymer chains in the deposited asymmetric unit do not contribute to the clash
proxy.
