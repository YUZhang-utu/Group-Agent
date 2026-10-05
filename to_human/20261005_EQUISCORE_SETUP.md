# EquiScore installation and saved-pose evaluation

User decision: use EquiScore as the public rescoring baseline. N-E remains a
separate user model awaiting its code. This delivery provides installation,
preflight, technical pilot, full saved-pose inference and Top-10 block analysis
through the terminal. Natural-language Chat dispatch for EquiScore is not yet
registered. Do not interpret an existing ChemPLP Chat report as model execution.

## Install on the Linux GPU workstation

From the existing repository and AIDD environment:

```bash
cd /mnt/medchem_taltio/wrk/yu_agent/Group-Agent
git pull --ff-only origin feature/structure-guided-chat
bash scripts/install_equiscore.sh
```

The default dedicated prefix is `/mnt/local/hand/yuzhang/aidd/tools/equiscore`.
An optional first argument selects another absolute, dedicated prefix. Conda,
Git, network access and an NVIDIA driver are prerequisites. No sudo is used and
the active AIDD environment is not modified. Existing unmanaged directories are
rejected. Repeating installation in the script's managed prefix is supported;
existing profile receptor overrides are preserved.

The recipe targets Linux x86_64: Python 3.9, PyTorch 1.11 CUDA 11.3, DGL 0.9.1,
RDKit 2022.9.5 and ProLIF 1.1.0. Versions deliberately match the older upstream
APIs. New GPU architectures may require a separately tested environment recipe;
a driver check alone does not establish CUDA kernel compatibility.

Outputs are `profile.json`, `installed-packages.txt` and `doctor.log` in that
prefix. Success requires strict checkpoint loading and CUDA tensor/graph kernels,
then prints `status: ready`. A failure stops installation; retain the error log.
The installation recipe has been syntax-checked locally, not installed or run
on the user's Linux GPU workstation yet.

## Pilot, full inference, then analysis

Use the existing AIDD Python 3.11+ environment to launch these commands. The worker
uses the separate Python 3.9 environment named in `profile.json`. If needed, set
`AIDD_PY` to the absolute existing AIDD interpreter; do not activate the EquiScore
environment for the coordinator.

```bash
bash scripts/run_e103_equiscore.sh pilot
```

The preset reuses the completed `plants-parallel20/report.json` in run
`PROMPT-1d4f82d50d294925`. It selects 32 CIDs independently of ChemPLP, across
all receptors (96 pairs for the user's three-receptor panel). It verifies the
original output seals, which may take substantial I/O before any GPU work.
Inspect `execution/equiscore-pilot/report.json` and `worker.log`. Require
`status: complete`, `scored: 96`, `failed_pairs: 0`. This is a compatibility gate,
not evidence of predictive accuracy or applicability to all macrocycles.

```bash
bash scripts/run_e103_equiscore.sh full
bash scripts/run_e103_equiscore.sh analyze
```

Full scoring requires the successful pilot under identical source, profile,
worker, model, environment and coordinator code. It consumes all saved poses,
not just the ChemPLP Top-10, and makes no new docking calls. One pose is processed
at a time to preserve explicit identity and bounded GPU memory. Runtime for
821,619 pairs is unmeasured; the pilot informs practical scheduling. This initial
worker is correctness-oriented and does not promise parallel screening speed.

Outputs sit beside `plants-parallel20`: `equiscore-pilot`, `equiscore-full`,
`equiscore-analysis`. Overrides: `PLANTS_REPORT`, `EQUISCORE_PROFILE`,
`EQUISCORE_OUTPUT` (parent of those three directories).

An interrupted run can resume with the same command and sealed inputs. Stored
successes and explicit failures are not recomputed; absent rows are resumed.
After correcting a failed-input or environment issue, use a new output parent
and repeat the pilot. A complete receipt is reused without launching inference.
Use a fresh analysis directory for a repeated analysis.

## Output interpretation

- `equiscore-full/scores.csv`: original CID/receptor/pose identity and ChemPLP,
  EquiScore value, explicit status and error. Failed values are null, never zero.
- `equiscore-analysis/block_rankings.csv`: mean of the highest-scoring ten
  distinct molecules, each represented by its best in-block conformer. Higher
  EquiScore is better. Top-5 summaries are also retained. Group by scheme/receptor.
- `block_top_candidates.csv`: ranked molecules and original pose references.
- `rank_comparison.csv`: independent ChemPLP and EquiScore block ranks, rank
  changes and overlap of their Top-10 molecules. No mixing of numerical scales.
- `chemplp_baseline/`: reproducible reference ranking on the same sample panel.

Partial full runs can be analyzed diagnostically; their primary block ranks are
withheld. Insufficient distinct molecules are always unranked. No experimental
EF1, affinity or superiority claim is produced without appropriate labels and
validation. The original ligand chemistry review remains outstanding.

PDB fallback: the original receptor's adjacent `protein-input.pdb` is discovered
from the sealed preparation profile. If unavailable, add `receptor_pdbs` to the
EquiScore profile with keys `1RV1_B`, `7NA2_A`, `5J7F_A` and paths to the original
aligned PDB files. Do not substitute unrelated RCSB coordinates: the worker
checks the heavy-atom inventory and coordinate frame against the prepared MOL2.
Conversion rejects chemistry/stereochemistry differences and does not silently
repair ligand structures or remove metals.

## Provenance and checks

Official repository: https://github.com/Intelligent-Drug-Discovery-Lab/EquiScore

Paper: https://www.nature.com/articles/s42256-024-00849-z

Frozen commit: `8b2a9289cf7d181fa49de6ac6712260e8c500c4a`.
Checkpoint: `workdir/official_weight/save_model_screen.pt`.
SHA256: `d4367bb73686b2363e238abb778fab55e2924458ec1ced561072bd82f711695d`.
Official screening uses softmax class 1 and descending ordering; this number is
not a calibrated probability of experimental activity or a binding affinity.

Checkpoint compatibility correction: the pinned screening checkpoint contains
`mu` and `dev`, which the current upstream model does not define or use. Upstream
filters all unknown checkpoint entries. Our adapter allows only those two unused
extra keys after verifying the checkpoint hash, then loads all model parameters
with `strict=True`. Other unexpected keys, missing parameters and shape mismatches
remain errors. This fixes the workstation-reported startup exception, but a
successful live pilot is still required.

After pulling this correction, preserve the failed run and select a fresh output
parent because the worker code seal changed. For the existing workstation run:

```bash
export EQUISCORE_OUTPUT=/mnt/local/hand/yuzhang/aidd/e097-chat-workspace/users/workstation/projects/prj-adf8a1b9f4f8-prompt-aidd/runs/PROMPT-1d4f82d50d294925/execution/equiscore-v2
bash scripts/run_e103_equiscore.sh pilot
```

Keep the same `EQUISCORE_OUTPUT` for later full/analyze commands. No environment
reinstallation or docking rerun is required for this checkpoint-key correction.

The DGL wheel is listed at https://data.dgl.ai/wheels/repo.html. ProLIF 1.1 API:
https://prolif.readthedocs.io/en/v1.1.0/source/modules/interaction-fingerprint.html.
Protocol: `experiments/E103-equiscore-rescoring.md`.
Local fixtures cover pilot/full gating, ranking reversal, molecule deduplication,
unknown/duplicate/missing predictions, resume, score bounds and sealed pose inputs.
No real EquiScore predictions are part of this repository delivery.
