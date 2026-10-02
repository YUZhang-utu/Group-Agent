# AIDD E097 continuity, 2026-10-02

Resume in `D:/agent/projects/aidd_structure_guided_chat`, branch
`feature/structure-guided-chat`, remote `YUZhang-utu/Group-Agent`.

User confirmed all three frozen regular partitions: E094 384, E095 1188, E096 1189.
Sample 100 conformers per block, 276,100 slots before CID deduplication. Same conformer
across partitions shares one docking result per receptor/settings; never deduplicate
by molecule. Special 132,885 and unsupported chemistry remain separate. No repartition.

E096 metadata and four matching ionizable controls passed on the workstation. Full
scientific recall remains unvalidated. E097 operates on source-backed memberships and
raw MOL2 records; it does not re-extract the 317-dimensional library features.

User additionally requires automated receptor preparation: give a PDB ID, identify
the target, assess single/multiple experimental receptor advice using existing pocket
analysis, ask human to accept representatives, then clean those chains and call SPORES.
User supplied `~/pdb/1/SPORES_64bit` and `--mode complete input.pdb output.mol2`.
The example `6Y4Q` is not an adopted MDM2 reference. Actual reference still unspecified.
Delete other chains, ligand, water and crystallization components; list nearby removed
non-water components. Native bound ligand heavy-atom centroid defines the initial site.
Verify SPORES preserves heavy-atom inventory and native coordinates before using it.

Implemented:

- `block_sampling.py`: per-block seeded reservoirs, shared source-verified MOL2 chunks,
  samples SQLite and CSV, source-file export resume, writer locks.
- `block_plants.py`: receptor-specific shared jobs, installed executable/config seals,
  rigid ligand default, first-job gate, failure/timeout handling and retry, strict
  score/pose identity mapping, block summaries and population weighting.
- `plants_receptors.py`: adopted representative extraction, SPORES, site geometry,
  output coordinate checks and generated `plants-profile.json`.
- `protein_from_pdb` and `receptor_assess` planner actions; existing pocket advice
  reused. Ambiguous target entities or native ligand instances need explicit resolution.
- Chat `block_evaluation` stages: sample, sample_prepare, adopt_receptors, receptors,
  prepare, run, dock_analyze, evaluate, analyze. Model intent, not fixed word matching.
  The same contract is available to the existing domain agent.
- Complete typed sample/prepare/dock/analyze chains and sealed prior-task follow-ups.
  Receptor selection/preparation is a separate explicitly adopted owned source.

Local verification: 85 related tests passed; fake SPORES/PLANTS only. Real binary
output conventions, full-library timing, live LLM routing quality and target docking
quality are not validated. N-E/binding-affinity models are not executed; analysis CSV
contains CID/receptor/pose-file/record-index keys for that later integration.

Workstation base:
`/mnt/local/hand/yuzhang/aidd/mc_runs/e081-001/audit-stage/blocks-20000`.
See [E097 protocol/setup](E097_BLOCK_PLANTS_EVALUATION.md). Configure paths once with
`scripts/configure_e097_blocks.py`; preserve the user's existing runtime using its
`--runtime` option or AIDD_RUNTIME_PROFILE. PLANTS is discovered from existing environment,
PATH/home or `--plants`. The generator supplies the known SPORES path/mode.

Use Chat sampling if subsequent conversational continuation is required. Independent
batch sampling is also available via `scripts/run_e097_block_samples.sh`, but its output
is not silently imported into the owned task database. Do not rerun an existing long
job; resume its stage using matching configuration/code. Do not claim a prepared job
has actually docked. No credentials, raw library, proprietary executables or sampled
MOL2 files should be uploaded to GitHub.
