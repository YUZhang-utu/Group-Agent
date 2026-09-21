# Workstation docking installation

The operator reported the Schrodinger installation at
`/site/app9/schrodinger/2025-1`. PLANTS is provisionally configured as
`~/PLANTS1.2/PLANTS1.2_64bit`; the workstation probe must confirm the filename.
Tilde expansion uses the account running the probe.

```bash
git pull --ff-only
bash scripts/setup_prompt_profiles.sh
bash scripts/probe_docking_tools.sh
```

Setup creates `/mnt/local/hand/yuzhang/aidd/config/docking.local.json` without
overwriting existing profiles. Edit this file if the binary has another location.
Set `AIDD_DOCKING_PROFILE` to use a different profile file. This profile is separate
from the prompt runtime; do not add docking fields to that runtime yet.

The probe does not execute binaries or consume licenses. `executable` means an
existing file with execution permission; it does not prove compatible libraries,
license availability or successful docking. An explicitly configured missing path
does not fall back to another installation. All output remains English.

## Remaining integration and validation

Docking execution and chat docking actions are not implemented by this change.
The existing chat supports protein, AF3 and calibrated WEE1 search workflows.

For Glide, provide a prepared receptor/grid ZIP and prepared ligand input, then
validate licensed execution and a fixed reference-ligand redocking protocol.
For PLANTS, provide prepared receptor and ligand MOL2 files and a verified binding
site definition. Record preparation, engine configuration, input hashes, timing,
pose exports and reference-pose comparisons separately from 3D retrieval metrics.
Do not infer these inputs from the installation path or accept docking quality
because a program exits successfully. E031 remains annotation-only.
