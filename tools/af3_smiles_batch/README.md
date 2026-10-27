# Generic AlphaFold 3 SMILES batch

Each input table row produces a separate protein-ligand complex prediction.
The example molecules are syntax examples, not recommended binders.
Use your previously verified protein JSON as the base: monomer, homodimer and
heterodimer assemblies are preserved, including inline MSA/templates. No protein
sequence is inferred or invented. A protein dimer plus one ligand is the default;
use --ligand-copies 2 only when two identical ligand copies are intended.

Alternatively, copy protein.template.json to myprotein.json and replace the
sequence placeholder. For a homodimer change protein.id from "A" to ["A","B"].
For a heterodimer use two separate protein entries with different sequences and
IDs "A" and "B". The unedited template is deliberately rejected by the generator.

## Prepare and run

```bash
conda activate /mnt/medchem_taltio/wrk/yu_agent/Group_Agent/envs/aidd_workstation
python prepare_batch.py --base /path/to/myprotein.json \
  --ligands ligands.tsv --seeds 1 2 3 4 5 --output myprotein_ligands
python run_batch.py myprotein_ligands --check
python run_batch.py myprotein_ligands --run
```

Run on the school GPU allocation with Apptainer available. Host Python needs
RDKit for preparation only. The AF3 installation stays inside the existing SIF.
The runner uses the previous school paths, configurable with AF3_INSTALL,
AF3_IMAGE, AF3_MODELS and AF3_DB_DIR. --check verifies host paths only, not AF3
schema compatibility, available GPU memory or successful inference.

The table has `id` and `smiles` columns, separated by tabs (.tsv) or commas (.csv).
SMILES are copied literally after strict RDKit parsing. Include desired charge
and stereochemistry explicitly. Disconnected salts are rejected; supply the
intended ligand and add ions as explicit fixed entities if appropriate.
IDs must use letters/digits/underscore/hyphen and be unique ignoring case.

Base JSON must use dialect alphafold3. Its version is preserved, not upgraded:
the installed container must support it. The usual simple protein/SMILES inputs
work with version 1. External *Path fields must be converted to inline content
before use; relative files are not silently lost during input relocation.
Existing ligands are rejected unless --keep-base-ligands is explicit. That flag
keeps fixed cofactors in every job; do not accidentally retain a previous test
ligand. Covalent attachment of the newly added SMILES ligand is not implemented.

Output inputs/ contains one JSON per ligand; manifest.json records hashes,
literal SMILES, seeds and added chain IDs. The batch runner executes sequentially
to avoid launching many independent AF3 jobs on the same GPU. It stops at the
first error and writes a separate timestamped run report. It does not auto-resume:
rerunning starts all jobs again, using fresh output directories.

Prediction output follows the existing runner convention:
AF3_INSTALL/tests/output/LIGAND_ID/TIMESTAMP_PID/. Each job's final *_model.cif
and confidence JSON are referenced by the batch report. Seed/sample outputs
remain in the AF3 job directory. A successful process without a final model is
reported as failed. Confidence scores are not binding affinities.

If using an ordinary sequence-only base, each job may repeat MSA/template search.
For many ligands, reuse an appropriate prior AF3 *_data.json containing inline
MSA/templates, retaining only the intended proteins and fixed cofactors. Do not
replace missing MSA with empty strings as a speed shortcut. This helper does not
promise automatic cross-job feature caching.

The result is a predicted complex, not an experimental co-crystal or proof that
the ligand binds. Compare poses and interface confidence across seeds; AF3 is
not guaranteed to enforce a desired pocket. This toolkit was locally validated
for input creation and execution handling; school GPU inference is not run here.

Official format: https://github.com/google-deepmind/alphafold3/blob/main/docs/input.md
