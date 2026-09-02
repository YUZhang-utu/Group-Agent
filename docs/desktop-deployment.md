# Desktop deployment and validation

The agent source can run on Windows, macOS, or Linux. The long-lived desktop is
the data authority: its Project folders, registry, molecular libraries, and
result metadata are not stored in Git. Laptops run the same client code and
reach that machine through an authenticated service or SSH bridge in a later
deployment stage.

## 1. Install Miniforge and create the workstation environment

From an Anaconda/Miniforge prompt in the repository root:

```powershell
conda env create -f environments/aidd-workstation.yml
conda activate aidd-workstation
python -m pip install --no-deps -e .
aidd-agent doctor
pytest -q
```

For the Linux chemistry workstation, update the same named environment with
the Reduce-enabled specification:

```bash
conda env update -n aidd-workstation -f environments/aidd-workstation-linux.yml --prune
conda activate aidd-workstation
python -m pip install --no-deps -e .
cp configs/tools/reduce.example.json configs/tools/reduce.local.json
reduce -help
```

Machine-specific HET dictionary paths belong in the ignored
`configs/tools/reduce.local.json`. Windows does not require Reduce; it can still
compute heavy-atom distance and SASA evidence.

Prepare and run the query-only 5-angstrom pocket on Linux:

```bash
aidd-agent export-reduce-pocket --mmcif /path/to/8BJU.cif --output query-pocket-5A.pdb --ccd-id QT9 --chain A --residue 601
aidd-agent run-reduce-pocket --profile configs/tools/reduce.local.json --pocket-pdb query-pocket-5A.pdb --output-dir reduce-QT9
aidd-agent validate-reduce-pocket --output-dir reduce-QT9 --ccd-id QT9 --chain A --residue 601
```

The output directory preserves the hydrogenated PDB, Reduce stderr/flip report,
input/output hashes, exact argument array, and scope. These data annotate query
anchors and never filter library candidates.

`het_dictionary: null` in the run manifest means Reduce used its default lookup;
it is acceptable only when validation confirms that QT9 and protein hydrogen
counts increased and no dictionary/connectivity warning was emitted.

PyMOL is included in the Conda environment. If its package is unavailable for a
particular platform, remove `pymol-open-source` from the environment file,
create the environment, and install a licensed/local PyMOL build separately.

## 2. Validate a real RCSB structure

```powershell
aidd-agent fetch-pdb --pdb 4YHJ --project-root D:\AIDD\users\alice\projects\PRJ-ID
```

The command writes `inputs/structures/4YHJ.cif` and a metadata JSON containing
the source URL, byte count, and SHA-256 digest.

## 3. Configure AlphaFold

AlphaFold is deliberately not installed into the general workstation
environment. AlphaFold 3 and its databases require Linux, substantial storage,
GPU resources, and separately obtained model parameters. Clone the official
repository on the desktop or cluster, follow its installation documentation,
request the model parameters under the applicable terms, and copy
`configs/models/alphafold3.example.json` to a private path outside Git. Update
all paths and set `license_acknowledged` to `true` only after reviewing those
terms.

Generate an input without running the model:

```powershell
aidd-agent prepare-alphafold3 --name wee1 --sequence MSEQUENCE... `
  --project-root D:\AIDD\users\alice\projects\PRJ-ID `
  --output D:\AIDD\users\alice\projects\PRJ-ID\runs\af3\input.json `
  --profile D:\AIDD\private\alphafold3.json
```

The returned command is an argument array suitable for a local executor or a
validated Slurm Bundle. It does not contain credentials. AlphaFold 2 is retained
as a fallback profile for sites that already maintain that installation.

The same reviewed sequence can be prepared for any supported prediction
backend through `prepare-structure-prediction`. Supported backends are
AlphaFold 2, AlphaFold 3, Boltz-2, and Chai-1; each stays in a separate
environment configured by a private profile copied from `configs/models/`.
The Agent emits a command array but does not accept license terms, store model
weights in Git, or run GPU work on a login node.

```powershell
aidd-agent fetch-alphafold-db --uniprot P30291 `
  --project-root D:\AIDD\users\alice\projects\PRJ-ID

aidd-agent prepare-structure-prediction --backend boltz2 --name wee1 `
  --sequence MSEQUENCE... --project-root D:\AIDD\users\alice\projects\PRJ-ID `
  --output D:\AIDD\users\alice\projects\PRJ-ID\runs\boltz2\input.yaml `
  --profile D:\AIDD\private\boltz2.json
```

## 4. Data and source-control boundary

Never commit molecular libraries, PDB/mmCIF files, SQLite registries, prediction
outputs, model parameters, SSH keys, API keys, or `.env` files. Git contains
only code, environment definitions, templates, tests, and documentation.
