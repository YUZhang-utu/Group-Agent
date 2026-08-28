# University storage and HPC quick start

This procedure recreates the tested workstation environment on a Linux login
node and keeps all scientific data outside Git. Replace `/GROUP/AIDD` with the
allocated university path.

## Install and validate

```bash
cd /GROUP/AIDD
mkdir -p software runtime
git clone https://github.com/YUZhang-utu/Group-Agent.git software/Group-Agent
cd software/Group-Agent

module load miniconda  # use the site's actual Conda/Miniforge module name
conda env create -f environments/aidd-workstation.yml
conda activate aidd-workstation
aidd-agent doctor
pytest -q
```

If the site has no Conda module, install Miniforge under the user's home or the
group software directory, subject to site policy. A Windows virtual environment
must not be copied to Linux; the YAML recreates equivalent packages for Linux.

## Initialize a test deployment

```bash
export AIDD_ROOT=/GROUP/AIDD/runtime
mkdir -p "$AIDD_ROOT"
aidd-agent init-storage --storage-root "$AIDD_ROOT"
aidd-agent init --db "$AIDD_ROOT/registry/aidd.sqlite3"
```

Create the first identity and Project, copying the returned IDs into subsequent
commands:

```bash
aidd-agent create-user --db "$AIDD_ROOT/registry/aidd.sqlite3" \
  --username test-user --display-name "Test User"

aidd-agent create-project --db "$AIDD_ROOT/registry/aidd.sqlite3" \
  --user USR-ID --name "Target Structure Pilot" \
  --objective "Validate target structure selection on university storage" \
  --storage-root "$AIDD_ROOT"

aidd-agent activate-project --db "$AIDD_ROOT/registry/aidd.sqlite3" \
  --user USR-ID --project PRJ-ID
```

## Exercise the target-structure workflow

```bash
aidd-agent create-campaign --db "$AIDD_ROOT/registry/aidd.sqlite3" \
  --user USR-ID --project PRJ-ID --name structure-pilot \
  --objective "Review experimental target structures"

aidd-agent set-target --db "$AIDD_ROOT/registry/aidd.sqlite3" \
  --user USR-ID --campaign CAM-ID --name "Target name" \
  --organism "Homo sapiens" --uniprot P12345

aidd-agent search-pdb --db "$AIDD_ROOT/registry/aidd.sqlite3" \
  --user USR-ID --campaign CAM-ID --uniprot P12345 --max-resolution 3.0

aidd-agent fetch-pdb --pdb 4YHJ \
  --project-root "$AIDD_ROOT/users/test-user/projects/PRJ-ID"
```

Use a real UniProt accession, PDB ID, Project path, and scientific rationale.
Do not run docking, prediction, or MD on a login node. Generate a Run/Slurm
Bundle and submit compute work through the scheduler after the smoke test.
