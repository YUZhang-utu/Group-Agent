# Campaign ligand comparison and similarity search

This workflow keeps three identities separate: the public CCD component, one
ligand instance in one crystal structure, and the standardized chemical query.
Only components retained by the PDB ligand classifier can be registered.

## 1. Register co-crystal ligand instances

Create `ligands.json` from reviewed structures. Each crystal instance has its
own chain/residue identity. `usrcat` is optional and, when present, must contain
60 numeric values from the actual 3D conformer.

```json
{
  "ligands": [
    {
      "pdb_id": "8BJU",
      "ccd_id": "LIG",
      "chain_id": "A",
      "residue_number": "901",
      "altloc": "",
      "standardized_smiles": "REPLACE_WITH_VALID_STANDARDIZED_SMILES",
      "source_path": "/absolute/path/to/extracted_ligand.sdf",
      "metadata": {"standardization": "neutral-parent-v1"}
    }
  ]
}
```

```bash
aidd-agent init --db "$AIDD_DB"
aidd-agent register-campaign-ligands --db "$AIDD_DB" --user USR-... \
  --campaign CAM-... --manifest ligands.json
aidd-agent campaign-ligands --db "$AIDD_DB" --user USR-... --campaign CAM-...
```

## 2. Compare and optionally request an LLM recommendation

The comparison reports Morgan/Tanimoto for every pair and USRCAT only when both
actual 3D descriptors exist.

```bash
aidd-agent compare-campaign-ligands --db "$AIDD_DB" --user USR-... \
  --campaign CAM-... > ligand-comparison.json
aidd-agent create-ligand-query-review --db "$AIDD_DB" --user USR-... \
  --campaign CAM-... --requirements-json query-requirements.json \
  --comparison-json ligand-comparison.json
```

An LLM recommendation remains advisory. After review, lock exactly one query:

```bash
aidd-agent select-query-ligand --db "$AIDD_DB" --user USR-... \
  --campaign CAM-... --ligand LIG-... \
  --rationale "Reviewed representative ATP-site chemotype"
```

## 3. Search the macrocycle library

The Morgan index must use molecule IDs. The USRCAT index may use conformer IDs;
in that case provide a JSON object mapping each conformer ID to its molecule ID.

```bash
aidd-agent search-similar-ligands --db "$AIDD_DB" --user USR-... \
  --campaign CAM-... --library LIB-... \
  --morgan-index /path/to/morgan.json.gz \
  --usrcat-index /path/to/usrcat.npz \
  --conformer-map-json /path/to/conformer-to-molecule.json \
  --two-d-pool 1000 --limit 100
```

The search record stores the locked query snapshot, index SHA-256 values,
parameters, best conformer, 2D and 3D scores, completed stages, and final rank.
It does not submit docking automatically.
# Production chemistry preparation

The registry and hierarchical search do not infer chemistry from component
names. On a workstation with `.[workstation]`, prepare retained crystal ligands
from downloaded Campaign mmCIF files and authoritative RCSB CCD definitions:

```bash
aidd-agent prepare-campaign-ligands --db "$AIDD_DB" --user USR-... \
  --campaign CAM-... --structures-dir "$PROJECT_ROOT/inputs/structures" \
  --ccd-cache "$RUNTIME_ROOT/cache/rcsb-ccd" \
  --output-dir "$RUNTIME_ROOT/campaigns/CAM-.../ligands"
aidd-agent register-campaign-ligands --db "$AIDD_DB" --user USR-... \
  --campaign CAM-... \
  --manifest "$RUNTIME_ROOT/campaigns/CAM-.../ligands/campaign-ligands.json"
```

Preparation joins crystal coordinates to CCD atoms by atom name and accepts
connectivity only from CCD bond tables. It does not guess bonds or generate
missing conformers. Missing heavy atoms and invalid chemistry are recorded in
the manifest's `errors` array.

After reviewing and locking a query ligand, build the registered MOL2 library
indices:

```bash
aidd-agent build-library-indices --db "$AIDD_DB" \
  --library LIB-AFA68EE6888C \
  --output-dir "$RUNTIME_ROOT/libraries/LIB-AFA68EE6888C/indices/v1"
```

The output includes the Morgan molecule index, conformer-level USRCAT index,
conformer-to-molecule map, and a checksummed version manifest. A production run
must archive its console output and verify manifest counts before search.
