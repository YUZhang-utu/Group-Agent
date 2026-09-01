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
