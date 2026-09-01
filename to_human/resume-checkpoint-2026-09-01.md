# AIDD Macrocycle Agent — Resume Checkpoint

Date: 2026-09-01  
Repository: `https://github.com/YUZhang-utu/Group-Agent.git`  
Checkpoint commit: `2859048`  
Validated offline suite: `49 passed`

## Project boundary

Continue only the Group-Agent AIDD macrocycle project. Do not mix it with the
separate KRAS necessity/enhancement project in the workspace root.

Desktop repository:

```text
D:\agent\projects\aidd_macrocycle_agent
```

University workstation repository:

```text
/mnt/medchem_taltio/wrk/yu_agent/Group-Agent
```

Scientific runtime/storage:

```text
/mnt/medchem_taltio/wrk/yu_agent/runtime
```

AlphaFold 3 is an independent external backend under:

```text
/mnt/medchem_taltio/wrk/yu_agent/alphafold3
```

## Active production context

- User: `USR-667167B3141D`
- Project: `PRJ-600CDDACB90A`
- Task: `TSK-A0C50CA9CF19`
- Campaign: `CAM-C0B597CE7ADE`
- Target: human WEE1, UniProt `P30291`
- Campaign state: `structures_review`
- Registered RCSB structures: 25
- Macrocycle library: `LIB-AFA68EE6888C`
- Library size: 100,000 molecules / 299,999 conformers

## Completed milestones

- Campaign/Project/Task ownership and provenance infrastructure.
- WEE1 RCSB search and conservative ligand-component filtering. GOL, CL, NA,
  EDO, PO4, MG and similar crystallization components are not candidate ligands.
- Multi-PDB and unified experimental/predicted receptor review in PyMOL.
- Audited receptor exclusions. 9TG7 remains Campaign provenance but is excluded
  from the WEE1 kinase receptor ensemble because WEE1 is only a degron peptide.
- AlphaFold 3 deployed via Apptainer on RTX 5090; WEE1 kinase prediction and
  production Campaign registration succeeded.
- Provider-neutral LLM receptor-eligibility interface. Recommendations are
  structured and mandatory-human-review; no direct Campaign mutation.
- Campaign co-crystal ligand registry with PDB/CCD/chain/residue/altloc,
  standardized chemistry, optional source checksum, and optional USRCAT.
- Pairwise Morgan/Tanimoto and available-USRCAT comparison.
- Provider-neutral LLM query-ligand recommendation interface.
- Immutable, rationale-required human query-ligand lock.
- Hierarchical Morgan 2D pool plus best-conformer USRCAT 3D reranking, explicit
  2D-only fallback, and immutable search provenance.

## Important current limitation

The 100,000-molecule MOL2 registry preserves molecule/conformer provenance but
does not yet contain standardized SMILES or production Morgan/USRCAT artifacts.
Do not claim that a real library similarity search has run. The current 2D/3D
search orchestration is implemented and tested, but production chemistry
preprocessing and workstation validation remain.

Campaign ligands also require automatic extraction from the downloaded mmCIF
files with authoritative RCSB CCD bond orders. Do not manually guess ligand
bonds or fabricate missing 3D conformers.

## Exact next development milestone

Implement automatic ligand preparation and production index building:

1. Read each eligible Campaign mmCIF and enumerate only retained ligand
   instances with model, chain, residue, insertion code, and altloc.
2. Fetch/cache authoritative RCSB CCD chemical definitions and bond orders.
3. Produce standardized parent SMILES plus crystal-coordinate SDF artifacts,
   with hashes and explicit standardization version.
4. Compute real USRCAT descriptors only from valid 3D conformers.
5. Register the WEE1 ligand manifest and generate the pairwise comparison.
6. Optionally create an LLM recommendation request; human reviews and locks one
   query ligand.
7. Parse the registered macrocycle MOL2 conformers with RDKit, build the real
   Morgan molecule index and USRCAT conformer index, plus conformer-to-molecule
   mapping and checksummed manifests.
8. Execute and validate the first `search-similar-ligands` run.

## Workstation resume commands

```bash
cd /mnt/medchem_taltio/wrk/yu_agent/Group-Agent
git pull --ff-only origin main
python -m pip install -e ".[workstation]"
aidd-agent init --db "$AIDD_DB"
python -m pytest -q
```

Expected test result at this checkpoint:

```text
49 passed
```

## Continuation prompt

Use this prompt next time:

```text
Continue only the Group-Agent AIDD macrocycle project from
to_human/resume-checkpoint-2026-09-01.md and research-state.yaml. Start from
commit 2859048. Implement the next milestone: automatic Campaign mmCIF/RCSB CCD
ligand extraction and standardized ligand artifacts, followed by production
Morgan/USRCAT index construction for LIB-AFA68EE6888C. Preserve the LLM as an
advisory interface, require human query-ligand locking, do not fabricate bond
orders or conformers, and do not mix in the KRAS necessity/enhancement project.
```
