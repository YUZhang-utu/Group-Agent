from __future__ import annotations

import argparse
import json
import multiprocessing as mp
from pathlib import Path
import time

from .mol2 import iter_mol2_records


_FACTORY = None


def _init_worker():
    global _FACTORY
    from rdkit import RDConfig
    from rdkit.Chem import ChemicalFeatures
    _FACTORY = ChemicalFeatures.BuildFeatureFactory(
        str(Path(RDConfig.RDDataDir) / "BaseFeatures.fdef"))


def _process_molecule(task):
    from rdkit import Chem
    from rdkit.Chem import Descriptors3D, rdMolDescriptors
    molecule_name, records = task
    template = None; conformers = features = heavy_atoms = 0
    for raw in records:
        mol = Chem.MolFromMol2Block(raw, sanitize=True, removeHs=False)
        if mol is None: return molecule_name, 0, 0, 0, False
        if template is None:
            template = tuple((f.GetFamily(), tuple(f.GetAtomIds()))
                             for f in _FACTORY.GetFeaturesForMol(mol))
        # Exercise the per-conformer production computations.
        rdMolDescriptors.GetUSRCAT(mol)
        Descriptors3D.PMI1(mol); Descriptors3D.PMI2(mol); Descriptors3D.PMI3(mol)
        conformers += 1; features += len(template)
        heavy_atoms += mol.GetNumHeavyAtoms()
    return molecule_name, conformers, features, heavy_atoms, True


def load_groups(source: Path, conformer_limit: int):
    groups=[]; current_name=None; current=[]; total=0
    for record in iter_mol2_records(source):
        if current_name is not None and record.molecule_name != current_name:
            groups.append((current_name, tuple(current))); current=[]
        current_name=record.molecule_name; current.append(record.raw_text); total += 1
        if total >= conformer_limit: break
    if current: groups.append((current_name, tuple(current)))
    return groups


def benchmark(source: Path, conformers: int, worker_counts: list[int]) -> dict:
    groups=load_groups(source, conformers); reports=[]
    for workers in worker_counts:
        started=time.perf_counter()
        if workers == 1:
            _init_worker(); results=map(_process_molecule, groups)
            result=list(results)
        else:
            context=mp.get_context("spawn")
            with context.Pool(workers, initializer=_init_worker) as pool:
                result=list(pool.imap(_process_molecule, groups, chunksize=16))
        seconds=time.perf_counter()-started
        completed=sum(x[1] for x in result); failures=sum(not x[4] for x in result)
        reports.append({"workers":workers,"molecules":len(groups),"conformers":completed,
                        "failures":failures,"seconds":seconds,
                        "conformers_per_second":completed/seconds,
                        "speedup_vs_one":None})
    baseline=reports[0]["conformers_per_second"]
    for row in reports: row["speedup_vs_one"]=row["conformers_per_second"]/baseline
    return {"format":"aidd-parallel-ingest-benchmark","version":1,
            "source":str(source.resolve()),"requested_conformers":conformers,
            "results":reports}


def main():
    p=argparse.ArgumentParser();p.add_argument("--source",type=Path,required=True)
    p.add_argument("--conformers",type=int,default=3000);p.add_argument("--workers",default="1,4,8,16")
    p.add_argument("--output",type=Path,required=True);a=p.parse_args()
    report=benchmark(a.source,a.conformers,[int(x) for x in a.workers.split(",")])
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(report,indent=2),encoding="utf-8")
    print(json.dumps(report,indent=2));return 0

if __name__=="__main__":raise SystemExit(main())
