"""Replay real mmCIF complexes offline or through an installed PyMOL runtime.

Usage: python scripts/check_pymol_interactions.py --structure 6Q9L.cif
       --ligand HTZ --chain A --resi 201 --output check-6Q9L [--offline]
"""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path

from aidd_agent.pymol_chemistry import component_metadata, apply_components
from aidd_agent.pymol_interactions import detect, display, COLORS


def offline(path, ligand_spec, metadata):
    import gemmi
    atoms = {}; groups = {}; ligand = set(); protein = set(); edges = []
    structure = gemmi.read_structure(str(path))
    for chain in structure[0]:
        for res in chain:
            names = {}
            for a in res:
                if a.altloc not in ('\x00', '', 'A') or a.occ <= 0:continue
                i = len(atoms) + 1
                atoms[i] = dict(name=a.name, elem=a.element.name, charge=a.charge,
                    chain=chain.name, resi=str(res.seqid), resn=res.name, segi='',
                    coord=[a.pos.x,a.pos.y,a.pos.z], acceptor=False)
                names[a.name] = i
                if dict(chain=chain.name,resi=str(res.seqid),resn=res.name)==ligand_spec:ligand.add(i)
                elif res.het_flag=='A':protein.add(i)
            groups[(chain.name,str(res.seqid),res.name)] = names
    for (_, _, comp), names in groups.items():
        for a,b in metadata.get('components',{}).get(comp,{}).get('bonds',[]):
            if a in names and b in names:edges.append((names[a],names[b]))
    if not ligand or not protein:raise ValueError('Ligand or protein selection is empty')
    import numpy as np
    ligand_xyz=np.array([atoms[i]['coord'] for i in ligand])
    near={tuple(atoms[i][k] for k in ('chain','resi','resn')) for i in protein
          if np.linalg.norm(ligand_xyz-np.array(atoms[i]['coord']),axis=1).min()<=7}
    protein={i for i in protein if tuple(atoms[i][k] for k in ('chain','resi','resn')) in near}
    aromatic, typing=apply_components(atoms,edges,metadata,set())
    raw,reps=detect(atoms,edges,ligand,protein,aromatic,[])
    return dict(mode='offline_real_coordinates', typing=typing,
        accepted_counts=dict(Counter(r['kind'] for r in raw)),
        displayed_counts=dict(Counter(r['kind'] for r in reps)), representatives=reps,
        limitations=['Not a PyMOL execution or rendering test.',
                    'Polar and halogen acceptor typing not evaluated offline.',
                    'Deposited charges only; protonation not inferred.'])


def live(path, ligand_spec, metadata, output):
    import pymol2
    with pymol2.PyMOL() as pm:
        cmd=pm.cmd
        cmd.load(str(path),'v001')
        catalog=[dict(id='v001',path=str(path),ligand=ligand_spec,chemical_components=metadata)]
        reports=[]
        for run in range(2):
            report=display(cmd,dict(id=f'replay-{run}',view=dict(operation='typed_interactions')),catalog,output)
            verified={}
            for kind,count in report['displayed_counts']['v001'].items():
                name='v001_ix_'+kind
                if name not in cmd.get_names('objects') or cmd.get_type(name)!='object:measurement':
                    raise AssertionError('Missing measurement object: '+name)
                actual=cmd.get('dash_color',name)
                # PyMOL versions may return either the color name or its index.
                if str(actual) not in {COLORS[kind],str(cmd.get_color_index(COLORS[kind]))}:
                    raise AssertionError(f'Wrong color for {name}: {actual}')
                verified[kind]=dict(count=count,color=actual)
            reports.append(dict(report=report,verified_objects=verified))
        if reports[0]['verified_objects']!=reports[1]['verified_objects']:
            raise AssertionError('Repeated invocation changed interaction objects')
        return dict(mode='real_pymol',runs=reports)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--structure',type=Path,required=True)
    parser.add_argument('--ligand',required=True)
    parser.add_argument('--chain',required=True)
    parser.add_argument('--resi',required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--offline',action='store_true')
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    path=args.structure.resolve();metadata=component_metadata(path)
    ligand_spec=dict(resn=args.ligand,chain=args.chain,resi=args.resi)
    result=offline(path,ligand_spec,metadata) if args.offline else live(path,ligand_spec,metadata,args.output)
    result.update(structure=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest())
    (args.output/'report.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k not in {'representatives','runs'}},indent=2))


if __name__=='__main__':main()
