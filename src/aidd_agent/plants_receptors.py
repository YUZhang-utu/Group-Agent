"""Prepare explicitly adopted experimental pocket representatives using SPORES."""
from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path
import re
import subprocess

import numpy as np
from scipy.spatial import cKDTree

from .block_sampling import check_outputs, exclusive_output
from .final_work_blocks import read, sha
from .joint_spatial_profiles import save
from .screening_selection import check_hashes
from .mol2 import iter_mol2_blocks


def site_geometry(points,padding=5.):
    xyz=np.asarray(points,dtype=float)
    if xyz.ndim!=2 or xyz.shape[1]!=3 or not len(xyz) or not np.isfinite(xyz).all(): raise ValueError('Invalid reference heavy-atom coordinates')
    center=xyz.mean(axis=0)
    return dict(center=center.tolist(),radius=float(np.ceil(np.linalg.norm(xyz-center,axis=1).max()+padding)),
                policy='Reference ligand heavy-atom geometric centroid, not whole-protein centroid or mass center',padding_angstrom=padding)


def verify_heavy_coordinates(input_atoms, mol2):
    records=list(iter_mol2_blocks(Path(mol2)))
    if len(records)!=1: raise ValueError('SPORES must produce one protein MOL2 record')
    output=[];inside=False
    for line in records[0][1].splitlines():
        if line.startswith('@<TRIPOS>'): inside=line.strip()=='@<TRIPOS>ATOM';continue
        if inside and line.strip():
            fields=line.split();element=fields[5].split('.')[0].upper()
            if element not in ('H','D'): output.append((element,[float(v) for v in fields[2:5]]))
    if Counter(e for e,_ in input_atoms)!=Counter(e for e,_ in output): raise ValueError('SPORES changed heavy-atom inventory')
    maximum=0.
    for element in {e for e,_ in input_atoms}:
        a=np.array([p for e,p in input_atoms if e==element]); b=np.array([p for e,p in output if e==element])
        if not np.isfinite(b).all(): raise ValueError('Nonfinite prepared coordinates')
        dist,idx=cKDTree(a).query(b)
        if len(set(idx.tolist()))!=len(a): raise ValueError('Non-bijective heavy atom correspondence')
        maximum=max(maximum,float(dist.max()))
    if maximum>.05: raise ValueError('SPORES moved heavy atoms; binding-site frame requires review')
    return dict(heavy_atoms=len(output),maximum_heavy_coordinate_deviation=maximum)


def extract(row, output, padding=5.):
    from Bio.PDB import MMCIFParser, PDBIO, Select
    from .chemistry_prep import enumerate_ligand_instances
    queries=row.get('queries',[])
    if not queries: raise ValueError('Adopted representative has no coordinate-backed ligand site')
    # Deterministic native ligand selection is recorded, never an LLM coordinate guess.
    query=sorted(queries,key=lambda q:q['query_id'])[0]
    pdb,ccd,chain,resnum=query['query_id'].split(':')
    if pdb.upper()!=row['pdb_id'].upper(): raise ValueError('Reference PDB mismatch')
    source=Path(row['structure_path']); output=Path(output);output.mkdir(parents=True,exist_ok=True)
    instances=[r for r in enumerate_ligand_instances(source,[ccd]) if r['model']=='1' and r['chain_id']==chain and r['residue_number']==resnum]
    if len(instances)!=1: raise ValueError('Reference ligand alternate locations/insertion are ambiguous')
    ligand=np.array([[a[k] for k in ('x','y','z')] for a in instances[0]['atoms'] if a['element'].upper() not in ('H','D')])
    geometry=site_geometry(ligand,padding)
    structure=MMCIFParser(QUIET=True).get_structure(pdb,str(source)); model=next(structure.get_models())
    if row['target_chain'] not in model: raise ValueError('Adopted target chain absent')
    target=model[row['target_chain']];selected=[];removed=Counter();pending=[]
    for protein_chain in model:
        for residue in protein_chain:
            points=[a.coord for a in residue.get_atoms() if a.element not in ('H','D')]
            distance=float(cKDTree(ligand).query(np.asarray(points))[0].min()) if points else float('inf')
            reference=protein_chain.id==chain and str(residue.id[1])==resnum and residue.resname==ccd
            if protein_chain.id==target.id and residue.id[0]==' ': continue
            removed[residue.resname]+=1
            if not reference and residue.resname not in ('HOH','DOD') and distance<5:
                pending.append(dict(chain=protein_chain.id,residue=str(residue.id),name=residue.resname,distance=distance))
    if pending:
        save(output/'removed-nearby-components.json',dict(reason='Removed under the explicit selected-chain-only preparation policy; inspect cofactors, ions and interfaces',components=pending))
    # Model/chain filtering precedes SPORES; input coordinates stay in the native frame.
    for residue in target:
        if residue.id[0]!=' ':continue
        if residue.is_disordered()==2: raise ValueError('Disordered residue identity requires review')
        for atom in residue:
            choices=atom.disordered_get_list() if atom.is_disordered() else [atom]
            best=sorted(choices,key=lambda a:(-(a.occupancy or 0),a.altloc))[0]
            if best.element not in ('H','D'): selected.append(best)
    if not selected: raise ValueError('Empty receptor chain')
    chosen={id(a) for a in selected}
    class ProteinOnly(Select):
        def accept_model(self,m): return int(m.id==model.id)
        def accept_chain(self,c): return int(c.id==target.id)
        def accept_residue(self,r): return int(r.id[0]==' ')
        def accept_atom(self,a): return int(id(a) in chosen)
    # Use the selected source model directly when possible; block ambiguous long-chain exports.
    if len(target.id)!=1: raise ValueError('Multi-character chain requires explicit PDB chain remapping')
    writer=PDBIO();writer.set_structure(structure);writer.save(str(output/'protein-input.pdb'),ProteinOnly())
    evidence=dict(pdb_id=pdb,target_chain=target.id,query_id=query['query_id'],geometry=geometry,
                  source=str(source),source_sha256=sha(source),removed_components=dict(removed),
                  alternate_location_policy='Highest occupancy per atom, lexical alternate ID tie break',
                  nearby_removed_components=pending,
                  preparation_policy='User-selected target polymer chain only; remove waters, bound ligand and other components; explicitly report nearby removed components; no missing-heavy-atom rebuilding')
    save(output/'site.json',evidence)
    return evidence,[(a.element.upper(),a.coord.tolist()) for a in selected]


@exclusive_output
def prepare(adoption_report, profile, output):
    adopted=read(adoption_report)
    if adopted.get('kind')!='pocket_adoption' or adopted.get('status')!='complete': raise ValueError('Explicitly adopted experimental receptors required')
    check_hashes(adopted['sources']); pockets=read(adopted['pocket_report']); check_hashes(pockets['sources'])
    cfg=read(profile);allowed={'spores_executable','spores_mode','plants_executable','padding_angstrom','ligand_chemistry_reviewed','workers','timeout_seconds'}
    if set(cfg)-allowed: raise ValueError('Unknown receptor tool setting')
    mode=cfg.get('spores_mode','complete')
    if mode not in ('completepdb','reprotpdb','complete'): raise ValueError('Use a documented installed SPORES PDB mode')
    def path(key):
        p=Path(cfg[key]).expanduser(); return (p if p.is_absolute() else Path(profile).resolve().parent/p).resolve()
    spores=path('spores_executable'); plants=path('plants_executable')
    if not spores.is_file(): raise ValueError('SPORES executable unavailable')
    padding=cfg.get('padding_angstrom',5.)
    if type(padding) not in (int,float) or not 0<padding<=30: raise ValueError('Invalid site padding')
    timeout=cfg.get('timeout_seconds',600)
    if type(timeout) is not int or not 1<=timeout<=86400: raise ValueError('Invalid preparation timeout')
    output=Path(output).resolve()
    signature=dict(adoption=str(Path(adoption_report).resolve()),adoption_sha256=sha(Path(adoption_report)),profile=cfg,
                   spores_sha256=sha(spores),implementation_sha256=sha(Path(__file__)))
    if output.exists():
        if not (output/'signature.json').exists() or read(output/'signature.json')!=signature: raise ValueError('Use a fresh receptor output')
        if (output/'report.json').exists():return check_outputs(output)
    else: output.mkdir(parents=True)
    save(output/'signature.json',signature);receptors=[];geometries=[]
    groups={g['id']:g for g in pockets['clusters']+pockets.get('receptor_options',[])}
    rows={r['id']:r for r in pockets['structures']};seen=set()
    for group_id in adopted['selected_cluster_ids']:
        group=groups[group_id];row=rows[group['representative']]
        if row['id'] in seen:continue
        seen.add(row['id']);rid=re.sub('[^A-Za-z0-9_-]','_',row['id']);directory=output/rid
        evidence,atoms=extract(row,directory,padding)
        command=[str(spores),'--mode',mode,'protein-input.pdb','protein.mol2']
        if (directory/'prepared.json').exists():
            receipt=read(directory/'prepared.json')
            if receipt['mol2_sha256']!=sha(directory/'protein.mol2'):raise ValueError('Prepared receptor changed')
        else:
            with (directory/'spores.log').open('w',encoding='utf-8') as log:
                result=subprocess.run(command,cwd=directory,stdout=log,stderr=subprocess.STDOUT,shell=False,timeout=timeout)
            if result.returncode or not (directory/'protein.mol2').is_file(): raise ValueError('SPORES preparation failed; inspect receptor log')
            receipt=dict(mol2_sha256=sha(directory/'protein.mol2'),coordinates=verify_heavy_coordinates(atoms,directory/'protein.mol2'),command=command)
            save(directory/'prepared.json',receipt)
        geometry=evidence['geometry'];geometries.append(evidence)
        receptors.append(dict(id=rid,mol2=str(directory/'protein.mol2'),center=geometry['center'],radius=geometry['radius'],reviewed=True,
                              evidence='Explicit pocket adoption '+str(Path(adoption_report).resolve())+'; source/site and preparation checks in '+str(directory/'site.json')))
    if not receptors: raise ValueError('No adopted representatives')
    plants_profile=dict(executable=str(plants),receptors=receptors,ligand_mode='rigid',search_speed='speed1',cluster_structures=1,
                        workers=cfg.get('workers',4),timeout_seconds=7200,ligand_chemistry_reviewed=cfg.get('ligand_chemistry_reviewed',False))
    save(output/'plants-profile.json',plants_profile)
    report=dict(status='complete',kind='plants_receptors',receptors=len(receptors),sites=geometries,plants_profile=str(output/'plants-profile.json'),
                readiness='prepared_receptors_not_docked',adoption_report=str(Path(adoption_report).resolve()),
                output_hashes={str(p.relative_to(output)):sha(p) for p in output.rglob('*') if p.is_file() and p.name!='report.json'},
                limitations=['SPORES default protonation/typing, not pKa enumeration or missing-heavy-atom reconstruction.',
                             'Native ligand centroid/radius is an initial box; large macrocycles may need a reviewed larger site.',
                             'Experimental receptor selection was explicitly adopted; docking and affinity remain unvalidated.'])
    save(output/'report.json',report);return report


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--adoption',required=True);p.add_argument('--profile',required=True);p.add_argument('--output',required=True)
    a=p.parse_args();print(json.dumps(prepare(a.adoption,a.profile,a.output),indent=2))


if __name__=='__main__':main()
