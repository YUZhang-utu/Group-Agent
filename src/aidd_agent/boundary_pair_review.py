"""Source-verified boundary/definite pair panels with continuous property and steric metrics."""
import argparse
from collections import Counter,defaultdict
import csv
import hashlib
import itertools
import json
from pathlib import Path
import sqlite3

import numpy as np
from rdkit import Chem

from .macrocycle_descriptors import describe_peptide
from .macrocycle_property_profiles import profile,property_difference
from .mol2 import iter_mol2_blocks,parse_mol2_block,load_rdkit_mol2


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def readonly(path):return sqlite3.connect(Path(path).resolve().as_uri()+'?mode=ro',uri=True)


def band(value):
    return 'up_to_35' if value<=35 else '35_to_45' if value<=45 else 'above_45_control'


def select_panel(db,limit):
    # One minimum-hash CID per original class per band, then a bounded class panel.
    chosen={};counts=Counter()
    for cid,group,dev,encoded in db.execute('SELECT cid,original_group,maximum_deviation,candidates FROM candidate'):
        matches=json.loads(encoded)
        if not matches or dev is None:continue
        category=band(dev);counts[category]+=1;k=(category,group)
        priority=hashlib.sha256(('E089-v1:'+cid).encode()).hexdigest()
        row=dict(cid=cid,original_group=group,maximum_deviation=dev,band=category,matches=matches,priority=priority)
        if k not in chosen or priority<chosen[k]['priority']:chosen[k]=row
    result=[]
    for category in sorted(counts):
        result.extend(sorted((r for r in chosen.values() if r['band']==category),key=lambda r:r['priority'])[:limit])
    return result,counts


def rmsd(first,second):
    a=np.asarray(first,float);b=np.asarray(second,float)
    if a.shape!=b.shape or a.ndim!=2 or a.shape[1]!=3 or not np.isfinite(a).all() or not np.isfinite(b).all():
        raise ValueError('Incompatible backbone coordinates')
    a=a-a.mean(0);b=b-b.mean(0);u,_,vh=np.linalg.svd(a.T@b)
    if np.linalg.det(u@vh)<0:u[:,-1]*=-1
    return float(np.sqrt(np.mean(np.sum((a@(u@vh)-b)**2,axis=1))))


def aligned_profile(value,stored_units):
    observed=[tuple(u) for u in value['units']];wanted=[tuple(u) for u in stored_units]
    if len(observed)!=len(wanted) or len(set(observed))!=len(observed) or set(observed)!=set(wanted):
        raise ValueError('Property/stored unit atom mapping mismatch')
    order=[observed.index(u) for u in wanted]
    return {**value,**{k:[value[k][i] for i in order] for k in
                      ('units','omega_states','properties','typed_spatial_moments','sterics')}}


def prepare_record(record,row,variant):
    p=row['provenance']
    if (record.content_sha256,record.name,record.molecule_name)!=(p['content_sha256'],p['source_record_name'],p['source_name']):
        raise ValueError('Raw source record identity/hash mismatch')
    mol,mode=load_rdkit_mol2(record.raw_text,record.name)
    if mode!='strict':raise ValueError('Strict chemistry required')
    name=None if 'cross_conformer_alignment_unverified' in p['flags'] else record.molecule_name
    description=describe_peptide(mol,name,variant)
    if description['hard_group']!=row['hard_group'] or description['units']!=p['units'] or not np.allclose(description['descriptor'],row['descriptor'],rtol=1e-9,atol=1e-9):
        raise ValueError('Sample does not reproduce stored source descriptor')
    heavy=Chem.RemoveHs(Chem.Mol(mol));xyz=np.asarray(heavy.GetConformer().GetPositions())
    value=aligned_profile(profile(mol,name),p['units'])
    value['backbone_coordinates']=[xyz[u[:3]].tolist() for u in p['units']]
    return value


def compare(first,second,rotation):
    # Candidate rotation belongs to second (boundary), bringing it to first's order.
    score=property_difference(first,second,rotation)
    a=np.asarray(first['backbone_coordinates']);b=np.roll(np.asarray(second['backbone_coordinates']),-rotation,axis=0)
    score.pop('note');score['backbone_rmsd_angstrom']=rmsd(a.reshape(-1,3),b.reshape(-1,3))
    return score


def run(build,candidates,output,per_band=100,references=3):
    if per_band<1 or references<2:raise ValueError('Positive panel size and at least two references required')
    build=Path(build).resolve();candidates=Path(candidates).resolve();output=Path(output).resolve()
    if any(output==p or p in output.parents for p in (build,candidates)):raise ValueError('Fresh output outside inputs required')
    receipt=json.loads((build/'report.json').read_text());upstream=json.loads((candidates/'report.json').read_text())
    if receipt['status']!='complete' or upstream['status']!='complete':raise ValueError('Completed inputs required')
    if upstream['sources'].get(str(build/'report.json'))!=sha(build/'report.json'):raise ValueError('Candidate/build provenance mismatch')
    if sha(candidates/'candidates.sqlite')!=upstream['output_hashes']['candidates.sqlite']:raise ValueError('Candidate database changed')
    diagnostic=next(Path(p) for p in upstream['sources'] if p!=str(build/'report.json') and Path(p).name=='report.json')
    if sha(diagnostic)!=upstream['sources'][str(diagnostic)]:raise ValueError('Diagnostic changed')
    diag=json.loads(diagnostic.read_text());model=Path(diag['model'])
    if sha(model/'report.json')!=diag['provenance']['model_report_sha256']:raise ValueError('Model receipt changed')
    with readonly(candidates/'candidates.sqlite') as db:panel,populations=select_panel(db,per_band)
    if not panel:raise ValueError('No angular candidates available for pair review')
    print(f'Selected {len(panel)} boundary panel records',flush=True)
    targets=sorted({m['hard_group'] for r in panel for m in r['matches']});refs={};selected={}
    with readonly(model/'blocks.sqlite') as blocks,readonly(build/'descriptors.sqlite') as source:
        for number,group in enumerate(targets,1):
            seen=set();ids=[]
            for cid,mid in blocks.execute('SELECT cid,mid FROM point WHERE group_id=?',(group,)):
                if mid in seen:continue
                seen.add(mid);ids.append(cid)
                if len(ids)>=references:break
            if not ids:raise ValueError('Candidate target has no reference records')
            refs[group]=ids
            if number%25==0:print(f'Selected references for {number}/{len(targets)} target classes',flush=True)
        for cid in sorted({r['cid'] for r in panel}|{c for ids in refs.values() for c in ids}):
            row=source.execute('SELECT payload FROM descriptor WHERE cid=?',(cid,)).fetchone()
            if row is None:raise ValueError('Missing source descriptor')
            selected[cid]=json.loads(row[0])
    for r in panel:
        if selected[r['cid']]['hard_group']!=r['original_group']:raise ValueError('Candidate class mismatch')
    for group,ids in refs.items():
        if any(selected[c]['hard_group']!=group for c in ids):raise ValueError('Reference class mismatch')
    output.mkdir(parents=True,exist_ok=False)
    (output/'panel.json').write_text(json.dumps(dict(panel=panel,target_references=refs),indent=2),encoding='utf-8')
    wanted=defaultdict(dict)
    for cid,row in selected.items():
        p=row['provenance'];path=str(Path(p['source_path']).resolve());index=p['source_record_index']
        if index in wanted[path]:raise ValueError('Duplicate selected source identity')
        wanted[path][index]=cid
    profiles={}
    with (output/'profiles.jsonl').open('w',encoding='utf-8') as features,(output/'selected.mol2').open('w',encoding='utf-8') as raw:
        for number,(path,records) in enumerate(sorted(wanted.items()),1):
            print(f'Reading source file {number}/{len(wanted)}; {len(records)} selected records',flush=True)
            remaining=set(records)
            for index,text in iter_mol2_blocks(Path(path)):
                if index and index%50000==0:print(f'Located source record {index}; {len(remaining)} selections remain in this file',flush=True)
                if index not in records:continue
                cid=records[index];record=parse_mol2_block(Path(path),index,text)
                value=prepare_record(record,selected[cid],receipt['schema']['variant']);profiles[cid]=value
                features.write(json.dumps(dict(cid=cid,source=selected[cid]['provenance'],**value))+'\n')
                raw.write(record.raw_text);raw.write('\n')
                remaining.remove(index)
                if not remaining:break
            if remaining:raise ValueError('Selected source records missing')
    scores=[]
    for r in panel:
        for m in r['matches']:
            for cid in refs[m['hard_group']]:
                scores.append(dict(kind='boundary_vs_reference',band=r['band'],query=r['cid'],reference=cid,
                    target_group=m['hard_group'],rotation=m['rotation'],maximum_deviation=r['maximum_deviation'],
                    same_molecule=selected[r['cid']]['molecule_id']==selected[cid]['molecule_id'],
                    **compare(profiles[cid],profiles[r['cid']],m['rotation'])))
    for group,ids in refs.items():
        for a,b in itertools.combinations(ids,2):
            scores.append(dict(kind='within_definite_control',band='control',query=a,reference=b,target_group=group,
                rotation=0,maximum_deviation='',same_molecule=selected[a]['molecule_id']==selected[b]['molecule_id'],
                **compare(profiles[b],profiles[a],0)))
    with (output/'pair_metrics.csv').open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=list(scores[0]));writer.writeheader();writer.writerows(scores)
    summaries={}
    for category in sorted({r['band'] for r in scores}):
        current=[r for r in scores if r['band']==category]
        summaries[category]=dict(pairs=len(current),metrics={k:dict(median=float(np.median([r[k] for r in current])),
            p90=float(np.quantile([r[k] for r in current],.9))) for k in
            ('property_rms','steric_rms','maximum_local_steric_difference','backbone_rmsd_angstrom')})
    result=dict(status='complete',readiness='pair_review_not_membership_adoption',selected_conformers=len(selected),
        boundary_panel_counts=dict(Counter(r['band'] for r in panel)),eligible_population_by_band=dict(populations),
        target_classes=len(refs),targets_without_within_class_controls=sum(len(v)<2 for v in refs.values()),
        pair_summaries=summaries,sources={str(build/'report.json'):sha(build/'report.json'),
            str(candidates/'report.json'):sha(candidates/'report.json')},
        output_hashes={n:sha(output/n) for n in ('panel.json','profiles.jsonl','selected.mol2','pair_metrics.csv')},
        code_hashes={name:sha(Path(__file__).with_name(name)) for name in
                    ('boundary_pair_review.py','macrocycle_property_profiles.py','macrocycle_descriptors.py')},
        limitations=['Class-balanced deterministic review panel, not population-unbiased coverage or recall estimation',
            'References use first indexed distinct molecules; few references cannot characterize an entire class',
            'All candidate rotations retained; do not pick minimum scores without disclosing selection',
            'Raw selected records verified; unselected blocks are read as text only to locate selected records',
            'Aggregate pair summaries have unequal class weights; inspect per-class controls in CSV',
            'No final reassignment, property cutoff, PLANTS or N-E execution'])
    (output/'report.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('build','candidates','output'):p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--per-band',type=int,default=100);p.add_argument('--references',type=int,default=3)
    a=p.parse_args();print(json.dumps(run(a.build,a.candidates,a.output,a.per_band,a.references),indent=2))


if __name__=='__main__':main()
