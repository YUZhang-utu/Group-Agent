"""Read-only diagnostics of saved pilot poses and original molecule records."""
import argparse
from collections import Counter, defaultdict
import json
import hashlib
from pathlib import Path
import sqlite3

import numpy as np

from .expanded_wee1 import fingerprint
from .screening_selection import check_hashes


def distribution(values):
    return dict(zip(('min','q25','median','q75','max'),
                    np.quantile(values,[0,.25,.5,.75,1]).tolist())) if values else None


def summarize(db, total):
    hits=defaultdict(set); masks=Counter(); scores=defaultdict(list); rows=[]
    for mid,template,payload in db.execute('SELECT mid,template,payload FROM poses'):
        row=json.loads(payload);hits[template].add(mid);rows.append((mid,row))
        masks[int(row['matched_anchor_count'])]+=1
        for key in ('composite_score','gaussian_same_pose','optional_score'):
            if key in row:scores[key].append(row[key])
    union=set().union(*hits.values()) if hits else set()
    return dict(sampled_molecules=total,matching_molecules=len(union),
        templates=[dict(template=t,hits=len(ids),exclusive_hits=len(ids-set().union(
            *(v for k,v in hits.items() if k!=t)))) for t,ids in sorted(hits.items())],
        retained_pose_contact_counts=dict(masks),retained_pose_scores={k:distribution(v) for k,v in scores.items()},
        grouping_status=dict(db.execute("SELECT CASE WHEN error='' THEN 'ok' ELSE error END,COUNT(*) FROM members GROUP BY error")),
        valid_scaffold_groups=db.execute("SELECT COUNT(DISTINCT cluster) FROM members WHERE error='' ").fetchone()[0],
        limitation='Distributions cover retained representatives only, not every evaluated or rejected pose')


def compare_topology(mol, chem):
    from .chemical_companion import _bond_order
    heavy=[a for a in mol.GetAtoms() if a.GetAtomicNum()>1]
    indices={a.GetIdx():i for i,a in enumerate(heavy)}
    bonds=sorted((min(indices[b.GetBeginAtomIdx()],indices[b.GetEndAtomIdx()]),
                  max(indices[b.GetBeginAtomIdx()],indices[b.GetEndAtomIdx()]),_bond_order(b))
                 for b in mol.GetBonds() if b.GetBeginAtomIdx() in indices and b.GetEndAtomIdx() in indices)
    saved=sorted((min(int(b['begin']),int(b['end'])),max(int(b['begin']),int(b['end'])),int(b['order'])) for b in chem.bonds)
    return dict(atomic_numbers_equal=[a.GetAtomicNum() for a in heavy]==chem.atomic_numbers.tolist(),
        charges_equal=[a.GetFormalCharge() for a in heavy]==chem.formal_charges.tolist(),
        aromatic_flags_equal=[a.GetIsAromatic() for a in heavy]==((chem.atom_flags&1)!=0).tolist(),
        bonds_equal=bonds==saved,
        aromatic_heteroatom_hydrogens=[dict(index=i,element=a.GetSymbol(),hydrogens=a.GetTotalNumHs(includeNeighbors=True))
            for i,a in enumerate(heavy) if a.GetIsAromatic() and a.GetAtomicNum()!=6],
        hydrogen_metadata_in_companion=False)


def recover(batch, db, limit):
    from .chemical_companion import ChemicalCompanionReader
    from .mol2 import iter_mol2_blocks, load_rdkit_mol2
    reader=ChemicalCompanionReader(batch/'chemical/catalog.json');wanted=defaultdict(list);results=[]
    with sqlite3.connect((batch/'registry.sqlite3').resolve().as_uri()+'?mode=ro',uri=True) as registry:
        rows=db.execute("SELECT m.mid,MIN(p.gid) FROM members m JOIN poses p ON p.mid=m.mid WHERE m.error<>'' GROUP BY m.mid ORDER BY m.mid LIMIT ?",(limit,)).fetchall()
        for mid,gid in rows:
            chem=reader.get(gid)
            if chem.molecule_id!=mid:raise ValueError('Chemical molecule identity mismatch')
            record=registry.execute('SELECT source_path,source_record_index,content_sha256 FROM conformer WHERE id=? AND molecule_id=?',(chem.conformer_id,mid)).fetchone()
            if record is None:raise ValueError('Missing source conformer')
            wanted[record[0]].append((record[1],record[2],chem))
    for source,entries in wanted.items():
        pending={i:(sha,chem) for i,sha,chem in entries}
        if not Path(source).is_file():
            results.extend(dict(molecule_id=c.molecule_id,status='source_missing',path=source) for _,_,c in entries);continue
        for record_index,raw_text in iter_mol2_blocks(Path(source)):
            if record_index not in pending:continue
            sha,chem=pending.pop(record_index)
            if hashlib.sha256(raw_text.encode('utf-8')).hexdigest()!=sha:raise ValueError('Original MOL2 content hash mismatch')
            mol,mode=load_rdkit_mol2(raw_text,chem.conformer_id)
            results.append(dict(molecule_id=chem.molecule_id,conformer_id=chem.conformer_id,
                source=source,record_index=record_index,content_sha256=sha,status='recovered',
                original_sanitization=mode,**compare_topology(mol,chem)))
            if not pending:break
        if pending:raise ValueError('Original MOL2 records missing')
    return results


def compare_runs(left, right):
    reports=[json.loads((p/'report.json').read_text()) for p in (left,right)]
    for r in reports:check_hashes(r['output_hashes'])
    if any(not r.get('sample_complete') for r in reports):raise ValueError('Complete runs required')
    samples=[json.loads((p/'sample.json').read_text()) for p in (left,right)]
    same_sample=samples[0]==samples[1]
    same_design=reports[0]['design']==reports[1]['design']
    def poses(path):
        with sqlite3.connect((path/'candidates.sqlite').resolve().as_uri()+'?mode=ro',uri=True) as db:
            return {(mid,t,mask,spatial):(score,gid,json.loads(payload))
                    for mid,t,mask,spatial,score,gid,payload in db.execute('SELECT * FROM poses')}
    a,b=poses(left),poses(right)
    exact=a==b
    return dict(same_sample=same_sample,same_design=same_design,exact_pose_payloads=exact,
        equivalence_passed=same_sample and same_design and exact,
        baseline_scan_seconds=reports[0]['scan_wall_seconds'],new_scan_seconds=reports[1]['scan_wall_seconds'],
        scan_speedup=reports[0]['scan_wall_seconds']/reports[1]['scan_wall_seconds'])


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--pilot',type=Path,required=True);p.add_argument('--batch',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--recover',type=int,default=5)
    p.add_argument('--compare',type=Path,help='Optional complete rerun for exact output equivalence')
    a=p.parse_args()
    if not 0<=a.recover<=20:p.error('--recover must be between 0 and 20')
    original=json.loads((a.pilot/'report.json').read_text());check_hashes(original['output_hashes'])
    if not original.get('sample_complete'):raise ValueError('Review requires a complete pilot')
    a.output.mkdir(parents=True,exist_ok=False)
    with sqlite3.connect((a.pilot/'candidates.sqlite').resolve().as_uri()+'?mode=ro',uri=True) as db:
        report=summarize(db,original['sampled_molecules'])
        report['source_recovery']=recover(a.batch,db,a.recover) if a.recover else []
    report['sources']=fingerprint([a.pilot/'report.json',a.pilot/'candidates.sqlite'])
    report['production_changed']=False
    if a.compare:report['rerun_comparison']=compare_runs(a.pilot,a.compare)
    (a.output/'report.json').write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
