"""Export immutable rank pages from original hash-verified MOL2 source records."""
import argparse
from collections import defaultdict
import csv
import hashlib
import json
from pathlib import Path
import re
import sqlite3

import numpy as np

from .budget_screen import readonly
from .expanded_wee1 import fingerprint
from .gaussian_batch import _atomic_json
from .library_acceptance import read
from .screening_selection import check_hashes
from .prompt_workflow import file_lock


def transformed_mol2(raw, transform):
    matrix=np.asarray(transform,dtype=float).reshape(4,4)
    if not np.isfinite(matrix).all():raise ValueError('Invalid pose transform')
    lines=[];atoms=False
    for line in raw.splitlines(keepends=True):
        if line.startswith('@<TRIPOS>'):atoms=line.strip()=='@<TRIPOS>ATOM'
        elif atoms and line.strip():
            fields=line.split()
            if len(fields)<6:raise ValueError('Malformed MOL2 atom row')
            xyz=np.array(list(map(float,fields[2:5])))@matrix[:3,:3].T+matrix[:3,3]
            fields[2:5]=[f'{v:.6f}' for v in xyz];line=' '.join(fields)+'\n'
        lines.append(line)
    return ''.join(lines)


def hydrogen_roundtrip(mol):
    """Heavy-graph reconstruction retains attached hydrogen counts and charges.

    Explicit atom hydrogens are collapsed into heavy-atom H counters only in this
    audit. Original MOL2 output retains their coordinates and representation.
    """
    from rdkit import Chem
    heavy=[a for a in mol.GetAtoms() if a.GetAtomicNum()>1]
    metadata=[dict(atomic_number=a.GetAtomicNum(),formal_charge=a.GetFormalCharge(),
        explicit_h_counter=a.GetNumExplicitHs(),explicit_h_neighbors=sum(n.GetAtomicNum()==1 for n in a.GetNeighbors()),
        total_h=a.GetTotalNumHs(includeNeighbors=True)) for a in heavy]
    lookup={a.GetIdx():i for i,a in enumerate(heavy)};copy=Chem.RWMol()
    for atom,record in zip(heavy,metadata):
        new=Chem.Atom(atom);new.SetNumExplicitHs(record['total_h']);new.SetNoImplicit(True);copy.AddAtom(new)
    for b in mol.GetBonds():
        if b.GetBeginAtomIdx() in lookup and b.GetEndAtomIdx() in lookup:
            copy.AddBond(lookup[b.GetBeginAtomIdx()],lookup[b.GetEndAtomIdx()],b.GetBondType())
    reconstructed=copy.GetMol()
    try:
        Chem.SanitizeMol(reconstructed)
        same=all(a.GetFormalCharge()==r['formal_charge'] and a.GetTotalNumHs()==r['total_h']
                 and a.GetNumExplicitHs()==r['total_h'] for a,r in zip(reconstructed.GetAtoms(),metadata))
        status='ok' if same else 'charge_or_hydrogen_mismatch'
    except Exception as exc:status=type(exc).__name__
    return dict(status=status,atoms=metadata)


def export(batch,run,output,start=1,count=100000):
    from .chemical_companion import ChemicalCompanionReader
    from .mol2 import iter_mol2_blocks,load_rdkit_mol2
    from .pilot_review import compare_topology
    batch=Path(batch).resolve();run=Path(run).resolve();out=Path(output).resolve()
    if start<1 or count<1:raise ValueError('Positive one-based rank range required')
    if out.is_relative_to(batch) or batch.is_relative_to(out) or out==run or run.is_relative_to(out):
        raise ValueError('Export overlaps protected input')
    report=read(run/'report.json')
    if report['status']!='complete':raise ValueError('Complete ranking required')
    check_hashes(report['outputs'])
    out.mkdir(parents=True,exist_ok=True)
    with file_lock(out/'export.lock'):
        protocol=dict(ranking_sha256=report['outputs'][str(run/'ranking.sqlite')],start=start,count=count,batch=str(batch))
        if (out/'protocol.json').exists() and read(out/'protocol.json')!=protocol:raise ValueError('Changed export page')
        _atomic_json(out/'protocol.json',protocol)
        if (out/'report.json').exists():
            result=read(out/'report.json');check_hashes(result['outputs']);return result
        chemical=ChemicalCompanionReader(batch/'chemical/catalog.json');wanted=defaultdict(dict);rows=[]
        with readonly(run/'ranking.sqlite') as db,readonly(batch/'registry.sqlite3') as registry:
            selected=db.execute('SELECT rank,mid,rrf,template_support FROM ranking WHERE rank>=? AND rank<? ORDER BY rank',(start,start+count))
            for rank,mid,rrf,support in selected:
                if not re.fullmatch(r'[A-Za-z0-9_-]+',mid):raise ValueError('Unsafe molecule ID')
                poses=[json.loads(db.execute('SELECT p.payload FROM poses p JOIN template_ranks t ON p.mid=t.mid AND p.template=t.template WHERE p.mid=? ORDER BY t.rank,p.template LIMIT 1',(mid,)).fetchone()[0])]
                pose=poses[0];chem=chemical.get(pose['global_id'])
                if chem.molecule_id!=mid:raise ValueError('Chemical identity mismatch')
                record=registry.execute('SELECT c.source_path,c.source_record_index,c.content_sha256,c.source_record_name,m.source_name FROM conformer c JOIN molecule m ON c.molecule_id=m.id WHERE c.id=? AND c.molecule_id=?',(chem.conformer_id,mid)).fetchone()
                if record is None:raise ValueError('Source record missing from registry')
                source,index,sha,record_name,name=record
                entry=dict(rank=rank,molecule_id=mid,name=name,record_name=record_name,rrf_score=rrf,template_support=support,
                    conformer_id=chem.conformer_id,source_path=source,source_record_index=index,source_sha256=sha,
                    representative_policy='Best within-template rank; template-ID tie',poses=poses,
                    original_mol2=f'original/{mid}.mol2',posed_mol2=f'posed/{mid}.mol2')
                if index in wanted[source]:raise ValueError('Duplicate original source record')
                wanted[source][index]=(entry,chem);rows.append(entry)
        (out/'original').mkdir(exist_ok=True);(out/'posed').mkdir(exist_ok=True)
        for source,pending in wanted.items():
            for index,raw in iter_mol2_blocks(Path(source)):
                if index not in pending:continue
                entry,chem=pending.pop(index)
                if hashlib.sha256(raw.encode('utf-8')).hexdigest()!=entry['source_sha256']:raise ValueError('MOL2 source hash mismatch')
                mol,mode=load_rdkit_mol2(raw,chem.conformer_id);topology=compare_topology(mol,chem)
                identity=all(topology[k] for k in ('atomic_numbers_equal','charges_equal','aromatic_flags_equal','bonds_equal'))
                audit=hydrogen_roundtrip(mol)
                entry['chemical_qc']=dict(status='ok' if identity and audit['status']=='ok' and mode=='strict' else 'review_required',
                    original_sanitization=mode,companion_topology=topology,hydrogen_roundtrip=audit)
                if not identity:raise ValueError('Source/companion topology mismatch')
                for key,text in [('original_mol2',raw),('posed_mol2',transformed_mol2(raw,entry['poses'][0]['transform']))]:
                    target=out/entry[key];temp=target.with_suffix('.partial');temp.write_text(text,encoding='utf-8',newline='');temp.replace(target)
                entry['mol2_hashes']=fingerprint([out/entry['original_mol2'],out/entry['posed_mol2']])
                if not pending:break
            if pending:raise ValueError('Missing source MOL2 records')
            print(f'Exported source: {source}',flush=True)
        with readonly(run/'ranking.sqlite') as db,(out/'molecules.jsonl').open('w',encoding='utf-8') as stream:
            for row in rows:
                full=dict(row,poses=[json.loads(p[0]) for p in db.execute('SELECT p.payload FROM poses p JOIN template_ranks t ON p.mid=t.mid AND p.template=t.template WHERE p.mid=? ORDER BY t.rank,p.template',(row['molecule_id'],))])
                stream.write(json.dumps(full,ensure_ascii=False)+'\n')
        fields=['rank','molecule_id','name','record_name','rrf_score','template_support','conformer_id','original_mol2','posed_mol2','qc_status']
        with (out/'molecules.csv').open('w',encoding='utf-8-sig',newline='') as stream:
            writer=csv.DictWriter(stream,fieldnames=fields);writer.writeheader()
            for row in rows:writer.writerow({**{k:row[k] for k in fields if k!='qc_status'},'qc_status':row['chemical_qc']['status']})
        outputs=fingerprint([out/'molecules.csv',out/'molecules.jsonl'])
        for row in rows:outputs.update(row['mol2_hashes'])
        result=dict(status='complete',requested_molecules=count,exported_molecules=len(rows),start_rank=start,
            end_rank=start+len(rows)-1,shortfall=count-len(rows),review_required=sum(r['chemical_qc']['status']!='ok' for r in rows),
            outputs=outputs,scope='Original input and transformed representative; no protonation or docking preparation applied')
        _atomic_json(out/'report.json',result);return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('batch','run','output'):p.add_argument('--'+key,type=Path,required=True)
    p.add_argument('--start-rank',type=int,default=1);p.add_argument('--count',type=int,default=100000)
    a=p.parse_args();r=export(a.batch,a.run,a.output,a.start_rank,a.count)
    print(json.dumps({k:v for k,v in r.items() if k!='outputs'},indent=2))


if __name__=='__main__':main()
