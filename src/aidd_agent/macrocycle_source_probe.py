"""Stream a bounded MOL2 prefix to validate macrocycle extraction on real source data."""
import argparse
from collections import Counter
import json
import hashlib
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from .macrocycle_blocks import describe,digest,partition,VERSION
from .mol2 import iter_mol2_records,load_rdkit_mol2
from .chemical_companion import BOND_DTYPE,_bond_order


def probe(sources,output,per_file=100,capacity=10000):
    import itertools
    if not 1<=per_file<=10000:raise ValueError('per_file must be 1..10000')
    if capacity<1:raise ValueError('capacity must be positive')
    output=Path(output)
    if output.exists():raise ValueError('Use a new output directory')
    output.mkdir(parents=True);rows=[];failures=[];sources_meta=[]
    with (output/'descriptors.jsonl').open('w',encoding='utf-8') as stream:
        for path in map(Path,sources):
            sources_meta.append(dict(path=str(path.resolve()),bytes=path.stat().st_size,mtime_ns=path.stat().st_mtime_ns))
            for record in itertools.islice(iter_mol2_records(path),per_file):
                base=dict(global_id=len(rows)+len(failures),molecule_id='LOCAL-M-'+digest(record.molecule_name)[:20],
                    conformer_id='LOCAL-C-'+digest([str(path.resolve()),record.record_index,record.content_sha256])[:20],
                    source_name=record.molecule_name,source_record_name=record.name,conformer_index=record.conformer_index,
                    source_path=str(path.resolve()),source_record_index=record.record_index,content_sha256=record.content_sha256)
                try:
                    mol,mode=load_rdkit_mol2(record.raw_text,record.name)
                    heavy=[a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum()>1];index={a:i for i,a in enumerate(heavy)}
                    bonds=[]
                    for b in mol.GetBonds():
                        a,c=b.GetBeginAtomIdx(),b.GetEndAtomIdx()
                        if a in index and c in index:bonds.append((index[a],index[c],_bond_order(b),0,0))
                    chem=SimpleNamespace(atomic_numbers=np.array([mol.GetAtomWithIdx(i).GetAtomicNum() for i in heavy]),
                        formal_charges=np.array([mol.GetAtomWithIdx(i).GetFormalCharge() for i in heavy]),bonds=np.array(bonds,dtype=BOND_DTYPE))
                    xyz=np.array([list(mol.GetConformer().GetAtomPosition(i)) for i in heavy])
                    descriptor=describe(chem,xyz)
                    # Preserve original atom IDs and names, rather than assuming ID=index+1.
                    atom_section=record.raw_text.split('@<TRIPOS>ATOM',1)[1].split('@<TRIPOS>',1)[0]
                    atoms=[line.split() for line in atom_section.splitlines() if line.strip()]
                    if len(atoms)!=mol.GetNumAtoms():raise ValueError('source_atom_mapping_length_mismatch')
                    base.update(sanitization=mode,heavy_atom_source_ids=[int(atoms[i][0]) for i in heavy],
                        heavy_atom_names=[atoms[i][1] for i in heavy],substructure_labels=sorted({tuple(a[6:8]) for a in atoms}),
                        annotated_backbone_bonds=sum('BACKBONE' in line for line in record.raw_text.splitlines()))
                    row=dict(base,**descriptor);rows.append(row);stream.write(json.dumps(row)+'\n')
                except ValueError as exc:failures.append(dict(base,reason=str(exc)))
    memberships,blocks=partition(rows,capacity)
    for name,value in [('memberships',memberships),('blocks',blocks),('unassigned',failures)]:
        with (output/(name+'.jsonl')).open('w',encoding='utf-8') as stream:
            for row in value:stream.write(json.dumps(row)+'\n')
    report=dict(status='complete',version=VERSION,scope='First records per source file, not a uniform sample or full library',
        sources=sources_meta,records=len(rows)+len(failures),assigned=len(rows),unassigned=len(failures),blocks=len(blocks),capacity=capacity,
        ring_sizes=dict(Counter(r['ring_size'] for r in rows)),amide_counts=dict(Counter(r['amide_count'] for r in rows)),
        peptide_unit_counts=dict(Counter(len(r['peptide_units']) for r in rows)),
        assigned_molecule_names=len({r['source_name'] for r in rows}),
        reasons=dict(Counter(r['reason'] for r in failures)),
        identity_scope='LOCAL IDs for this source probe only; production registry IDs are neither inferred nor changed',
        limitations=['Name tokens are preserved but not assigned to building-block atoms',
                    'Only backbone torsion partitioning; no performance/enrichment or energy-barrier validation'])
    def sha(path):
        return hashlib.sha256(path.read_bytes()).hexdigest()
    report['output_hashes']={p.name:sha(p) for p in output.glob('*.jsonl')}
    report['implementation_hashes']={p.name:sha(p) for p in
        (Path(__file__),Path(__file__).with_name('macrocycle_blocks.py'))}
    (output/'report.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    return report


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--mol2',type=Path,nargs='+',required=True)
    parser.add_argument('--output',type=Path,required=True);parser.add_argument('--per-file',type=int,default=100)
    parser.add_argument('--capacity',type=int,default=10000);a=parser.parse_args()
    print(json.dumps(probe(a.mol2,a.output,a.per_file,a.capacity),indent=2))


if __name__=='__main__':main()
