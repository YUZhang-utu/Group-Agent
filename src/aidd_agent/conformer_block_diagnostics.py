"""Read-only attribution of capacity splitting and boundary hard-group fragmentation."""
import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
import sqlite3


def readonly(path):
    return sqlite3.connect(Path(path).resolve().as_uri()+'?mode=ro', uri=True)


def bucket(n):
    for bound in (100, 500, 2000, 10000, 20000):
        if n < bound:
            return f'less_than_{bound}'
    return 'at_least_20000'


def run(model, descriptors, output):
    model=Path(model).resolve();descriptors=Path(descriptors).resolve();output=Path(output).resolve()
    if output==model or model in output.parents or output==descriptors.parent:
        raise ValueError('Use a separate fresh diagnostic output directory')
    report_path=model/'report.json'
    receipt=json.loads(report_path.read_text(encoding='utf-8'))
    if receipt.get('status')!='complete':
        raise ValueError('Completed block model required')
    groups=[];hist=Counter();summary_by_boundary={}
    with readonly(model/'blocks.sqlite') as db, readonly(descriptors) as source:
        roots=list(db.execute("SELECT node,group_id,n FROM tree WHERE path='' ORDER BY group_id"))
        leaves=list(db.execute('SELECT group_id,n FROM tree WHERE axis IS NULL'))
        leaf_counts=Counter(g for g,n in leaves)
        leaf_totals=Counter()
        for g,n in leaves:
            leaf_totals[g]+=n;hist[n]+=1
        for number,(node,group,n) in enumerate(roots,1):
            row=db.execute('SELECT cid FROM point WHERE group_id=? LIMIT 1',(group,)).fetchone()
            if row is None: raise ValueError('Root without a source conformer')
            payload=source.execute('SELECT payload FROM descriptor WHERE cid=?',row).fetchone()
            if payload is None: raise ValueError('Group representative missing from descriptor database')
            descriptor=json.loads(payload[0]);states=descriptor['provenance']['omega_states']
            if descriptor['hard_group']!=group or set(states)-{'cis','trans','boundary'}:
                raise ValueError('Group identity or omega states mismatch')
            if leaf_totals[group]!=n: raise ValueError('Leaf/root populations disagree')
            boundary=states.count('boundary');category='with_boundary' if boundary else 'without_boundary'
            entry=summary_by_boundary.setdefault(category,dict(hard_groups=0,conformers=0,blocks=0,capacity_extra_blocks=0))
            for key,value in dict(hard_groups=1,conformers=n,blocks=leaf_counts[group],capacity_extra_blocks=leaf_counts[group]-1).items():
                entry[key]+=value
            groups.append(dict(hard_group=group,conformers=n,blocks=leaf_counts[group],
                capacity_extra_blocks=leaf_counts[group]-1,boundary_positions=boundary,
                cis_positions=states.count('cis'),trans_positions=states.count('trans'),
                omega_pattern=';'.join(states),representative_conformer=row[0]))
            if number%1000==0: print(f'Inspected {number}/{len(roots)} hard groups',flush=True)
    total=sum(n for g,n in leaves)
    if total!=receipt['counts']['conformers'] or len(leaves)!=receipt['counts']['blocks']:
        raise ValueError('Tree counts disagree with completed model report')
    if set(leaf_counts)!={g for node,g,n in roots}: raise ValueError('Root/leaf group mismatch')
    bins={}
    for n,count in sorted(hist.items()):
        entry=bins.setdefault(bucket(n),dict(blocks=0,conformers=0))
        entry['blocks']+=count;entry['conformers']+=n*count
    result=dict(status='complete',scope='Read-only metadata diagnostic; no regrouping or raw-coordinate recomputation',
        model=str(model),descriptors=str(descriptors),model_id=receipt['model_id'],capacity=receipt['capacity'],
        conformers=total,blocks=len(leaves),hard_groups=len(roots),
        blocks_without_capacity_limit=len(roots),capacity_extra_blocks=len(leaves)-len(roots),
        mean_block_size=total/len(leaves) if leaves else 0,
        mean_group_size_without_capacity_limit=total/len(roots) if roots else 0,
        block_size_bins=bins,boundary_attribution=summary_by_boundary,
        maximum_uncapped_group_size=max((n for node,g,n in roots),default=0),
        representative_descriptors_inspected=len(groups),
        provenance=dict(model_report_sha256=hashlib.sha256(report_path.read_bytes()).hexdigest(),
                        implementation_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),
        limitations=['Requires the existing successful full integrity validation; large input hashes are not rechecked here',
            'One descriptor per hard group supplies its omega pattern; population counts come from the fitted tree',
            'Boundary-associated groups are not necessarily caused solely by boundary: size and chirality also stratify',
            'No boundary-to-cis/trans reassignment; ambiguous angular states must remain traceable',
            'Removing capacity bounds changes logical group size, not chemical identity; task batching can remain separate'])
    output.mkdir(parents=True,exist_ok=False)
    for name,records in [('hard_groups.csv',groups),('block_size_histogram.csv',
                          [dict(size=n,blocks=c) for n,c in sorted(hist.items())])]:
        with (output/name).open('w',newline='',encoding='utf-8') as f:
            writer=csv.DictWriter(f,fieldnames=list(records[0]) if records else ['empty'])
            writer.writeheader();writer.writerows(records)
    result['output_hashes']={name:hashlib.sha256((output/name).read_bytes()).hexdigest()
                             for name in ('hard_groups.csv','block_size_histogram.csv')}
    (output/'report.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--model',type=Path,required=True)
    p.add_argument('--descriptors',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();print(json.dumps(run(a.model,a.descriptors,a.output),indent=2))


if __name__=='__main__':main()
