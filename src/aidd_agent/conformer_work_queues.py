"""Uncapped identity classes and bounded evaluation queues from an E086 audit."""
import argparse
from collections import defaultdict
import csv
import hashlib
import json
from pathlib import Path
import sqlite3


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def key(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True).encode()).hexdigest()[:24]


def table(path, rows, fields):
    with path.open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader();writer.writerows(rows)


def plan(diagnostic, output, rare_threshold=100, budgets=(100,200)):
    diagnostic=Path(diagnostic).resolve();output=Path(output).resolve()
    if rare_threshold<1 or not budgets or any(n<1 for n in budgets) or len(set(budgets))!=len(budgets):
        raise ValueError('Positive threshold and unique positive budgets required')
    report=json.loads((diagnostic/'report.json').read_text())
    model=Path(report['model']).resolve()
    if any(output==p or p in output.parents for p in (model,diagnostic)):
        raise ValueError('Use a fresh output outside model and diagnostic trees')
    if report['status']!='complete' or sha(model/'report.json')!=report['provenance']['model_report_sha256']:
        raise ValueError('Diagnostic model evidence changed')
    for name,digest in report['output_hashes'].items():
        if sha(diagnostic/name)!=digest: raise ValueError('Diagnostic artifact changed: '+name)
    with (diagnostic/'hard_groups.csv').open(newline='',encoding='utf-8') as f:
        original=list(csv.DictReader(f))
    classes=[];queues=defaultdict(list);seen=set()
    for row in original:
        group=row['hard_group'];n=int(row['conformers']);states=row['omega_pattern'].split(';')
        if group in seen or n<1 or set(states)-{'cis','trans','boundary'}:
            raise ValueError('Invalid class identity, population or omega pattern')
        seen.add(group);boundary=states.count('boundary')
        if boundary!=int(row['boundary_positions']): raise ValueError('Boundary count mismatch')
        if boundary:
            # Execution pool only: no inferred chirality equivalence or cis/trans assignment.
            kind='boundary_review_pool';signature=[kind,len(states)]
        elif n<rare_threshold:
            # Definite cis/trans patterns stay separated even in execution queues.
            kind='rare_definite_pool';signature=[kind,len(states),states]
        else:
            kind='definite_class';signature=[kind,group]
        queue='queue-'+key(signature)
        item=dict(hard_group=group,logical_class=group,queue_id=queue,queue_kind=kind,
                  conformers=n,old_blocks=int(row['blocks']),unit_count=len(states),
                  omega_pattern=row['omega_pattern'],boundary_positions=boundary)
        classes.append(item);queues[queue].append(item)
    if len(classes)!=report['hard_groups'] or sum(r['conformers'] for r in classes)!=report['conformers']:
        raise ValueError('Diagnostic population mismatch')
    by_group={r['hard_group']:r for r in classes};mapping=[];totals=defaultdict(int);counts=defaultdict(int)
    # Validate every old leaf population against the diagnostic, without rewriting its points.
    with sqlite3.connect((model/'blocks.sqlite').as_uri()+'?mode=ro',uri=True) as db:
        for block,group,n in db.execute('SELECT node,group_id,n FROM tree WHERE axis IS NULL ORDER BY node'):
            if group not in by_group or n<1: raise ValueError('Unknown or invalid old block')
            r=by_group[group];totals[group]+=n;counts[group]+=1
            mapping.append(dict(old_block=block,logical_class=group,queue_id=r['queue_id'],conformers=n))
    if len(mapping)!=report['blocks'] or any(totals[g]!=r['conformers'] or counts[g]!=r['old_blocks'] for g,r in by_group.items()):
        raise ValueError('Existing block populations changed')
    queue_rows=[];allocations=[]
    for queue,members in sorted(queues.items()):
        members=sorted(members,key=lambda r:r['hard_group']);kind=members[0]['queue_kind']
        population=sum(r['conformers'] for r in members)
        queue_rows.append(dict(queue_id=queue,kind=kind,logical_classes=len(members),conformers=population,
                               unit_count=members[0]['unit_count']))
        for budget in sorted(budgets):
            # Nested deterministic balanced class allocation, never an activity ranking.
            ordered=sorted(members,key=lambda r:key([queue,r['hard_group'],'allocation-v1']))
            quota={r['hard_group']:0 for r in ordered};remaining=min(budget,population)
            while remaining:
                for r in ordered:
                    if quota[r['hard_group']]<r['conformers']:
                        quota[r['hard_group']]+=1;remaining-=1
                    if not remaining: break
            allocations.extend(dict(queue_id=queue,budget=budget,logical_class=r['hard_group'],
                planned_conformer_slots=quota[r['hard_group']],
                coverage='planned' if quota[r['hard_group']] else 'unassessed_do_not_reject') for r in members)
    output.mkdir(parents=True,exist_ok=False)
    for name,rows in [('classes.csv',classes),('queues.csv',queue_rows),('old_block_mapping.csv',mapping),
                      ('sampling_plan.csv',allocations)]:
        table(output/name,rows,list(rows[0]) if rows else ['empty'])
    kinds={kind:dict(queues=sum(r['kind']==kind for r in queue_rows),
                     conformers=sum(r['conformers'] for r in queue_rows if r['kind']==kind),
                     logical_classes=sum(r['logical_classes'] for r in queue_rows if r['kind']==kind))
           for kind in sorted({r['kind'] for r in queue_rows})}
    summary=dict(status='complete',readiness='offline_queue_plan_not_screening_dispatch',
        source_model=str(model),source_diagnostic=str(diagnostic),model_id=report['model_id'],
        source_conformers=report['conformers'],original_blocks=report['blocks'],
        logical_classes=len(classes),evaluation_queues=len(queues),queue_kinds=kinds,
        logical_capacity_limit=None,rare_threshold=rare_threshold,
        sampling=[dict(budget_per_queue=b,planned_conformer_slots=sum(r['planned_conformer_slots'] for r in allocations if r['budget']==b),
                       classes_without_slots=sum(r['planned_conformer_slots']==0 for r in allocations if r['budget']==b)) for b in sorted(budgets)],
        integrity=dict(old_leaf_mapping_complete=True,class_populations_preserved=True,
                       original_cis_trans_classes_unchanged=True,original_databases_modified=False),
        sources={str(diagnostic/'report.json'):sha(diagnostic/'report.json'),
                 str(diagnostic/'hard_groups.csv'):sha(diagnostic/'hard_groups.csv'),
                 str(model/'report.json'):sha(model/'report.json')},
        output_hashes={name:sha(output/name) for name in ('classes.csv','queues.csv','old_block_mapping.csv','sampling_plan.csv')},
        implementation_sha256=sha(__file__),
        limitations=['Execution pools are not merged chemical or conformational states',
            'Sampling plan is conformer slots, not exported ligands or unique-molecule counts',
            'No score, geometry, docking, affinity or recall calculation; unassessed classes cannot be rejected',
            'Boundary angles and chirality are not reassigned; all original hard-group identities are retained',
            'Existing successful full integrity validation required; large SQLite hashes are not rechecked',
            'Old frozen model remains valid; this sidecar does not alter existing search or docking dispatch'])
    (output/'report.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    return summary


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--diagnostic',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--rare-threshold',type=int,default=100);p.add_argument('--budgets',type=int,nargs='+',default=[100,200])
    a=p.parse_args();print(json.dumps(plan(a.diagnostic,a.output,a.rare_threshold,tuple(a.budgets)),indent=2))


if __name__=='__main__': main()
