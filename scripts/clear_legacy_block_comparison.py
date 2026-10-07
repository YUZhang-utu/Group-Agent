"""Offline, narrowly scoped removal of one twelve-arm legacy comparison."""
import argparse
import hashlib
import itertools
import json
import os
from pathlib import Path
import re
import shutil
import sqlite3
import time


def inside(path, root):
    path, root = Path(path).resolve(), Path(root).resolve()
    if path == root or root not in path.parents:
        raise ValueError('Path is outside its expected root: '+str(path))
    return path


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def inventory(storage, anchor):
    storage = Path(storage).resolve()
    dbpath = storage/'chat/conversations.sqlite3'
    with sqlite3.connect(dbpath.as_uri()+'?mode=ro',uri=True) as db:
        db.row_factory = sqlite3.Row
        jobs = [dict(r) for r in db.execute('SELECT * FROM jobs')]
    first = next((j for j in jobs if j['id']==anchor),None)
    if first is None:
        raise ValueError('Anchor task does not exist')
    first_plan = inside(first['plan'],storage)
    project = first_plan.parent.parent.parent
    if first_plan.name!='plan.json' or first_plan.parent.parent.name!='runs':
        raise ValueError('Unexpected task layout')

    def load(job):
        plan = inside(job['plan'],project/'runs')
        if plan.parent.parent != project/'runs' or plan.name!='plan.json':
            raise ValueError('Unexpected run path')
        steps = read(plan)['plan']['steps']
        if len(steps)!=1 or steps[0]['action']!='block_search_dock':
            return None
        params=steps[0]['params'];rid=params['request_id']
        if not re.fullmatch('[a-f0-9]+',rid):
            raise ValueError('Invalid request ID')
        request_path=inside(project/'block-campaigns'/(rid+'.json'),project/'block-campaigns')
        if hashlib.sha256(request_path.read_bytes()).hexdigest()!=params['request_sha256']:
            raise ValueError('Request seal changed')
        request=read(request_path)
        return plan,request_path,request

    loaded=load(first)
    if loaded is None:
        raise ValueError('Anchor is not a block campaign')
    base=loaded[2]
    def eligible(r):
        return (r.get('search_engine','legacy-budget-v1')=='legacy-budget-v1'
                and r.get('search_policy')=='ranked_blocks_until_budget' and r.get('dock') is False)
    if not eligible(base):
        raise ValueError('Only the legacy search-only comparison can be removed')
    keys=('session','analysis','analysis_sha256','query','query_sha256','conformers')
    rows=[]
    for job in jobs:
        if job['session']!=first['session'] or not job.get('plan'):
            continue
        if not job['request'].startswith('Search leading blocks for conformers '):
            continue
        item=load(job)
        if item is None:
            continue
        plan,request_path,r=item
        if not eligible(r) or any(r.get(k)!=base.get(k) for k in keys):
            continue
        if r['selection']['receptor']!=base['selection']['receptor']:
            continue
        rows.append(dict(job=job,run=str(plan.parent),request=str(request_path),
                         arm=[r['selection'][k] for k in ('scheme','method','blocks_per_method')]))
    expected=set(itertools.product(('E094','E095','E096'),('chemplp','equiscore'),(5,10)))
    if len(rows)!=12 or {tuple(r['arm']) for r in rows}!=expected:
        raise ValueError('Expected exactly one complete twelve-arm cohort; refusing ambiguous removal')
    selected={r['job']['id'] for r in rows}
    runs=[Path(r['run']) for r in rows]
    if len(set(runs))!=12:
        raise ValueError('Shared run directories cannot be deleted')
    for key in ('analysis','query'):
        source=Path(base[key]).resolve()
        if any(source==p or p in source.parents for p in runs):
            raise ValueError('Shared scientific input would be removed')
    for job in jobs:
        if job['id'] in selected or not job.get('plan'):
            continue
        plan=Path(job['plan'])
        if plan.exists() and any(p.name in plan.read_text(encoding='utf-8') for p in runs):
            raise ValueError('Another task references a selected run; inspect dependencies')
    for run in runs:
        inside(run,project/'runs')
        if not re.fullmatch('PROMPT-[a-f0-9]+',run.name):
            raise ValueError('Unexpected run directory name')
        if any(p.is_symlink() for p in run.rglob('*')):
            raise ValueError('Symlink in run directory; inspect manually')
    return dict(storage=str(storage),project=str(project),anchor=anchor,
                preserved_inputs={k:base[k] for k in ('analysis','query')},tasks=rows)


def require_stopped():
    if os.name!='posix' or not Path('/proc').is_dir():
        raise ValueError('Apply is supported only on the Linux workstation')
    for entry in Path('/proc').iterdir():
        if not entry.name.isdigit() or int(entry.name)==os.getpid():
            continue
        try:
            command=(entry/'cmdline').read_bytes().replace(b'\0',b' ').decode(errors='replace')
        except (FileNotFoundError,ProcessLookupError,PermissionError):
            continue
        if any(s in command for s in ('aidd_agent.chat_agent','aidd_agent.chat_web','aidd_agent.prompt_workflow','multiprocessing.spawn')):
            raise ValueError('Stop Chat and its computation workers first; active PID '+entry.name)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--storage-root',type=Path,required=True)
    parser.add_argument('--anchor-task',required=True)
    parser.add_argument('--apply',action='store_true',help='Delete exactly the validated legacy cohort after stopping Chat')
    args=parser.parse_args()
    report=inventory(args.storage_root,args.anchor_task)
    print(json.dumps(report,indent=2))
    if not args.apply:
        print('Preview only. No files or task records changed.')
        return
    require_stopped()
    storage=Path(report['storage'])
    audit=storage/'maintenance'/('removed-legacy-'+args.anchor_task+'-'+str(time.time_ns()))
    audit.mkdir(parents=True)
    (audit/'inventory.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    # Record cancellation before removing outputs so accidental later startup cannot launch them.
    dbpath=storage/'chat/conversations.sqlite3'
    ids=[r['job']['id'] for r in report['tasks']]
    with sqlite3.connect(dbpath) as db:
        db.executemany("UPDATE jobs SET status='cancelled',cancel=1,error='Legacy comparison removal requested' WHERE id=?",[(i,) for i in ids])
    for row in report['tasks']:
        shutil.rmtree(inside(row['run'],Path(report['project'])/'runs'))
        inside(row['request'],Path(report['project'])/'block-campaigns').unlink()
    with sqlite3.connect(dbpath) as db:
        db.executemany('DELETE FROM jobs WHERE id=?',[(i,) for i in ids])
    print('Removed twelve legacy tasks and their run directories. Audit: '+str(audit))
    print('Scoring inputs, adopted query, completed union campaign and conversation messages retained.')


if __name__=='__main__':
    main()
