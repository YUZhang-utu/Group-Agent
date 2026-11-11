"""Immutable project candidate collections and reusable MOL2 delivery."""
import json
from pathlib import Path
import re
import sqlite3
import time

from .expanded_wee1 import fingerprint
from .library_acceptance import read,sha
from .project_context import ensure_within
from .screening_selection import check_hashes
from .budget_screen import readonly
from .gaussian_batch import _atomic_json


def get(root,identifier,verify=True):
    if not isinstance(identifier,str) or not re.fullmatch('[a-f0-9]{24}',identifier):raise ValueError('Invalid candidate set ID')
    path=ensure_within(Path(root)/'candidate-sets'/identifier/'manifest.json',Path(root))
    if not path.is_file():raise ValueError('Candidate set is unavailable in this project')
    result=read(path)
    if result.get('set_id')!=identifier:raise ValueError('Candidate set identity changed')
    members=ensure_within(Path(result['members']),Path(root))
    if str(members) not in result['outputs']:raise ValueError('Candidate membership is not sealed')
    for output in result['outputs']:ensure_within(Path(output),Path(root))
    if verify:check_hashes(result['outputs'])
    return result


def handle(app,sid,args):
    from .research_tasks import queue
    ctx=app.context_for(sid);root=Path(ctx['project_root']);operation=args.get('operation')
    if operation=='list':
        return dict(sets=[{k:r.get(k) for k in ('set_id','name','count','operation','parents','created')} for p in sorted((root/'candidate-sets').glob('*/manifest.json')) for r in [read(p)]])
    if operation=='show':
        item=get(root,args.get('set_id'));offset=args.get('offset',0);limit=args.get('limit',20)
        if type(offset) is not int or offset<0 or type(limit) is not int or not 1<=limit<=100:raise ValueError('Invalid page')
        with readonly(item['members']) as db:
            rows=[json.loads(r[0]) for r in db.execute('SELECT payload FROM members ORDER BY rank LIMIT ? OFFSET ?',(limit,offset))]
        return dict(set={k:v for k,v in item.items() if k!='outputs'},members=rows,offset=offset)
    if operation=='label':
        item=get(root,args.get('set_id'));ids=args.get('molecule_ids')
        if not isinstance(ids,list) or not 1<=len(ids)<=1000 or any(not isinstance(v,str) for v in ids):raise ValueError('Supply 1..1000 molecule IDs')
        label=args.get('label')
        if not isinstance(label,str) or not 1<=len(label)<=80:raise ValueError('A label is required')
        with readonly(item['members']) as db:
            if any(db.execute('SELECT 1 FROM members WHERE mid=?',(mid,)).fetchone() is None for mid in ids):raise ValueError('Label only members of this set')
        with sqlite3.connect(root/'candidate-labels.sqlite') as db:
            db.execute('CREATE TABLE IF NOT EXISTS labels(mid TEXT,label TEXT,note TEXT,created REAL,PRIMARY KEY(mid,label))')
            db.executemany('INSERT OR REPLACE INTO labels VALUES(?,?,?,?)',[(mid,label,str(args.get('note',''))[:2000],time.time()) for mid in ids])
        return dict(label=label,molecules=len(ids),scope='Project annotation; does not imply experimentally verified activity')
    allowed={'operation','name','task_id','set_id','other_set_id','limit','score_min','exclude_labels','max_per_scaffold','start_rank','count'}
    if operation not in {'import','union','intersection','difference','filter','diverse','export'} or set(args)-allowed:
        raise ValueError('Unsupported candidate operation or fields')
    values=dict(args);inputs={}
    if operation=='import':
        job=app.task(sid,args.get('task_id'))
        if job['status']!='complete' or not job.get('report'):raise ValueError('Import requires a completed ranking task')
        report=ensure_within(Path(job['report']),root)
        values['task_report']=str(report);inputs.update(fingerprint([report]))
    else:
        item=get(root,args.get('set_id'));inputs.update(item['outputs'])
        inputs.update(fingerprint([root/'candidate-sets'/item['set_id']/'manifest.json']))
        if operation in {'union','intersection','difference'}:
            other=get(root,args.get('other_set_id'));inputs.update(other['outputs'])
            inputs.update(fingerprint([root/'candidate-sets'/other['set_id']/'manifest.json']))
            if (item['catalog_sha256'],item['batch'])!=(other['catalog_sha256'],other['batch']):raise ValueError('Candidate sets use different library versions')
    if operation!='export' and (not isinstance(args.get('name'),str) or not 1<=len(args['name'])<=160):raise ValueError('Name the new candidate set')
    for key in ('limit','count','start_rank','max_per_scaffold'):
        if key in args and (type(args[key]) is not int or not 1<=args[key]<=1000000):raise ValueError('Invalid '+key)
    if 'score_min' in args and (type(args['score_min']) not in (float,int) or not __import__('math').isfinite(args['score_min'])):raise ValueError('Invalid score_min')
    if 'exclude_labels' in args and (not isinstance(args['exclude_labels'],list) or any(not isinstance(x,str) for x in args['exclude_labels'])):raise ValueError('exclude_labels must be strings')
    if operation=='export' and 'count' not in args:raise ValueError('Specify the MOL2 delivery count')
    # Freeze label membership at submission so later annotations cannot alter this task.
    if args.get('exclude_labels'):
        path=root/'candidate-labels.sqlite';excluded=[]
        if path.exists():
            with readonly(path) as db:
                for label in args['exclude_labels']:excluded.extend(r[0] for r in db.execute('SELECT mid FROM labels WHERE label=?',(label,)))
        values['excluded_molecule_ids']=sorted(set(excluded))
    return queue(app,sid,'candidate_sets',values,inputs)


def imported(args,root,cfg):
    report=read(ensure_within(Path(args['task_report']),root))
    options=[]
    for step in report.get('steps',{}).values():
        if step.get('status')!='complete':continue
        value=step.get('result',{})
        if value.get('report'):
            nested=read(ensure_within(Path(value['report']),root));value={**value,**nested}
        ranking=value.get('outputs',{}).get('ranking') or value.get('ranking')
        if ranking and isinstance(ranking,str) and ranking.endswith('ranking.sqlite'):
            options.append(ensure_within(Path(ranking),root))
    if len(options)!=1:raise ValueError('Choose a task with one completed molecule ranking')
    ranking=options[0];source=read(ranking.parent/'report.json');check_hashes(source['outputs'])
    if str(ranking) not in source['outputs']:raise ValueError('Ranking is not sealed by its report')
    batch=Path(cfg.get('search',{}).get('batch','')).resolve();catalog=batch/'artifacts/catalog.json'
    digest=sha(catalog)
    proof=source.get('sources',{})
    protocol=ranking.parent/'protocol.json'
    if protocol.exists():proof={**proof,**read(protocol).get('sources',{})}
    if proof.get(str(catalog))!=digest:raise ValueError('Ranking/library provenance does not match the configured library')
    return ranking,batch,digest


def execute(args,directory,root,cfg):
    operation=args['operation'];root=Path(root);directory=Path(directory)
    if operation=='export':return export(args,directory,root)
    identifier=args['output_id'];folder=root/'candidate-sets'/identifier;folder.mkdir(parents=True,exist_ok=True)
    manifest=folder/'manifest.json'
    if manifest.exists():return get(root,identifier)
    members=folder/'members.sqlite';members.unlink(missing_ok=True)
    db=sqlite3.connect(members);db.execute('CREATE TABLE members(rank INTEGER UNIQUE,mid TEXT PRIMARY KEY,payload TEXT)')
    parents=[];scoring='source_ranking';count=0
    try:
        if operation=='import':
            ranking,batch,digest=imported(args,root,cfg)
            with readonly(ranking) as source:
                limit=args.get('limit',1000000000)
                for rank,mid in source.execute('SELECT rank,mid FROM ranking ORDER BY rank LIMIT ?',(limit,)):
                    value=source.execute('SELECT p.payload FROM poses p JOIN template_ranks t ON p.mid=t.mid AND p.template=t.template WHERE p.mid=? ORDER BY t.rank,p.template LIMIT 1',(mid,)).fetchone()
                    if not value:raise ValueError('Ranking member has no representative pose')
                    row=json.loads(value[0]);row.update(molecule_id=mid,source_rank=rank,ranking_source=str(ranking))
                    row['score']=row.get('gaussian_score',row.get('optional_score'))
                    row['score_kind']='gaussian_shape_color' if 'gaussian_score' in row else 'source_optional'
                    count+=1
                    db.execute('INSERT INTO members VALUES(?,?,?)',(count,mid,json.dumps(row)))
        else:
            item=get(root,args['set_id']);batch=Path(item['batch']);digest=item['catalog_sha256'];parents=[item['set_id']]
            other=get(root,args['other_set_id']) if operation in {'union','intersection','difference'} else None
            if other:parents.append(other['set_id'])
            excluded=set(args.get('excluded_molecule_ids',[]));scaffolds={};chemistry=None
            if operation=='diverse':
                from .chemical_companion import ChemicalCompanionReader
                from .preselection_full import scaffold_key
                if sha(batch/'artifacts/catalog.json')!=digest:raise ValueError('Library changed before scaffold grouping')
                chemistry=ChemicalCompanionReader(batch/'chemical/catalog.json')
            secondary=readonly(other['members']) if other else None
            try:
                for source_item in ([item,other] if operation=='union' else [item]):
                    with readonly(source_item['members']) as source:
                        for mid,payload in source.execute('SELECT mid,payload FROM members ORDER BY rank'):
                            if count>=args.get('limit',1000000000):break
                            if db.execute('SELECT 1 FROM members WHERE mid=?',(mid,)).fetchone():continue
                            in_other=secondary.execute('SELECT 1 FROM members WHERE mid=?',(mid,)).fetchone() if secondary else None
                            if operation=='intersection' and not in_other:continue
                            if operation=='difference' and in_other:continue
                            row=json.loads(payload)
                            if mid in excluded:continue
                            if 'score_min' in args:
                                if row.get('score') is None or row.get('score_kind')!='gaussian_shape_color':raise ValueError('score_min requires a Gaussian shape/color source; preserve source rank for other methods')
                                if row['score']<args['score_min']:continue
                            if chemistry:
                                chem=chemistry.get(row['global_id'])
                                if chem.molecule_id!=mid:raise ValueError('Scaffold molecule identity mismatch')
                                try:key,_=scaffold_key(chem)
                                except (ValueError,RuntimeError):key='unresolved-'+mid
                                if scaffolds.get(key,0)>=args.get('max_per_scaffold',1):continue
                                scaffolds[key]=scaffolds.get(key,0)+1;row['scaffold_group']=key
                            count+=1;db.execute('INSERT INTO members VALUES(?,?,?)',(count,mid,json.dumps(row)))
            finally:
                if secondary:secondary.close()
        db.commit()
    finally:db.close()
    result=dict(kind='candidate_set',status='complete',set_id=identifier,name=args['name'],count=count,
        members=str(members),batch=str(batch),catalog_sha256=digest,operation=operation,parents=parents,
        created=time.time(),ordering='Preserve left/source rank; union appends unseen right members. No cross-method score fusion.',
        outputs=fingerprint([members]))
    _atomic_json(manifest,result);return result


def export(args,directory,root):
    from .budget_export import export as deliver
    item=get(root,args['set_id']);batch=Path(item['batch'])
    if sha(batch/'artifacts/catalog.json')!=item['catalog_sha256']:raise ValueError('Library version changed')
    run=directory/'ranking';run.mkdir(exist_ok=True);ranking=run/'ranking.sqlite'
    if not (run/'report.json').exists():
        ranking.unlink(missing_ok=True)
        with sqlite3.connect(ranking) as db,readonly(item['members']) as source:
            db.executescript('CREATE TABLE ranking(rank,mid,rrf,template_support); CREATE TABLE poses(mid,template,payload); CREATE TABLE template_ranks(mid,template,rank);')
            for rank,mid,payload in source.execute('SELECT rank,mid,payload FROM members ORDER BY rank'):
                db.execute('INSERT INTO ranking VALUES(?,?,?,1)',(rank,mid,1/(60+rank)))
                db.execute('INSERT INTO poses VALUES(?,?,?)',(mid,'candidate_set',payload))
                db.execute('INSERT INTO template_ranks VALUES(?,?,?)',(mid,'candidate_set',rank))
            db.execute('CREATE INDEX pose_lookup ON poses(mid,template)');db.execute('CREATE INDEX rank_lookup ON template_ranks(mid,template)')
        _atomic_json(run/'report.json',dict(status='complete',outputs=fingerprint([ranking])))
    result=deliver(batch,run,directory/'mol2',args.get('start_rank',1),args['count'])
    return dict(kind='candidate_export',status='complete',set_id=item['set_id'],exported_molecules=result['exported_molecules'],
        shortfall=result['shortfall'],original_mol2_directory=str(directory/'mol2/original'),posed_mol2_directory=str(directory/'mol2/posed'),
        names=str(directory/'mol2/molecules.csv'),export_manifest=str(directory/'mol2/report.json'),
        outputs=fingerprint([directory/'mol2/report.json',directory/'mol2/molecules.csv',directory/'mol2/molecules.jsonl']))
