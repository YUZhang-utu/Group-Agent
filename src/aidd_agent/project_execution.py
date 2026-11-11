"""Concrete-plan mission checks and content-based task reuse before scientific launch."""
import hashlib
import json
from pathlib import Path

from .library_acceptance import read,sha
from .project_context import ensure_within
from .project_guidance import encoded,mission,render


def expanded_steps(root,plan):
    folders={'research_task':'research-requests','ligand_library_search':'ligand-library','target_query_prepare':'target-search','target_library_search':'target-search'}
    result=[]
    for step in plan['steps']:
        p=step['params'];extra=None
        if step['action'] in folders:
            # Use the adapters' actual request resolver; their locations can differ.
            execution=root/'runs'/'placeholder'/'execution'
            module={'research_task':'research_tasks','ligand_library_search':'ligand_library'}.get(step['action'])
            if module:
                from importlib import import_module
                path=import_module('aidd_agent.'+module).request_file(execution,p)
            else:
                from .target_search_workflow import request_file
                path=request_file(execution,p)
            extra=read(path)
        result.append(dict(step,request=extra))
    return result


def validate(root,plan,current,cfg):
    rules=current['rules'];actions={s['action'] for s in plan['steps']}
    if 'allowed_actions' in rules and actions-set(rules['allowed_actions']):raise ValueError('Mission action allowlist rejects: '+', '.join(sorted(actions-set(rules['allowed_actions']))))
    forbidden=actions & set(rules.get('forbidden_actions',[]))
    if forbidden:raise ValueError('Mission forbids: '+', '.join(sorted(forbidden)))
    if rules.get('full_library') and actions & {'consensus_budget','search_3d'}:raise ValueError('Mission requires full-catalog search; a budget/retrieval route cannot silently replace it')
    if rules.get('library_catalog_sha256') and actions & {'ligand_library_search','target_library_search','consensus_budget','search_3d','guided_funnel','full_library_screen'}:
        batch=cfg.get('search',{}).get('batch')
        if not batch or sha(Path(batch)/'artifacts/catalog.json')!=rules['library_catalog_sha256']:raise ValueError('Configured library differs from the pinned project version')
    from .af3_inputs import seeds
    for step in expanded_steps(root,plan):
        p=step['params']
        if step['action']=='af3_prepare':
            for field in ('protein_sequences','ligand_smiles'):
                if field in rules and p.get(field)!=rules[field]:raise ValueError('Mission requires exact '+field+'; use the pinned input without substitution')
            if 'protein_copies' in rules and p.get('protein_copies',1)!=rules['protein_copies']:raise ValueError('Protein assembly/copy count differs from mission')
            if 'seeds' in rules and seeds(p)!=rules['seeds']:raise ValueError('AF3 seeds differ from mission')
            if rules.get('protein_sequences') and ('start' in p or 'end' in p):
                if p.get('start',1)!=1 or any(p.get('end',len(s))!=len(s) for s in rules['protein_sequences']):raise ValueError('Pinned mission sequence cannot be trimmed')
    return dict(status='typed_rules_passed',actions=sorted(actions),scope='Configured structured constraints enforced; free-text goal relevance requires explicit requirement linkage and review')


def normalized(value,root):
    if isinstance(value,list):return [normalized(v,root) for v in value]
    if isinstance(value,dict):
        out={}
        for k,v in value.items():
            if k=='output_id':continue
            if k in {'input_hashes','sources'} and isinstance(v,dict):
                for path,digest in v.items():
                    if sha(Path(path))!=digest:raise ValueError('Input changed before reuse check: '+path)
                out[k]=sorted(v.values());continue
            out[k]=normalized(v,root)
        return out
    if isinstance(value,str):
        try:
            p=Path(value)
            if p.is_absolute() and p.is_file():return dict(content_sha256=sha(p),extension=p.suffix)
        except (OSError,ValueError):pass
    return value


def signature(root,plan,cfg,config_inputs,search_only=False):
    rows=expanded_steps(root,plan);lookup={s['id']:'step-'+str(i) for i,s in enumerate(rows)}
    normalized_rows=[]
    for s in rows:
        p=dict(s['params'])
        request=s['request']
        if search_only and s['action']=='ligand_library_search':
            request={k:v for k,v in request.items() if k!='export_molecules'}
        if s['request'] is not None:
            p.pop('request_id',None);p.pop('request_sha256',None)
        if s['action']=='af3_prepare':p.pop('name',None)
        for field in ('input_step','protein_step','protein_steps','receptor_step'):
            if field in p:p[field]=[lookup.get(v,v) for v in p[field]] if isinstance(p[field],list) else lookup.get(p[field],p[field])
        for field in ('source_run','receptor_run'):
            if field in p:
                directory=ensure_within(root/'runs'/p[field]/'execution',root)
                if not directory.exists():raise ValueError('Referenced upstream run is unavailable')
                p[field]=dict(run_id=p[field],artifacts={str(v.relative_to(directory)):sha(v) for v in directory.rglob('*') if v.is_file() and not v.name.endswith(('.lock','.log','.tmp'))})
        normalized_rows.append(dict(action=s['action'],params=normalized(p,root),request=normalized(request,root)))
    # Include runtime file contents, all Python code, library catalog and referenced inputs.
    stamps={str(p):sha(Path(p)) for p in config_inputs}
    batch=cfg.get('search',{}).get('batch')
    if batch:
        for name in ('artifacts/catalog.json','chemical/catalog.json','registry.sqlite3'):
            path=Path(batch)/name
            if path.exists():stamps[name]=sha(path)
    code={p.name:sha(p) for p in sorted(Path(__file__).parent.glob('*.py'))}
    import sys
    import platform
    from importlib.metadata import version,PackageNotFoundError
    environment=dict(python=sys.version,platform=platform.platform())
    for package in ('numpy','rdkit'):
        try:environment[package]=version(package)
        except PackageNotFoundError:environment[package]='unavailable'
    return hashlib.sha256(encoded(dict(steps=normalized_rows,configuration=normalized(cfg,root),stamps=stamps,code=code,environment=environment)).encode()).hexdigest()


def reusable(root,plan):
    """Only automatically reuse lanes whose complete inputs are fingerprinted here."""
    rows=expanded_steps(root,plan)
    if not rows:return False
    for s in rows:
        if s['action']=='research_task':
            if s['request']['kind']!='molecule_prepare':return False
        elif s['action']=='af3_prepare':
            if not s['params'].get('protein_sequences'):return False
        elif s['action'] not in {'ligand_library_search','target_query_prepare','target_library_search'}:
            # Live databases, mutable candidate labels and external compute runtime
            # data require their existing explicit reuse/resume adapters.
            return False
    return True


def outputs(job,root):
    """Seal all execution artifacts, including model/pose files, not just a summary."""
    if not job.get('report') or not job.get('plan'):raise ValueError('No completed artifact inventory available')
    directory=ensure_within(Path(job['plan']).parent/'execution',root)
    if not directory.exists():raise ValueError('Execution artifacts unavailable')
    result={str(p):sha(ensure_within(p,root)) for p in directory.rglob('*') if p.is_file() and not p.name.endswith(('.lock','.log','.tmp'))}
    if not result or str(Path(job['report'])) not in result:raise ValueError('Completed report missing')
    report=read(Path(job['report']))
    if report.get('status')!='complete':raise ValueError('Scientific report is not complete')
    # Set members can live outside the task directory, but must stay in this project.
    for s in report.get('steps',{}).values():
        r=s.get('result',{})
        for path,digest in r.get('outputs',{}).items():
            p=Path(path)
            if len(str(digest))==64 and p.is_absolute():
                p=ensure_within(p,root)
                if sha(p)!=digest:raise ValueError('Recorded output changed')
                result[str(p)]=digest
    return result


def check(app,job):
    from .prompt_workflow import plan_runtime
    sid=job['session'];ctx=app.context_for(sid);root=Path(ctx['project_root']);current=mission(app,sid)
    path=ensure_within(Path(job['plan']),root);envelope=read(path)
    seal=read(path.parent/'plan-seal.json')
    if seal.get('plan_sha256')!=sha(path):raise ValueError('Plan changed before project alignment check')
    plan=envelope['plan'];cfg,inputs=plan_runtime(plan,app.runtime)
    alignment=validate(root,plan,current,cfg);fingerprint=signature(root,plan,cfg,inputs)
    with app.connect() as db:
        old=db.execute('SELECT payload FROM project_task_cards WHERE task_id=?',(job['id'],)).fetchone()
        prior=json.loads(old['payload']) if old else {}
        replica=db.execute('SELECT reason FROM project_replicates WHERE task_id=?',(job['id'],)).fetchone()
        others=[json.loads(r['payload']) for r in db.execute('SELECT payload FROM project_task_cards WHERE project_id=? AND task_id!=?',(ctx['project_id'],job['id']))]
    requirement_ids=prior.get('requirement_ids',[])
    actions=set(alignment['actions']);requirements=current['requirements']
    if not requirement_ids:
        requirement_ids=[r['id'] for r in requirements if set(r['actions']) & actions]
        if not requirement_ids and len(requirements)==1 and not requirements[0]['actions']:requirement_ids=[requirements[0]['id']]
    if not requirement_ids or set(requirement_ids)-{r['id'] for r in requirements}:raise ValueError('Link this task to a current mission requirement before execution')
    selected=[r for r in requirements if r['id'] in requirement_ids]
    if all(r['actions'] for r in selected) and actions-set(a for r in selected for a in r['actions']):raise ValueError('Plan includes actions outside its linked mission requirements')
    if prior.get('executed_revision') and prior['executed_revision']!=current['revision']:
        raise ValueError('Mission changed since this task began; keep its original snapshot and explicitly create a revised task')
    card=dict(prior,task_id=job['id'],mission_revision=current['revision'],mission_snapshot=current,
        requirement_ids=requirement_ids,contribution=prior.get('contribution',plan['summary']),alignment=alignment,signature=fingerprint,
        plan_sha256=sha(path),reuse_task_id=None)
    card['reuse_eligible']=reusable(root,plan)
    if len(plan['steps'])==1 and plan['steps'][0]['action']=='ligand_library_search':
        card['search_signature']=signature(root,plan,cfg,inputs,search_only=True)
    if replica:card['replicate_reason']=replica['reason']
    reuse=None
    for prior_card in others if not replica and card['reuse_eligible'] else []:
        if not prior_card.get('reuse_eligible'):continue
        same_search=card.get('search_signature') and card['search_signature']==prior_card.get('search_signature')
        if same_search and prior_card.get('signature')!=fingerprint:
            oldjob=app.task(sid,prior_card['task_id'])
            if oldjob['status']=='complete' and prior_card.get('output_hashes'):
                from .screening_selection import check_hashes
                check_hashes(prior_card['output_hashes'])
                raise ValueError('Completed full-library ranking already exists in task '+oldjob['id']+
                    '. Only the export count changed. Use candidate_sets import from that task, then export the requested count; no new full-library scan is needed.')
        if prior_card.get('signature')!=fingerprint:continue
        oldjob=app.task(sid,prior_card['task_id'])
        if oldjob['status']=='complete':
            inventory=prior_card.get('output_hashes')
            if not inventory:continue
            from .screening_selection import check_hashes
            check_hashes(inventory)
            card['reuse_task_id']=oldjob['id'];card['output_hashes']=inventory;reuse=oldjob;break
        if oldjob['status'] in {'queued','planning','running','failed','blocked','interrupted'}:
            raise ValueError('Equivalent project task '+oldjob['id']+' is '+oldjob['status']+'. Inspect/resume that task; use an explicit independent-replicate control only if a separate calculation is intended.')
    with app.connect() as db:db.execute('INSERT OR REPLACE INTO project_task_cards VALUES(?,?,?)',(job['id'],ctx['project_id'],encoded(card)))
    render(app,sid)
    return card,reuse


def started(app,job):
    with app.connect() as db:
        row=db.execute('SELECT payload FROM project_task_cards WHERE task_id=?',(job['id'],)).fetchone()
        if not row:return
        card=json.loads(row['payload']);card['executed_revision']=card['mission_revision']
        db.execute('UPDATE project_task_cards SET payload=? WHERE task_id=?',(encoded(card),job['id']))


def finished(app,job):
    root=Path(app.context_for(job['session'])['project_root'])
    with app.connect() as db:
        row=db.execute('SELECT payload FROM project_task_cards WHERE task_id=?',(job['id'],)).fetchone()
        if not row:return
        card=json.loads(row['payload'])
    try:card['output_hashes']=outputs(job,root)
    except (OSError,ValueError,KeyError) as exc:card['reuse_unavailable']=str(exc)
    with app.connect() as db:db.execute('UPDATE project_task_cards SET payload=? WHERE task_id=?',(encoded(card),job['id']))
