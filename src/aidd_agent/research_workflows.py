"""Durable project dependency graphs. Uncertain dispatch is never retried."""
import json
import re
import time
import uuid

TOOLS = {'molecule_prepare','ligand_library','candidate_sets','workflow','target_search','af3_batch','pose_qc','screening_benchmark'}
REFERENCE = re.compile(r'\{\{(task|query|set):([A-Za-z][A-Za-z0-9_-]{0,39})\}\}')


def initialize(app):
    with app.connect() as db:
        db.execute('CREATE TABLE IF NOT EXISTS research_workflows(id TEXT PRIMARY KEY,project_id TEXT,session TEXT,payload TEXT)')


def save(app, row):
    row['updated'] = time.time()
    with app.connect() as db:
        db.execute('INSERT OR REPLACE INTO research_workflows VALUES(?,?,?,?)',
                   (row['id'],row['project_id'],row['session'],json.dumps(row)))


def rows(app, sid):
    pid = app.context_for(sid)['project_id']
    with app.connect() as db:
        return [json.loads(r['payload']) for r in db.execute('SELECT payload FROM research_workflows WHERE project_id=? ORDER BY id',(pid,))]


def validate_nodes(nodes):
    from .domain_tools import validate_arguments
    if not isinstance(nodes,list) or not 1 <= len(nodes) <= 20: raise ValueError('Use 1..20 workflow nodes')
    seen = set()
    for node in nodes:
        if not isinstance(node,dict) or set(node) != {'id','tool','arguments','depends_on'}: raise ValueError('Node requires id, tool, arguments and depends_on')
        nid = node['id']; deps = node['depends_on']; args = node['arguments']
        if not isinstance(nid,str) or not re.fullmatch('[A-Za-z][A-Za-z0-9_-]{0,39}',nid) or nid in seen: raise ValueError('Node IDs must be unique simple identifiers')
        if not isinstance(deps,list) or any(not isinstance(v,str) or v not in seen for v in deps): raise ValueError('Dependencies must name preceding nodes')
        if node['tool'] not in TOOLS: raise ValueError('Unsupported workflow tool')
        validate_arguments(node['tool'],args)
        if node['tool']=='ligand_library' and args.get('operation')!='start': raise ValueError('Workflow ligand node must start a search')
        if node['tool']=='candidate_sets' and args.get('operation') not in {'import','union','intersection','difference','filter','diverse','export'}: raise ValueError('Workflow candidate node must queue work')
        if node['tool']=='workflow' and args.get('decision',{}).get('intent') in {'resume','cancel','pymol','pymol_agent','status','results','capabilities'}: raise ValueError('Workflow nodes must create scientific tasks, not control existing tasks')
        for _,ref in REFERENCE.findall(json.dumps(args)):
            if ref not in deps: raise ValueError('Each bound reference must be an explicit dependency')
        seen.add(nid)


def handle(app,sid,args,provider='deepseek',request=''):
    operation=args.get('operation','status')
    with app.lock:
        if operation=='create':
            objective=args.get('objective')
            if not isinstance(objective,str) or not 1<=len(objective)<=4000: raise ValueError('A concrete objective is required')
            validate_nodes(args.get('nodes'))
            row=dict(id=uuid.uuid4().hex[:24],project_id=app.context_for(sid)['project_id'],session=sid,
                objective=objective,request=request or objective,provider=provider,status='active',created=time.time(),
                nodes=[dict(n,status='pending') for n in args['nodes']])
            save(app,row)
            return dict(status='active',workflow_id=row['id'],objective=objective,message='Saved durable workflow; the background worker will dispatch its ready nodes.')
        owned=rows(app,sid)
        if operation=='status':
            return dict(workflows=[r for r in owned if not args.get('workflow_id') or r['id']==args['workflow_id']])
        found=[r for r in owned if r['id']==args.get('workflow_id')]
        if len(found)!=1: raise ValueError('Workflow is unavailable in this project')
        row=found[0]
        if operation in {'pause','cancel'}:
            if row['status']=='complete': raise ValueError('Workflow already completed')
            row['status']='paused' if operation=='pause' else 'cancelled'
        elif operation=='resume':
            if row['status'] in {'cancelled','complete'}: raise ValueError('This workflow cannot be resumed')
            if any(n['status']=='dispatching' for n in row['nodes']): raise ValueError('Uncertain dispatch: reconcile the recorded task before resuming')
            row['status']='active';row.pop('error',None)
        elif operation in {'resolve','reconcile'}:
            matches=[n for n in row['nodes'] if n['id']==args.get('node_id')]
            if len(matches)!=1 or row['status'] in {'cancelled','complete'}: raise ValueError('Node cannot be changed')
            node=matches[0]
            if operation=='resolve':
                if node['status']!='pending': raise ValueError('Only an unsubmitted node can be resolved')
                candidate=dict(node,arguments=args.get('arguments'))
                proposed=[{k:(candidate if n is node else n)[k] for k in ('id','tool','arguments','depends_on')} for n in row['nodes']]
                validate_nodes(proposed);node['arguments']=args['arguments']
                row['request']+='\nExplicit workflow resolution:\n'+request
            else:
                if node['status']!='dispatching': raise ValueError('Reconcile only an uncertain dispatch')
                task=app.task(sid,args.get('task_id'))
                candidates=node.get('candidate_task_ids',[])
                if task['id'] not in candidates: raise ValueError('Task is not a recorded candidate for this dispatch; inspect its execution log')
                node.update(task_id=task['id'],status='submitted')
            row.update(status='active');row.pop('error',None)
        else: raise ValueError('Unsupported workflow operation')
        save(app,row)
        return dict(workflow_id=row['id'],status=row['status'],message='Future dispatch updated; current scientific tasks are controlled separately.')


def resolve(app,row,value):
    if isinstance(value,list): return [resolve(app,row,v) for v in value]
    if isinstance(value,dict): return {k:resolve(app,row,v) for k,v in value.items()}
    if not isinstance(value,str): return value
    def substitute(match):
        kind,nid=match.groups();node=next(n for n in row['nodes'] if n['id']==nid)
        if node['status']!='complete': raise ValueError('Referenced node is not complete')
        task_id=node['task_id']
        if kind=='task': return task_id
        if kind=='query':
            from .ligand_library import queries
            choices=queries(app,row['session'],task_id)
            if len(choices)!=1: raise ValueError('Choose an explicit query_id for '+nid+': '+json.dumps([q['query_id'] for q in choices[:100]]))
            return choices[0]['query_id']
        from .library_acceptance import read
        from .project_context import ensure_within
        from pathlib import Path
        job=app.task(row['session'],task_id);root=Path(app.context_for(row['session'])['project_root'])
        report=read(ensure_within(Path(job['report']),root))
        choices=[s['result']['set_id'] for s in report.get('steps',{}).values() if s.get('status')=='complete' and s.get('result',{}).get('kind')=='candidate_set']
        if len(choices)!=1: raise ValueError('Dependency must produce one candidate set')
        return choices[0]
    return REFERENCE.sub(substitute,value)


def advance(app,row):
    from .domain_tools import DomainTools
    sid=row['session']
    for node in row['nodes']:
        if node['status']=='dispatching':
            # A crash may occur between dispatch and receipt persistence. Never guess.
            before=set(node.get('before_task_ids',[]))
            node['candidate_task_ids']=sorted(set(node.get('candidate_task_ids',[])) | {j['id'] for j in app.jobs(sid) if j['id'] not in before and j['created']>=node.get('dispatch_started',0)})
            row.update(status='blocked',error='Uncertain dispatch for '+node['id']+'; inspect candidate tasks and reconcile. No automatic resubmission.')
            save(app,row);return
        if node['status']=='submitted':
            task=app.task(sid,node['task_id'])
            if task['status']=='complete': node['status']='complete'
            elif task['status'] in {'failed','blocked','cancelled','interrupted'}:
                row.update(status='blocked',error='Task '+task['id']+': '+task['status']+'. '+str(task.get('error') or 'Inspect its report; resume the existing task if appropriate.'))
                save(app,row);return
            else: save(app,row);return
        if node['status']=='pending':
            if any(next(n for n in row['nodes'] if n['id']==dep)['status']!='complete' for dep in node['depends_on']): continue
            try: arguments=resolve(app,row,node['arguments'])
            except (ValueError,OSError,KeyError) as exc:
                row.update(status='blocked',error=str(exc));save(app,row);return
            node.update(status='dispatching',dispatch_started=time.time(),before_task_ids=[j['id'] for j in app.jobs(sid)],resolved_arguments=arguments)
            save(app,row)
            try:
                result=DomainTools(app,sid,row['request'],row['provider']).call(node['tool'],arguments)
                ids=[result['task_id']] if result.get('task_id') else [t['id'] for t in result.get('tasks',[])]
                if len(ids)!=1: raise ValueError('Expected exactly one queued scientific task')
                app.task(sid,ids[0]);node.update(status='submitted',task_id=ids[0]);save(app,row);return
            except Exception as exc:
                candidates=[j['id'] for j in app.jobs(sid) if j['id'] not in node['before_task_ids']]
                node['candidate_task_ids']=candidates
                if not candidates: node['status']='pending'
                row.update(status='blocked',error=type(exc).__name__+': '+str(exc));save(app,row);return
    if all(n['status']=='complete' for n in row['nodes']):
        row['status']='complete'
        app.message(sid,'assistant','Project workflow '+row['id']+' completed. '+row['objective'])
    save(app,row)


def tick(app):
    with app.lock:
        with app.connect() as db:
            active=[json.loads(r['payload']) for r in db.execute('SELECT payload FROM research_workflows')]
        for row in active:
            if row['status']!='active': continue
            try: advance(app,row)
            except Exception as exc:
                row.update(status='blocked',error=type(exc).__name__+': '+str(exc));save(app,row)
