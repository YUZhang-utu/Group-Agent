"""Sealed project-owned requests for preparation and candidate operations."""
import json
from pathlib import Path
import time
import uuid

from .gaussian_batch import _atomic_json
from .library_acceptance import read, sha
from .project_context import ensure_within


def queue(app,sid,kind,args,inputs=None):
    from .prompt_workflow import create_plan
    ctx=app.context_for(sid); root=Path(ctx['project_root'])
    identifier=uuid.uuid4().hex[:24]
    if kind=='candidate_sets':args=dict(args,output_id=identifier)
    folder=root/'research-requests'/identifier; folder.mkdir(parents=True)
    request=dict(kind=kind,arguments=args,input_hashes=inputs or {})
    path=folder/'request.json';_atomic_json(path,request)
    title=('Prepare molecules' if kind=='molecule_prepare' else 'Candidates: '+args['operation'])
    title+=' / '+str(args.get('name') or args.get('set_id') or args.get('attachment_id') or 'literal SMILES')[:160]
    plan=dict(version=1,summary=title,clarifications=[],steps=[dict(id='work',action='research_task',
        params=dict(request_id=identifier,request_sha256=sha(path)))])
    plan_path=create_plan(Path(ctx['db']),ctx['user_id'],ctx['project_id'],json.dumps(args),local_plan=plan)
    jid=uuid.uuid4().hex[:16]
    with app.connect() as db:
        db.execute('INSERT INTO jobs(id,session,request,provider,profile,status,plan,report,log,created) VALUES(?,?,?,?,?,?,?,?,?,?)',
            (jid,sid,plan['summary'],'local-coordinator','','queued',str(plan_path),str(plan_path.parent/'execution/report.json'),str(app.root/(jid+'.log')),time.time()))
    return dict(status='queued',task_id=jid,project_id=ctx['project_id'])


def request_file(execution,params):
    root=execution.parent.parent.parent
    path=ensure_within(root/'research-requests'/params['request_id']/'request.json',root)
    if sha(path)!=params['request_sha256']:raise ValueError('Research request changed')
    return path


def execute(params,directory,execution,cfg,compute):
    from .screening_selection import check_hashes
    request=read(request_file(execution,params));check_hashes(request['input_hashes'])
    if request['kind']=='molecule_prepare':
        from .molecule_preparation import execute as prepare
        result=prepare(request['arguments'],directory)
    elif request['kind']=='candidate_sets':
        from .candidate_sets import execute as candidates
        result=candidates(request['arguments'],directory,execution.parent.parent.parent,cfg)
    else:raise ValueError('Unknown research task')
    check_hashes(request['input_hashes'])
    result.update(report=str(directory/'report.json'))
    _atomic_json(directory/'report.json',result)
    return result
