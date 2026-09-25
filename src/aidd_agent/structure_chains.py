"""Durable, conservative continuation of explicitly requested structure workflows."""
import json
import time
import uuid


def initialize(app):
    with app.connect() as db:
        db.execute('CREATE TABLE IF NOT EXISTS structure_chains(id TEXT PRIMARY KEY,session TEXT,status TEXT,payload TEXT,updated REAL)')
        db.execute("UPDATE structure_chains SET status='uncertain' WHERE status='dispatching'")


def rows(app,sid):
    with app.connect() as db:
        return [dict(json.loads(r['payload']),id=r['id'],status=r['status']) for r in
                db.execute('SELECT * FROM structure_chains WHERE session=? ORDER BY updated',(sid,))]


def save(app,chain,status):
    payload={k:v for k,v in chain.items() if k not in {'id','status'}}
    with app.connect() as db:
        db.execute('UPDATE structure_chains SET status=?,payload=?,updated=? WHERE id=?',
                   (status,json.dumps(payload),time.time(),chain['id']))
    chain['status']=status


def control(app,sid,provider,request,args):
    app.session(sid)
    operation=args.get('operation','status')
    with app.lock:
        if operation=='status':return dict(chains=rows(app,sid))
        if operation=='start':
            if set(args)-{'operation','task_id','goal','reference','budget'}:raise ValueError('Unknown chain fields')
            job=app.task(sid,args['task_id']);goal=args.get('goal','recommendation')
            if goal not in {'recommendation','budget'}:raise ValueError('Goal must be recommendation or budget')
            if goal=='budget' and not app.allow_compute:raise ValueError('Budget continuation requires --allow-compute')
            from .chat_agent import validate_route
            reference=args.get('reference',{});budget=args.get('budget',{})
            validate_route(dict(intent='consensus',task_id=job['id'],message='',request='',reference=reference))
            from .budget_workflow import validate_budget
            validate_budget(budget)
            if any(c['initial_task']==job['id'] and c['status'] not in {'complete','cancelled'} for c in rows(app,sid)):
                raise ValueError('A chain already tracks this source task; inspect its status')
            chain=dict(id=uuid.uuid4().hex,initial_task=job['id'],current_task=job['id'],goal=goal,
                       reference=reference,budget=budget,provider=provider,request=request,history=[],error=None)
            with app.connect() as db:
                db.execute('INSERT INTO structure_chains VALUES(?,?,?,?,?)',
                           (chain['id'],sid,'active',json.dumps(chain),time.time()))
            return dict(status='queued',chain_id=chain['id'],current_task=job['id'],goal=goal)
        if set(args)-{'operation','chain_id'}:raise ValueError('Only chain_id is accepted for chain controls')
        chain=next((c for c in rows(app,sid) if c['id']==args.get('chain_id')),None)
        if chain is None:raise ValueError('Unknown chain in this session')
        if operation not in {'pause','resume','cancel'}:raise ValueError('Unknown chain operation')
        if operation=='resume' and chain['status'] in {'uncertain','dispatching','complete','cancelled'}:
            raise ValueError('This chain cannot auto-resume. Inspect linked jobs; uncertain dispatch must be reconciled manually.')
        chain['error']=None
        save(app,chain,{'pause':'paused','resume':'active','cancel':'cancelled'}[operation])
        return dict(chain=chain,scope='Controls future continuation only; an existing scientific task is not cancelled or resumed')


def tick(app):
    """Called by the existing single worker, before choosing its next queued job."""
    from .domain_tools import DomainTools
    with app.lock:
        with app.connect() as db:
            active=[dict(r) for r in db.execute("SELECT * FROM structure_chains WHERE status='active'")]
        for record in active:
            chain=dict(json.loads(record['payload']),id=record['id']);sid=record['session']
            try:
                job=app.task(sid,chain['current_task'])
                if job['status'] in {'queued','planning','running'}:continue
                if job['status']!='complete':raise ValueError('Resolve current task '+job['id']+': '+job['status'])
                tools=DomainTools(app,sid,chain['request'],chain['provider'])
                context=tools.call('structure_workflow',dict(task_id=job['id']))['tasks'][0]
                if context.get('error'):raise ValueError(context['error'])
                stages=context['stages']
                if len(stages)!=1:raise ValueError('Choose a task with exactly one supported structure workflow stage')
                stage=stages[0];action=stage['action']
                final={'recommendation':{'consensus_recommend'},'budget':{'consensus_budget'}}[chain['goal']]
                if action in final:
                    if stage['status']!='complete' or stage['evidence'].get('status')!='complete' or stage['evidence'].get('readiness','').startswith('needs_'):
                        raise ValueError('Final stage still requires review')
                    chain['error']=None;save(app,chain,'complete')
                    app.message(sid,'assistant',f"Structure chain {chain['id']} complete. Final task: {job['id']}.")
                    continue
                intent={'structure_diversity':'consensus','structure_consensus':'recommend',
                        'consensus_recommend':'budget','consensus_design':'budget'}.get(action)
                if intent not in stage['next_intents']:raise ValueError('Stage requires input or has no supported continuation: '+action)
                if intent=='budget' and (chain['goal']!='budget' or not app.allow_compute):raise ValueError('Budget execution is not enabled for this chain')
                if intent=='consensus' and not all(chain['reference'].get(k) for k in ('reference_query','target_chain')):
                    raise ValueError('Supply reference_query and target_chain when starting the chain; no reference identity is guessed')
                decision=dict(intent=intent,task_id=job['id'],message='',request='')
                if intent=='consensus':decision['reference']=chain['reference']
                if intent=='budget':decision['budget']=chain['budget']
                chain['pending_decision']=decision;save(app,chain,'dispatching')
                try:
                    result=tools.call('workflow',dict(decision=decision))
                    tasks=result.get('tasks',[])
                    if len(tasks)!=1:raise ValueError('Dispatch returned no unique job receipt')
                    next_task=tasks[0]['id']
                except Exception as exc:
                    chain['dispatch_error']=str(exc) if isinstance(exc,ValueError) else type(exc).__name__
                    chain['error']='Dispatch outcome uncertain; inspect task list before any retry'
                    save(app,chain,'uncertain')
                    app.message(sid,'assistant',f"Structure chain {chain['id']}: {chain['error']}.")
                    continue
                chain['history'].append(dict(source_task=job['id'],intent=intent,task_id=next_task))
                chain['current_task']=next_task;chain.pop('pending_decision',None);chain['error']=None
                save(app,chain,'active')
                app.message(sid,'assistant',f"Structure chain {chain['id']}: queued {intent}, task {next_task}.")
            except (ValueError,KeyError,TypeError,OSError) as exc:
                chain['error']=str(exc);save(app,chain,'blocked')
                app.message(sid,'assistant',f"Structure chain {chain['id']} blocked: {exc}")
