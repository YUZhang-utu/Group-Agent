"""Versioned mission, evidence-linked feedback and an inspectable project dossier."""
import hashlib
import json
from pathlib import Path
import time
import uuid

from .project_context import ensure_within


def encoded(value):
    return json.dumps(value,sort_keys=True,ensure_ascii=False,separators=(',',':'))


def initialize(app):
    with app.connect() as db:
        db.executescript('''CREATE TABLE IF NOT EXISTS project_missions(project_id TEXT,revision INTEGER,payload TEXT,PRIMARY KEY(project_id,revision));
        CREATE TABLE IF NOT EXISTS project_learning(id TEXT PRIMARY KEY,project_id TEXT,kind TEXT,payload TEXT,created REAL);
        CREATE TABLE IF NOT EXISTS project_task_cards(task_id TEXT PRIMARY KEY,project_id TEXT,payload TEXT);
        CREATE TABLE IF NOT EXISTS project_proposals(id TEXT PRIMARY KEY,project_id TEXT,base_revision INTEGER,payload TEXT,status TEXT);
        CREATE TABLE IF NOT EXISTS project_dossier(project_id TEXT PRIMARY KEY,sha256 TEXT);
        CREATE TABLE IF NOT EXISTS project_replicates(task_id TEXT PRIMARY KEY,reason TEXT);''')


def mission(app,sid):
    ctx=app.context_for(sid)
    with app.connect() as db:
        row=db.execute('SELECT payload FROM project_missions WHERE project_id=? ORDER BY revision DESC LIMIT 1',(ctx['project_id'],)).fetchone()
    if row:return json.loads(row['payload'])
    from .research_projects import projects
    meta=next(p for p in projects(app) if p['project_id']==ctx['project_id'])
    original=dict(objective=meta['objective'],target=meta.get('target','unspecified'),
        requirements=[dict(id='R1',text=meta['objective'],actions=[],acceptance='User review of the requested deliverable')],
        rules={},excluded_work='',constraints=[])
    record=dict(original,revision=1,project_id=ctx['project_id'],created=time.time(),source='Project creation metadata',baseline=original)
    with app.connect() as db:
        db.execute('INSERT OR IGNORE INTO project_missions VALUES(?,?,?)',(ctx['project_id'],1,encoded(record)))
    return mission(app,sid)


def validate_mission(value):
    from .prompt_plan import ACTION_FIELDS
    fields={'objective','target','requirements','rules','excluded_work','constraints'}
    if not isinstance(value,dict) or set(value)!=fields:raise ValueError('Mission requires '+', '.join(sorted(fields)))
    for key in ('objective','target','excluded_work'):
        if not isinstance(value[key],str) or len(value[key])>2000 or (key!='excluded_work' and not value[key].strip()):raise ValueError('Invalid mission '+key)
    constraints=value['constraints']
    if not isinstance(constraints,list) or len(constraints)>30 or any(not isinstance(v,str) or not 1<=len(v)<=300 for v in constraints):raise ValueError('Use at most 30 constraints of 1..300 characters')
    requirements=value['requirements'];seen=set()
    if not isinstance(requirements,list) or not 1<=len(requirements)<=20:raise ValueError('Use 1..20 mission requirements')
    for r in requirements:
        if not isinstance(r,dict) or set(r)!={'id','text','actions','acceptance'}:raise ValueError('Requirement requires id,text,actions,acceptance')
        if not isinstance(r['id'],str) or not __import__('re').fullmatch('[A-Za-z][A-Za-z0-9_-]{0,39}',r['id']) or r['id'] in seen:raise ValueError('Unique requirement IDs required')
        seen.add(r['id'])
        if any(not isinstance(r[k],str) or not 1<=len(r[k])<=800 for k in ('text','acceptance')):raise ValueError('Requirement text/acceptance must be 1..800 characters')
        if not isinstance(r['actions'],list) or any(a not in ACTION_FIELDS for a in r['actions']):raise ValueError('Unknown required action')
    rules=value['rules'];allowed={'allowed_actions','forbidden_actions','protein_sequences','protein_copies','seeds','ligand_smiles','full_library','library_catalog_sha256'}
    if not isinstance(rules,dict) or set(rules)-allowed:raise ValueError('Unknown typed mission rule')
    for key in ('allowed_actions','forbidden_actions'):
        if key in rules and (not isinstance(rules[key],list) or any(a not in ACTION_FIELDS for a in rules[key])):raise ValueError('Invalid action allowlist')
    for key in ('protein_sequences','ligand_smiles'):
        if key in rules and (not isinstance(rules[key],list) or not rules[key] or any(not isinstance(v,str) or not v for v in rules[key])):raise ValueError('Invalid '+key)
    if 'protein_sequences' in rules:
        from .af3_inputs import sequence
        rules['protein_sequences']=[sequence(s) for s in rules['protein_sequences']]
    if 'protein_copies' in rules and (type(rules['protein_copies']) is not int or not 1<=rules['protein_copies']<=26):raise ValueError('Invalid protein_copies')
    if 'seeds' in rules:
        from .af3_inputs import seeds
        seeds(dict(seeds=rules['seeds']))
    if 'full_library' in rules and type(rules['full_library']) is not bool:raise ValueError('full_library must be boolean')
    if 'library_catalog_sha256' in rules and not __import__('re').fullmatch('[a-f0-9]{64}',str(rules['library_catalog_sha256'])):raise ValueError('Invalid library hash')
    if len(encoded(value))>16000:raise ValueError('Mission exceeds 16000 characters; use concise constraints and requirements')
    return value


def revise(app,sid,value,base_revision,source):
    """Only direct user controls call this; model tools can propose, not apply."""
    with app.lock:
        old=mission(app,sid)
        if type(base_revision) is not int or base_revision!=old['revision']:raise ValueError('Mission changed; reload before applying your edits')
        validated=validate_mission(json.loads(encoded(value)))
        row=dict(validated,revision=old['revision']+1,project_id=old['project_id'],baseline=old['baseline'],source=source,created=time.time())
        with app.connect() as db:
            db.execute('INSERT INTO project_missions VALUES(?,?,?)',(row['project_id'],row['revision'],encoded(row)))
        render(app,sid)
        return row


def entries(app,sid):
    pid=app.context_for(sid)['project_id']
    with app.connect() as db:
        return [dict(json.loads(r['payload']),id=r['id'],kind=r['kind'],created=r['created']) for r in db.execute('SELECT * FROM project_learning WHERE project_id=? ORDER BY created,id',(pid,))]


def add_entry(app,sid,kind,payload,identifier=None):
    identifier=identifier or uuid.uuid4().hex[:24]
    with app.connect() as db:
        db.execute('INSERT OR IGNORE INTO project_learning VALUES(?,?,?,?,?)',(identifier,app.context_for(sid)['project_id'],kind,encoded(payload),time.time()))
    return identifier


def evidence_value(app,sid,task_id,pointer):
    task=app.task(sid,task_id)
    if task['status']!='complete' or not task.get('report'):raise ValueError('Evidence requires a completed task with a saved report')
    path=ensure_within(Path(task['report']),Path(app.context_for(sid)['project_root']))
    raw=path.read_bytes();value=json.loads(raw)
    if not isinstance(pointer,str) or not pointer.startswith('/'):raise ValueError('Use an exact JSON pointer into the task report')
    for part in pointer[1:].split('/'):
        key=part.replace('~1','/').replace('~0','~');value=value[int(key)] if isinstance(value,list) else value[key]
    if len(encoded(value))>4000:raise ValueError('Select a smaller evidence field')
    return dict(task_id=task_id,report=str(path),report_sha256=hashlib.sha256(raw).hexdigest(),pointer=pointer,value=value,
        verification='Recorded report field; not independent scientific validation')


def handle(app,sid,args,request=''):
    op=args.get('operation','read');current=mission(app,sid)
    if op=='read':
        with app.connect() as db:
            proposals=[dict(id=r['id'],base_revision=r['base_revision'],mission=json.loads(r['payload'])) for r in db.execute("SELECT * FROM project_proposals WHERE project_id=? AND status='proposed'",(current['project_id'],))]
            cards=[json.loads(r['payload']) for r in db.execute('SELECT payload FROM project_task_cards WHERE project_id=?',(current['project_id'],))]
        return dict(mission=current,learning=entries(app,sid),task_cards=cards,proposals=proposals,dossier=render(app,sid))
    if op=='propose':
        value=validate_mission(args.get('mission'));identifier=uuid.uuid4().hex[:24]
        with app.connect() as db:db.execute('INSERT INTO project_proposals VALUES(?,?,?,?,?)',(identifier,current['project_id'],current['revision'],encoded(value),'proposed'))
        render(app,sid)
        return dict(proposal_id=identifier,status='proposed',message='Review the concrete mission changes in Project direction; only a user action can apply them')
    if op=='feedback':
        kind=args.get('kind');quote=args.get('quote')
        if kind not in {'constraint','preference','correction','observation','hypothesis','next_step'}:raise ValueError('Invalid feedback kind')
        if not isinstance(quote,str) or not 1<=len(quote)<=1200 or quote not in request:raise ValueError('Quote exact user feedback from the current request (1..1200 characters)')
        with app.connect() as db:
            source=db.execute("SELECT id,text FROM messages WHERE session=? AND role='user' ORDER BY id DESC LIMIT 1",(sid,)).fetchone()
        if source is None or quote not in source['text']:raise ValueError('Feedback must come from an actual user message, not generated tool context')
        active=[r for r in entries(app,sid) if r['kind'] in {'constraint','preference','correction'} and not r.get('resolved')]
        if kind in {'constraint','preference','correction'} and len(active)>=30:raise ValueError('Review/archive existing feedback in Project direction before adding more pinned entries')
        payload=dict(text=quote,source_message_id=source['id'],source_session=sid,verification='User statement; not independently verified',resolved=False)
        if args.get('task_id'):payload['task_id']=app.task(sid,args['task_id'])['id']
        identifier=hashlib.sha256(encoded([current['project_id'],kind,quote,source['id'],payload.get('task_id')]).encode()).hexdigest()[:24]
        add_entry(app,sid,kind,payload,identifier);render(app,sid)
        return dict(entry_id=identifier,status='recorded',scope='Project-specific feedback. Typed execution rules change only through mission controls.')
    if op=='evidence':
        payload=evidence_value(app,sid,args.get('task_id'),args.get('pointer'))
        identifier=hashlib.sha256(encoded([current['project_id'],payload]).encode()).hexdigest()[:24]
        add_entry(app,sid,'evidence',payload,identifier);render(app,sid)
        return dict(entry_id=identifier,**payload)
    if op=='link':
        job=app.task(sid,args.get('task_id'));ids=args.get('requirement_ids');reason=args.get('contribution')
        if not isinstance(ids,list) or not ids or set(ids)-{r['id'] for r in current['requirements']}:raise ValueError('Use actual mission requirement IDs')
        if not isinstance(reason,str) or not 1<=len(reason)<=1000:raise ValueError('Describe the concrete contribution')
        with app.connect() as db:
            found=db.execute('SELECT payload FROM project_task_cards WHERE task_id=?',(job['id'],)).fetchone()
            card=json.loads(found['payload']) if found else dict(task_id=job['id'])
            if card.get('executed_revision'):raise ValueError('Executed task links are historical records; link a revised task instead')
            card.update(requirement_ids=ids,contribution=reason,link_source='agent_proposal',linked_revision=current['revision'])
            db.execute('INSERT OR REPLACE INTO project_task_cards VALUES(?,?,?)',(job['id'],current['project_id'],encoded(card)))
        render(app,sid);return card
    raise ValueError('Unsupported guidance operation; mission application/entry resolution are direct user controls')


def user_control(app,sid,args):
    with app.lock:return _user_control(app,sid,args)


def _user_control(app,sid,args):
    op=args.get('operation')
    if op=='apply':return dict(mission=revise(app,sid,args['mission'],args['base_revision'],'Direct user project control'))
    if op=='adopt':
        m=mission(app,sid)
        with app.connect() as db:r=db.execute("SELECT * FROM project_proposals WHERE id=? AND project_id=? AND status='proposed'",(args['proposal_id'],m['project_id'])).fetchone()
        if not r:raise ValueError('Unknown pending proposal')
        result=revise(app,sid,json.loads(r['payload']),r['base_revision'],'User adopted proposal '+r['id'])
        with app.connect() as db:db.execute("UPDATE project_proposals SET status='adopted' WHERE id=?",(r['id'],))
        return dict(mission=result)
    if op=='resolve_feedback':
        m=mission(app,sid)
        with app.connect() as db:
            r=db.execute('SELECT kind,payload FROM project_learning WHERE id=? AND project_id=?',(args['entry_id'],m['project_id'])).fetchone()
            if not r:raise ValueError('Unknown project entry')
            if r['kind'] not in {'constraint','preference','correction','observation','hypothesis','next_step'}:raise ValueError('Only user feedback can be archived through this control')
            value=json.loads(r['payload']);value.update(resolved=True,resolution='Archived by user; historical record retained')
            db.execute('UPDATE project_learning SET payload=? WHERE id=?',(encoded(value),args['entry_id']))
        add_entry(app,sid,'resolution',dict(entry_id=args['entry_id'],source='Direct user action'));render(app,sid);return dict(status='resolved')
    if op=='replicate':
        job=app.task(sid,args['task_id']);reason=args.get('reason')
        if job['status'] not in {'queued','blocked','failed','interrupted'}:raise ValueError('Select an unstarted/blocked task for independent execution')
        if not isinstance(reason,str) or not 1<=len(reason)<=1000:raise ValueError('Record why an independent calculation is requested')
        with app.connect() as db:db.execute('INSERT OR REPLACE INTO project_replicates VALUES(?,?)',(job['id'],reason))
        return dict(task_id=job['id'],status='replicate_authorized',message='Resume this task explicitly; existing scientific outputs remain sealed')
    if op=='accept_requirement':
        m=mission(app,sid);rid=args.get('requirement_id');reason=args.get('reason');ids=args.get('evidence_ids')
        if args.get('base_revision')!=m['revision']:raise ValueError('Mission changed; review current acceptance criteria')
        if rid not in {r['id'] for r in m['requirements']}:raise ValueError('Unknown requirement')
        if not isinstance(reason,str) or not 1<=len(reason)<=1000:raise ValueError('Record the acceptance rationale')
        recorded={e['id']:e for e in entries(app,sid) if e['kind']=='evidence'}
        if not isinstance(ids,list) or not ids or any(i not in recorded for i in ids):raise ValueError('Select project evidence IDs for acceptance')
        for identifier in ids:
            e=recorded[identifier];path=ensure_within(Path(e['report']),Path(app.context_for(sid)['project_root']))
            if hashlib.sha256(path.read_bytes()).hexdigest()!=e['report_sha256']:raise ValueError('Acceptance evidence changed; inspect the current report')
        identifier=add_entry(app,sid,'acceptance',dict(requirement_id=rid,mission_revision=m['revision'],evidence_ids=ids,
            reason=reason,verification='User acceptance of this requirement; not independent experimental validation'))
        render(app,sid);return dict(entry_id=identifier,status='user_accepted',mission_revision=m['revision'])
    raise ValueError('Unknown direct project control')


def task_event(app,jid,status):
    with app.connect() as db:
        job=db.execute('SELECT * FROM jobs WHERE id=?',(jid,)).fetchone()
        row=db.execute('SELECT payload FROM project_task_cards WHERE task_id=?',(jid,)).fetchone()
    if not job:return
    sid=job['session'];m=mission(app,sid)
    payload=dict(task_id=jid,status=status,error=(job['error'] or '')[:1600],report=job['report'],verification='Execution state, not biological validation')
    card=json.loads(row['payload']) if row else {}
    if card.get('reuse_task_id'):payload['reused_from_task']=card['reuse_task_id']
    # Stable signature makes repeated identical notifications idempotent.
    key=hashlib.sha256(encoded([m['project_id'],payload]).encode()).hexdigest()[:24]
    add_entry(app,sid,'task_event',payload,key)
    if status=='complete' and job['report']:
        try:
            root=Path(app.context_for(sid)['project_root']);path=ensure_within(Path(job['report']),root)
            raw=path.read_bytes();report=json.loads(raw)
            facts={}
            fields={'kind','status','prepared_molecules','failed_molecules','ranked_molecules','processed_conformers',
                'exported_molecules','full_coverage','shortfall','set_id','count','model','input_json','confidence_analysis','names',
                'pdb_id','report','report_sha256','isoform_assignment','polymer_entity_count','assembly_annotation_count'}
            for step_id,step in report.get('steps',{}).items():
                if step.get('status')!='complete':continue
                facts[step_id]={k:v for k,v in step.get('result',{}).items() if k in fields and isinstance(v,(str,int,float,bool))}
            observed=dict(task_id=jid,report=str(path),report_sha256=hashlib.sha256(raw).hexdigest(),pointer='/steps',value=facts,
                verification='Allowlisted fields from completed execution report; not independently validated science')
            if card.get('reuse_task_id'):observed['reused_from_task']=card['reuse_task_id']
            add_entry(app,sid,'evidence',observed,hashlib.sha256(encoded([m['project_id'],observed]).encode()).hexdigest()[:24])
        except (OSError,ValueError,KeyError):pass
    render(app,sid)


def render(app,sid):
    """Render from records. Never silently overwrite a manually edited dossier."""
    with app.lock:
        m=mission(app,sid);ctx=app.context_for(sid);path=Path(ctx['project_root'])/'PROJECT.md'
        with app.connect() as db:
            previous=db.execute('SELECT sha256 FROM project_dossier WHERE project_id=?',(m['project_id'],)).fetchone()
            cards=[json.loads(r['payload']) for r in db.execute('SELECT payload FROM project_task_cards WHERE project_id=?',(m['project_id'],))]
        if path.exists() and (not previous or hashlib.sha256(path.read_bytes()).hexdigest()!=previous['sha256']):
            return dict(path=str(path),status='manual_edit_conflict',message='Existing PROJECT.md preserved. Review edits against canonical project records; automatic overwrite stopped.')
        lines=['# Project direction and work ledger','',f"Project: {m['project_id']} | Mission revision: {m['revision']}",'',
            '## Original mission','',m['baseline']['objective'],'', 'Original target: '+m['baseline']['target'],
            '', 'Original requirements and limits:', '```json',json.dumps(m['baseline'],indent=2),'```',
            '','## Current authorized mission','',m['objective'],
            '', 'Target: '+m['target'],'','Excluded work: '+(m['excluded_work'] or 'Not specified'),'',
            '## Requirements and definition of done','']
        for r in m['requirements']:lines.extend([f"- {r['id']}: {r['text']}",'  Acceptance: '+r['acceptance']])
        lines.extend(['','Requirement acceptance is not inferred from task completion. Review the linked evidence.','',
            '## Active constraints','',*['- '+v for v in m['constraints']],'','Typed execution rules:', '```json',json.dumps(m['rules'],indent=2),'```',
            '', '## Project feedback and evidence',''])
        with app.connect() as db:
            notes=[dict(r) for r in db.execute('SELECT id,kind,text,task_id FROM research_notes WHERE project_id=? ORDER BY created',(m['project_id'],))]
        for note in notes:lines.append('- Legacy project statement (not independently verified): '+encoded(note))
        for e in entries(app,sid):
            if e['kind']=='task_event':continue
            lines.append('- '+e['id']+' / '+e['kind']+' / '+('archived' if e.get('resolved') else 'recorded')+': '+encoded({k:v for k,v in e.items() if k not in {'id','kind','created'}}))
        lines.extend(['','## Work ledger',''])
        cardmap={c['task_id']:c for c in cards}
        for job in app.project_jobs(sid):
            c=cardmap.get(job['id'],{})
            lines.extend(['- '+job['id']+' / '+job['status']+': '+job['request'][:180].replace('\n',' '),
                '  Requirements: '+', '.join(c.get('requirement_ids',[])),
                '  Alignment/reuse: '+encoded({k:c[k] for k in ('alignment','reuse_task_id','signature','mission_revision') if k in c}),
                '  Error: '+str(job.get('error') or 'none'),
                '  Report: '+str(job.get('report') or 'not available')])
        from .research_workflows import rows
        lines.extend(['','## Dependency workflows and next work',''])
        for flow in rows(app,sid):
            lines.append('- '+flow['id']+' / '+flow['status']+': '+flow['objective'])
            for n in flow['nodes']:lines.append('  - '+n['id']+' / '+n['status']+' / dependencies: '+', '.join(n['depends_on'])+' / task: '+str(n.get('task_id','not submitted')))
            if flow.get('error'):lines.append('  Blocker: '+flow['error'])
        lines.extend(['','This view is generated from project records. User feedback is not independently verified evidence.',
            'Task completion does not imply mission acceptance or experimental validity.',''])
        raw='\n'.join(lines).encode('utf-8');temp=path.with_suffix('.md.tmp');temp.write_bytes(raw);temp.replace(path)
        with app.connect() as db:db.execute('INSERT OR REPLACE INTO project_dossier VALUES(?,?)',(m['project_id'],hashlib.sha256(raw).hexdigest()))
        return dict(path=str(path),status='current',revision=m['revision'])


def brief(app,sid,request='',budget=22000):
    """Pinned mission and user corrections first; task-specific evidence fills remaining space."""
    m=mission(app,sid);all_entries=entries(app,sid)
    active=[e for e in all_entries if e['kind'] in {'constraint','preference','correction'} and not e.get('resolved')]
    compact_m={k:m[k] for k in ('project_id','revision','objective','target','requirements','rules','constraints','excluded_work')}
    # Literal molecules stay in the mission file/executor; context carries identity hashes.
    if 'protein_sequences' in compact_m['rules']:
        compact_m['rules']=dict(compact_m['rules'],protein_sequences=[dict(length=len(s),sha256=hashlib.sha256(s.encode()).hexdigest()) for s in m['rules']['protein_sequences']])
    with app.connect() as db:
        legacy=[dict(r) for r in db.execute("SELECT id,kind,text,task_id FROM research_notes WHERE project_id=? AND kind IN ('constraint','decision') ORDER BY created",(m['project_id'],))]
    payload=dict(mission=compact_m,pinned_feedback=active,notes=legacy,evidence=[],recent_tasks=[],next_work=[],omitted_entries=0,context_characters=0,
        accepted_requirements=list({e['requirement_id']:dict(requirement_id=e['requirement_id'],entry_id=e['id']) for e in all_entries
            if e['kind']=='acceptance' and e['mission_revision']==m['revision']}.values()),
        use_of_experience='Use completed artifacts and sourced feedback to choose the next relevant step. Reused tasks are not independent evidence. Failures diagnose execution, not scientific invalidity. Conflicting feedback remains unresolved until clarified; never silently rewrite the mission.',
        source='Canonical project records, not a recursively rewritten summary',project_file=str(Path(app.context_for(sid)['project_root'])/'PROJECT.md'))
    if len(encoded(payload))>budget:raise ValueError('Pinned project constraints exceed context budget. Review/archive resolved feedback; constraints were not silently dropped.')
    ids=set(__import__('re').findall(r'\b[a-f0-9]{16,32}\b',request))
    candidates=[e for e in all_entries if e not in active and not e.get('resolved') and
        (e['kind']!='task_event' or e.get('status') in {'complete','blocked','failed','interrupted'})]
    candidates.sort(key=lambda e:(e.get('task_id') in ids or e['id'] in ids,e['kind'] in {'evidence','acceptance','next_step'},e['created']),reverse=True)
    for e in candidates:
        entry=dict(e)
        if e['kind']=='evidence':
            path=Path(e['report'])
            try:unchanged=hashlib.sha256(path.read_bytes()).hexdigest()==e['report_sha256']
            except OSError:unchanged=False
            entry['freshness']='unchanged' if unchanged else 'stale_or_unavailable'
        if len(payload['evidence'])>=12 or len(encoded(payload))+len(encoded(entry))+10>budget-3500:payload['omitted_entries']+=1;continue
        payload['evidence'].append(entry)
    jobs=app.project_jobs(sid)
    jobs.sort(key=lambda j:(j['id'] in ids,j['created']),reverse=True)
    with app.connect() as db:
        cards={r['task_id']:json.loads(r['payload']) for r in db.execute('SELECT task_id,payload FROM project_task_cards WHERE project_id=?',(m['project_id'],))}
    for j in jobs[:6]:
        item={k:j[k] for k in ('id','status','error','report')};item['request']=j['request'][:250]
        c=cards.get(j['id'],{})
        item.update({k:c[k] for k in ('reuse_task_id','requirement_ids','mission_revision') if k in c})
        if len(encoded(payload))+len(encoded(item))+10<=budget:payload['recent_tasks'].append(item)
    from .research_workflows import rows
    ready=[]
    for flow in sorted(rows(app,sid),key=lambda f:f['created'],reverse=True):
        if flow['status'] in {'complete','cancelled'}:continue
        for node in flow['nodes']:
            if node['status']=='pending' and all(next(n for n in flow['nodes'] if n['id']==d)['status']=='complete' for d in node['depends_on']):
                ready.append(dict(workflow_id=flow['id'],node_id=node['id'],tool=node['tool'],workflow_status=flow['status'],blocker=flow.get('error')));break
        if len(ready)>=3:break
    if len(encoded(payload))+len(encoded(ready))+10<=budget:payload['next_work']=ready
    from .literature_reviews import brief as literature_brief
    literature=literature_brief(app,sid,request)
    if literature['reviews'] and len(encoded(payload))+len(encoded(literature))+40<=budget:
        payload['literature_reviews']=literature
    # Account for the length field itself without evicting mandatory constraints.
    for _ in range(3):payload['context_characters']=len(encoded(payload))
    if len(encoded(payload))>budget:raise ValueError('Pinned project constraints exceed context budget; constraints were not silently dropped.')
    return payload
