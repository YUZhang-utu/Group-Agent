"""Owned project workspaces, cross-conversation evidence and durable notes."""
import json
from pathlib import Path
import time
import uuid

from .registry import connect
from .project_context import create_scientific_project, list_owned_projects, require_project_owner


def initialize(app):
    with app.connect() as db:
        columns = {row['name'] for row in db.execute('PRAGMA table_info(sessions)')}
        if 'project_id' not in columns: db.execute('ALTER TABLE sessions ADD COLUMN project_id TEXT')
        db.execute('UPDATE sessions SET project_id=? WHERE project_id IS NULL', (app.context['project_id'],))
        db.executescript('CREATE TABLE IF NOT EXISTS research_project_metadata(project_id TEXT PRIMARY KEY,payload TEXT); '
            'CREATE TABLE IF NOT EXISTS research_notes(id TEXT PRIMARY KEY,project_id TEXT,kind TEXT,text TEXT,task_id TEXT,created REAL); '
            'CREATE TABLE IF NOT EXISTS research_events(id INTEGER PRIMARY KEY,project_id TEXT,task_id TEXT,status TEXT,created REAL);')


def context(app, sid):
    app.session(sid)
    with app.connect() as db:
        row = db.execute('SELECT project_id FROM sessions WHERE id=?', (sid,)).fetchone()
    with connect(Path(app.context['db'])) as db:
        project = require_project_owner(db, app.context['user_id'], row['project_id'])
    return dict(app.context, project_id=row['project_id'], project_root=project['workspace_path'])


def projects(app):
    with connect(Path(app.context['db'])) as db:
        rows = list_owned_projects(db, app.context['user_id'])
    with app.connect() as db:
        metadata = {r['project_id']: json.loads(r['payload']) for r in db.execute('SELECT * FROM research_project_metadata')}
    return [dict(row, **metadata.get(row['project_id'], {})) for row in rows]


def project_jobs(app, sid):
    pid = context(app, sid)['project_id']
    with app.connect() as db:
        return [dict(r) for r in db.execute('SELECT j.* FROM jobs j JOIN sessions s ON j.session=s.id WHERE s.project_id=? ORDER BY j.created,j.id', (pid,))]


def remember_status(app, jid, status):
    with app.connect() as db:
        row = db.execute('SELECT s.project_id FROM jobs j JOIN sessions s ON j.session=s.id WHERE j.id=?', (jid,)).fetchone()
        if row:
            previous = db.execute('SELECT status FROM research_events WHERE task_id=? ORDER BY id DESC LIMIT 1', (jid,)).fetchone()
            if not previous or previous['status'] != status:
                db.execute('INSERT INTO research_events(project_id,task_id,status,created) VALUES(?,?,?,?)',
                           (row['project_id'], jid, status, time.time()))


def memory(app, sid, limit=20):
    ctx = context(app, sid); pid = ctx['project_id']
    project = next(r for r in projects(app) if r['project_id'] == pid)
    jobs = project_jobs(app, sid); counts = {}
    tasks = []
    for job in jobs:
        counts[job['status']] = counts.get(job['status'], 0) + 1
    for job in jobs[-limit:]:
        row = {k:job[k] for k in ('id','session','status','report','error','created')}
        row['request'] = job['request'][:350]
        row['actions'] = []
        if job.get('report'):
            from .project_context import ensure_within
            try:
                report = json.loads(ensure_within(Path(job['report']), Path(ctx['project_root'])).read_text(encoding='utf-8'))
                row['actions'] = [dict(id=k,action=v.get('action'),status=v.get('status'),
                    result={key:value for key,value in v.get('result',{}).items() if key in {
                        'kind','report','input_json','model','confidence_analysis','set_id','count','prepared_molecules',
                        'failed_molecules','ranked_molecules','exported_molecules','names','ranking','original_mol2_directory'}})
                    for k,v in report.get('steps',{}).items()]
            except (OSError, ValueError): pass
        tasks.append(row)
    with app.connect() as db:
        notes = [dict(r) for r in db.execute('SELECT id,kind,text,task_id,created FROM research_notes WHERE project_id=? ORDER BY created DESC LIMIT ?', (pid, limit))]
    from .research_workflows import rows
    from .candidate_sets import handle as candidates
    workflows=[dict(id=r['id'],objective=r['objective'],status=r['status'],error=r.get('error'),
        nodes=[{k:n[k] for k in ('id','status','task_id') if k in n} for n in r['nodes']]) for r in rows(app,sid)]
    return dict(project=project, total_tasks=len(jobs), status_counts=counts, recent_tasks=tasks,
                workflows=workflows,candidate_sets=candidates(app,sid,dict(operation='list'))['sets'],
                notes=notes, scope='This owned project across conversations; notes are recorded statements, task statuses are execution evidence')


def handle(app, sid, args):
    operation = args.get('operation','memory')
    if operation == 'list': return dict(projects=projects(app))
    ctx = context(app, sid)
    if operation == 'memory': return memory(app, sid)
    if operation == 'search':
        query=args.get('query');limit=args.get('limit',20);offset=args.get('offset',0)
        if not isinstance(query,str) or not 1<=len(query)<=200 or type(limit) is not int or not 1<=limit<=100 or type(offset) is not int or offset<0:
            raise ValueError('Search requires query, optional offset >= 0 and limit 1..100')
        pattern='%'+query.replace('\\','\\\\').replace('%','\\%').replace('_','\\_')+'%'
        with app.connect() as db:
            hits=[dict(r) for r in db.execute("""SELECT 'task' AS kind,j.id,j.request AS text,j.status,j.created
                FROM jobs j JOIN sessions s ON j.session=s.id WHERE s.project_id=? AND j.request LIKE ? ESCAPE '\\'
                UNION ALL SELECT n.kind,n.id,n.text,'recorded_statement',n.created FROM research_notes n
                WHERE n.project_id=? AND n.text LIKE ? ESCAPE '\\'
                UNION ALL SELECT 'message',CAST(m.id AS TEXT),m.text,m.role,m.created FROM messages m
                JOIN sessions s ON m.session=s.id WHERE s.project_id=? AND m.text LIKE ? ESCAPE '\\'
                ORDER BY created DESC LIMIT ? OFFSET ?""",(ctx['project_id'],pattern)*3+(limit+1,offset))]
        return dict(matches=[dict(r,text=r['text'][:4000],truncated=len(r['text'])>4000) for r in hits[:limit]],
                    offset=offset,has_more=len(hits)>limit,scope='Project tasks, notes and conversation records; statements are not verified findings')
    if operation in {'create','open'}:
        if operation == 'create':
            allowed = {'operation','name','objective','target','organism'}
            if set(args)-allowed: raise ValueError('Unknown project fields')
            for key in ('name','objective','target'):
                if not isinstance(args.get(key),str) or not 1 <= len(args[key].strip()) <= 2000:
                    raise ValueError('Project name, objective and target are required')
            with connect(Path(ctx['db'])) as db:
                existing = db.execute('SELECT id FROM project WHERE owner_user_id=? AND name=?', (ctx['user_id'],args['name'].strip())).fetchone()
                if existing: raise ValueError('A project already has this name; use open to continue it')
                pid = create_scientific_project(db, ctx['user_id'], args['name'], args['objective'], app.root.parent)
            with app.connect() as db:
                db.execute('INSERT INTO research_project_metadata VALUES(?,?)', (pid, json.dumps(dict(target=args['target'],organism=args.get('organism','unspecified')))))
        else:
            pid = args.get('project_id')
            with connect(Path(ctx['db'])) as db: require_project_owner(db,ctx['user_id'],pid)
        # Never reassign a conversation with sealed tasks or uploaded assets.
        new_sid = app.new_session(project_id=pid)
        from .project_guidance import render
        render(app,new_sid)
        return dict(project_id=pid, session=new_sid, message='Opened a conversation in the selected project. Existing conversations and tasks keep their original project.')
    if operation == 'note':
        if set(args)-{'operation','kind','text','task_id'}: raise ValueError('Unknown project note fields')
        if args.get('kind') not in {'decision','constraint','finding','next_step'}: raise ValueError('Use decision, constraint, finding or next_step')
        if not isinstance(args.get('text'),str) or not 1 <= len(args['text']) <= 4000: raise ValueError('Note text must be 1..4000 characters')
        if args.get('task_id'): app.task(sid,args['task_id'])
        identifier = uuid.uuid4().hex
        with app.connect() as db:
            db.execute('INSERT INTO research_notes VALUES(?,?,?,?,?,?)', (identifier,ctx['project_id'],args['kind'],args['text'],args.get('task_id'),time.time()))
        from .project_guidance import render
        render(app,sid)
        return dict(note_id=identifier, project_id=ctx['project_id'], provenance='Recorded project statement; not independently verified scientific fact')
    if operation == 'assets':
        rows=[]
        from .chat_uploads import attachment
        for path in sorted((Path(ctx['project_root'])/'uploads').glob('*/*/attachment.json')):
            record=json.loads(path.read_text(encoding='utf-8'))
            try: rows.append(attachment(app,sid,record['id']))
            except ValueError: continue
        return dict(assets=rows[-100:], total=len(rows))
    raise ValueError('Unknown project operation')
