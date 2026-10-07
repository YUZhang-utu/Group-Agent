"""Owned adoption and dispatch of existing scores to block-restricted 3D search."""
import hashlib
import json
from pathlib import Path
import re
import shutil
import time
import uuid

from .final_work_blocks import read, sha
from .joint_spatial_profiles import save
from .project_context import ensure_within

OPERATIONS = {
    'list': (set(), set()),
    'comparison': (set(), set()),
    'discover': (set(), set()),
    'retry_config': ({'task_id'}, set()),
    'import': (set(), {'report','candidate_id'}),
    'import_query': (set(), {'report','candidate_id'}),
    'ranks': ({'task_id'}, {'scheme', 'receptor', 'limit', 'score'}),
    'top': ({'task_id', 'scheme', 'receptor', 'block_id'}, {'limit', 'score'}),
    'preview': ({'task_id', 'scheme', 'receptor'}, {'method', 'blocks'}),
    'start': ({'task_id', 'query_task_id', 'scheme', 'receptor'},
              {'method', 'blocks', 'conformers', 'dock', 'search_policy'}),
    'compare': ({'task_id', 'query_task_id', 'receptor'}, {'conformers'}),
}


def validate(args):
    if not isinstance(args, dict):
        raise ValueError('Block campaign arguments must be an object')
    op = args.get('operation', 'list')
    if not isinstance(op, str) or op not in OPERATIONS:
        raise ValueError('Unknown block campaign operation')
    required, optional = OPERATIONS[op]
    if not required <= args.keys() or set(args) - required - optional - {'operation'}:
        raise ValueError('Invalid ' + op + ' arguments; required: ' + ', '.join(sorted(required)) +
                         '; optional: ' + ', '.join(sorted(optional)))
    if op in {'import','import_query'}:
        if ('report' in args)+('candidate_id' in args) != 1:
            raise ValueError('Supply either the exact report path or a discovered candidate_id')
        if 'report' in args and (not isinstance(args['report'],str) or not args['report'].strip()):
            raise ValueError('Supply a report path')
        if 'candidate_id' in args and (not isinstance(args['candidate_id'],str) or not re.fullmatch('[a-f0-9]{24}',args['candidate_id'])):
            raise ValueError('Use a discovered candidate_id')
    for key in ('task_id', 'query_task_id'):
        if key in args and (not isinstance(args[key], str) or not re.fullmatch('[a-f0-9]{16}', args[key])):
            raise ValueError('Use an existing conversation task ID for ' + key)
    if args.get('method', 'union') not in {'union', 'intersection', 'equiscore', 'chemplp'}:
        raise ValueError('Choose union, intersection, equiscore or chemplp')
    if type(args.get('blocks', 5)) is not int or args.get('blocks', 5) not in (5, 10):
        raise ValueError('Choose 5 or 10 leading blocks per scoring method')
    for key, default, hi in [('conformers', 100000, 1000000), ('limit', 5, 100)]:
        if type(args.get(key, default)) is not int or not 1 <= args.get(key, default) <= hi:
            raise ValueError('Invalid ' + key)
    if type(args.get('dock', True)) is not bool:
        raise ValueError('dock must be boolean')
    if args.get('search_policy', 'all_selected') not in {'all_selected', 'ranked_blocks_until_budget'}:
        raise ValueError('Unknown search_policy')
    if args.get('search_policy') == 'ranked_blocks_until_budget' and args.get('method', 'union') not in {'equiscore','chemplp'}:
        raise ValueError('Sequential block search requires one scoring method')
    if args.get('score', 'EquiScore') not in {'EquiScore', 'ChemPLP'}:
        raise ValueError('Choose EquiScore or ChemPLP')
    if 'scheme' in args and args['scheme'] not in {'E094', 'E095', 'E096'}:
        raise ValueError('Choose E094, E095 or E096')
    if 'receptor' in args and (not isinstance(args['receptor'], str) or not re.fullmatch('[A-Za-z0-9_-]{1,64}', args['receptor'])):
        raise ValueError('Invalid receptor ID')
    return op


def owned_report(app, sid, task_id, project, kinds):
    job = app.task(sid, task_id)
    if job['status'] != 'complete':
        raise ValueError('A completed task in this conversation is required')
    parent = ensure_within(Path(job['report']), project)
    data = read(parent)
    paths = [parent] if data.get('kind') in kinds else []
    for step in data.get('steps', {}).values():
        result = step.get('result', {})
        if step.get('status') == 'complete' and result.get('report'):
            p = ensure_within(Path(result['report']), parent.parent)
            if read(p).get('kind') in kinds:
                paths.append(p)
    if len(paths) != 1:
        raise ValueError('Choose a task containing exactly one compatible completed result')
    return paths[0], job


def analysis_header(path, project):
    path = ensure_within(Path(path).resolve(), project)
    r = read(path)
    ranking = r.get('ranking', {})
    if (path.name != 'report.json' or r.get('kind') != 'block_analyze' or r.get('status') != 'complete'
            or ranking.get('score') != 'EquiScore' or ranking.get('better') != 'higher'
            or ranking.get('ranking_unit') != 'molecule' or ranking.get('top_n') != 10):
        raise ValueError('Select the completed full EquiScore analysis report (Top-10 distinct molecules)')
    ensure_within(Path(r['source_report']), project)
    return path, r


def select_blocks(path, scheme, receptor, method='union', blocks=5):
    from .block_ranking import inspect
    rows = {}
    for score, report in [('EquiScore', path), ('ChemPLP', path.parent / 'chemplp_baseline/report.json')]:
        inspect(report, operation='ranks', scheme=scheme, receptor=receptor, top_n=10, limit=blocks)
        # Read the verified database directly, preserving full rank provenance.
        from .block_sampling import readonly
        with readonly(report.parent / 'block_ranking.sqlite') as db:
            db.row_factory = __import__('sqlite3').Row
            rows[score] = [dict(r) for r in db.execute(
                'SELECT * FROM rankings WHERE scheme=? AND receptor=? AND top_n=10 AND status=? AND rank<=? ORDER BY rank,block_id',
                (scheme, receptor, 'eligible', blocks))]
    eq = {r['block_id'] for r in rows['EquiScore']}
    ch = {r['block_id'] for r in rows['ChemPLP']}
    chosen = {'union': eq | ch, 'intersection': eq & ch, 'equiscore': eq, 'chemplp': ch}[method]
    populations = {r['block_id']:r['population'] for values in rows.values() for r in values}
    return dict(scheme=scheme, receptor=receptor, method=method, blocks_per_method=blocks,
                block_ids=sorted(chosen), selected_blocks=len(chosen), rankings=rows,
                block_population_conformers=sum(populations[b] for b in chosen),
                block_statistic='Mean of best 10 distinct molecules; best in-block conformer',
                candidate_unit='conformer', multiple_conformers_per_molecule=True)


def queue(app, sid, project, request, action):
    from .prompt_workflow import create_plan
    ctx = app.context
    key = hashlib.sha256(json.dumps(request, sort_keys=True).encode()).hexdigest()[:24]
    root = project / 'block-campaigns'
    root.mkdir(exist_ok=True)
    entry = root / (key + '.json')
    if not entry.exists():
        save(entry, request)
    title = ('Import saved multi-cocrystal query ' if action == 'block_adopt_query' else
             'Import saved EquiScore analysis ' if action == 'block_adopt_scores' else
             'Search leading blocks and dock conformers ' if request.get('dock') else
             'Search leading blocks for conformers ') + key
    matches = [j for j in app.jobs(sid) if j['request'] == title]
    if matches:
        job = matches[-1]
        return dict(status=job['status'], task_id=job['id'], reused=True,
                    message='Existing receipt retained. Resume this task explicitly if interrupted; no duplicate submitted.')
    plan = dict(version=1, summary=title, clarifications=[], steps=[dict(id='campaign', action=action,
                params=dict(request_id=key, request_sha256=sha(entry)))])
    path = create_plan(Path(ctx['db']), ctx['user_id'], ctx['project_id'], title, local_plan=plan)
    jid = uuid.uuid4().hex[:16]
    with app.connect() as db:
        db.execute('INSERT INTO jobs(id,session,request,provider,profile,status,plan,report,log,created) VALUES(?,?,?,?,?,?,?,?,?,?)',
                   (jid, sid, title, 'local', '', 'queued', str(path), str(path.parent/'execution/report.json'),
                    str(path.parent/'execution.log'), time.time()))
    return dict(status='queued', task_id=jid, reused=False)


def retry_config(app, sid, project, task_id):
    """Create one fresh plan for an explicit retry before scientific work began."""
    job = app.task(sid, task_id)
    if job['status'] != 'blocked':
        raise ValueError('retry_config requires a configuration-blocked campaign')
    plan_path = ensure_within(Path(job['plan']).resolve(), project)
    report_path = ensure_within(Path(job['report']).resolve(), plan_path.parent/'execution')
    if read(plan_path.parent/'plan-seal.json')['plan_sha256'] != sha(plan_path):
        raise ValueError('Original plan changed')
    steps = read(plan_path)['plan']['steps']
    if len(steps) != 1 or steps[0]['action'] != 'block_search_dock' or steps[0]['id'] != 'campaign':
        raise ValueError('retry_config requires a block-search campaign')
    report = read(report_path)
    step = report.get('steps',{}).get('campaign',{})
    error = step.get('error','')
    if (report.get('status') != 'blocked' or step.get('status') != 'blocked' or
        not (error == 'Configure search.batch, block_evaluation.sampling_profile and plants_profile' or
             error.startswith('Missing block-search runtime configuration: '))):
        raise ValueError('Only missing runtime configuration can use retry_config')
    directory = report_path.parent/'campaign'
    # execute_step creates blocks/ before run() checks runtime configuration.
    # Only that known empty directory is scaffolding, never a scientific receipt.
    entries = list(directory.iterdir()) if directory.exists() else []
    artifacts = [p for p in entries if not (
        p.name == 'blocks' and not p.is_symlink() and p.is_dir() and not any(p.iterdir()))]
    if artifacts or (report_path.parent/'campaign.stage.json').exists():
        raise ValueError('Scientific artifacts already exist; inspect the original task instead')
    params = steps[0]['params']
    entry = ensure_within(project/'block-campaigns'/(params['request_id']+'.json'), project)
    if sha(entry) != params['request_sha256']:
        raise ValueError('Original campaign request changed')
    request = read(entry)
    if request.get('session') != sid:
        raise ValueError('Campaign belongs to another conversation')
    request['retry_of'] = task_id
    return queue(app, sid, project, request, 'block_search_dock')


def handle(app, sid, args):
    from .prompt_workflow import project_root
    op = validate(args)
    app.session(sid)
    ctx = app.context
    project = project_root(Path(ctx['db']), ctx['user_id'], ctx['project_id'])
    if op == 'comparison':
        rows = []
        for job in app.jobs(sid):
            if not job['request'].startswith('Search leading blocks for conformers '):
                continue
            plan = ensure_within(Path(job['plan']), project)
            steps = read(plan)['plan']['steps']
            if len(steps)!=1 or steps[0]['action']!='block_search_dock':
                continue
            params = steps[0]['params']
            entry = ensure_within(project/'block-campaigns'/(params['request_id']+'.json'),project)
            if sha(entry)!=params['request_sha256']:
                raise ValueError('Campaign request changed')
            request = read(entry)
            if request.get('session')!=sid or request.get('search_policy')!='ranked_blocks_until_budget':
                continue
            selection = request['selection']
            row = dict(task_id=job['id'],status=job['status'],scheme=selection['scheme'],
                       method=selection['method'],blocks=selection['blocks_per_method'],receptor=selection['receptor'],
                       query_sha256=request['query_sha256'],analysis_sha256=request['analysis_sha256'],
                       requested_conformers=request['conformers'],report=job.get('report'))
            if job['status']=='complete':
                path,_ = owned_report(app,sid,job['id'],project,{'block_search_dock'})
                report = read(path)
                row.update({k:report[k] for k in ('searched_conformers','ranked_conformers','exported_conformers',
                    'unique_molecules','shortfall','stop_reason','elapsed_seconds','resumed_execution',
                    'mol2_directory','searched_blocks') if k in report})
            rows.append(row)
        return dict(arms=rows,docking_requested=False,
                    limitations='Compare matching query/analysis hashes and budgets only. Resumed times are not fresh latency benchmarks; no affinity or enrichment measured.')
    if op == 'retry_config':
        with app.lock:
            return retry_config(app, sid, project, args['task_id'])
    if op == 'discover':
        from .block_saved_inputs import discover
        return discover(app,sid,project)
    if op == 'list':
        tasks = app.jobs(sid)
        queries = []
        for job in tasks:
            if job['status'] != 'complete':
                continue
            try:
                path, _ = owned_report(app, sid, job['id'], project, {'consensus_design', 'consensus_recommendation'})
                r = read(path)
                queries.append(dict(task_id=job['id'], kind=r['kind'], target=r.get('target'), report=str(path)))
            except (ValueError, OSError, KeyError):
                continue
        return dict(tasks=[{k:j.get(k) for k in ('id','status','request','report','error')}
                           for j in tasks if 'leading blocks' in j['request'] or j['request'].startswith('Import saved ')],
                    query_tasks=queries, candidate_unit='conformer', default_conformers=100000,
                    recovery='If inputs are missing here, call discover to locate original CLI analysis and multi-cocrystal query packages before asking for a task ID.')
    if op == 'import_query':
        from .block_saved_inputs import resolve, query_header
        path, _ = query_header(resolve(args,sid,project,'query'))
        with app.lock:
            return queue(app,sid,project,dict(session=sid,saved_query=str(path),saved_query_sha256=sha(path)), 'block_adopt_query')
    if op == 'import':
        from .block_saved_inputs import resolve
        path, _ = analysis_header(resolve(args,sid,project,'scores'), project)
        with app.lock:
            return queue(app, sid, project, dict(session=sid, report=str(path), report_sha256=sha(path)), 'block_adopt_scores')
    path, _ = owned_report(app, sid, args['task_id'], project, {'block_analyze'})
    analysis_header(path, project)
    if op == 'compare':
        # Validate all input selections before queuing any of the twelve arms.
        owned_report(app, sid, args['query_task_id'], project, {'consensus_design', 'consensus_recommendation'})
        arms = [dict(operation='start', task_id=args['task_id'], query_task_id=args['query_task_id'],
                     scheme=scheme, receptor=args['receptor'], method=method, blocks=blocks,
                     conformers=args.get('conformers',100000), dock=False,
                     search_policy='ranked_blocks_until_budget')
                for scheme in ('E094','E095','E096') for method in ('chemplp','equiscore') for blocks in (5,10)]
        for arm in arms:
            if not select_blocks(path,arm['scheme'],arm['receptor'],arm['method'],arm['blocks'])['block_ids']:
                raise ValueError('Missing eligible blocks for '+arm['scheme']+'/'+arm['method'])
        tasks = [dict(scheme=a['scheme'],method=a['method'],blocks=a['blocks'],**handle(app,sid,a)) for a in arms]
        return dict(status='queued' if any(t['status'] in {'queued','running'} for t in tasks) else 'recorded',
                    tasks=tasks, docking_requested=False, candidate_unit='conformer',
                    message='Independent arms; inspect existing task IDs. Top-5 and Top-10 may stop at the same prefix.')
    if op in {'ranks', 'top'}:
        from .block_ranking import inspect
        if args.get('score') == 'ChemPLP':
            path = path.parent/'chemplp_baseline/report.json'
        return inspect(path, operation=op, top_n=10, **{k:v for k,v in args.items() if k in {'scheme','receptor','block_id','limit'}})
    selection = select_blocks(path, args['scheme'], args['receptor'], args.get('method','union'), args.get('blocks',5))
    if op == 'preview':
        return selection
    if not selection['block_ids']:
        raise ValueError('No blocks selected; inspect the model rankings or choose union')
    query, job = owned_report(app, sid, args['query_task_id'], project, {'consensus_design', 'consensus_recommendation'})
    run = Path(job['plan']).parent.name
    request = dict(session=sid, analysis=str(path), analysis_sha256=sha(path), selection=selection,
                   query=str(query), query_sha256=sha(query), query_run=run,
                   conformers=args.get('conformers',100000), dock=args.get('dock',True))
    if 'search_policy' in args:
        request['search_policy'] = args['search_policy']
    with app.lock:
        return queue(app, sid, project, request, 'block_search_dock')


def execute(action, params, output, execution, cfg, allow_compute):
    project = execution.parent.parent.parent
    entry = ensure_within(project/'block-campaigns'/(params['request_id']+'.json'), project)
    if sha(entry) != params['request_sha256']:
        raise ValueError('Campaign request changed')
    request = read(entry)
    output.mkdir(parents=True, exist_ok=True)
    if action == 'block_adopt_scores':
        return adopt(request, output, project)
    if action == 'block_adopt_query':
        from .block_saved_inputs import adopt_query
        return adopt_query(request,output)
    from .prompt_workflow import Blocked
    if not allow_compute:
        raise Blocked('Enable --allow-compute for the requested block search and docking')
    from .screening_selection import upstream
    query, _ = upstream(execution, request['query_run'], ('consensus_design','consensus_recommend','block_adopt_query'))
    if query.resolve() != Path(request['query']).resolve() or sha(query) != request['query_sha256']:
        raise ValueError('Selected query receipt changed')
    path, _ = analysis_header(request['analysis'], project)
    if sha(path) != request['analysis_sha256']:
        raise ValueError('Selected analysis changed')
    selection = request['selection']
    if select_blocks(path, selection['scheme'], selection['receptor'], selection['method'], selection['blocks_per_method']) != selection:
        raise ValueError('Block selection changed')
    from .block_conformer_search import run
    return run(request, output, cfg)


def adopt(request, output, project):
    from .block_sampling import check_outputs
    path, report = analysis_header(request['report'], project)
    if sha(path) != request['report_sha256']:
        raise ValueError('Analysis report changed after submission')
    print('Verifying saved EquiScore analysis and full-run receipts; no inference or docking.', flush=True)
    check_outputs(path.parent)
    source = ensure_within(Path(report['source_report']), project)
    full = check_outputs(source.parent)
    if (full.get('kind') != 'block_equiscore_run' or full.get('scope') != 'full' or
        full.get('model_executed') is not True or
        full.get('pairs',0) <= 0 or full.get('pairs') != full.get('scored') or full.get('failed_pairs') != 0 or
        full.get('worker_exit_code') != 0 or report['ranking']['source_report_sha256'] != sha(source)):
        raise ValueError('Complete, failure-free full EquiScore run and matching analysis required')
    plants = ensure_within(Path(full['source_report']), project)
    if sha(plants) != full['identity']['source_sha256']:
        raise ValueError('Original PLANTS report changed')
    check_outputs(path.parent/'chemplp_baseline')
    baseline = read(path.parent/'chemplp_baseline/report.json')
    if baseline['ranking']['source_report_sha256'] != sha(plants):
        raise ValueError('ChemPLP baseline has a different source')
    if (output/'report.json').exists():
        return check_outputs(output)
    # Copy compact inspection artifacts, keeping the original full panel immutable.
    adopted = None
    for parent, dest, original in [(path.parent, output, report),
                                  (path.parent/'chemplp_baseline', output/'chemplp_baseline', baseline)]:
        dest.mkdir(exist_ok=True)
        hashes = {}
        for name in ('block_ranking.sqlite','block_rankings.csv','block_top_candidates.csv','rank_comparison.csv','ranking_review.md'):
            if name in original['output_hashes']:
                shutil.copyfile(parent/name, dest/name)
                hashes[name] = sha(dest/name)
        copied = dict(original, output_hashes=hashes,
                      outputs={name:str(dest/name) for name in hashes},
                      adoption=dict(source=str(parent/'report.json'), sha256=sha(parent/'report.json'),
                                    full_pairs=full['pairs'], failed_pairs=0, computation_launched=False))
        if dest == output:
            adopted = copied
        else:
            save(dest/'report.json', copied)
    adopted['output_hashes'].update({p.relative_to(output).as_posix():sha(p) for p in (output/'chemplp_baseline').iterdir() if p.is_file()})
    # Publish completion only after both scoring views are available.
    save(output/'report.json', adopted)
    return adopted
