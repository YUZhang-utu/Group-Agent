"""Session-owned adoption of existing CLI docking reports, without engine dispatch."""
import hashlib
import csv
import json
from pathlib import Path
import re
import time
import uuid

from .project_context import ensure_within
from .final_work_blocks import sha
from .joint_spatial_profiles import save


def read_report(path, project):
    path = ensure_within(Path(path).resolve(), project)
    if path.name != 'report.json' or not path.is_file():
        raise ValueError('Select the existing PLANTS report.json inside the active Project')
    if path.stat().st_size > 512 * 1024 * 1024:
        raise ValueError('Report exceeds 512 MiB')
    raw = path.read_bytes()
    report = json.loads(raw)
    if report.get('kind') != 'block_plants_run' or report.get('status') not in {'complete', 'partial'}:
        raise ValueError('A final complete or partial block_plants_run report is required')
    for key in ('jobs', 'attempted_jobs', 'failed_jobs', 'scored'):
        if type(report.get(key)) is not int or report[key] < 0:
            raise ValueError('Invalid docking counter: ' + key)
    if not 0 <= report['failed_jobs'] <= report['attempted_jobs'] <= report['jobs'] or report['jobs'] == 0:
        raise ValueError('Inconsistent docking counters')
    if report['status'] == 'complete' and (report['attempted_jobs'] != report['jobs'] or report['failed_jobs'] or report.get('first_job_gate') != 'passed'):
        raise ValueError('Complete report has incomplete or failed jobs')
    hashes = report.get('output_hashes', {})
    if not isinstance(hashes, dict) or not {'scores.csv', 'signature.json'} <= hashes.keys():
        raise ValueError('Missing docking output seals')
    for name, digest in hashes.items():
        if Path(name).is_absolute() or '..' in Path(name).parts:
            raise ValueError('Output path escapes docking directory')
        if not isinstance(digest, str) or not re.fullmatch('[a-f0-9]{64}', digest):
            raise ValueError('Invalid output digest')
    for key in ('sample_report', 'prepared_report'):
        ensure_within(Path(report[key]).resolve(), project)
    return path, report, hashlib.sha256(raw).hexdigest()


def summary(report):
    return {k: report[k] for k in ('status', 'jobs', 'attempted_jobs', 'failed_jobs', 'scored',
                                  'first_job_gate', 'score', 'ligand_mode', 'limitations') if k in report}


def handle(app, sid, args):
    from .prompt_workflow import project_root, create_plan
    app.session(sid)
    ctx = app.context
    project = project_root(Path(ctx['db']), ctx['user_id'], ctx['project_id'])
    root = project / 'block-results'
    operation = args.get('operation', 'list')
    if operation in {'ranks', 'top'}:
        from .block_ranking import inspect
        job = app.task(sid, args['task_id'])
        if job['status'] != 'complete':
            return dict(status=job['status'], message='Wait for the analysis task to complete')
        report_path = ensure_within(Path(job['report']), project)
        parent = json.loads(report_path.read_text(encoding='utf-8'))
        reports = []
        if parent.get('kind') == 'block_analyze':
            reports.append(report_path)
        for step in parent.get('steps', {}).values():
            result = step.get('result', {})
            if result.get('status') == 'block_analyze' and result.get('report'):
                reports.append(ensure_within(Path(result['report']), report_path.parent))
        if len(reports) != 1:
            raise ValueError('Select a task containing exactly one completed block analysis')
        return inspect(reports[0], **{k: v for k, v in args.items() if k in
            {'operation', 'scheme', 'receptor', 'block_id', 'top_n', 'limit'}})
    if operation not in {'attach', 'list', 'analyze'}:
        raise ValueError('Choose attach, list, analyze, ranks or top')
    if operation == 'list':
        rows = [json.loads(p.read_text(encoding='utf-8')) for p in sorted(root.glob('*.json'))]
        return dict(attachments=[r for r in rows if r.get('session') == sid])
    with app.lock:
        root.mkdir(exist_ok=True)
        if operation == 'attach':
            path, report, digest = read_report(args['report'], project)
            identifier = hashlib.sha256((sid + str(path) + digest).encode()).hexdigest()[:24]
            dest = root / (identifier + '.json')
            if not dest.exists():
                save(dest, dict(id=identifier, session=sid, report=str(path), report_sha256=digest,
                                summary=summary(report), verification='Report counters only; all output seals are checked by analysis.'))
            return dict(status='attached', attachment=json.loads(dest.read_text(encoding='utf-8')), computation_launched=False)
        rows = [json.loads(p.read_text(encoding='utf-8')) for p in sorted(root.glob('*.json'))]
        rows = [r for r in rows if r.get('session') == sid]
        selected = [r for r in rows if r['id'] == args.get('attachment_id')]
        if len(selected) != 1:
            raise ValueError('Choose an attachment ID from this conversation')
        row = selected[0]
        if sha(Path(row['report'])) != row['report_sha256']:
            raise ValueError('Source report changed; attach its final version before analysis')
        from .block_ranking import options
        ranking_options = options(args.get('top_n', 10), args.get('ranking_unit', 'molecule'))
        title = 'Top-N v1 analyze attached docking ' + row['id'] + ' ' + json.dumps(ranking_options, sort_keys=True)
        existing = [j for j in app.jobs(sid) if j['request'] == title]
        if existing and existing[-1]['status'] in {'queued', 'running', 'complete'}:
            job = existing[-1]
            return dict(status=job['status'], tasks=[dict(id=job['id'], status=job['status'])], reused=True)
        plan = dict(version=1, summary=title, clarifications=[], steps=[dict(id='analysis',
                    action='block_import_analysis', params=dict(attachment_id=row['id'], attachment_sha256=sha(root / (row['id'] + '.json')), **ranking_options))])
        path = create_plan(Path(ctx['db']), ctx['user_id'], ctx['project_id'], title, local_plan=plan)
        jid = uuid.uuid4().hex[:16]
        with app.connect() as db:
            db.execute('INSERT INTO jobs(id,session,request,provider,profile,status,plan,report,log,created) VALUES(?,?,?,?,?,?,?,?,?,?)',
                       (jid, sid, title, 'local', '', 'queued', str(path), str(path.parent / 'execution/report.json'), str(path.parent / 'execution.log'), time.time()))
        return dict(status='queued', tasks=[dict(id=jid, status='queued')], engine_execution=False)


def analyze_attachment(project, params, output):
    from .block_plants import analyze
    from .prompt_workflow import file_lock
    path = ensure_within(project / 'block-results' / (params['attachment_id'] + '.json'), project)
    if sha(path) != params['attachment_sha256']:
        raise ValueError('Attachment changed')
    attachment = json.loads(path.read_text(encoding='utf-8'))
    source, report, digest = read_report(attachment['report'], project)
    if digest != attachment['report_sha256']:
        raise ValueError('Docking report changed')
    # The CLI uses the same adjacent lock. Never analyze an active scheduler.
    with file_lock(source.parent.parent / ('.' + source.parent.name + '.e097.lock')):
        print('Inspecting attached docking report and prepared/sample bindings; no engine will run.', flush=True)
        if sha(source) != digest:
            raise ValueError('Docking report changed while acquiring its lock')
        for name in report['output_hashes']:
            ensure_within(source.parent / name, source.parent)
        from .block_sampling import check_outputs
        prepared = Path(report['prepared_report'])
        signature = json.loads((source.parent / 'signature.json').read_text(encoding='utf-8'))
        if signature.get('prepared_report') != str(prepared) or signature.get('prepared_sha256') != sha(prepared):
            raise ValueError('Prepared report does not match the docking signature')
        prep = check_outputs(prepared.parent)
        if prep.get('jobs') != report['jobs']:
            raise ValueError('Docking job count differs from preparation')
        prep_signature = json.loads((prepared.parent / 'signature.json').read_text(encoding='utf-8'))
        sample = Path(report['sample_report'])
        if prep.get('sample_report') != str(sample) or prep_signature.get('sample_sha256') != sha(sample):
            raise ValueError('Sampling report does not match preparation')
        with (source.parent / 'scores.csv').open(encoding='utf-8', newline='') as stream:
            scored = 0
            incomplete = 0
            for row in csv.DictReader(stream):
                scored += row['status'] == 'ok'
                incomplete += row['status'] != 'ok'
        if scored != report['scored'] or (report['status'] == 'complete' and incomplete):
            raise ValueError('Score table does not match report counters')
        print(f"Checking {len(report['output_hashes']):,} original output seals and analyzing saved scores. Large panels may require substantial file I/O.", flush=True)
        result = analyze(source, output)
        from .block_ranking import add_rankings
        result = add_rankings(source, output, top_n=params.get('top_n', 10), ranking_unit=params.get('ranking_unit', 'molecule'))
        if sha(source) != digest:
            raise ValueError('Docking report changed during analysis')
        result.setdefault('outputs', {}).update({name: str(output / name) for name in ('block_summary.csv', 'block_scores.csv', 'review.md')})
        result['source_report'] = str(source)
        result['source_report_sha256'] = digest
        result['verification'] = 'Original docking output seals, preparation and sampling bindings checked; engine not invoked.'
        result['ligand_chemistry_reviewed'] = prep.get('ligand_chemistry_reviewed', False)
        save(output / 'report.json', result)
        print('Saved block analysis and pose handoff; original docking outputs preserved.', flush=True)
    return result
