"""Local, bounded presentation of saved PDB inspection evidence."""
import hashlib
import json
from pathlib import Path
import re
import sqlite3
import time
import uuid

from .project_context import ensure_within


def affirmative(request):
    text = request.casefold()
    text = re.sub(r"\b(?:do not|don't|never)\b[^.\n;]*", '', text)
    return re.sub(r'(?:\u4e0d\u8981|\u7981\u6b62)[^\u3002\n;]*', '', text)


def note_requested(request):
    text = affirmative(request)
    return bool(re.search(r'\b(?:record|save|add)\b[^.\n;]{0,100}\b(?:project note|project finding)\b|(?:\u4fdd\u5b58|\u8bb0\u5f55)[^\u3002\n;]{0,40}(?:\u9879\u76ee\u7b14\u8bb0|\u9879\u76ee\u8bb0\u5f55)', text))


def request_task(request):
    ids = set(re.findall(r'\b[a-f0-9]{16}\b', request, re.I))
    text = affirmative(request)
    if len(ids) != 1 or not re.search(r'\b(?:read|summari[sz]e|review)\b|\u8bfb\u53d6|\u603b\u7ed3|\u67e5\u770b', text):
        return None
    if not re.search(r'\b(?:inspection|pdb|report)\b|\u62a5\u544a|\u7ed3\u6784\u68c0\u67e5', text):
        return None
    # This narrow route must not swallow new work or unrelated knowledge requests.
    if re.search(r'\b(?:run|launch|predict|dock|screen|execute|submit|resume|cancel|delete|archive|change|compare|literature|paper|search|download|fetch|create|prepare|export|analyze|af3|qc)\b|\u8fd0\u884c|\u542f\u52a8|\u9884\u6d4b|\u5bf9\u63a5|\u5220\u9664|\u6587\u732e', text):
        return None
    return next(iter(ids)).lower()


def read_json(path):
    if path.stat().st_size > 20 * 1024 * 1024:
        raise ValueError('Inspection JSON exceeds the 20 MiB read limit')
    raw = path.read_bytes()
    value = json.loads(raw)
    if not isinstance(value, dict):
        raise ValueError('Saved inspection/task report must be a JSON object')
    return value, hashlib.sha256(raw).hexdigest()


def task_steps(report):
    steps = report.get('steps', {})
    if not isinstance(steps, dict) or any(not isinstance(row, dict) for row in steps.values()):
        raise ValueError('Invalid saved task step records')
    return steps.values()


def digest_file(path):
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def render(report, task_id, report_path):
    lines = [f"Saved PDB {report['pdb_id']} inspection for task {task_id}.",
             'Read locally; no new inspection or scientific calculation. Report and structure hashes verified.',
             'Evidence mode: ' + str(report.get('evidence_mode', 'not recorded')) + '.']
    sections = [
        ('Deposited polymers', 'polymers', ('entity_id', 'description', 'type', 'sequence_length', 'deposited_author_chain_ids')),
        ('Observed polymer chains and coordinate coverage', 'observed_polymer_chains', ('model', 'label_asym_id', 'auth_asym_id', 'entity_id', 'observed_residue_count', 'coordinate_coverage')),
        ('Asymmetric-unit chain annotations', 'asymmetric_unit_chains', ('id', 'entity_id')),
        ('Biological assembly annotations', 'assemblies', ('id', 'details', 'method_details', 'oligomeric_details', 'oligomeric_count')),
        ('Assembly generators (label chain identifiers)', 'assembly_generators', ('assembly_id', 'oper_expression', 'asym_id_list')),
        ('Nonwater nonpolymer components; roles not automatically assigned', 'nonwater_nonpolymer_instances', ('model', 'comp_id', 'name', 'label_asym_id', 'auth_asym_id', 'auth_seq_id', 'insertion_code')),
        ('Deposited database associations', 'database_references', ('id', 'entity_id', 'db_name', 'pdbx_db_accession', 'db_code')),
        ('Deposited database alignment boundaries', 'database_alignments', ('align_id', 'ref_id', 'pdbx_strand_id', 'seq_align_beg', 'seq_align_end', 'db_align_beg', 'db_align_end')),
        ('Deposited sequence differences (not independent mutation verification)', 'sequence_differences', ('align_id', 'seq_num', 'mon_id', 'db_mon_id', 'pdbx_seq_db_seq_num', 'details'))]
    coverage = {}
    for title, key, fields in sections:
        rows = report.get(key, [])
        if not isinstance(rows, list) or any(not isinstance(row, dict) for row in rows):
            raise ValueError('Invalid inspection field: ' + key)
        selected = [{k: row[k] for k in fields if row.get(k) is not None} for row in rows[:20]]
        # Bound individual strings as well as row counts; disclose each shortened field.
        shortened = 0
        for row in selected:
            for k, value in row.items():
                if isinstance(value, str) and len(value) > 250:
                    row[k] = value[:250] + ' [field shortened]'
                    shortened += 1
        coverage[key] = dict(total=len(rows), displayed=len(selected), shortened_fields=shortened)
        lines.extend(['', title + f' ({len(selected)}/{len(rows)} rows):'])
        if selected:
            labels = dict(entity_id='entity', description='protein/polymer', type='polymer type',
                sequence_length='length', deposited_author_chain_ids='deposited author chains',
                label_asym_id='label chain', auth_asym_id='author chain',
                observed_residue_count='observed residues', coordinate_coverage='construct coverage',
                oligomeric_details='assembly', oligomeric_count='subunits',
                oper_expression='operators', asym_id_list='label chains', comp_id='component',
                auth_seq_id='author residue', insertion_code='insertion', pdbx_db_accession='accession',
                pdbx_seq_db_seq_num='database position', mon_id='deposited residue', db_mon_id='reference residue')
            def display(value):
                if isinstance(value, float):
                    return f'{value:.4f}'
                return re.sub(r'\s+', ' ', str(value)).replace('<', '&lt;').replace('>', '&gt;')
            lines.extend('- ' + '; '.join(labels.get(k, k.replace('_', ' ')) + ': ' + display(v)
                                         for k, v in row.items()) for row in selected)
        else:
            lines.append('No rows recorded in the saved report; not evidence of absence.')
    lines.extend(['', 'Interpretation and coverage:',
        '- Coordinate coverage is per model against the deposited construct, not full-length UniProt coverage.',
        '- Asymmetric-unit chain counts do not establish biological oligomerization.',
        '- GDP/GTP and ions are not automatically adopted as inhibitor screening references.',
        '- Isoform assignment: ' + str(report.get('isoform_assignment', 'not reported')) + '.',
        '- Sequence differences are deposited annotations. No fresh reference sequence alignment was performed.',
        '- At most 20 rows per section are displayed. Complete sequences, residue lists, assembly operators and other fields remain in the saved report.',
        '', 'Saved files:', '- Structure: ' + str(report['structure']), '- Inspection report: ' + str(report_path)])
    return '\n'.join(lines), coverage


def review(app, sid, args, request):
    task_id = args.get('task_id')
    if not isinstance(task_id, str) or not re.fullmatch('[a-f0-9]{16}', task_id):
        raise ValueError('Use an exact recorded task ID')
    save_note = args.get('save_note', False)
    if type(save_note) is not bool:
        raise ValueError('save_note must be boolean')
    if save_note and not note_requested(request):
        raise ValueError('Saving a project note requires an explicit current user request')
    job = app.task(sid, task_id)
    if job['status'] != 'complete':
        return dict(status='not_reviewed', task_id=task_id, task_status=job['status'],
            answer=f"Task {task_id}: {job['status']}. No completed inspection was reviewed or resubmitted.")
    ctx = app.context_for(sid)
    root = Path(ctx['project_root'])
    task_path = ensure_within(Path(job['report']), root)
    task_report, task_hash = read_json(task_path)
    steps = [row for row in task_steps(task_report)
             if row.get('action') == 'pdb_inspect' and row.get('status') == 'complete']
    if len(steps) != 1:
        raise ValueError('Select a task with exactly one completed pdb_inspect step')
    result = steps[0]['result']
    path = ensure_within(Path(result['report']), root)
    report, report_hash = read_json(path)
    if not result.get('report_sha256') or report_hash != result['report_sha256']:
        raise ValueError('Saved inspection report hash mismatch; no summary or note produced')
    if report.get('kind') != 'pdb_inspection' or report.get('pdb_id') != result.get('pdb_id'):
        raise ValueError('Inspection report identity mismatch')
    structure = ensure_within(Path(report['structure']), root)
    if structure != ensure_within(Path(result['path']), root):
        raise ValueError('Inspection structure path mismatch')
    structure_hash = digest_file(structure)
    if structure_hash != result.get('sha256') or structure_hash != report.get('structure_sha256'):
        raise ValueError('Saved structure hash mismatch; no summary or note produced')
    answer, coverage = render(report, task_id, path)
    note = dict(status='not_requested')
    if save_note:
        # A stable source-linked receipt is idempotent across repeated requests/chats.
        identifier = hashlib.sha256(f"pdb-review:{ctx['project_id']}:{task_id}:{report_hash}".encode()).hexdigest()[:32]
        text = (f"Reviewed saved PDB {report['pdb_id']} inspection from task {task_id}. "
                f"Deposited polymer entities: {len(report.get('polymers', []))}; "
                f"assembly annotations: {len(report.get('assemblies', []))}; "
                f"nonwater nonpolymer instances: {len(report.get('nonwater_nonpolymer_instances', []))}. "
                f"Isoform assignment: {report.get('isoform_assignment', 'not reported')}. "
                f"Report: {path}; SHA-256: {report_hash}. "
                'Deposited evidence, not independent biological verification. No new calculation or receptor adoption.')
        try:
            with app.connect() as db:
                inserted = db.execute('INSERT OR IGNORE INTO research_notes VALUES(?,?,?,?,?,?)',
                    (identifier, ctx['project_id'], 'finding', text, task_id, time.time())).rowcount
            note = dict(status='saved' if inserted else 'already_saved', note_id=identifier,
                        project_id=ctx['project_id'], task_id=task_id)
            from .project_guidance import render as render_project
            try:
                render_project(app, sid)
            except (ValueError, OSError, sqlite3.Error) as exc:
                note['document_refresh_error'] = type(exc).__name__
        except (ValueError, OSError, sqlite3.Error) as exc:
            note = dict(status='not_saved', error=type(exc).__name__)
    answer += '\n\nProject note: ' + note['status'] + ('. ID: ' + note['note_id'] if note.get('note_id') else '.')
    if note.get('document_refresh_error'):
        answer += ' The note is stored, but the generated project document could not be refreshed.'
    return dict(status='reviewed', task_id=task_id, answer=answer, source=str(path),
        report_sha256=report_hash, task_report_sha256=task_hash, structure_sha256=structure_hash,
        display_coverage=coverage, project_note=note, operation='local_saved_report_review_no_inference')


def coordinate(app, sid, request, on_event=None):
    task_id = request_task(request)
    if task_id is None:
        return None
    # Only take over requests for this explicit artifact type; other task reviews
    # retain their existing tools and scientific semantics.
    try:
        job = app.task(sid, task_id)
        if job['status'] == 'complete':
            root = Path(app.context_for(sid)['project_root'])
            report, _ = read_json(ensure_within(Path(job['report']), root))
            if not any(s.get('action') == 'pdb_inspect' for s in task_steps(report)):
                return None
    except (ValueError, OSError, KeyError, TypeError):
        pass
    from .pymol_bridge import write_json
    trace_id = uuid.uuid4().hex
    path = app.root / 'agents' / sid / (trace_id + '.json')
    path.parent.mkdir(parents=True, exist_ok=True)
    args = dict(task_id=task_id, save_note=note_requested(request))
    step = dict(step=dict(kind='tool', tool='pdb_review', arguments=args,
                         purpose='Read and summarize the saved inspection locally'), status='started')
    audit = dict(id=trace_id, request=request, status='running', started=time.time(),
                 steps=[step], answer_mode='local_pdb_report_review_no_model_call')
    write_json(path, audit)
    emit = on_event or (lambda event: None)
    emit(dict(type='tool_start', tool='pdb_review', text='Reading the saved inspection report locally.'))
    try:
        result = review(app, sid, args, request)
        answer = result['answer']
        step.update(status='complete', result=result, evidence_id='E1')
        audit.update(status='complete' if result['status'] == 'reviewed' else 'incomplete', evidence_ids=['E1'])
    except (ValueError, OSError, KeyError, TypeError, sqlite3.Error) as exc:
        error = app.safe_error(exc)
        step.update(status='failed', error=error)
        audit.update(status='incomplete', failure_stage='saved_report_reading')
        answer = 'Saved inspection review could not be completed: ' + error + '. No scientific task was submitted or resumed.'
    audit.update(answer=answer, finished=time.time())
    write_json(path, audit)
    emit(dict(type='tool_end', tool='pdb_review', status=step['status'], text='Saved inspection review: ' + step['status']))
    return answer + '\n\nAgent trace: ' + str(path)
