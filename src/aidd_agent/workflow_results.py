"""Read-only, session-owned snapshots of external scientific workflow reports."""
import hashlib
import json
from pathlib import Path
import re
import time

MAX_BYTES = 16 * 1024 * 1024
KINDS = {'kras_export', 'af3_batch', 'af3_prediction'}


def count(value, label):
    if type(value) is not int or value < 0:
        raise ValueError(label + ' must be a nonnegative integer')
    return value


def digest(value):
    if not isinstance(value, str) or not re.fullmatch(r'[a-fA-F0-9]{64}', value):
        raise ValueError('Missing or invalid reported SHA256')
    return value


def read(path, max_bytes=MAX_BYTES):
    from .prompt_plan import strict_json
    with path.open('rb') as stream:
        raw = stream.read(max_bytes + 1)
    if len(raw) > max_bytes:
        raise ValueError('Report exceeds the read size limit')
    report = strict_json(raw.decode('utf-8-sig'))
    if not isinstance(report, dict):
        raise ValueError('Report must be an object')
    return report, hashlib.sha256(raw).hexdigest()


def summarize(kind, report):
    """Validate counters, not scientific truth or referenced payload contents."""
    if kind == 'kras_export':
        if report.get('scope') != 'existing_pool_grid_docking_handoff' or report.get('status') != 'complete':
            raise ValueError('Expected a completed E123 existing-pool export report')
        names = {'KRAS_G12D', 'KRAS_G12V'}
        if not isinstance(report.get('targets'), dict) or not isinstance(report.get('exports'), dict):
            raise ValueError('Expected target and export objects')
        if set(report.get('targets', {})) != names or set(report.get('exports', {})) != names:
            raise ValueError('Expected both KRAS target counters and export receipts')
        targets = {}
        for name in sorted(names):
            row, export = report['targets'][name], report['exports'][name]
            if not isinstance(row, dict) or not isinstance(export, dict):
                raise ValueError('Expected per-target counters and export objects')
            for key in ('pool', 'processed', 'strict', 'balanced', 'lenient', 'unresolved'):
                count(row.get(key), name + '.' + key)
            if not (row['pool'] == row['processed'] and
                    row['strict'] <= row['balanced'] <= row['lenient'] <= row['processed'] and
                    row['unresolved'] + row['lenient'] == row['processed']):
                raise ValueError('Inconsistent pool/profile counters for ' + name)
            for key in ('selected', 'exported', 'quarantined', 'requested_cap'):
                count(export.get(key), name + '.' + key)
            if (export.get('status') != 'complete' or export.get('target') != name or
                    export['selected'] != export['exported'] + export['quarantined'] or
                    export['selected'] != min(row['balanced'], export['requested_cap'])):
                raise ValueError('Inconsistent export counters for ' + name)
            batches = export.get('batches')
            if not isinstance(batches, list) or len(batches) > 1000:
                raise ValueError('Invalid export batch list')
            paths = set()
            for batch in batches:
                if not isinstance(batch, dict):
                    raise ValueError('Expected a batch object')
                path = batch.get('path')
                if not isinstance(path, str) or not path or path in paths:
                    raise ValueError('Missing or duplicate batch path')
                paths.add(path)
                count(batch.get('molecules'), 'batch.molecules')
                digest(batch.get('sha256'))
            if sum(b['molecules'] for b in batches) != export['exported']:
                raise ValueError('Batch totals do not match exported molecules')
            targets[name] = dict(**row, selected=export['selected'], exported=export['exported'],
                quarantined=export['quarantined'], requested_cap=export['requested_cap'],
                folder=export.get('folder'), batch_count=len(batches), batches=batches)
        return dict(status='complete', scope=report['scope'], targets=targets,
            note='Per-target molecule IDs; cross-target overlap is not resolved. Original MOL2 coordinates are not docking poses.')
    if kind == 'af3_batch':
        mode, status, jobs = report.get('mode'), report.get('status'), report.get('jobs')
        if mode not in {'prediction', 'host_check'} or status not in {
                'running', 'failed', 'complete', 'host_checks_passed_not_inference'} or not isinstance(jobs, list):
            raise ValueError('Expected an AF3 SMILES batch run report')
        if (status == 'complete' and mode != 'prediction') or (
                status == 'host_checks_passed_not_inference' and mode != 'host_check'):
            raise ValueError('AF3 mode and completion status disagree')
        ids, rows = set(), []
        for job in jobs:
            if not isinstance(job, dict):
                raise ValueError('Expected an AF3 job object')
            identifier = job.get('id')
            if not isinstance(identifier, str) or not identifier or identifier in ids:
                raise ValueError('Missing or duplicate AF3 job ID')
            ids.add(identifier)
            state = job.get('status')
            allowed = {'model_generated', 'failed'} if mode == 'prediction' else {'host_checks_passed', 'failed'}
            if state not in allowed:
                raise ValueError('AF3 job status does not match batch mode')
            digest(job.get('input_sha256'))
            if state == 'model_generated':
                if not isinstance(job.get('model'), str) or not job['model']:
                    raise ValueError('AF3 generated model has no path')
                digest(job.get('model_sha256'))
            rows.append({key: job[key] for key in ('id', 'status', 'model', 'model_sha256', 'confidence', 'error') if key in job})
        expected = 'model_generated' if mode == 'prediction' else 'host_checks_passed'
        total = report.get('total_jobs')
        if total is not None:
            count(total, 'total_jobs')
            if total < len(jobs) or (status in {'complete', 'host_checks_passed_not_inference'} and total != len(jobs)):
                raise ValueError('AF3 total_jobs does not match recorded completion')
        if status in {'complete', 'host_checks_passed_not_inference'} and (
                not jobs or any(j['status'] != expected for j in jobs)):
            raise ValueError('AF3 completed batch contains missing or failed results')
        return dict(status=status, mode=mode, recorded_jobs=len(jobs), total_jobs=total,
            pending_jobs=None if total is None else total-len(jobs),
            generated_models=sum(j['status'] == 'model_generated' for j in jobs),
            failed_jobs=sum(j['status'] == 'failed' for j in jobs), jobs=rows[:100],
            jobs_truncated=len(rows) > 100,
            note='Recorded jobs may exclude pending inputs. Confidence is not affinity or experimental binding evidence.')
    if kind == 'af3_prediction':
        if report.get('backend') != 'alphafold3' or not report.get('structure_path'):
            raise ValueError('Expected an AF3 prediction manifest')
        digest(report.get('structure_sha256'))
        return dict(status='model_reported', **{key: report.get(key) for key in (
            'construct_name', 'chain_ids', 'structure_path', 'structure_sha256', 'confidence', 'ranking_score')},
            note='Prediction evidence only; confidence is not affinity or experimental binding evidence.')
    raise ValueError('Unsupported workflow result kind')


def source_path(app, supplied):
    from .prompt_workflow import project_root
    if not isinstance(supplied, str) or not supplied or not Path(supplied).is_absolute():
        raise ValueError('Provide an exact absolute report path on the Chat host')
    ctx = app.context
    project = project_root(Path(ctx['db']), ctx['user_id'], ctx['project_id']).resolve()
    path = Path(supplied).resolve()
    roots = [project, *getattr(app, 'result_import_roots', ())]
    if not any(path.is_relative_to(root) for root in roots):
        raise ValueError('Report is outside the active Project and configured --result-import-root directories')
    if path.suffix.lower() != '.json':
        raise ValueError('Expected a JSON report')
    return path


def handle(app, sid, args):
    from .prompt_workflow import _atomic_json
    app.session(sid)
    if not isinstance(args, dict):
        raise ValueError('workflow_results arguments must be an object')
    operation = args.get('operation', 'list')
    fields = {'list': set(), 'attach': {'kind', 'report'}, 'check': {'attachment_id'},
              'read': {'attachment_id', 'pointer', 'offset', 'limit'}}
    if not isinstance(operation, str) or operation not in fields or set(args) - fields[operation] - {'operation'}:
        raise ValueError('Invalid workflow_results operation or fields')
    root = app.root / 'workflow-results' / sid
    with app.lock:
        if operation == 'list':
            attachments = []
            for path in sorted(root.glob('*.json')):
                item, _ = read(path, MAX_BYTES * 4)
                attachments.append({key: value for key, value in item.items() if key != 'report_snapshot'})
            return dict(attachments=attachments, computation_launched=False)
        if operation == 'attach':
            kind = args.get('kind')
            if kind not in KINDS:
                raise ValueError('Choose kras_export, af3_batch or af3_prediction')
            path = source_path(app, args.get('report'))
            report, sha = read(path)
            summary = summarize(kind, report)
            identifier = hashlib.sha256((sid + kind + str(path) + sha).encode()).hexdigest()[:24]
            row = dict(id=identifier, session=sid, kind=kind, source=str(path), report_sha256=sha,
                imported_at=time.time(), verification='report_schema_and_counters_only',
                payload_verification='not_performed', source_status='unchanged_at_import',
                summary=summary, report_snapshot=report)
            root.mkdir(parents=True, exist_ok=True)
            dest = root / (identifier + '.json')
            if not dest.exists():
                _atomic_json(dest, row)
            return dict(status='attached', attachment_id=identifier, summary=summary,
                        payload_verification='not_performed', computation_launched=False)
        identifier = args.get('attachment_id')
        if not isinstance(identifier, str) or not re.fullmatch('[a-f0-9]{24}', identifier):
            raise ValueError('Invalid attachment ID')
        saved = root / (identifier + '.json')
        if not saved.is_file():
            raise ValueError('Attachment does not belong to this session')
        row, _ = read(saved, MAX_BYTES * 4)
        if operation == 'read':
            from .domain_tools import page
            return dict(attachment_id=identifier, report_sha256=row['report_sha256'],
                        content=page(row['report_snapshot'], args), computation_launched=False)
        try:
            _, sha = read(source_path(app, row['source']))
            state = 'unchanged' if sha == row['report_sha256'] else 'changed_reimport_to_update'
        except (OSError, ValueError):
            state = 'unavailable_or_invalid'
        return dict(attachment_id=identifier, source_status=state, checked_at=time.time(),
                    snapshot_retained=True, computation_launched=False)
