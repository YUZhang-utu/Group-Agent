"""Sealed, resumable EquiScore evaluation of existing PLANTS poses."""
import argparse
from contextlib import closing
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import sqlite3
import subprocess
import time

from .block_sampling import check_outputs, exclusive_output, readonly
from .final_work_blocks import sha
from .joint_spatial_profiles import save
from .project_context import ensure_within

COMMIT = '8b2a9289cf7d181fa49de6ac6712260e8c500c4a'
WEIGHT_SHA = 'd4367bb73686b2363e238abb778fab55e2924458ec1ced561072bd82f711695d'


def load_profile(path):
    path = Path(path).resolve()
    config = json.loads(path.read_text(encoding='utf-8'))
    if set(config) - {'python', 'repository', 'receptor_pdbs'} or not {'python', 'repository'} <= config.keys():
        raise ValueError('EquiScore profile requires python and repository; optional receptor_pdbs')
    for key in ('python', 'repository'):
        item = Path(config[key])
        item = item if item.is_absolute() else path.parent / item
        # A venv Python symlink must retain its venv path to load overlay packages.
        config[key] = os.path.abspath(item) if key == 'python' else str(item.resolve())
    if not isinstance(config.get('receptor_pdbs', {}), dict):
        raise ValueError('receptor_pdbs must map receptor IDs to PDB paths')
    config['receptor_pdbs'] = {key: str((Path(value) if Path(value).is_absolute() else path.parent / value).resolve())
                               for key, value in config.get('receptor_pdbs', {}).items()}
    repo = Path(config['repository'])
    result = subprocess.run(['git', '-C', str(repo), 'rev-parse', 'HEAD'], capture_output=True, text=True, check=True)
    if result.stdout.strip() != COMMIT: raise ValueError('Install the pinned EquiScore revision')
    if subprocess.run(['git', '-C', str(repo), 'diff', '--exit-code', 'HEAD', '--', '*.py'], capture_output=True).returncode:
        raise ValueError('Official EquiScore Python sources have local changes')
    weights = repo / 'workdir/official_weight/save_model_screen.pt'
    if sha(weights) != WEIGHT_SHA: raise ValueError('Official screening checkpoint differs')
    worker = Path(__file__).resolve().parents[2] / 'scripts/equiscore_worker.py'
    inputs = {str(p): sha(p) for p in [path, worker, weights, *sorted(repo.rglob('*.py'))]}
    return config, worker, inputs


def pilot_cids(scores, count=32):
    with Path(scores).open(encoding='utf-8', newline='') as stream:
        ids = {row['cid'] for row in csv.DictReader(stream)}
    return sorted(ids, key=lambda cid: hashlib.sha256(('20261005:' + cid).encode()).hexdigest())[:count]


def collect(source, predictions, selected=None):
    """Identity join with one explicit result for every requested docking pair."""
    with readonly(predictions) as db:
        predicted = {}
        for cid, receptor, status, score, error in db.execute('SELECT cid,receptor,status,score,error FROM scores'):
            if (cid, receptor) in predicted: raise ValueError('Duplicate prediction identity')
            predicted[cid, receptor] = status, score, error
    with Path(source).open(encoding='utf-8', newline='') as stream:
        rows = []
        seen = set()
        for row in csv.DictReader(stream):
            if selected is not None and row['cid'] not in selected: continue
            key = row['cid'], row['receptor']
            if key in seen: raise ValueError('Duplicate source conformer/receptor')
            seen.add(key)
            status, value, error = predicted.pop(key, ('not_run', None, 'Worker did not return this pair'))
            if status not in {'ok', 'unsupported_input', 'model_failed', 'not_run'}:
                raise ValueError('Unknown prediction status')
            if status == 'ok':
                if row['status'] != 'ok' or type(value) not in (int, float) or not math.isfinite(value) or not 0 <= value <= 1:
                    raise ValueError('Invalid EquiScore probability or unscored docking source')
            elif value is not None:
                raise ValueError('Failed EquiScore results must be null')
            rows.append(dict(row, chemplp=row['score'], status=status, score=value, error=error))
    if predicted: raise ValueError('Unknown/unrequested prediction identities')
    if not rows: raise ValueError('No requested source pairs')
    return rows


def verified_result(directory):
    """Verify completed or partial receipts without treating partial as success."""
    directory = Path(directory).resolve()
    report = json.loads((directory / 'report.json').read_text(encoding='utf-8'))
    if report.get('kind') != 'block_equiscore_run' or report.get('status') not in {'complete', 'partial'}:
        raise ValueError('Unsupported EquiScore receipt')
    for name, digest in report['output_hashes'].items():
        path = ensure_within(directory / name, directory)
        if sha(path) != digest: raise ValueError('EquiScore artifact changed: ' + name)
    return report


def write_csv(path, rows, fields=None):
    with Path(path).open('w', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields or list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def worker_failure(output, reason, exit_code, log_start=0):
    """Expose a bounded tail of this invocation, without masking the worker error."""
    log = Path(output) / 'worker.log'
    with log.open('rb') as stream:
        stream.seek(0, 2)
        stream.seek(max(log_start, stream.tell() - 12000))
        tail = '\n'.join(stream.read().decode('utf-8', errors='replace').splitlines()[-60:])
    return RuntimeError(f'{reason} (worker exit code {exit_code}). Log: {log}\n'
                        + (tail or 'No output was recorded by this worker invocation.'))


@exclusive_output
def run(source, profile, output, *, scope='pilot', pilot_report=None):
    from .prompt_workflow import file_lock
    source, output = Path(source).resolve(), Path(output).resolve()
    if scope not in {'pilot', 'full'}: raise ValueError('Choose pilot or full')
    config, worker, inputs = load_profile(profile)
    with file_lock(source.parent.parent / ('.' + source.parent.name + '.e097.lock')):
        original = check_outputs(source.parent)
        if original.get('kind') != 'block_plants_run' or original['status'] != 'complete':
            raise ValueError('A complete, sealed PLANTS run is required')
        prepared = Path(original['prepared_report'])
        prep = check_outputs(prepared.parent)
        run_signature = json.loads((source.parent / 'signature.json').read_text())
        if run_signature.get('prepared_sha256') != sha(prepared): raise ValueError('Preparation binding changed')
        sample = Path(original['sample_report'])
        check_outputs(sample.parent)
        prep_signature = json.loads((prepared.parent / 'signature.json').read_text())
        if prep_signature.get('sample_sha256') != sha(sample) or prep.get('sample_report') != str(sample):
            raise ValueError('Sample binding changed')
        receptors = {}
        for row in json.loads((prepared.parent / 'profile.json').read_text())['receptors']:
            pdb = Path(config.get('receptor_pdbs', {}).get(row['id'], str(Path(row['mol2']).parent / 'protein-input.pdb'))).resolve()
            mol2 = prepared.parent / 'receptors' / (row['id'] + '.mol2')
            if not pdb.is_file(): raise ValueError('Configure receptor_pdbs.' + row['id'] + ' with the original aligned protein-input.pdb')
            receptors[row['id']] = dict(pdb=str(pdb), mol2=str(mol2))
            inputs.update({str(pdb): sha(pdb), str(mol2): sha(mol2)})
        inputs.update({str(p): sha(p) for p in [source, source.parent / 'scores.csv', prepared, sample]})
        identity = dict(source_sha256=sha(source), profile_inputs=inputs, implementation_sha256=sha(Path(__file__)))
        if scope == 'full':
            if not pilot_report: raise ValueError('A successful technical pilot report is required before full scoring')
            pilot_report = Path(pilot_report).resolve()
            pilot = json.loads(pilot_report.read_text())
            check_outputs(pilot_report.parent)
            if pilot.get('scope') != 'pilot' or pilot.get('status') != 'complete' or pilot.get('identity') != identity or pilot.get('failed_pairs') != 0:
                raise ValueError('Pilot is incomplete, failed or bound to different inputs/code')
        selected = pilot_cids(source.parent / 'scores.csv') if scope == 'pilot' else None
        job = dict(profile=config, input_hashes=inputs, scores=str(source.parent / 'scores.csv'),
                   source_root=str(source.parent), source_hashes=original['output_hashes'],
                   receptors=receptors, selected_cids=selected)
        if scope == 'full': job['pilot_environment'] = pilot['environment']
        signature = dict(identity=identity, scope=scope, selected_cids=selected)
        if output.exists():
            if not (output / 'signature.json').is_file() or json.loads((output / 'signature.json').read_text()) != signature:
                raise ValueError('Use a fresh output or matching resume inputs')
            if (output / 'report.json').is_file():
                previous = verified_result(output)
                if previous['status'] == 'complete': return previous
                # Preserve the sealed snapshot before resuming the mutable database.
                # An interruption must not leave a stale receipt blocking recovery.
                history = output / 'history'
                history.mkdir(exist_ok=True)
                (output / 'report.json').replace(history / ('report-' + sha(output / 'report.json') + '.json'))
        else: output.mkdir(parents=True)
        save(output / 'signature.json', signature); save(output / 'job.json', job)
        started = time.monotonic()
        worker_env = os.environ.copy()
        worker_env['DGLBACKEND'] = 'pytorch'
        if os.name != 'nt':
            lib = str(Path(config['python']).parent.parent / 'lib')
            worker_env['LD_LIBRARY_PATH'] = lib + (':' + worker_env['LD_LIBRARY_PATH'] if worker_env.get('LD_LIBRARY_PATH') else '')
        with (output / 'worker.log').open('a', encoding='utf-8') as log:
            log_start = log.tell()
            result = subprocess.run([config['python'], str(worker), '--job', str(output / 'job.json')], stdout=log, stderr=subprocess.STDOUT, shell=False, env=worker_env)
        for name, digest in inputs.items():
            if sha(name) != digest: raise ValueError('Source/configuration changed during EquiScore execution')
        database = output / 'predictions.sqlite'
        if not database.is_file():
            raise worker_failure(output, 'EquiScore did not produce predictions', result.returncode, log_start)
        rows = collect(source.parent / 'scores.csv', database, set(selected) if selected is not None else None)
        write_csv(output / 'scores.csv', rows)
        good = sum(row['status'] == 'ok' for row in rows)
        receipt = output / 'environment.json'
        if not receipt.is_file(): raise worker_failure(output, 'Worker environment receipt is missing', result.returncode, log_start)
        report = dict(kind='block_equiscore_run', status='complete' if good == len(rows) and result.returncode == 0 else 'partial',
            scope=scope, identity=identity, score='EquiScore', better='higher', model='EquiScore screening',
            model_commit=COMMIT, weights_sha256=WEIGHT_SHA, model_executed=True,
            pairs=len(rows), scored=good, failed_pairs=len(rows)-good, worker_exit_code=result.returncode,
            elapsed_seconds=time.monotonic()-started, source_report=str(source), sample_report=str(sample),
            ligand_chemistry_reviewed=prep.get('ligand_chemistry_reviewed', False),
            environment=json.loads(receipt.read_text()),
            output_hashes={name: sha(output / name) for name in ['scores.csv', 'predictions.sqlite', 'signature.json', 'job.json', 'environment.json']},
            limitations=['Exploratory screening score, not measured affinity or experimental enrichment.',
                         'Pilot establishes technical compatibility only; original ligand chemistry review remains separate.'])
        save(output / 'report.json', report)
        return report


@exclusive_output
def analyze(source, output):
    """Create model-native block tables, then reuse the Top-N ranking interface."""
    from collections import defaultdict
    from .block_ranking import add_rankings
    source, output = Path(source).resolve(), Path(output).resolve()
    r = verified_result(source.parent)
    if r.get('kind') != 'block_equiscore_run' or r.get('scope') != 'full':
        raise ValueError('Only a full EquiScore panel can produce block rankings')
    sample_path = Path(r['sample_report'])
    if sha(sample_path) != r['identity']['profile_inputs'][str(sample_path)]:
        raise ValueError('Sample report changed after inference')
    check_outputs(sample_path.parent)
    if output.exists(): raise ValueError('Use a fresh EquiScore analysis directory')
    output.mkdir(parents=True)
    with (source.parent / 'scores.csv').open(encoding='utf-8', newline='') as stream:
        scores = {(row['cid'], row['receptor']): row for row in csv.DictReader(stream)}
    receptors = sorted({rec for _, rec in scores})
    counts = defaultdict(lambda: [0, 0])
    fields = ['scheme', 'block_id', 'cid', 'receptor', 'status', 'score', 'pose_name', 'pose_file', 'pose_index']
    with readonly(Path(r['sample_report']).parent / 'samples.sqlite') as db, (output / 'block_scores.csv').open('w', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields); writer.writeheader()
        populations = {(s, b): n for s, b, n in db.execute('SELECT scheme,block_id,n FROM population')}
        for scheme, bid, cid in db.execute('SELECT scheme,block_id,cid FROM sample'):
            for rec in receptors:
                row = scores.get((cid, rec))
                if row is None: raise ValueError('Full score table lost a conformer/receptor')
                writer.writerow({key: row[key] for key in fields if key not in {'scheme', 'block_id'}} | dict(scheme=scheme, block_id=bid))
                counts[scheme, bid, rec][0] += 1
                counts[scheme, bid, rec][1] += row['status'] == 'ok'
    summary = [dict(scheme=s, block_id=b, receptor=rec, population=populations[s, b],
                    sampled=n, scored=k, missing=n-k) for (s, b, rec), (n, k) in sorted(counts.items())]
    write_csv(output / 'block_summary.csv', summary)
    save(output / 'signature.json', dict(sha256=sha(source)))
    result = dict(kind='block_analyze', status='complete', comparisons=[],
        model='EquiScore', model_executed=True, source_report=str(source),
        output_hashes={name: sha(output / name) for name in ['block_scores.csv', 'block_summary.csv', 'signature.json']},
        limitations=r['limitations'])
    save(output / 'report.json', result)
    result = add_rankings(source, output)
    # Recompute the comparison baseline from exactly the same block memberships.
    baseline = output / 'chemplp_baseline'
    baseline.mkdir()
    with (output / 'block_scores.csv').open(encoding='utf-8', newline='') as stream, (baseline / 'block_scores.csv').open('w', encoding='utf-8', newline='') as target:
        writer = csv.DictWriter(target, fieldnames=fields); writer.writeheader()
        for row in csv.DictReader(stream):
            raw = scores[row['cid'], row['receptor']]
            row['score'], row['status'] = raw['chemplp'], 'ok'
            if not math.isfinite(float(row['score'])): raise ValueError('Invalid original ChemPLP score')
            writer.writerow(row)
    complete_summary = [dict(row, scored=row['sampled'], missing=0) for row in summary]
    write_csv(baseline / 'block_summary.csv', complete_summary)
    original_source = Path(r['source_report'])
    if sha(original_source) != r['identity']['source_sha256']: raise ValueError('Original docking report changed')
    save(baseline / 'signature.json', dict(sha256=sha(original_source)))
    save(baseline / 'report.json', dict(kind='block_analyze', status='complete', output_hashes={
        name: sha(baseline / name) for name in ['block_scores.csv', 'block_summary.csv', 'signature.json']}))
    add_rankings(original_source, baseline)
    comparison = []
    with readonly(baseline / 'block_ranking.sqlite') as old, readonly(output / 'block_ranking.sqlite') as new:
        old_rows = {(s, b, rec): (rank, mean) for s, b, rec, rank, mean in old.execute(
            'SELECT scheme,block_id,receptor,rank,top_n_mean FROM rankings WHERE top_n=10')}
        for s, b, rec, rank, mean in new.execute('SELECT scheme,block_id,receptor,rank,top_n_mean FROM rankings WHERE top_n=10'):
            previous, old_mean = old_rows[s, b, rec]
            query = 'SELECT molecule_id FROM candidates WHERE scheme=? AND block_id=? AND receptor=? AND candidate_rank<=10'
            a = {row[0] for row in old.execute(query, (s, b, rec))}
            z = {row[0] for row in new.execute(query, (s, b, rec))}
            comparison.append(dict(scheme=s, block_id=b, receptor=rec, chemplp_rank=previous,
                equiscore_rank=rank, rank_improvement=previous-rank if previous is not None and rank is not None else None,
                chemplp_top10_mean=old_mean, equiscore_top10_mean=mean,
                shared_top10_molecules=len(a & z), chemplp_candidates=len(a), equiscore_candidates=len(z)))
    write_csv(output / 'rank_comparison.csv', comparison)
    result['output_hashes']['rank_comparison.csv'] = sha(output / 'rank_comparison.csv')
    result['output_hashes'].update({p.relative_to(output).as_posix(): sha(p) for p in baseline.iterdir() if p.is_file()})
    result['outputs']['rank_comparison.csv'] = str(output / 'rank_comparison.csv')
    save(output / 'report.json', result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='action', required=True)
    run_parser = sub.add_parser('run')
    run_parser.add_argument('--source', required=True); run_parser.add_argument('--profile', required=True)
    run_parser.add_argument('--output', required=True); run_parser.add_argument('--scope', choices=['pilot', 'full'], default='pilot')
    run_parser.add_argument('--pilot-report')
    analysis = sub.add_parser('analyze')
    analysis.add_argument('--source', required=True); analysis.add_argument('--output', required=True)
    args = parser.parse_args(); options = vars(args).copy(); action = options.pop('action')
    result = run(**options) if action == 'run' else analyze(**options)
    print(json.dumps(result, indent=2))
    if result.get('status') != 'complete': raise SystemExit(2)


if __name__ == '__main__': main()
