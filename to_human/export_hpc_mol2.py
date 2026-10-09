"""Export the twelve completed E106 arms into a deduplicated HPC MOL2 bundle.

Standard-library only. No chemistry conversion, search, scoring or job submission.
Only molecule title lines change; original coordinates and atom/bond data remain.
"""
import argparse
import csv
import hashlib
import json
import tarfile
from pathlib import Path

RUNS = ('5f7c28e4bd704616 3a90516120f74140 4f0fad82c9714b1e '
        '34bf4cc4a192449e 7a7fddb92b744a38 9bd13fe58ad14df0 '
        'dcc64417cf534000 317795cf58ad4f04 803b293dd57c4aec '
        '82f400cb6c1c46df c71019343cb04b62 bfa50ed8eb9a489e').split()


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def verified(path, expected):
    raw = path.read_bytes()
    if not expected or sha(raw) != expected:
        raise ValueError(f'Missing hash or changed input: {path}')
    return raw


def records(raw):
    text = raw.decode('utf-8').replace('\r\n', '\n')
    marker = '@<TRIPOS>MOLECULE'
    parts = text.split(marker)
    if parts[0].strip():
        raise ValueError('Unexpected content before first MOL2 record')
    for part in parts[1:]:
        lines = (marker + part).rstrip('\n').splitlines(keepends=True)
        if len(lines) < 4 or '@<TRIPOS>ATOM' not in part or '@<TRIPOS>BOND' not in part:
            raise ValueError('Incomplete MOL2 record')
        yield lines


def export(runs_root, output, run_ids, batch_size=1000):
    if batch_size < 1:
        raise ValueError('Batch size must be positive')
    if output.exists():
        raise FileExistsError('Use a new output directory; existing files are never overwritten')
    output.mkdir(parents=True)
    (output/'ligands').mkdir()
    (output/'provenance').mkdir()
    seen, batches, arms, pending = {}, [], [], []
    manifest = dict(status='incomplete', coordinates='original_source',
                    deduplication='cid_only', batch_size=batch_size)
    (output/'manifest.json').write_text(json.dumps(manifest, indent=2))

    def flush():
        if not pending:
            return
        name = f'ligands/batch_{len(batches)+1:06d}.mol2'
        raw = ''.join(pending).encode('utf-8')
        (output/name).write_bytes(raw)
        batches.append(dict(batch_id=len(batches)+1, path=name,
                            conformers=len(pending), sha256=sha(raw)))
        pending.clear()

    with (output/'candidates.csv').open('w', encoding='utf-8', newline='') as cf, \
            (output/'memberships.csv').open('w', encoding='utf-8', newline='') as mf:
        cw, mw = csv.writer(cf), csv.writer(mf)
        cw.writerow(['ligand_name','cid','molecule_id','source_sha256','batch_id','record_index_1based'])
        mw.writerow(['arm','run','ligand_name','cid','molecule_id','block_id','export_order_1based'])
        for run in run_ids:
            folder = runs_root/('PROMPT-'+run)/'execution/campaign/blocks'
            report = read(folder/'report.json')
            if report.get('status') != 'complete' or report.get('kind') != 'block_search_dock':
                raise ValueError(f'Incomplete or wrong campaign: {run}')
            if report.get('docking_requested') is not False:
                raise ValueError('Expected search-only campaign')
            if report.get('mol2_coordinates') != 'original_source':
                raise ValueError('Unexpected coordinate provenance')
            selection = report['selection']
            arm = f"{selection['scheme']}/{selection['method']}/Top-{selection['blocks_per_method']}"
            if arm in {a['arm'] for a in arms}:
                raise ValueError('Duplicate comparison arm')
            selected = folder/'selected'
            sr = read(selected/'report.json')
            hashes = sr['output_hashes']
            raw = verified(selected/'candidates.jsonl', hashes.get('candidates.jsonl'))
            candidates = {}
            for order, line in enumerate(raw.decode('utf-8').splitlines(), 1):
                row = json.loads(line)
                if row['cid'] in candidates:
                    raise ValueError('Duplicate CID within arm')
                candidates[row['cid']] = (order, row)
            mids = {r['molecule_id'] for _, r in candidates.values()}
            if len(candidates) != report['exported_conformers'] or len(mids) != report['unique_molecules']:
                raise ValueError('Report and candidate counts disagree')
            chunks = json.loads(verified(selected/'exports.json', hashes.get('exports.json')))['chunks']
            visited = set()
            for chunk in chunks:
                path = (selected/chunk['path']).resolve()
                if not path.is_relative_to(selected.resolve()):
                    raise ValueError('Chunk outside selected directory')
                if chunk['sha256'] != hashes.get(chunk['path']):
                    raise ValueError('Inconsistent chunk hashes')
                mols = list(records(verified(path, chunk['sha256'])))
                if len(mols) != chunk['count'] or len(mols) != len(chunk['records']):
                    raise ValueError('Chunk record count mismatch')
                for lines, meta in zip(mols, chunk['records']):
                    cid = meta['cid']
                    if cid in visited or cid not in candidates:
                        raise ValueError('Unexpected or repeated exported CID')
                    visited.add(cid)
                    order, row = candidates[cid]
                    if (meta['mid'], meta['alias'], meta['source_sha256']) != (
                            row['molecule_id'], row['alias'], row['source_sha256']):
                        raise ValueError('Export metadata mismatch')
                    if lines[1].strip() != row['alias']:
                        raise ValueError('MOL2 title does not match export metadata')
                    lines[1] = 'CANONICAL_TITLE\n'
                    body_hash = sha((''.join(lines).rstrip('\n')+'\n').encode('utf-8'))
                    identity = (row['molecule_id'], row['source_sha256'], body_hash)
                    if cid in seen:
                        name, previous = seen[cid]
                        if previous != identity:
                            raise ValueError(f'Conflicting chemistry/coordinates for CID {cid}')
                    else:
                        name = f'conf_{len(seen)+1:09d}'
                        seen[cid] = (name, identity)
                        lines[1] = name+'\n'
                        cw.writerow([name,cid,row['molecule_id'],row['source_sha256'],len(batches)+1,len(pending)+1])
                        pending.append(''.join(lines).rstrip('\n')+'\n')
                        if len(pending) >= batch_size:
                            flush()
                    mw.writerow([arm,run,name,cid,row['molecule_id'],row['block_id'],order])
            if visited != set(candidates):
                raise ValueError('Some candidates have no exported MOL2')
            for name in ['report.json', 'protocol.json']:
                if (folder/name).exists():
                    (output/'provenance'/f'{run}-{name}').write_bytes((folder/name).read_bytes())
            arms.append(dict(arm=arm,run=run,conformers=len(candidates),molecules=len(mids),
                             shortfall=report['shortfall'],report_sha256=sha((folder/'report.json').read_bytes())))
            print(arm, 'verified;', len(seen), 'union conformers', flush=True)
        flush()
    with (output/'batches.tsv').open('w', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['batch_id','path','conformers','sha256'], delimiter='\t')
        writer.writeheader(); writer.writerows(batches)
    (output/'README.txt').write_text(
        'HPC docking input, original-source coordinates, not docked poses.\n'
        'Each ligands/batch_*.mol2 contains up to the configured number of records.\n'
        'Dispatch one independent batch per array task, with a separate output directory.\n'
        'All arms share globally unique conf_* ligand names. Preserve these in docking results.\n'
        'candidates.csv maps names to conformers/molecules; memberships.csv maps them to every arm.\n'
        'No molecular conversion, protonation change or geometry change was performed.\n'
        'Only MOL2 title lines and newline normalization changed. No molecule-level deduplication.\n'
        'Receptor, pocket, docking software and scheduler configuration are not included.\n'
        'These must be fixed identically across batches before submission.\n'
        'Do not submit an incomplete bundle; check manifest.json status first.\n', encoding='utf-8')
    manifest.update(status='complete', unique_conformers=len(seen), batches=len(batches), arms=arms,
                    outputs={p.relative_to(output).as_posix():sha(p.read_bytes())
                             for p in sorted(output.rglob('*')) if p.is_file() and p.name != 'manifest.json'})
    (output/'manifest.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    return manifest


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--runs-root', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--batch-size', type=int, default=1000)
    p.add_argument('--archive', action='store_true')
    a = p.parse_args()
    archive = a.output.with_name(a.output.name+'.tar.gz')
    if a.archive and archive.exists():
        raise FileExistsError(archive)
    result = export(a.runs_root,a.output,RUNS,a.batch_size)
    if a.archive:
        with tarfile.open(archive, 'x:gz') as tf:
            tf.add(a.output, arcname=a.output.name)
        print('Archive:', archive)
    print('COMPLETE:', result['unique_conformers'], 'conformers,', result['batches'], 'batches')


if __name__ == '__main__':
    main()
