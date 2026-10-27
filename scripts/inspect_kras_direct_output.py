"""Inspect completed direct-contact outputs without changing or rescanning them."""
import argparse
import hashlib
import json
from pathlib import Path
import sqlite3
import tarfile


def inspect(root):
    root = Path(root).resolve()
    report_path = root / 'report.json'
    raw = report_path.read_bytes()
    report = json.loads(raw)
    expected = {'KRAS_G12D': 5939689, 'KRAS_G12V': 64431310}
    checks = {
        'complete': report.get('status') == 'complete',
        'full_catalog': report.get('full_catalog') is True,
        'scope': report.get('scope') == 'direct_contacts_only_not_final_screening_hits',
        'processed': report.get('processed') == report.get('scheduled') == 25813808,
    }
    for target, poses in expected.items():
        checks[target + '_pose_count'] = report.get('targets', {}).get(target, {}).get('direct_poses') == poses
    targets = report.get('targets', {})
    counts = [targets.get(t, {}).get('distinct_molecules') for t in expected]
    union, shared = report.get('union_distinct_molecules'), report.get('shared_distinct_molecules')
    checks['union_identity'] = (all(isinstance(x, int) for x in counts + [union, shared])
                                and sum(counts) - shared == union
                                and 0 <= shared <= min(counts))
    result = dict(report=report, checks=checks, report_sha256=hashlib.sha256(raw).hexdigest(),
                  output_directory=str(root), files=[], databases=[],
                  limitations=['Database rows are samples, not an independent full recount.',
                               'Direct contacts do not establish full steric acceptance or binding.',
                               'No selection, export, docking or pose search is performed.'])
    for path in sorted(root.rglob('*')):
        if not path.is_file():
            continue
        relative = str(path.relative_to(root))
        result['files'].append(dict(path=relative, bytes=path.stat().st_size))
        if path.suffix.lower() not in ('.db', '.sqlite', '.sqlite3'):
            continue
        db = dict(path=relative, tables=[])
        connection = sqlite3.connect(path.as_uri() + '?mode=ro', uri=True, timeout=10)
        try:
            connection.execute('PRAGMA query_only=ON')
            for name, sql in connection.execute("SELECT name, sql FROM sqlite_master WHERE type='table' AND name NOT LIKE 'sqlite_%'").fetchall():
                quoted = '"' + name.replace('"', '""') + '"'
                columns = connection.execute('PRAGMA table_info(' + quoted + ')').fetchall()
                cursor = connection.execute('SELECT * FROM ' + quoted + ' LIMIT 3')
                rows = []
                for row in cursor:
                    rows.append([{'blob_bytes': len(v)} if isinstance(v, bytes)
                                 else v[:3000] if isinstance(v, str) else v for v in row])
                db['tables'].append(dict(name=name, sql=sql, columns=columns, sample_rows=rows))
        finally:
            connection.close()
        result['databases'].append(db)
    result['status'] = 'inspection_complete' if all(checks.values()) else 'report_checks_failed'
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-dir', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--count-script', type=Path)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('Use a new output directory; existing output is preserved.')
    if args.output.resolve().is_relative_to(args.run_dir.resolve()):
        parser.error('Write the inspection outside the completed run directory.')
    result = inspect(args.run_dir)
    args.output.mkdir(parents=True)
    destination = args.output / 'inspection.json'
    destination.write_text(json.dumps(result, indent=2), encoding='utf-8')
    bundle = args.output / 'kras-direct-inspection.tar.gz'
    with tarfile.open(bundle, 'w:gz') as archive:
        archive.add(destination, arcname='inspection.json')
        for name in ('report.json', 'protocol.json'):
            path = args.run_dir / name
            if path.is_file():
                archive.add(path, arcname=name)
        if args.count_script:
            archive.add(args.count_script, arcname='KRAS_count_direct_molecules.py')
    print(json.dumps(dict(status=result['status'], checks=result['checks'],
                          databases=len(result['databases']), bundle=str(bundle)), indent=2))
    if not all(result['checks'].values()):
        raise SystemExit(2)


if __name__ == '__main__':
    main()
