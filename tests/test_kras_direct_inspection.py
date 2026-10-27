import importlib.util
import json
from pathlib import Path
import sqlite3


spec = importlib.util.spec_from_file_location('direct_inspection', Path(__file__).resolve().parents[1] / 'scripts/inspect_kras_direct_output.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_report_and_readonly_database(tmp_path):
    report = dict(status='complete', full_catalog=True,
                  scope='direct_contacts_only_not_final_screening_hits',
                  processed=25813808, scheduled=25813808,
                  targets={'KRAS_G12D': dict(direct_poses=5939689, distinct_molecules=2632299),
                           'KRAS_G12V': dict(direct_poses=64431310, distinct_molecules=4080670)},
                  union_distinct_molecules=4144784, shared_distinct_molecules=2568185)
    path = tmp_path / 'report.json'
    path.write_text(json.dumps(report))
    db = tmp_path / 'candidates.sqlite'
    with sqlite3.connect(db) as connection:
        connection.execute('CREATE TABLE candidates(target TEXT, mid TEXT, gid INTEGER)')
        connection.execute('INSERT INTO candidates VALUES (?,?,?)', ('KRAS_G12D', 'MOL-test', 12))
    connection.close()
    before = db.read_bytes()
    result = module.inspect(tmp_path)
    assert all(result['checks'].values())
    assert result['databases'][0]['tables'][0]['sample_rows'] == [['KRAS_G12D', 'MOL-test', 12]]
    assert db.read_bytes() == before
    report['processed'] = 1
    path.write_text(json.dumps(report))
    assert module.inspect(tmp_path)['status'] == 'report_checks_failed'


if __name__ == '__main__':
    import tempfile
    with tempfile.TemporaryDirectory() as directory:
        test_report_and_readonly_database(Path(directory))
    print('PASS: completed report, read-only SQLite samples, incomplete report rejection')
