import importlib.util
import json
from pathlib import Path
import sqlite3
import sys

import numpy as np
import pytest

from test_export_selected_conformers import fixture, module as exporter

sys.path.insert(0, str(Path(__file__).parents[1] / 'scripts'))
spec = importlib.util.spec_from_file_location('top_conformers', Path(__file__).parents[1] / 'scripts/select_top_conformers.py')
select = importlib.util.module_from_spec(spec)
spec.loader.exec_module(select)


def setup(tmp):
    batch, search, exported, source, records = fixture(tmp)
    report = exporter.read(search / 'report.json')
    report['retrieval']['templates'] = [{'template': 'T0'}, {'template': 'T1'}]
    exporter.save(search / 'report.json', report)
    (search / 'adopted-design').mkdir()
    exporter.save(search / 'adopted-design/report.json', {'templates': [{'query_id': 'T0'}, {'query_id': 'T1'}]})
    exporter.save(search / 'protocol.json', {'template_quota': 2, 'rrf_k': 60})
    for ti in range(2):
        root = search / 'chunks' / f'{ti:02d}'
        root.mkdir(parents=True)
        p = root / '0000000000.npz'
        np.savez(p, global_ids=np.array([0, 1, 2]), molecule_ids=np.array(['m0', 'm0', 'm1']), stage_levels=np.array([7, 7, 7]))
        rows = [dict(global_id=g, molecule_id=m, conformer_id=f'c{g}', optional_score=score, gaussian_same_pose=.5)
                for g, m, score in [(0, 'm0', 10), (1, 'm0', 9), (2, 'm1', 1)]]
        p.with_suffix('.poses.jsonl').write_text(''.join(json.dumps(r)+'\n' for r in rows))
        exporter.save(p.with_suffix('.receipt.json'), {'files': {str(f): exporter.sha(f) for f in (p, p.with_suffix('.poses.jsonl'))}})
    exporter.run(batch, search, exported)
    return search, exported, tmp / 'top', records


def test_conformer_ranking_not_molecule_ranking_and_resume(tmp_path):
    search, exported, output, records = setup(tmp_path)
    r = select.run(search, exported, output, count=2, part_size=1)
    assert r['exported_conformers'] == 2 and r['distinct_molecules'] == 1
    assert r['scored_conformers'] == 3
    assert (output / 'top-00001.mol2').read_text() == records[0]
    assert (output / 'top-00002.mol2').read_text() == records[1]
    assert select.run(search, exported, output, count=2, part_size=1, resume=True)['exported_conformers'] == 2


def test_missing_template_and_resume(tmp_path):
    search, exported, output, _ = setup(tmp_path)
    root = search / 'chunks/01'
    moved = search / 'saved-template'
    root.rename(moved)
    with pytest.raises(ValueError, match='Missing per-conformer'):
        select.run(search, exported, output, count=2)
    assert not (output / 'report.json').exists()
    with sqlite3.connect(output / 'selection.sqlite') as db:
        assert db.execute('SELECT count(*) FROM completed').fetchone()[0] == 1
    moved.rename(root)
    assert select.run(search, exported, output, count=2, resume=True)['exported_conformers'] == 2
    with sqlite3.connect(output / 'selection.sqlite') as db:
        assert db.execute('SELECT min(support),max(support) FROM aggregate').fetchone() == (2, 2)


def test_hash_tampering_fails(tmp_path):
    search, exported, output, _ = setup(tmp_path)
    p = search / 'chunks/00/0000000000.poses.jsonl'
    p.write_text(p.read_text().replace('10', '11'))
    with pytest.raises(ValueError, match='Hash mismatch'):
        select.run(search, exported, output, count=2)


def test_no_unscored_padding(tmp_path):
    search, exported, output, _ = setup(tmp_path)
    for p in search.glob('chunks/*/*.poses.jsonl'):
        lines = p.read_text().splitlines(keepends=True)
        p.write_text(lines[0])
        receipt = p.with_name(p.name.replace('.poses.jsonl', '.receipt.json'))
        r = exporter.read(receipt)
        r['files'][str(p)] = exporter.sha(p)
        npz = p.with_name(p.name.replace('.poses.jsonl', '.npz'))
        np.savez(npz, global_ids=np.array([0, 1, 2]), molecule_ids=np.array(['m0', 'm0', 'm1']), stage_levels=np.array([7, 4, 4]))
        r['files'][str(npz)] = exporter.sha(npz)
        exporter.save(receipt, r)
    with pytest.raises(ValueError, match='No unscored padding'):
        select.run(search, exported, output, count=2)
