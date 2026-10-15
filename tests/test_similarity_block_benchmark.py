"""E110 synthetic tests reuse E109's mocked scientific evaluator."""
import argparse
import csv
import importlib.util
from pathlib import Path
import sqlite3

import pytest
from test_matched_block_benchmark import prepared, engine, arguments
from aidd_agent import matched_block_benchmark as b

spec = importlib.util.spec_from_file_location('similarity_benchmark', Path(__file__).parents[1] / 'scripts/benchmark_similarity_blocks.py')
e = importlib.util.module_from_spec(spec)
spec.loader.exec_module(e)


def setup_plan(prepared, tmp_path, fraction=0.5, size=1):
    reference = tmp_path / 'reference'
    b.evaluate(arguments(prepared, reference, fraction=fraction))
    args = argparse.Namespace(prepared=prepared, reference=reference, output=tmp_path / 'plan',
                              schemes=['E094'], sample_size=size, seed=17)
    e.plan(args)
    return args


def fit_args(tmp_path):
    return argparse.Namespace(plan=tmp_path / 'plan', output=tmp_path / 'fit', workers=1,
                              chunk_size=3, hardware='synthetic', resume=False)


def replay_args(tmp_path):
    return argparse.Namespace(plan=tmp_path / 'plan', fit=tmp_path / 'fit', output=tmp_path / 'replay',
                              random_repeats=3, seed=42, work_fractions=[.2,.5,1], targets=[1,999])


def test_exact_distinct_sampling_global_exclusion_and_missing_cache(prepared, engine, tmp_path):
    args = setup_plan(prepared, tmp_path)
    with b.readonly(args.output / 'schedule.sqlite') as db:
        assert list(db.execute('SELECT count(*) FROM draws GROUP BY block')) == [(1,), (1,), (1,)]
        assert db.execute('SELECT count(*) FROM validation JOIN fit_mid USING(mid)').fetchone()[0] == 0
        assert db.execute('SELECT count(*) FROM selected WHERE mid=?', ('mid-0',)).fetchone()[0] == 0
        assert db.execute('SELECT count(*) FROM selected').fetchone()[0] == 6
        b.attach(db, args.reference / 'schedule.sqlite', 'ref')
        assert db.execute('SELECT count(*) FROM missing JOIN ref.selected USING(gid)').fetchone()[0] == 0
    report = e.fit(fit_args(tmp_path))
    assert report['new_conformers'] + report['reused_conformers'] == 6
    assert report['new_conformers'] > 0
    with b.readonly(tmp_path / 'fit/fit.sqlite') as db:
        assert all(gid % 3 != 2 for gid, in db.execute('SELECT gid FROM positive'))
    result = e.replay(replay_args(tmp_path))
    assert result['status'] == 'complete'
    rows = list(csv.DictReader((tmp_path / 'replay/yield_at_work.csv').open()))
    final = [r for r in rows if float(r['budget_fraction']) == 1]
    assert len({r['qualified_molecules'] for r in final}) == 1
    assert {r['method'] for r in rows} == {'similarity', 'chemplp', 'equiscore', 'random'}
    assert all(int(r['cold_conformer_work']) == int(r['conformers'])+6 for r in rows if r['method']=='similarity')
    targets = list(csv.DictReader((tmp_path / 'replay/work_to_candidates.csv').open()))
    assert all(r['reached']=='False' and r['conformers']=='' for r in targets if r['target']=='999')


def test_all_cached_negative_outcomes_and_resume(prepared, engine, tmp_path):
    setup_plan(prepared, tmp_path, fraction=1)
    args = fit_args(tmp_path)
    report = e.fit(args)
    assert report['new_conformers'] == 0
    assert report['reused_conformers'] == 6
    with b.readonly(tmp_path / 'plan/schedule.sqlite') as db:
        assert db.execute('SELECT count(*) FROM positive').fetchone()[0] < 6
    with pytest.raises(ValueError, match='fresh output'):
        e.fit(args)
    args.resume = True
    assert e.fit(args) == report
    args.workers = 2
    with pytest.raises(ValueError, match='identical'):
        e.fit(args)


def test_insufficient_distinct_molecules_fails(prepared, engine, tmp_path):
    with pytest.raises(ValueError, match='distinct eligible'):
        setup_plan(prepared, tmp_path, size=500)


def test_deterministic_sampling_and_in_block_positive_counts(prepared, engine, tmp_path):
    args = setup_plan(prepared, tmp_path, fraction=1)
    args.output = tmp_path / 'plan2'
    e.plan(args)
    with b.readonly(tmp_path / 'plan/schedule.sqlite') as a, b.readonly(args.output / 'schedule.sqlite') as c:
        assert list(a.execute('SELECT * FROM draws')) == list(c.execute('SELECT * FROM draws'))
    e.fit(fit_args(tmp_path))
    rows = list(csv.DictReader((tmp_path / 'fit/priorities.csv').open()))
    assert all(r['sampled_molecules']=='1' for r in rows)
    assert all(int(r['qualified_molecules']) <= 1 for r in rows)
    assert [r['block_id'] for r in rows] == sorted(r['block_id'] for r in rows)


def test_complete_blocks_and_molecule_deduplication():
    assert e.curves(['a','empty','b'], {'a':2,'b':3}, {'a':{'m'},'b':{'m','n'}}) == [
        (0,0,0), (1,2,1), (2,2,1), (3,5,2)]


def test_tampered_fit_and_changed_script_rejected(prepared, engine, tmp_path, monkeypatch):
    setup_plan(prepared, tmp_path)
    args = fit_args(tmp_path)
    original = e.own_hash()
    monkeypatch.setattr(e, 'own_hash', lambda: 'changed')
    with pytest.raises(ValueError, match='script changed'):
        e.fit(args)
    monkeypatch.setattr(e, 'own_hash', lambda: original)
    e.fit(args)
    with (tmp_path / 'fit/priorities.csv').open('a') as stream:
        stream.write('tampered')
    with pytest.raises(ValueError, match='changed file'):
        e.replay(replay_args(tmp_path))


def test_in_block_membership_does_not_borrow_other_conformer_positives(prepared, engine, tmp_path):
    reference = tmp_path / 'reference'
    b.evaluate(arguments(prepared, reference))
    args = argparse.Namespace(prepared=prepared, reference=reference, output=tmp_path / 'plan',
                              schemes=['E095'], sample_size=1, seed=17)
    e.plan(args)
    e.fit(fit_args(tmp_path))
    with b.readonly(tmp_path / 'plan/schedule.sqlite') as db:
        expected = dict(db.execute('SELECT d.block,count(DISTINCT CASE WHEN s.gid % 3 != 2 THEN d.mid END) '
            'FROM draws d JOIN selected s ON d.mid=s.mid AND d.block=s.E095 GROUP BY d.block'))
    rows = list(csv.DictReader((tmp_path / 'fit/priorities.csv').open()))
    assert {r['block_id']: int(r['qualified_molecules']) for r in rows} == expected
    assert expected['b2'] == 0


def test_interrupted_fit_resume_and_sealed_schedule(prepared, engine, tmp_path):
    setup_plan(prepared, tmp_path)
    args = fit_args(tmp_path)
    first = e.fit(args)
    b.save(args.output / 'report.json', {'status':'running'})
    args.resume = True
    second = e.fit(args)
    assert second['resumed_execution']
    assert first['output_hashes']['priorities.csv'] == second['output_hashes']['priorities.csv']
