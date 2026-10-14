"""Synthetic matched-pipeline tests; no actual MDM2 screening is executed."""
import argparse
import json
from pathlib import Path
import sqlite3
import sys
from types import SimpleNamespace

import numpy as np
import pytest

from aidd_agent import matched_block_benchmark as benchmark


def test_prepare_common_mapping_and_molecule_training_exclusion(tmp_path, monkeypatch):
    campaign, analysis, sample, batch = [tmp_path / name for name in ('campaign', 'analysis', 'sample', 'batch')]
    for folder in (campaign / 'query', analysis / 'chemplp_baseline', sample, batch):
        folder.mkdir(parents=True)
    paths = {name: tmp_path / name for name in ('profiles', 'descriptors', *benchmark.SCHEMES)}
    for folder in paths.values():
        folder.mkdir()
    for scheme in ('E095', 'E096'):
        with sqlite3.connect(paths[scheme] / 'property_blocks.sqlite') as database:
            database.execute('CREATE TABLE membership(cid TEXT PRIMARY KEY,block_id TEXT)')
            database.executemany('INSERT INTO membership VALUES(?,?)', [('s0', 'b0'), ('s1', 'b0'), ('s2', 'b1')])
    with sqlite3.connect(batch / 'registry.sqlite3') as database:
        database.execute('CREATE TABLE conformer(id TEXT PRIMARY KEY,molecule_id TEXT)')
        database.executemany('INSERT INTO conformer VALUES(?,?)', [('a0', 'm0'), ('a1', 'm0'), ('a2', 'm1')])
    with sqlite3.connect(sample / 'samples.sqlite') as database:
        database.execute('CREATE TABLE population(scheme TEXT,block_id TEXT,n INTEGER)')
        database.executemany('INSERT INTO population VALUES(?,?,?)',
                             [(scheme, block, count) for scheme in benchmark.SCHEMES for block, count in [('b0', 2), ('b1', 1)]])
        database.execute('CREATE TABLE selected(mid TEXT)')
        database.execute("INSERT INTO selected VALUES('source-m0')")
    def snapshot(path):
        path = Path(path).resolve()
        return dict(path=str(path), size=path.stat().st_size, mtime_ns=path.stat().st_mtime_ns)
    benchmark.save(sample / 'signature.json', dict(database_snapshots=[snapshot(paths[scheme] / 'property_blocks.sqlite')
                                                                      for scheme in ('E095', 'E096')]))
    benchmark.save(sample / 'report.json', dict(status='complete'))
    benchmark.save(analysis / 'source.json', dict(sample_report=str(sample / 'report.json')))
    for folder in (analysis, analysis / 'chemplp_baseline'):
        with sqlite3.connect(folder / 'block_ranking.sqlite') as database:
            database.execute('CREATE TABLE rankings(scheme TEXT,block_id TEXT,receptor TEXT,top_n INTEGER,'
                             'rank INTEGER,status TEXT,population INTEGER)')
            database.executemany('INSERT INTO rankings VALUES(?,?,?,?,?,?,?)',
                                 [(scheme, block, '7NA2_A', 10, rank, 'eligible', count)
                                  for scheme in benchmark.SCHEMES for block, rank, count in [('b0', 1, 2), ('b1', 2, 1)]])
        benchmark.save(folder / 'report.json', dict(source_report=str(analysis / 'source.json'),
                       output_hashes={'block_ranking.sqlite': benchmark.sha(folder / 'block_ranking.sqlite')}))
    profile = tmp_path / 'sampling.json'
    benchmark.save(profile, {name: str(path) for name, path in paths.items()})
    benchmark.save(campaign / 'query/report.json', dict(readiness='ready_for_consensus_funnel'))
    benchmark.save(campaign / 'protocol.json', dict(request={'analysis': str(analysis / 'report.json')},
                   inputs={str(profile): benchmark.sha(profile), str(analysis / 'report.json'): benchmark.sha(analysis / 'report.json')},
                   library_registry=snapshot(batch / 'registry.sqlite3')))
    benchmark.seal(campaign, dict(status='complete', search_engine='consensus-threshold-v1'), ['protocol.json', 'query/report.json'])
    def members(profile, selection, report, output):
        assert selection == dict(scheme='E094', block_ids=['b0', 'b1'])
        with sqlite3.connect(output) as database:
            database.execute('CREATE TABLE member(cid TEXT,mid TEXT,block_id TEXT,artifact_cid TEXT,gid INTEGER)')
            database.executemany('INSERT INTO member VALUES(?,?,?,?,?)',
                                 [('s0', 'source-m0', 'b0', 'a0', 0), ('s1', 'source-m0', 'b0', 'a1', 1),
                                  ('s2', 'source-m1', 'b1', 'a2', 2)])
    monkeypatch.setitem(sys.modules, 'aidd_agent.block_sampling', SimpleNamespace(snapshot=snapshot))
    monkeypatch.setitem(sys.modules, 'aidd_agent.block_conformer_search',
                        SimpleNamespace(members=members, map_artifacts=lambda *args: None, profile_paths=lambda profile: paths))
    result = benchmark.prepare(argparse.Namespace(campaign=campaign, output=tmp_path / 'prepared', seed=42, receptor='7NA2_A'))
    assert result['regular_conformers'] == 3
    assert result['priority_fit_conformers'] == 2
    assert result['heldout_conformers'] == 1
    assert len(benchmark.load_priorities(tmp_path / 'prepared')) == 6
    with benchmark.readonly(tmp_path / 'prepared/universe.sqlite') as database:
        assert list(database.execute('SELECT excluded FROM member ORDER BY gid')) == [(1,), (1,), (0,)]


def arguments(prepared, output, **changes):
    values = dict(prepared=prepared, output=output, fraction=1.0, include_training=False,
                  mode='full', scheme=None, scorer=None, top_blocks='20', workers=1,
                  chunk_size=3, output_budget=4, hardware='synthetic-test', repeat=1, resume=False)
    values.update(changes)
    return argparse.Namespace(**values)


@pytest.fixture
def prepared(tmp_path):
    folder = tmp_path / 'prepared'
    folder.mkdir()
    with sqlite3.connect(folder / 'universe.sqlite') as database:
        database.execute('CREATE TABLE member(gid INTEGER,artifact_cid TEXT,canonical_mid TEXT,block_id TEXT,'
                         'E095 TEXT,E096 TEXT,sample_key INTEGER,excluded INTEGER)')
        database.executemany('INSERT INTO member VALUES(?,?,?,?,?,?,?,?)',
                             [(index, f'cid-{index}', f'mid-{index//2}', f'b{index//4}',
                               f'b{index%3}', f'b{(index//2)%3}', index//2*(1 << 60), int(index < 2))
                              for index in range(12)])
    priorities = [dict(scheme=scheme, scorer=scorer, rank=index+1, block_id=f'b{index}')
                  for scheme in benchmark.SCHEMES for scorer in benchmark.SCORERS for index in range(3)]
    benchmark.table(folder / 'priorities.csv', priorities)
    benchmark.save(folder / 'query.json', dict(readiness='ready_for_consensus_funnel', sources={},
                   design={}, anchors=[], templates=[dict(query_id='q1'), dict(query_id='q2')]))
    benchmark.seal(folder, dict(status='complete', kind='matched_benchmark_preparation', batch=str(tmp_path),
                   query_sha256=benchmark.sha(folder / 'query.json'), code=benchmark.code_identity(),
                   registry_snapshot={'fixture': True}), ['universe.sqlite', 'query.json', 'priorities.csv'])
    return folder


@pytest.fixture
def engine(monkeypatch):
    def refine(batch, output, query, ids, workers, chunk_size):
        assert query['block_search_mode'] == 'threshold'
        assert query['assignment_backend'] == 'numpy'
        assert '_block_cache' not in query
        for ordinal, template in enumerate(query['templates']):
            folder = output / 'chunks' / f'{ordinal:02d}'
            folder.mkdir(parents=True, exist_ok=True)
            for start in range(0, len(ids), chunk_size):
                path = folder / f'{start:010d}.npz'
                subset = ids[start:start+chunk_size]
                stages = np.array([7 if index % 3 != 2 else 2 for index in subset])
                np.savez(path, global_ids=subset, stage_levels=stages)
                rows = [dict(global_id=int(index), optional_score=float(12-index), gaussian_same_pose=0.8,
                             composite_score=0.6, conformer_id=f'cid-{index}', molecule_id=f'mid-{index//2}')
                        for index, stage in zip(subset, stages) if stage == 7]
                poses = path.with_suffix('.poses.jsonl')
                poses.write_text(''.join(json.dumps(row)+'\n' for row in rows))
                benchmark.save(path.with_suffix('.receipt.json'),
                               dict(files={str(item): benchmark.sha(item) for item in (path, poses)}))
        return {'fixture': True}
    monkeypatch.setitem(sys.modules, 'aidd_agent.budget_screen', SimpleNamespace(run_refinement=refine))
    monkeypatch.setitem(sys.modules, 'aidd_agent.screening_selection', SimpleNamespace(check_hashes=lambda values: None))
    monkeypatch.setitem(sys.modules, 'aidd_agent.block_sampling', SimpleNamespace(snapshot=lambda path: {'fixture': True}))
    monkeypatch.setitem(sys.modules, 'aidd_agent.contact_policy', SimpleNamespace(require_protein_contacts=lambda *args: None))
    return refine


def test_nested_molecule_scales_and_training_exclusion(prepared, tmp_path):
    sets = []
    for index, fraction in enumerate((0.25, 0.5, 1.0)):
        with sqlite3.connect((tmp_path / f'scale{index}.sqlite').as_uri(), uri=True) as connection:
            benchmark.select_members(connection, prepared, fraction, False)
            rows = list(connection.execute('SELECT gid,mid FROM selected'))
            assert all(mid != 'mid-0' for _, mid in rows)
            assert all(sum(candidate == mid for _, candidate in rows) == 2 for _, mid in rows)
            sets.append({gid for gid, _ in rows})
    assert sets[0] < sets[1] < sets[2]


def test_full_and_all_blocks_exact_equivalence(prepared, engine, tmp_path):
    full, blocked = tmp_path / 'full', tmp_path / 'blocks'
    result = benchmark.evaluate(arguments(prepared, full))
    second = benchmark.evaluate(arguments(prepared, blocked, mode='blocks', scheme='E096', scorer='chemplp', top_blocks='all'))
    assert result['input_conformers'] == 10
    assert result['qualified_conformers'] == 6
    assert result['semantic_pose_sha256'] == second['semantic_pose_sha256']
    assert result['timing_eligible']
    verified = benchmark.verify_equivalence(argparse.Namespace(full=full, blocks=blocked, output=tmp_path / 'verified'))
    assert verified['equivalence'] == 'passed'
    prefix = benchmark.evaluate(arguments(prepared, tmp_path / 'prefix', mode='blocks', scheme='E094', scorer='chemplp', top_blocks='1'))
    assert prefix['input_conformers'] == 2
    assert prefix['shortfall'] == 3


def test_no_overwrite_and_changed_resume_settings(prepared, engine, tmp_path):
    output = tmp_path / 'run'
    args = arguments(prepared, output)
    first = benchmark.evaluate(args)
    with pytest.raises(ValueError, match='fresh output'):
        benchmark.evaluate(args)
    args.resume = True
    assert benchmark.evaluate(args)['semantic_pose_sha256'] == first['semantic_pose_sha256']
    args.workers = 2
    with pytest.raises(ValueError, match='identical'):
        benchmark.evaluate(args)


def test_resume_execution_not_fresh_timing(prepared, engine, tmp_path):
    output = tmp_path / 'run'
    args = arguments(prepared, output)
    benchmark.evaluate(args)
    benchmark.save(output / 'report.json', dict(status='running'))
    args.resume = True
    report = benchmark.evaluate(args)
    assert report['resumed_execution'] and not report['timing_eligible']


def test_changed_prepared_input_fails(prepared, engine, tmp_path):
    with (prepared / 'priorities.csv').open('a') as stream:
        stream.write('\n')
    with pytest.raises(ValueError, match='changed file'):
        benchmark.evaluate(arguments(prepared, tmp_path / 'run'))


def test_collect_rejects_missing_pass_pose(prepared, engine, tmp_path):
    output = tmp_path / 'run'
    benchmark.evaluate(arguments(prepared, output))
    path = output / 'chunks/00/0000000000.npz'
    poses = path.with_suffix('.poses.jsonl')
    poses.write_text('')
    benchmark.save(path.with_suffix('.receipt.json'), dict(files={str(item): benchmark.sha(item) for item in (path, poses)}))
    with pytest.raises(ValueError, match='missing pose'):
        benchmark.collect_scores(output, np.arange(2, 12), benchmark.read(prepared / 'query.json'), 3)


def replay_args(prepared, reference, output):
    return argparse.Namespace(prepared=prepared, reference=reference, output=output,
                              fractions=[0.5, 1.0], top_k=[1, 2], top_blocks=[1, 2, 20],
                              target_recall=[0.9, 0.95], random_repeats=3, seed=42,
                              rerank_exports=True, output_budget=4)


def test_replay_all_schemes_random_scales_and_local_reranking(prepared, engine, tmp_path):
    reference = tmp_path / 'full'
    benchmark.evaluate(arguments(prepared, reference))
    args = replay_args(prepared, reference, tmp_path / 'replay')
    report = benchmark.replay(args)
    assert report['status'] == 'complete'
    with (args.output / 'routing_curves.csv').open() as stream:
        rows = list(__import__('csv').DictReader(stream))
    assert {row['scorer'] for row in rows} == {'chemplp', 'equiscore', 'random'}
    assert all(float(row['recall']) == 1 for row in rows if row['blocks'] == '3')
    original = (args.output / 'routing_curves.csv').read_bytes()
    args.output = tmp_path / 'replay-again'
    benchmark.replay(args)
    assert original == (args.output / 'routing_curves.csv').read_bytes()
    with (args.output / 'delivered_recall.csv').open() as stream:
        delivered = list(__import__('csv').DictReader(stream))
    assert all(float(row['recall']) == 1 for row in delivered if row['blocks'] == '3')


def test_replay_rejects_unmeasured_larger_scale(prepared, engine, tmp_path):
    reference = tmp_path / 'half'
    benchmark.evaluate(arguments(prepared, reference, fraction=0.5))
    with pytest.raises(ValueError, match='exceeds'):
        benchmark.replay(replay_args(prepared, reference, tmp_path / 'replay'))


def test_ranking_positive_tier_and_local_rrf():
    with sqlite3.connect(':memory:') as connection:
        connection.execute('CREATE TABLE poses(gid INTEGER,template TEXT,contact REAL,gaussian REAL)')
        connection.executemany('INSERT INTO poses VALUES(?,?,?,?)', [(1,'a',1,0.5), (2,'a',0,1),
                               (2,'b',0,1), (3,'a',1,0.7)])
        benchmark.aggregate(connection)
        assert [gid for gid, in connection.execute('SELECT gid FROM ranking ORDER BY rank')] == [3,1,2]


def test_timing_gates_recall_output_budget_and_repeats(prepared, engine, tmp_path):
    folders = []
    for repeat in (1,2,3):
        for mode, top in [('full','20'), ('blocks','all'), ('blocks','1')]:
            folder = tmp_path / f'{mode}-{top}-{repeat}'
            report = benchmark.evaluate(arguments(prepared, folder, repeat=repeat, mode=mode,
                                        scheme='E094' if mode=='blocks' else None,
                                        scorer='chemplp' if mode=='blocks' else None, top_blocks=top))
            report['total_wall_seconds'] = 20 if mode == 'full' else 10
            benchmark.save(folder / 'report.json', report)
            folders.append(folder)
    args = argparse.Namespace(runs=folders, output=tmp_path / 'timings', reference_k=2, target_recall=0.95)
    benchmark.timing(args)
    with (args.output / 'timing_comparison.csv').open() as stream:
        rows = list(__import__('csv').DictReader(stream))
    assert next(row for row in rows if row['method'].endswith('all'))['accepted_speedup'] == '2.0'
    assert next(row for row in rows if row['method'].endswith('1'))['accepted_speedup'] == ''
    args.output = tmp_path / 'unpaired'
    args.runs = folders[:-1]
    with pytest.raises(ValueError, match='matched repeats'):
        benchmark.timing(args)


def test_priority_duplicates_rejected(prepared):
    with (prepared / 'priorities.csv').open('a') as stream:
        stream.write('E094,chemplp,1,b0\n')
    with pytest.raises(ValueError, match='Invalid priority'):
        benchmark.load_priorities(prepared)


def test_fraction_and_seed_key():
    assert benchmark.molecule_key('same', 42) == benchmark.molecule_key('same', 42)
    assert benchmark.molecule_key('same', 42) != benchmark.molecule_key('same', 43)
    assert benchmark.fraction_rule(1.0) == ('1', ())
    for fraction in (0, -1, 1.1, float('nan')):
        with pytest.raises(ValueError):
            benchmark.fraction_rule(fraction)
