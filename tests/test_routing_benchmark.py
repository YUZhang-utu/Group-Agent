"""Synthetic benchmark acceptance; no measured MDM2 scientific result."""
import argparse
import json
from pathlib import Path
import sqlite3
import sys

import pytest

from aidd_agent import routing_benchmark as benchmark


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    benchmark.save(path, value)


@pytest.fixture
def historical_inputs(tmp_path):
    registry = tmp_path / 'registry.sqlite3'
    with sqlite3.connect(registry) as connection:
        connection.execute('CREATE TABLE conformer(id TEXT PRIMARY KEY,molecule_id TEXT,content_sha256 TEXT)')
        connection.executemany('INSERT INTO conformer VALUES(?,?,?)',
                               [(f'art-{index}', f'mol-{index//2}', f'hash-{index}') for index in range(8)])
    reference = tmp_path / 'reference'
    reference.mkdir()
    benchmark.table(reference / 'selected-conformers.csv', [dict(rank=index+1, conformer_id=f'art-{index}',
                    molecule_id=f'mol-{index//2}', global_id=index, source_sha256=f'hash-{index}') for index in range(6)])
    write_json(reference / 'report.json', dict(status='complete', exported_conformers=6, distinct_molecules=3,
               outputs={'selected-conformers.csv': benchmark.sha(reference / 'selected-conformers.csv')}))
    runs = tmp_path / 'runs'
    for run, (scheme, method, limit) in zip(benchmark.RUNS,
            [(scheme, method, limit) for scheme in benchmark.SCHEMES for method in benchmark.METHODS for limit in (5, 10)]):
        folder = runs / ('PROMPT-' + run) / 'execution/campaign/blocks'
        selected = folder / 'selected'
        selected.mkdir(parents=True)
        indices = [0, 2, 6] if limit == 5 else [0, 1, 2, 6]
        records = [dict(rank=rank, artifact_cid=f'art-{index}', cid=f'structural-{index}',
                        molecule_id=f'local-mol-{index//2}', gid=index, source_sha256=f'hash-{index}')
                   for rank, index in enumerate(indices, 1)]
        (selected / 'candidates.jsonl').write_text(''.join(json.dumps(row)+'\n' for row in records))
        write_json(selected / 'report.json', dict(status='complete', unique_conformers=len(indices), unique_molecules=3,
                   output_hashes={'candidates.jsonl': benchmark.sha(selected / 'candidates.jsonl')}))
        write_json(folder / 'protocol.json', dict(library_registry=dict(path=str(registry))))
        write_json(folder / 'query/report.json', dict(templates=['query']))
        with sqlite3.connect(folder / 'members.sqlite') as connection:
            connection.execute('CREATE TABLE member(artifact_cid TEXT,gid INTEGER,block_id TEXT)')
            connection.executemany('INSERT INTO member VALUES(?,?,?)',
                                   [(f'art-{index}', index, 'block') for index in range(8)])
        with sqlite3.connect(folder / 'ranking.sqlite') as connection:
            connection.execute('CREATE TABLE ranking(gid INTEGER)')
            connection.executemany('INSERT INTO ranking VALUES(?)', [(index,) for index in indices])
        write_json(folder / 'report.json', dict(status='complete', kind='block_search_dock',
                   search_engine='consensus-threshold-v1', docking_requested=False,
                   selection=dict(scheme=scheme, method=method, blocks_per_method=limit),
                   requested_conformers=4, exported_conformers=len(indices), unique_molecules=3,
                   shortfall=4-len(indices), ranked_conformers=len(indices), searched_conformers=8,
                   searched_blocks=[dict(block_id='block', searched_conformers=8, cache_reused=limit == 10)],
                   output_hashes={name: benchmark.sha(folder / name)
                                  for name in ('protocol.json', 'query/report.json', 'selected/report.json', 'members.sqlite', 'ranking.sqlite')}))
    return argparse.Namespace(reference=reference, runs_root=runs, registry=None, run_ids=None,
                              top_k=[1, 2, 6, 1000], budgets=[1, 4], skip_stages=False, output=tmp_path / 'analysis')


def test_historical_namespaces_and_shortfalls(historical_inputs):
    result = benchmark.historical(historical_inputs)
    assert result['status'] == 'complete'
    assert result['common_export_prefix'] == 3
    rows = list(benchmark.csv_rows(historical_inputs.output / 'historical_recovery.csv'))
    match = next(row for row in rows if row['arm'] == 'E094/chemplp/Top-5'
                 and row['unit'] == 'conformer' and row['reference_k'] == '2' and row['full_output'] == 'True')
    assert float(match['reference_recall']) == 0.5
    assert not any(row['reference_k'] == '1000' for row in rows)
    summaries = list(benchmark.csv_rows(historical_inputs.output / 'arm_summary.csv'))
    assert len(summaries) == 12
    assert all(row['outside_historical_reference'] == '1' for row in summaries)
    assert (historical_inputs.output / 'recovery_work.svg').exists()
    stages = list(benchmark.csv_rows(historical_inputs.output / 'reference_stage_coverage.csv'))
    assert stages[0]['outside_selected_blocks'] == '0'
    with pytest.raises(FileExistsError):
        benchmark.historical(historical_inputs)


def test_reference_hash_corruption(historical_inputs):
    with (historical_inputs.reference / 'selected-conformers.csv').open('a') as stream:
        stream.write('\n')
    with pytest.raises(ValueError, match='changed file'):
        benchmark.historical(historical_inputs)
    assert benchmark.read(historical_inputs.output / 'audit.json')['status'] == 'incomplete'


def test_stage_database_hash_corruption(historical_inputs):
    folder = historical_inputs.runs_root / ('PROMPT-'+benchmark.RUNS[0]) / 'execution/campaign/blocks'
    with sqlite3.connect(folder / 'ranking.sqlite') as connection:
        connection.execute('DELETE FROM ranking WHERE gid=0')
    with pytest.raises(ValueError, match='changed file'):
        benchmark.historical(historical_inputs)


def test_skip_stage_analysis_explicit(historical_inputs):
    historical_inputs.skip_stages = True
    benchmark.historical(historical_inputs)
    assert not (historical_inputs.output / 'reference_stage_coverage.csv').exists()


def test_registry_hash_mismatch(historical_inputs):
    registry_path = benchmark.read(historical_inputs.runs_root / ('PROMPT-'+benchmark.RUNS[0]) /
                                   'execution/campaign/blocks/protocol.json')['library_registry']['path']
    with sqlite3.connect(registry_path) as connection:
        connection.execute("UPDATE conformer SET content_sha256='changed' WHERE id='art-0'")
    with pytest.raises(ValueError, match='Source hash mismatch'):
        benchmark.historical(historical_inputs)


def test_bad_reference_rank(historical_inputs):
    path = historical_inputs.reference / 'selected-conformers.csv'
    rows = list(benchmark.csv_rows(path))
    rows[0]['rank'] = '2'
    benchmark.table(path, rows)
    report = benchmark.read(historical_inputs.reference / 'report.json')
    report['outputs'][path.name] = benchmark.sha(path)
    benchmark.save(historical_inputs.reference / 'report.json', report)
    with pytest.raises(ValueError, match='contiguous'):
        benchmark.historical(historical_inputs)


@pytest.fixture
def replay_inputs(tmp_path):
    universe, priorities = tmp_path / 'universe.csv', tmp_path / 'priorities.csv'
    rows = [dict(conformer_id=f'cid-{index}', molecule_id=f'mid-{index//2}', qualified=int(index < 4),
                 reference_rank=index+1 if index < 4 else '',
                 **{scheme: f'block-{index//4}' for scheme in benchmark.SCHEMES}) for index in range(8)]
    benchmark.table(universe, rows)
    benchmark.table(priorities, [dict(scheme=scheme, scorer=method, rank=index+1, block_id=f'block-{index}')
                                 for scheme in benchmark.SCHEMES for method in benchmark.METHODS for index in range(2)])
    manifest = tmp_path / 'manifest.json'
    benchmark.save(manifest, dict(scope='complete_fixed_sublibrary', conformers=8,
                   library_identity='fixture', query_identity='fixture-query', evaluation_policy='frozen',
                   priority_fit_exclusion='fixture disjoint molecules', universe_sha256=benchmark.sha(universe),
                   priorities_sha256=benchmark.sha(priorities)))
    return argparse.Namespace(output=tmp_path / 'replay', manifest=manifest, universe=universe,
                              priorities=priorities, top_k=[2], seed=42, repeats=3, uniform_chunk=2, target_recall=0.9)


def test_replay_random_controls_and_whole_blocks(replay_inputs):
    benchmark.replay(replay_inputs)
    summaries = list(benchmark.csv_rows(replay_inputs.output / 'work_to_recall.csv'))
    assert len(summaries) == (6 + 4*3)*2
    assert all(row['work_to_target'] == '4' for row in summaries if row['method'] == 'E094/chemplp')
    original = (replay_inputs.output / 'routing_curves.csv').read_bytes()
    replay_inputs.output = replay_inputs.output.parent / 'replay-again'
    benchmark.replay(replay_inputs)
    assert original == (replay_inputs.output / 'routing_curves.csv').read_bytes()


def test_replay_rejects_selected_union(replay_inputs):
    manifest = benchmark.read(replay_inputs.manifest)
    manifest['scope'] = 'selected_union'
    benchmark.save(replay_inputs.manifest, manifest)
    with pytest.raises(ValueError, match='complete evaluated'):
        benchmark.replay(replay_inputs)


def test_replay_rejects_missing_priority(replay_inputs):
    rows = list(benchmark.csv_rows(replay_inputs.priorities))[:-1]
    benchmark.table(replay_inputs.priorities, rows)
    manifest = benchmark.read(replay_inputs.manifest)
    manifest['priorities_sha256'] = benchmark.sha(replay_inputs.priorities)
    benchmark.save(replay_inputs.manifest, manifest)
    with pytest.raises(ValueError, match='Unranked'):
        benchmark.replay(replay_inputs)


def test_universe_rejected_rank_invalid(tmp_path):
    path = tmp_path / 'universe.csv'
    benchmark.table(path, [dict(conformer_id='cid', molecule_id='mid', qualified=0, reference_rank=1,
                               E094='a', E095='b', E096='c')])
    with pytest.raises(ValueError, match='qualified candidate'):
        benchmark.load_universe(path)


@pytest.mark.parametrize('cached', [False, True])
def test_measure_actual_command_and_cache_rejection(tmp_path, cached):
    receipt = tmp_path / 'result.json'
    spec = tmp_path / 'spec.json'
    payload = dict(status='complete', cache_reused=cached, resumed_execution=False)
    command = [sys.executable, '-c', f'from pathlib import Path; Path({str(receipt)!r}).write_text({json.dumps(payload)!r})']
    benchmark.save(spec, dict(method='test', comparison_id='fixture', hardware='local', workers=1,
                             protocol_identity='fixture', repeat=1, command=command, cwd=str(tmp_path),
                             result_receipt=str(receipt), phase='online', cache_policy='fresh_application_cache'))
    args = argparse.Namespace(spec=spec, output=tmp_path / 'measurement')
    if cached:
        with pytest.raises(ValueError, match='ineligible'):
            benchmark.measure(args)
    else:
        assert benchmark.measure(args)['wall_seconds'] > 0
    assert benchmark.read(args.output / 'audit.json')['eligible_fresh_execution'] is not cached


def test_timing_summary_requires_matched_repeats(tmp_path):
    receipts = []
    for method, seconds in [('baseline', 20), ('blocked', 10)]:
        for repeat in range(3):
            path = tmp_path / f'{method}-{repeat}.json'
            benchmark.save(path, dict(status='complete', eligible_fresh_execution=True, wall_seconds=seconds,
                           spec=dict(method=method, repeat=repeat, comparison_id='fixture', hardware='same',
                                     workers=1, protocol_identity='policy', phase='online', cache_policy='fresh_application_cache')))
            receipts.append(path)
    args = argparse.Namespace(output=tmp_path / 'summary', receipts=receipts, baseline='baseline')
    benchmark.timing_summary(args)
    rows = list(benchmark.csv_rows(args.output / 'timing_summary.csv'))
    assert rows[1]['median_paired_time_ratio'] == '2.0'
    args.output = tmp_path / 'bad-summary'
    args.receipts = receipts[:-1]
    with pytest.raises(ValueError, match='matched repeats'):
        benchmark.timing_summary(args)


def test_cli_negative_cutoff():
    with pytest.raises(ValueError, match='positive'):
        benchmark.main(['historical', '--reference', '.', '--runs-root', '.', '--output', 'unused', '--top-k', '-1'])
