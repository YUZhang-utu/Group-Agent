import csv
import json
import sqlite3

import numpy as np
import pytest

from aidd_agent.final_work_blocks import allocate, build_blocks, members, sha, validate
from aidd_agent.work_block_properties import partition, summarize, refine
from aidd_agent.work_block_delivery import property_members


def row(group, n, pattern='trans;trans'):
    return dict(hard_group=group, conformers=n, unit_count=len(pattern.split(';')), omega_pattern=pattern)


def test_no_small_tail_blocks_and_chirality_labels_preserved():
    classes = [row('R', 12000), row('S', 200), row('cis', 300, 'cis;trans'), row('long', 15, 'trans;trans;trans')]
    blocks = allocate(classes, 5000)
    assert len(blocks) == 1 and blocks[0]['conformers'] == 12515
    assert blocks[0]['groups'] == ['R', 'S', 'cis', 'long']
    assert blocks[0]['kind'] == 'mixed_execution_pool'
    assert not blocks[0]['below_minimum']


def test_small_classes_pool_before_attaching_residuals():
    classes = [row('a', 3000), row('b', 2500), row('c', 9000, 'cis;trans'), row('d', 30)]
    blocks = allocate(classes, 5000)
    assert len(blocks) == 2
    assert sorted(b['conformers'] for b in blocks) == [5530, 9000]
    assert blocks == allocate(list(reversed(classes)), 5000)


def test_small_library_is_one_explicit_exception():
    blocks = allocate([row('a', 30), row('b', 40, 'cis;trans')], 5000)
    assert len(blocks) == 1 and blocks[0]['below_minimum']
    assert allocate([], 5000) == []


def fixture(tmp_path):
    model, diagnostic, candidates, routing = [tmp_path / n for n in ('model', 'diagnostic', 'candidates', 'routing')]
    for path in (model, diagnostic, candidates, routing):
        path.mkdir()
    with sqlite3.connect(model / 'blocks.sqlite') as db:
        db.executescript('CREATE TABLE point(cid TEXT PRIMARY KEY,mid TEXT,group_id TEXT,node TEXT); CREATE INDEX grp ON point(group_id);')
        db.executemany('INSERT INTO point VALUES(?,?,?,?)', [(f'c{i}', f'm{i}', 'definite', 'leaf') for i in range(7)] +
                       [('b0', 'mb0', 'boundary', 'old0'), ('b1', 'mb1', 'boundary', 'old1')])
    def receipt(path, value):
        (path / 'report.json').write_text(json.dumps(value), encoding='utf-8')
    receipt(model, dict(status='complete', counts=dict(conformers=9, blocks=3), output_hashes={'blocks.sqlite': sha(model / 'blocks.sqlite')}))
    with (diagnostic / 'hard_groups.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=['hard_group', 'conformers', 'omega_pattern']); writer.writeheader()
        writer.writerows([dict(hard_group='definite', conformers=7, omega_pattern='trans;trans'),
                         dict(hard_group='boundary', conformers=2, omega_pattern='boundary;trans')])
    receipt(diagnostic, dict(status='complete', provenance={'model_report_sha256': sha(model / 'report.json')},
        output_hashes={'hard_groups.csv': sha(diagnostic / 'hard_groups.csv')}))
    with sqlite3.connect(candidates / 'candidates.sqlite') as db:
        db.execute('CREATE TABLE candidate(cid TEXT PRIMARY KEY,original_group TEXT,maximum_deviation REAL,candidates TEXT)')
        db.executemany('INSERT INTO candidate VALUES(?,?,?,?)',
            [('b0', 'boundary', 35, json.dumps([dict(hard_group='definite', rotation=0)])), ('b1', 'boundary', 60, '[]')])
    receipt(candidates, dict(status='complete', sources={str(diagnostic / 'report.json'): sha(diagnostic / 'report.json')},
        output_hashes={'candidates.sqlite': sha(candidates / 'candidates.sqlite')}))
    with sqlite3.connect(routing / 'boundary_assignments.sqlite') as db:
        db.execute('CREATE TABLE assignment(cid TEXT PRIMARY KEY,original_group TEXT,state TEXT,target_group TEXT,rotation INTEGER,score REAL,evidence TEXT)')
        evidence = json.dumps(dict(accepted=[dict(hard_group='definite', rotation=0, score=.5, prototype_id='c0')]))
        db.executemany('INSERT INTO assignment VALUES(?,?,?,?,?,?,?)', [
            ('b0', 'boundary', 'assigned_provisional', 'definite', 0, .5, evidence),
            ('b1', 'boundary', 'special', None, None, None, '{}')])
    (routing / 'routing_models.json').write_text(json.dumps({'definite': dict(status='routable', prototype_ids=['c0'])}))
    receipt(routing, dict(status='complete', total_conformers=9, provisional_regular_conformers=8,
        boundary_counts=dict(assigned_provisional=1, special=1),
        sources={str(candidates / 'report.json'): sha(candidates / 'report.json')},
        output_hashes={n: sha(routing / n) for n in ['boundary_assignments.sqlite', 'routing_models.json']}))
    return model, diagnostic, candidates, routing


def test_sidecar_full_enumeration_and_input_preservation(tmp_path):
    inputs = fixture(tmp_path); before = sha(inputs[0] / 'blocks.sqlite')
    output = tmp_path / 'result'; report = build_blocks(*inputs, output, 5)
    items = list(members(output))
    assert len(items) == len({r['cid'] for r in items}) == 9
    assert report['regular_conformers'] == 8 and report['special_conformers'] == 1
    adopted = next(r for r in items if r['cid'] == 'b0')
    assert adopted['original_group'] == 'boundary' and adopted['logical_class'] == 'definite'
    assert list(members(output, 'special-exhaustive'))[0]['cid'] == 'b1'
    assert validate(output)['structural_gate'] == 'passed'
    assert before == sha(inputs[0] / 'blocks.sqlite')


@pytest.mark.parametrize('change', ['cid', 'rotation', 'missing', 'unknown_target', 'score', 'prototype'])
def test_invalid_assignment_rejected_even_with_updated_receipt_hash(tmp_path, change):
    inputs = fixture(tmp_path); routing = inputs[-1]
    with sqlite3.connect(routing / 'boundary_assignments.sqlite') as db:
        sql = {'cid': "UPDATE assignment SET cid='unknown' WHERE cid='b0'",
               'rotation': "UPDATE assignment SET rotation=1 WHERE cid='b0'",
               'missing': "DELETE FROM assignment WHERE cid='b1'",
               'unknown_target': "UPDATE assignment SET target_group='unknown' WHERE cid='b0'",
               'score': "UPDATE assignment SET score=1.5 WHERE cid='b0'",
               'prototype': "UPDATE assignment SET evidence=replace(evidence,'c0','bad') WHERE cid='b0'"}[change]
        db.execute(sql)
    report = json.loads((routing / 'report.json').read_text())
    report['output_hashes']['boundary_assignments.sqlite'] = sha(routing / 'boundary_assignments.sqlite')
    (routing / 'report.json').write_text(json.dumps(report))
    with pytest.raises(ValueError):
        build_blocks(*inputs, tmp_path / 'result', 5)
    assert not (tmp_path / 'result/report.json').exists()


def test_hash_change_is_rejected(tmp_path):
    inputs = fixture(tmp_path)
    with (inputs[0] / 'blocks.sqlite').open('ab') as stream:
        stream.write(b'changed')
    with pytest.raises(ValueError, match='hash changed'):
        build_blocks(*inputs, tmp_path / 'result', 5)


def test_property_refinement_does_not_force_small_or_identical_splits():
    assert len(partition(np.zeros((120, 69)), 20)[0]) == 1
    assert len(partition(np.vstack([np.zeros((19, 69)), np.ones((81, 69)) * 4]), 20)[0]) == 1
    leaves, evidence = partition(np.vstack([np.zeros((40, 69)), np.ones((60, 69)) * 4]), 20)
    assert sorted(map(len, leaves)) == [40, 60] and evidence[0]['centroid_contrast'] == pytest.approx(4)
    assert sorted(np.concatenate(leaves).tolist()) == list(range(100))
    leaves, _ = partition(np.arange(1200).reshape(200, 6), 30, 3, .1)
    assert len(leaves) <= 3 and min(map(len, leaves)) >= 30


def test_property_summary_retains_local_bulk_and_not_sequence_labels():
    x = np.zeros((4, 29)); x[0, 12] = 2
    y = summarize(x)
    assert len(y) == 69 and y[12] == .5 and y[46 + 12] == 2
    np.testing.assert_array_equal(y, summarize(np.roll(x, 1, axis=0)))
    with pytest.raises(ValueError):
        summarize(np.zeros((4, 10)))


def test_property_delivery_covers_all_regular_and_retains_special(tmp_path):
    inputs = fixture(tmp_path); blocks = tmp_path / 'work'
    build_blocks(*inputs, blocks, 3)
    profiles = tmp_path / 'profiles'; profiles.mkdir()
    with sqlite3.connect(profiles / 'profiles.sqlite') as db:
        db.execute('CREATE TABLE item(cid TEXT PRIMARY KEY,parent TEXT,vector BLOB)')
        for index, item in enumerate(r for r in members(blocks) if r['state'] != 'special'):
            db.execute('INSERT INTO item VALUES(?,?,?)', (item['cid'], item['block_id'],
                np.full(69, 0 if index < 4 else 4, dtype='<f4').tobytes()))
    profile_report = dict(status='complete', blocks=str(blocks), coverage=1,
        signature=dict(blocks_report_sha256=sha(blocks / 'report.json')),
        output_hashes={'profiles.sqlite': sha(profiles / 'profiles.sqlite')})
    (profiles / 'report.json').write_text(json.dumps(profile_report))
    result = refine(profiles, tmp_path / 'properties')
    assert result['regular_work_blocks'] == 2 and result['smallest_regular_block'] == 4
    assert validate(blocks, tmp_path / 'properties')['property_gate'] == 'passed'
    with sqlite3.connect(tmp_path / 'properties/property_blocks.sqlite') as db:
        ids = [r[0] for r in db.execute('SELECT block_id FROM block')]
    exported = [r for bid in ids for r in property_members(blocks, tmp_path / 'properties', bid)]
    assert len(exported) == len({r['cid'] for r in exported}) == 8
    assert list(property_members(blocks, tmp_path / 'properties', 'special-exhaustive'))[0]['cid'] == 'b1'
    assert next(r for r in exported if r['cid'] == 'b0')['original_group'] == 'boundary'


def test_property_incomplete_coverage_cannot_refine(tmp_path):
    profiles = tmp_path / 'profiles'; profiles.mkdir()
    with sqlite3.connect(profiles / 'profiles.sqlite') as db:
        db.execute('CREATE TABLE item(cid TEXT PRIMARY KEY)')
    inputs = fixture(tmp_path); blocks = tmp_path / 'work'; build_blocks(*inputs, blocks, 3)
    report = dict(status='complete', blocks=str(blocks), coverage=.9,
        signature=dict(blocks_report_sha256=sha(blocks / 'report.json')),
        output_hashes={'profiles.sqlite': sha(profiles / 'profiles.sqlite')})
    (profiles / 'report.json').write_text(json.dumps(report))
    with pytest.raises(ValueError, match='Complete matching'):
        refine(profiles, tmp_path / 'bad')
