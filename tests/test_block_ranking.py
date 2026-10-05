"""Top-tail ranking, molecule deduplication and saved-result inspection."""
import csv
import json
import sqlite3

import pytest

from aidd_agent.block_ranking import add_rankings, inspect, options
from aidd_agent.final_work_blocks import sha
from aidd_agent.joint_spatial_profiles import save


def panel(tmp_path, groups, status='complete'):
    sample = tmp_path / 'sample'
    sample.mkdir()
    source = tmp_path / 'run.json'
    save(source, dict(score='ChemPLP TOTAL_SCORE (lower is better)', status=status,
                     sample_report=str(sample / 'report.json')))
    out = tmp_path / 'analysis'
    out.mkdir()
    rows, summaries = [], []
    with sqlite3.connect(sample / 'samples.sqlite') as db:
        db.execute('CREATE TABLE selected(cid TEXT PRIMARY KEY,mid TEXT)')
        for bid, values in groups.items():
            for i, (mid, score) in enumerate(values):
                cid = f'{bid}-{i}'
                db.execute('INSERT INTO selected VALUES(?,?)', (cid, mid))
                rows.append(dict(scheme='E094', block_id=bid, receptor='R1', cid=cid,
                    status='ok' if score is not None else 'failed', score=score,
                    pose_name=cid, pose_file=f'{bid}.mol2', pose_index=i))
            scored = sum(s is not None for _, s in values)
            summaries.append(dict(scheme='E094', block_id=bid, receptor='R1',
                sampled=len(values), scored=scored, missing=len(values)-scored, population=1000))
    for name, values in [('block_scores.csv', rows), ('block_summary.csv', summaries)]:
        with (out / name).open('w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=list(values[0]))
            writer.writeheader(); writer.writerows(values)
    save(out / 'signature.json', dict(sha256=sha(source)))
    save(out / 'report.json', dict(kind='block_analyze', output_hashes={
        name: sha(out / name) for name in ['block_scores.csv', 'block_summary.csv']}))
    return source, out


def test_top_ten_overrides_whole_block_mean_and_exposes_top_five(tmp_path):
    # A has a much worse whole-block mean but a better leading tail than B.
    a = [(f'a{i}', -100.0 if i < 10 else 100.0) for i in range(100)]
    b = [(f'b{i}', -20.0) for i in range(100)]
    source, out = panel(tmp_path, {'A': a, 'B': b})
    result = add_rankings(source, out)
    assert result['ranking']['leaders'][0]['block_id'] == 'A'
    ranks = inspect(out / 'report.json', scheme='E094', receptor='R1')
    assert [r['top_n_mean'] for r in ranks['blocks']] == [-100, -20]
    top = inspect(out / 'report.json', operation='top', scheme='E094', receptor='R1', block_id='A')
    assert len(top['candidates']) == 5
    assert len({r['molecule_id'] for r in top['candidates']}) == 5
    assert all(r['pose_file'] == 'A.mol2' for r in top['candidates'])
    assert inspect(out / 'report.json')['status'] == 'select_group'


def test_deduplicate_molecules_use_best_conformer_and_do_not_pad(tmp_path):
    values = [('same', -100.0 - i) for i in range(20)] + [(f'm{i}', -float(i)) for i in range(9)]
    source, out = panel(tmp_path, {'A': values, 'small': [('x', -999.0)]})
    add_rankings(source, out)
    ranks = inspect(out / 'report.json', scheme='E094', receptor='R1')
    assert len(ranks['blocks']) == 1
    assert ranks['blocks'][0]['top_n_mean'] == pytest.approx((-119 - 36) / 10)
    assert ranks['blocks'][0]['scored_units'] == 10
    assert ranks['excluded'] == [{'status': 'insufficient_candidates', 'blocks': 1}]
    top = inspect(out / 'report.json', operation='top', scheme='E094', receptor='R1', block_id='A')
    assert top['candidates'][0]['cid'] == 'A-19'
    assert top['candidates'][1]['molecule_id'] == 'm8'
    add_rankings(source, out, ranking_unit='conformer')
    ranks = inspect(out / 'report.json', scheme='E094', receptor='R1')
    assert ranks['blocks'][0]['top_n_mean'] == -114.5


def test_evicted_molecule_can_reenter_with_a_better_conformer(tmp_path):
    values = [(f'm{i}', -float(i)) for i in range(30)] + [('m0', -999)]
    source, out = panel(tmp_path, {'A': values})
    add_rankings(source, out)
    top = inspect(out / 'report.json', operation='top', scheme='E094', receptor='R1', block_id='A')
    assert top['candidates'][0]['molecule_id'] == 'm0'
    assert top['candidates'][0]['score'] == -999


@pytest.mark.parametrize('partial', [True, False])
def test_incomplete_scores_never_compete_with_complete_blocks(tmp_path, partial):
    values = [(f'm{i}', -float(i)) for i in range(12)]
    if not partial: values.append(('missing', None))
    source, out = panel(tmp_path, {'A': values}, status='partial' if partial else 'complete')
    add_rankings(source, out)
    ranks = inspect(out / 'report.json', scheme='E094', receptor='R1')
    assert not ranks['blocks']
    assert ranks['excluded'] == [{'status': 'incomplete_panel', 'blocks': 1}]


def test_custom_n_ties_and_artifact_integrity(tmp_path):
    values = [(f'm{i}', -1.0) for i in range(15)]
    source, out = panel(tmp_path, {'B': values, 'A': values})
    add_rankings(source, out, top_n=12)
    ranks = inspect(out / 'report.json', scheme='E094', receptor='R1')
    assert [r['block_id'] for r in ranks['blocks']] == ['A', 'B']
    assert ranks['top_n'] == 12
    assert inspect(out / 'report.json', scheme='E094', receptor='R1', top_n=5)['top_n'] == 5
    with pytest.raises(ValueError, match='Queue analysis'):
        inspect(out / 'report.json', scheme='E094', receptor='R1', top_n=7)
    with (out / 'block_ranking.sqlite').open('ab') as stream: stream.write(b'changed')
    with pytest.raises(ValueError, match='changed'):
        inspect(out / 'report.json')


@pytest.mark.parametrize('value', [0, -1, 101, True, 5.5, '10'])
def test_invalid_n(value):
    with pytest.raises(ValueError): options(value)
