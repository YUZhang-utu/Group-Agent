"""Synthetic configuration audit checks; no live screening."""
import importlib.util
from pathlib import Path

import pytest


SPEC = importlib.util.spec_from_file_location('screening_audit', Path(__file__).resolve().parents[1] /
                                            'scripts/audit_screening_configuration.py')
audit = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(audit)


def query_fixture(tmp_path):
    receptor = tmp_path / 'receptor.npz'
    receptor.write_bytes(b'fixture')
    query = dict(kind='consensus_design', target={'accession': 'fixture'}, reference={'query_id': 'fixture'},
                 design=dict(receptor_npz=str(receptor), mandatory_anchors=[], alternative_groups=[],
                             optional_weights={'anchor': 1}, gaussian_weight=0.5),
                 templates=[dict(query_id='fixture', query_npz=str(receptor))],
                 sources={str(receptor): audit.sha(receptor)}, anchors=[], consensus_npz=str(receptor))
    return query


def test_declared_settings_not_assumed_defaults(tmp_path):
    summary = audit.summarize_query(query_fixture(tmp_path), tmp_path)
    assert summary['enabled_declarations']['mandatory_anchor_count'] == 0
    assert summary['enabled_declarations']['spatial_reward_weight'] == 0
    assert summary['design']['minimum_pose_score'] is None
    assert summary['receptor']['status'] == 'verified'


def test_corrupt_resource_rejected(tmp_path):
    query = query_fixture(tmp_path)
    Path(query['design']['receptor_npz']).write_bytes(b'changed')
    with pytest.raises(ValueError, match='resource changed'):
        audit.summarize_query(query, tmp_path)


def test_missing_resources_are_not_equal_evidence(tmp_path):
    query = query_fixture(tmp_path)
    Path(query['design']['receptor_npz']).unlink()
    summary = audit.summarize_query(query, tmp_path)
    comparisons = {row['field']: row for row in audit.compare_queries(summary, summary)}
    assert comparisons['receptor_sha256']['comparison'] == 'unknown_or_absent'
    assert comparisons['template_sha256']['comparison'] == 'unknown_or_absent'
    assert comparisons['gaussian_weight']['comparison'] == 'equal'


def test_complete_audit_and_no_overwrite(tmp_path):
    query = query_fixture(tmp_path)
    runs = tmp_path / 'runs'
    for run in audit.RUNS:
        folder = runs / ('PROMPT-' + run) / 'execution/campaign/blocks'
        (folder / 'query').mkdir(parents=True)
        audit.save(folder / 'query/report.json', query)
        audit.save(folder / 'protocol.json', {})
        audit.save(folder / 'report.json', dict(status='complete', selection={'receptor': 'fixture'},
                   output_hashes={name: audit.sha(folder / name) for name in ('query/report.json', 'protocol.json')}))
    historical = tmp_path / 'historical'
    (historical / 'adopted-design').mkdir(parents=True)
    audit.save(historical / 'adopted-design/report.json', query)
    audit.save(historical / 'report.json', {'status': 'complete'})
    audit.save(historical / 'protocol.json', {})
    output = tmp_path / 'output'
    assert audit.audit(runs, historical, output)['status'] == 'complete'
    with pytest.raises(FileExistsError):
        audit.audit(runs, historical, output)
