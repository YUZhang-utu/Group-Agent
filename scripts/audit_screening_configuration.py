"""Export sealed scientific query settings without executing search or modifying inputs."""
import argparse
import hashlib
import json
from pathlib import Path


RUNS = ('5f7c28e4bd704616 3a90516120f74140 4f0fad82c9714b1e '
        '34bf4cc4a192449e 7a7fddb92b744a38 9bd13fe58ad14df0 '
        'dcc64417cf534000 317795cf58ad4f04 803b293dd57c4aec '
        '82f400cb6c1c46df c71019343cb04b62 bfa50ed8eb9a489e').split()
RUNS_ROOT = '/mnt/local/hand/yuzhang/aidd/e097-chat-workspace/users/workstation/projects/prj-adf8a1b9f4f8-prompt-aidd/runs'
OLD_SEARCH = '/mnt/local/hand/yuzhang/aidd/e054-chat-20260923-104806/users/workstation/projects/prj-42c81b1322ab-prompt-aidd/runs/PROMPT-59addb9ffd6444c3/execution/guided/guided/search'
DESIGN_FIELDS = ('mandatory_anchors', 'alternative_groups', 'optional_weights', 'optional_groups',
                 'optional_normalization', 'optional_budget', 'spatial_groups', 'resolved_spatial_groups',
                 'occupancy_weight', 'gaussian_weight', 'optional_weight', 'minimum_score',
                 'minimum_pose_score', 'pose_threshold_derivation', 'adaptive_coarse', 'coarse_constraints',
                 'pocket', 'exclusions', 'permissiveness', 'template_ids')


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1048576), b''):
            digest.update(chunk)
    return digest.hexdigest()


def save(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False), encoding='utf-8')


def require(condition, message):
    if not condition:
        raise ValueError(message)


def resource(path_value, sources, base):
    if not path_value:
        return dict(status='not_declared')
    path = Path(path_value)
    if not path.is_absolute():
        path = base / path
    if not path.is_file():
        return dict(path=str(path), status='missing')
    digest = sha(path)
    expected = sources.get(path_value) or sources.get(str(path.resolve()))
    require(not expected or expected == digest, f'Sealed query resource changed: {path}')
    return dict(path=str(path), sha256=digest, status='verified' if expected else 'present_unsealed_in_query_sources')


def summarize_query(query, base):
    require(query.get('kind') == 'consensus_design', 'Expected an adopted consensus_design, not a proposal')
    design = query['design']
    sources = query.get('sources', {})
    templates = []
    for template in query['templates']:
        templates.append(dict(query_id=template['query_id'],
                              resource=resource(template.get('query_npz'), sources, base),
                              metadata=template))
    require(len({row['query_id'] for row in templates}) == len(templates), 'Duplicate query template')
    anchors = query.get('anchors', [])
    receptor = resource(design.get('receptor_npz'), sources, base)
    return dict(target=query.get('target'), reference=query.get('reference'),
                readiness=query.get('readiness'), template_count=len(templates), templates=templates,
                receptor=receptor, consensus=resource(query.get('consensus_npz'), sources, base),
                design={field: design.get(field) for field in DESIGN_FIELDS}, anchors=anchors,
                scoring_contract=query.get('scoring_contract'),
                crystal_self_controls=query.get('crystal_self_controls'),
                failed_template_self_controls=query.get('failed_template_self_controls'),
                limitations=query.get('limitations'),
                enabled_declarations=dict(
                    mandatory_anchor_count=len(design.get('mandatory_anchors', [])),
                    alternative_group_count=len(design.get('alternative_groups', [])),
                    optional_anchor_count=len(design.get('optional_weights', {})),
                    optional_family_count=len(design.get('optional_groups', [])),
                    spatial_group_count=len(design.get('spatial_groups', [])),
                    spatial_reward_weight=design.get('occupancy_weight', 0),
                    gaussian_composite_weight=design.get('gaussian_weight'),
                    direct_receptor_resource_count=int(bool(design.get('receptor_npz'))),
                    water_metal_anchor_ids=[row['anchor_id'] for row in anchors
                                           if row.get('feature_class') in ('water_bridge', 'metal_coordination')]))


def compare_queries(old, current):
    rows = []
    left = {**old['design'], 'reference': old['reference'],
            'template_ids_ordered': [row['query_id'] for row in old['templates']],
            'receptor_sha256': old['receptor'].get('sha256'),
            'consensus_sha256': old['consensus'].get('sha256'),
            'template_sha256': [row['resource'].get('sha256') for row in old['templates']]}
    right = {**current['design'], 'reference': current['reference'],
             'template_ids_ordered': [row['query_id'] for row in current['templates']],
             'receptor_sha256': current['receptor'].get('sha256'),
             'consensus_sha256': current['consensus'].get('sha256'),
             'template_sha256': [row['resource'].get('sha256') for row in current['templates']]}
    for field in left:
        missing = left[field] is None or right[field] is None
        if field == 'template_sha256':
            missing = any(value is None for value in left[field] + right[field])
        rows.append(dict(field=field, comparison='unknown_or_absent' if missing else
                         'equal' if left[field] == right[field] else 'different',
                         historical=left[field], current=right[field]))
    return rows


def audit(runs_root, historical_search, output):
    output.mkdir(parents=True, exist_ok=False)
    save(output / 'audit.json', dict(status='incomplete'))
    arms, query_hashes, query_summary = [], set(), None
    for run in RUNS:
        folder = runs_root / ('PROMPT-' + run) / 'execution/campaign/blocks'
        report = read(folder / 'report.json')
        require(report.get('status') == 'complete', f'Incomplete arm: {folder}')
        for name in ('query/report.json', 'protocol.json'):
            require(sha(folder / name) == report.get('output_hashes', {}).get(name), f'Changed/unsealed {folder / name}')
        query = read(folder / 'query/report.json')
        protocol = read(folder / 'protocol.json')
        digest = sha(folder / 'query/report.json')
        if digest not in query_hashes:
            summary = summarize_query(query, folder / 'query')
            save(output / f'current-query-{digest[:16]}.json', summary)
            if query_summary is None:
                query_summary = summary
            query_hashes.add(digest)
        selection = report['selection']
        arms.append(dict(run=run, selection=selection, query_sha256=digest,
                         report_sha256=sha(folder / 'report.json'), protocol_sha256=sha(folder / 'protocol.json'),
                         engine=report.get('search_engine'), filter_policy=report.get('filter_policy'),
                         ranking=protocol.get('ranking'), library_registry=protocol.get('library_registry'),
                         code=protocol.get('code'), input_hashes=protocol.get('inputs'),
                         requested_conformers=report.get('requested_conformers'),
                         searched_blocks=report.get('searched_blocks'),
                         search_policy=report.get('search_policy'), stop_reason=report.get('stop_reason')))
    historical = None
    missing = []
    paths = [historical_search / name for name in ('report.json', 'protocol.json', 'adopted-design/report.json')]
    if all(path.is_file() for path in paths):
        report, protocol, query = [read(path) for path in paths]
        require(report.get('status') == 'complete', 'Historical search incomplete')
        historical = dict(report=report, protocol=protocol, query=summarize_query(query, paths[2].parent),
                          files={str(path): sha(path) for path in paths},
                          provenance='Recorded current hashes; old report/protocol lack an external parent seal here')
        save(output / 'historical-query-and-protocol.json', historical)
    else:
        missing.extend(str(path) for path in paths if not path.is_file())
    require(len(query_hashes) == 1, 'Current arms have different adopted queries; inspect exported files')
    differences = compare_queries(historical['query'], query_summary) if historical else []
    save(output / 'query-comparison.json', differences)
    save(output / 'current-arms.json', arms)
    findings = [
        'This is a configuration audit, not an independent validation of molecular physics or candidate quality.',
        'Equal query settings do not imply equal execution: historical budget and current threshold branches differ.',
        'Absent hard-anchor or spatial-group declarations are not proof that these are scientifically required.',
        'A single receptor resource is not an ensemble, even when multiple ligand templates are present.',
        'Gaussian composite weight describes the pose threshold; contact-first final RRF has separate semantics.',
        'Pocket-state recommendation, protonation, ligand chemistry, and receptor compatibility need their own evidence.']
    resource_gaps = []
    for label, summary in [('current', query_summary)] + ([('historical', historical['query'])] if historical else []):
        resources = [('receptor', summary['receptor']), ('consensus', summary['consensus'])]
        resources.extend((row['query_id'], row['resource']) for row in summary['templates'])
        resource_gaps.extend(dict(scope=label, resource=name, **record) for name, record in resources
                             if record['status'] != 'verified')
    result = dict(status='complete_with_gaps' if missing or resource_gaps else 'complete', current_arms=12,
                  current_query_variants=len(query_hashes), historical_files_missing=missing,
                  resource_verification_gaps=resource_gaps,
                  observed_current_receptors=sorted({str(row['selection'].get('receptor')) for row in arms}),
                  limitations=findings, source_sha256=sha(Path(__file__)), scientific_execution=False)
    save(output / 'audit.json', result)
    (output / 'report.md').write_text('# Screening configuration audit\n\n'
        'Inspect current-query-*.json, current-arms.json, historical-query-and-protocol.json and query-comparison.json.\n\n'
        + '\n'.join('- ' + finding for finding in findings)
        + '\n\nMissing historical files: ' + json.dumps(missing) + '\n', encoding='utf-8')
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--runs-root', type=Path, default=Path(RUNS_ROOT))
    parser.add_argument('--historical-search', type=Path, default=Path(OLD_SEARCH))
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = audit(args.runs_root, args.historical_search, args.output)
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
