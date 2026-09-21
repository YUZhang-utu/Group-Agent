"""Diagnostics for explicit rules; synthetic controls are not library positives."""
from types import SimpleNamespace

import numpy as np

from . import full_library_screen as full
from .interaction_fast import match_batch


def anchor_diagnostics(scores, assignments, anchors, threshold):
    hits = (scores >= threshold) & (assignments >= 0)
    cumulative = np.ones(len(scores), dtype=bool)
    rows = []
    for i, anchor in enumerate(anchors):
        cumulative &= hits[:, i]
        rows.append(dict(anchor_id=anchor['anchor_id'], individual_passes=int(hits[:, i].sum()),
                         cumulative_all_passes=int(cumulative.sum())))
    return rows


def self_control(original, expanded, definition):
    from .gaussian_batch import _score_ids
    from .joint_coarse import JointCoarse
    from .necessary_conditions import NecessaryConditions
    from .pose_feasibility import possible_seed_mask
    from .gaussian_batch import prepare_seeds
    policy = definition['condition_policy']
    anchors = [a for a in definition['anchors'] if a['anchor_id'] in policy['required_anchors']]
    indices = expanded['anchor_feature_indices']
    columns = [a['score_column'] for a in anchors]
    candidate = SimpleNamespace(**{k:original[k] for k in ('shape_points','feature_points','feature_types','feature_directions')},
        feature_kinds=original['feature_direction_kinds'], molecule_id='synthetic-query', conformer_id='synthetic-query')
    gate = NecessaryConditions(expanded,[a['feature_index'] for a in anchors],policy['match_mode'],policy['minimum_score'])
    joint = JointCoarse(original,policy['coarse_constraints']).check(candidate)[0]
    reader = SimpleNamespace(get=lambda gid:candidate)
    rigid = _score_ids(reader,original,np.array([0]),backend='reference',**full.POSE_PARAMETERS)
    seeds, _ = prepare_seeds(candidate,original,**full.POSE_PARAMETERS)
    transforms = np.concatenate([np.eye(4)[None],rigid[full.OBJECTIVE+'__transform'].reshape(1,4,4)])
    query = tuple(expanded[k][indices] for k in ('feature_points','feature_types','feature_directions','feature_direction_kinds','anchored_weights'))
    _, assignments, scores = match_batch(query,
        np.repeat(candidate.feature_points[None],2,axis=0), np.repeat(candidate.feature_types[None],2,axis=0),
        np.repeat(candidate.feature_directions[None],2,axis=0),np.repeat(candidate.feature_kinds[None],2,axis=0),transforms)
    hits = (scores[:,columns]>=policy['minimum_score']) & (assignments[:,columns]>=0)
    passed = hits.all(axis=1) if policy['match_mode']=='all' else hits.any(axis=1)
    return dict(scope='Synthetic crystal-query features; not an independently prepared ligand or library positive',
        joint_eligible=bool(joint), anchor_bound_passed=bool(gate.check(candidate)[0]),
        seed_feasibility_passed=bool(possible_seed_mask(gate,candidate,seeds).any()),
        identity_passed=bool(passed[0]), gaussian_selected_pose_passed=bool(passed[1]),
        anchor_ids=[a['anchor_id'] for a in anchors], identity_scores=scores[0,columns].tolist(),
        gaussian_selected_pose_scores=scores[1,columns].tolist())


def validate_survivors(definition, retained, rejected, output):
    from .funnel_benchmark import work, compare_filtered
    _, original, expanded, _ = full._STATE
    control = self_control(original,expanded,definition)
    rejected = np.asarray(rejected,dtype=np.int64)
    rejected = rejected[np.linspace(0,len(rejected)-1,min(64,len(rejected)),dtype=int)] if len(rejected) else rejected
    ids = np.unique(np.concatenate([np.asarray(retained,dtype=np.int64),rejected]))
    anchors = [a for a in definition['anchors'] if a['anchor_id'] in definition['condition_policy']['required_anchors']]
    columns = [a['score_column'] for a in anchors]
    checks, positive_ids, scores, assignments = [], [], [], []
    old_backend, old_flag = definition.get('pose_backend'), definition.get('pose_feasibility')
    try:
        for start in range(0,len(ids),32):
            batch = ids[start:start+32]
            definition.update(pose_backend='reference',pose_feasibility=False)
            reference, _ = work(batch)
            definition.update(pose_backend='numpy',pose_feasibility=True)
            actual, _ = work(batch)
            check = compare_filtered(reference,actual)
            checks.append(check)
            positive_ids.extend(reference['global_ids'][reference['condition_passed']].tolist())
            retained_rows = np.isin(reference['global_ids'],retained)
            scores.append(reference[full.OBJECTIVE+'__anchor_scores'][retained_rows][:,columns])
            assignments.append(reference[full.OBJECTIVE+'__anchor_assignments'][retained_rows][:,columns])
            full._atomic_savez(output/f'reference-{start:08d}.npz',reference)
            full._atomic_savez(output/f'optimized-{start:08d}.npz',actual)
            full.ev.log(f'Rule validation: checked {min(start+32,len(ids))}/{len(ids)}; reference positives {len(positive_ids)}')
    finally:
        definition.update(pose_backend=old_backend,pose_feasibility=old_flag)
    controls_passed = all(control[k] for k in ('joint_eligible','anchor_bound_passed','seed_feasibility_passed','identity_passed','gaussian_selected_pose_passed'))
    return dict(status='passed' if controls_passed and all(c['passed'] for c in checks) else 'failed',
        synthetic_control=control, checked_ids=ids.tolist(), retained_ids=list(map(int,retained)),
        rejected_audit_ids=rejected.tolist(), reference_positive_ids=positive_ids,
        reference_positive_count=len(positive_ids), checks=checks,
        positive_membership_validation='covered_on_panel' if positive_ids else 'not_covered_no_library_positives',
        anchor_diagnostics=anchor_diagnostics(np.concatenate(scores) if scores else np.empty((0,len(columns))),
            np.concatenate(assignments) if assignments else np.empty((0,len(columns))),anchors,definition['condition_policy']['minimum_score']),
        diagnostic_scope='Individual and cumulative ALL counts over coarse survivors; order is reported anchor_ids, thresholds unchanged',
        full_library_acceptance=False)
