"""E113: bounded contact-feature coverage comparison with native-state controls."""
import argparse
import multiprocessing
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from types import SimpleNamespace
import time
import numpy as np
import audit_native_state_poses as native
base = native.base
LEVELS = (0.5, 0.7, 0.9)
_WORK = None


def match_contact_features(state, features, transform):
    from aidd_agent.interaction_fast import match_batch
    indices = sorted({a['feature_index'] for a in state['anchors']})
    if not indices:
        return False, [], 0.0, []
    data = state['aligned']
    query = tuple(data[k][indices] for k in ('feature_points', 'feature_types',
                                            'feature_directions', 'feature_direction_kinds'))
    average, assignments, scores = match_batch((*query, np.ones(len(indices))),
        features.feature_points[None], features.feature_types[None],
        features.feature_directions[None], features.feature_kinds[None],
        native.rigid(transform)[None])
    lookup = {idx: j for j, idx in enumerate(indices)}
    hits = [a for a in state['anchors'] if assignments[0, lookup[a['feature_index']]] >= 0
            and scores[0, lookup[a['feature_index']]] >= base.POLICY['anchor_threshold']]
    residues = sorted({a['target_residue'] for a in hits})
    ile = any(a['target_residue'] == 61 and a['protein_atom'] == 'CD1'
              and a['feature_class'] == 'hydrophobic' for a in hits)
    detail = [dict(reference_feature=i, candidate_feature=int(assignments[0,j]),
                   score=float(scores[0,j])) for j,i in enumerate(indices)]
    return bool(ile and len(residues) >= 3), residues, float(average[0]), detail


def controlled_candidate(state, key):
    data = state['aligned']
    matrix = np.eye(4)
    matrix[:3,:3] = [[0,-1,0],[0,0,-1],[1,0,0]]
    matrix[:3,3] = [17,-11,6]
    points = native.move(data['feature_points'],matrix)
    candidate = SimpleNamespace(shape_points=native.move(data['shape_points'],matrix),
        feature_points=points, feature_types=data['feature_types'], conformer_id=key, molecule_id=key)
    features = SimpleNamespace(feature_points=points, feature_types=data['feature_types'],
        feature_directions=data['feature_directions'] @ matrix[:3,:3].T,
        feature_kinds=data['feature_direction_kinds'])
    return candidate, features


def search(candidate, features, states, deadline):
    from aidd_agent.gaussian_batch import prepare_seeds
    from aidd_agent.gaussian_overlay import score_overlay
    from aidd_agent.full_library_screen import POSE_PARAMETERS
    result = dict(complete=False, tested_seeds=0, contact=False, original_pass=False,
                  geometry_without_shape_color=False, best={},
                  **{'replacement_'+str(t):False for t in LEVELS})
    supported = {k:s for k,s in states.items() if s['record']['supported_for_requested_contact_rule']}
    result['supported_targets'] = len(supported)
    if not supported:
        result['reason'] = 'no_supported_target_state'
        return result
    for key,state in supported.items():
        if time.monotonic() > deadline: return result
        seeds,_ = prepare_seeds(candidate,state['aligned'],**POSE_PARAMETERS)
        for seed in seeds:
            if time.monotonic() > deadline: return result
            result['tested_seeds'] += 1
            matrix = native.rigid(seed.transform_matrix)
            ok,residues,coverage,assignments = match_contact_features(state,features,matrix)
            if not ok: continue
            result['contact'] = True
            data = state['aligned']
            score = score_overlay(data['shape_points'],candidate.shape_points,
                data['feature_points'],data['feature_types'],candidate.feature_points,
                candidate.feature_types,matrix.ravel(),cutoff=None)
            shape,color = score['shape']['tanimoto'],score['color']['tanimoto']
            metrics = base.geometry_metrics(native.move(candidate.shape_points,matrix),
                                           state['tree'],state['contact'],state['regions'])
            flags = base.gates(ok,shape,color,metrics)
            geometry = flags['steric'] and flags['envelope'] and flags['regions']
            result['original_pass'] |= base.accepted(flags)
            result['geometry_without_shape_color'] |= geometry
            for level in LEVELS:
                result['replacement_'+str(level)] |= geometry and flags['shape'] and coverage >= level
            pose = dict(template=key, shape=shape, color=color, coverage=coverage,
                        residues=residues, feature_assignments=assignments, **metrics,
                        original_pass=base.accepted(flags), geometry_pass=geometry,
                        transform=matrix.ravel().tolist())
            for objective in ('shape','color','coverage','coverage_with_geometry'):
                if objective == 'coverage_with_geometry' and not geometry: continue
                score_key = 'coverage' if objective == 'coverage_with_geometry' else objective
                old = result['best'].get(objective)
                if old is None or pose[score_key] > old[score_key]:
                    result['best'][objective] = pose
    result['complete'] = True
    # Export only the few selected poses, avoiding coordinate work inside seed loops.
    for pose in result['best'].values():
        pose['aligned_heavy_atom_coordinates'] = native.move(candidate.shape_points,
            np.asarray(pose['transform']).reshape(4,4)).tolist()
    return result


def initialize(batch, states, coarse, deadline):
    global _WORK
    native.fresh_init(batch,states,coarse,deadline)
    _WORK = native._FRESH


def evaluate(task):
    from aidd_agent.guided_filters import adaptive_stage
    reader,feature_reader,states,coarse,deadline = _WORK
    metadata = dict(kind=task['kind'],id=task['id'],mid=task['mid'])
    if time.monotonic() > deadline:
        return dict(metadata,complete=False,reason='deadline_before_task')
    if task['kind'] != 'library':
        key = task['id']
        candidate,features = controlled_candidate(states[key],key)
        targets = ({key:states[key]} if task['kind']=='self_control' else
                   {k:s for k,s in states.items() if k!=key})
    else:
        candidate = reader.get(int(task['id']))
        features = feature_reader.get(int(task['id']))
        base.require(candidate.conformer_id==features.conformer_id==task['cid'] and
                     candidate.molecule_id==features.molecule_id==task['mid'],'Candidate identity mismatch')
        if adaptive_stage(candidate,coarse)<3:
            return dict(metadata,complete=True,coarse_pass=False,tested_seeds=0,contact=False,
                        original_pass=False,geometry_without_shape_color=False,best={},
                        **{'replacement_'+str(t):False for t in LEVELS})
        targets = states
    return dict(metadata,coarse_pass=True,**search(candidate,features,targets,deadline))


def summarize(rows):
    completed = [r for r in rows if r['complete']]
    return dict(selected=len(rows),complete=len(completed),incomplete=len(rows)-len(completed),
        counts={k:sum(bool(r.get(k,False)) for r in completed) for k in
                ['contact','original_pass','geometry_without_shape_color']+
                ['replacement_'+str(t) for t in LEVELS]},
        original_pass_fraction=(sum(r['original_pass'] for r in completed)/len(rows)
                                if rows and len(completed)==len(rows) else None))


def pilot(args, selected, states, query, batch, out, controls, started):
    base.require(0<args.sample<=len(selected),'Sample exceeds available E111 sample')
    base.require(1<=args.workers<=32,'Invalid worker count')
    selected = selected[:args.sample]
    base.write_csv(out/'sample.csv',selected)
    tasks = [dict(kind=kind,id=key,mid=key) for kind in ('self_control','leave_self_out') for key in states]
    tasks += [dict(kind='library',id=r['gid'],mid=r['mid'],cid=r['cid']) for r in selected]
    protocol = base.read(out/'protocol.json')
    protocol.update(mode='E113_contact_feature_diagnostic',
                    diagnostic_script_sha256=base.sha(__file__),
                    coverage_levels=list(LEVELS),control_tasks=22,
                    coverage_definition='Mean directional one-to-one score over unique adopted native contact feature indices',
                    no_exposure_or_buried_polar_model=True)
    base.save(out/'protocol.json',protocol)
    setup = time.monotonic()-started
    tick = time.monotonic()
    with ProcessPoolExecutor(max_workers=args.workers,mp_context=multiprocessing.get_context('spawn'),
        initializer=initialize,initargs=(str(batch),states,query['design']['adaptive_coarse'],tick+args.seconds)) as pool:
        outcomes = []
        for row in pool.map(evaluate,tasks,chunksize=1):
            outcomes.append(row)
            if len(outcomes)%20==0: print('Returned '+str(len(outcomes))+'/'+str(len(tasks)),flush=True)
    base.save(out/'pose_diagnostics.json',outcomes)
    groups = {kind:summarize([r for r in outcomes if r['kind']==kind])
              for kind in ('library','self_control','leave_self_out')}
    # Unsupported self states are explicit; they are not timeouts or binding failures.
    timed_out = any(not r['complete'] and r.get('reason')!='no_supported_target_state' for r in outcomes)
    report = dict(status='time_budget_exhausted' if timed_out else 'complete',
        kind='E113_contact_feature_diagnostic',groups=groups,
        unsupported_controls=[r['id'] for r in outcomes if r.get('reason')=='no_supported_target_state'],
        tested_seeds=sum(r.get('tested_seeds',0) for r in outcomes),workers=args.workers,
        setup_seconds=setup,search_seconds=time.monotonic()-tick,total_seconds=time.monotonic()-started,
        limitations=['Exploratory sensitivity levels, not production thresholds.',
                     'Adopted contact subset only; geometric matching is not binding-energy evidence.',
                     'Additional features excluded from coverage denominator; no solvent exposure is inferred.',
                     'No independent buried-polar or direct protein interaction validation.',
                     'Controls bypass library coarse eligibility; leave-self-out references can be related.',
                     'Shape remains 0.50; changing color alone may not recover candidates.',
                     'Partial task maxima are diagnostic only and excluded from complete counts.'])
    report['output_hashes']={p.name:base.sha(p) for p in out.iterdir() if p.is_file() and p.name!='report.json'}
    base.save(out/'report.json',report)
    print(base.json.dumps(report,indent=2))
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('previous','prepared','reference','output'):
        parser.add_argument('--'+name,type=Path,required=True)
    parser.add_argument('--sample',type=int,default=200)
    parser.add_argument('--workers',type=int,default=20)
    parser.add_argument('--seconds',type=float,default=120)
    args=parser.parse_args()
    args.fresh=True
    native.fresh_pilot=pilot
    native.run(args)


if __name__=='__main__': main()
