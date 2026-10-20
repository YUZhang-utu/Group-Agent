"""E114: recheck 65 saved rigid proposals against original workstation receptor data."""
import argparse
from pathlib import Path
import time
import numpy as np
import audit_native_state_poses as native
from diagnose_contact_features import controlled_candidate, match_contact_features
base=native.base


def assess(candidate, features, state, transform):
    from aidd_agent.gaussian_overlay import score_overlay
    matrix=native.rigid(transform);data=state['aligned']
    ok,residues,coverage,assignments=match_contact_features(state,features,matrix)
    score=score_overlay(data['shape_points'],candidate.shape_points,data['feature_points'],data['feature_types'],
                        candidate.feature_points,candidate.feature_types,matrix.ravel(),cutoff=None)
    shape,color=score['shape']['tanimoto'],score['color']['tanimoto']
    metrics=base.geometry_metrics(native.move(candidate.shape_points,matrix),state['tree'],state['contact'],state['regions'])
    flags=base.gates(ok,shape,color,metrics)
    return dict(shape=shape,color=color,coverage=coverage,residues=residues,**metrics,flags=flags,
                original_policy_pass=base.accepted(flags),
                contacts_shape_geometry=bool(ok and flags['shape'] and flags['steric'] and flags['envelope'] and flags['regions']),
                feature_assignments=assignments)


def pilot(args, selected, states, query, batch, out, controls, started):
    from aidd_agent.gaussian_batch import ArtifactCatalogReader
    from aidd_agent.interaction_fast import InteractionFeatureReader
    receipt=args.receipt.resolve();bundle=base.read(args.proposals)
    base.verify(receipt/'report.json',bundle['source_report_sha256'])
    base.verify(receipt/'protocol.json',bundle['source_protocol_sha256'])
    report=base.read(receipt/'report.json');protocol=base.read(receipt/'protocol.json')
    base.require(report['status']=='complete','Incomplete E113 source')
    for name,digest in report['output_hashes'].items():base.verify(receipt/name,digest)
    base.require(protocol['query_sha256']==base.read(out/'protocol.json')['query_sha256'],'Query lineage mismatch')
    base.require(protocol['previous_report_sha256']==base.sha(args.previous/'report.json'),'Sample lineage mismatch')
    evidence={(r['kind'],r['id']):r for r in base.read(receipt/'pose_diagnostics.json')}
    sample={r['gid']:r for r in selected}
    reader=ArtifactCatalogReader(batch/'artifacts/catalog.json')
    features=InteractionFeatureReader(batch/'artifacts/catalog.json',batch/'chemical/catalog.json')
    output=[];tick=time.monotonic()
    try:
        for p in bundle['proposals']:
            if time.monotonic()-tick>args.seconds:break
            source=evidence[(p['kind'],p['id'])]['best'][p['objective']]
            base.require(source['template']==p['template'] and
                         np.allclose(source['transform'],p['original_transform'],atol=1e-12),'Proposal source mismatch')
            if p['kind']=='library':
                item=sample[p['id']];candidate=reader.get(int(p['id']));cf=features.get(int(p['id']))
                base.require(candidate.conformer_id==cf.conformer_id==item['cid'] and
                             candidate.molecule_id==cf.molecule_id==p['mid']==item['mid'],'Library identity mismatch')
            else:
                base.require(p['kind']=='leave_self_out','Unsupported proposal kind')
                candidate,cf=controlled_candidate(states[p['id']],p['id'])
            state=states[p['template']]
            before=assess(candidate,cf,state,p['original_transform'])
            after=assess(candidate,cf,state,p['transform'])
            base.require(abs(before['shape']-p['original_shape'])<1e-6 and
                         abs(after['shape']-p['optimized_shape'])<1e-6,'Shape replay mismatch')
            output.append(dict(kind=p['kind'],id=p['id'],mid=p['mid'],template=p['template'],category=p['category'],
                               before=before,after=after,transform=p['transform']))
    finally:
        features.close()
        if hasattr(reader,'close'):reader.close()
    base.save(out/'checked_poses.json',output)
    settings=base.read(out/'protocol.json')
    settings.update(mode='E114_rigid_proposal_recheck_no_search',checker_sha256=base.sha(__file__),
                    proposals_sha256=base.sha(args.proposals),e113_report_sha256=bundle['source_report_sha256'])
    base.save(out/'protocol.json',settings)
    result=dict(status='complete' if len(output)==len(bundle['proposals']) else 'time_budget_exhausted',
                selected=len(bundle['proposals']),processed=len(output),total_seconds=time.monotonic()-started,
                evaluation_seconds=time.monotonic()-tick,
                groups={kind:dict(poses=sum(r['kind']==kind for r in output),
                    shape_pass_after=sum(r['kind']==kind and r['after']['flags']['shape'] for r in output),
                    original_policy_pass_after=sum(r['kind']==kind and r['after']['original_policy_pass'] for r in output),
                    contacts_shape_geometry_after=sum(r['kind']==kind and r['after']['contacts_shape_geometry'] for r in output))
                    for kind in ['library','leave_self_out']},
                limitations=['Outcome-selected proposals; counts are poses, not library-wide pass rates.',
                             'Rigid shape-only optimization may lose contacts or introduce clashes.',
                             'No threshold change or new library seed search.'])
    result['output_hashes']={p.name:base.sha(p) for p in out.iterdir() if p.name!='report.json'}
    base.save(out/'report.json',result);print(base.json.dumps(result,indent=2));return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('previous','prepared','reference','output','receipt'):
        parser.add_argument('--'+name,type=Path,required=True)
    parser.add_argument('--proposals',type=Path,default=Path(__file__).resolve().parents[1]/'to_human/E114_SHAPE_PROPOSALS.json')
    parser.add_argument('--seconds',type=float,default=120)
    args=parser.parse_args();args.fresh=True;args.sample=10000;args.workers=1
    native.fresh_pilot=pilot;native.run(args)


if __name__=='__main__':main()
