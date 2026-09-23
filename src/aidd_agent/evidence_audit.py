"""Three-layer crystal matching and source-evidence audit without adoption."""
import argparse
import copy
import json
from pathlib import Path

import numpy as np

from . import contact_diagnostics as diag
from .consensus_design import validate
from .contact_groups import selected_anchors
from .expanded_wee1 import fingerprint
from .gaussian_batch import _load_query, _atomic_json
from .screening_selection import check_hashes


def run(recommendation, output, duplicate_distance=.5):
    source=Path(recommendation)
    proposal=json.loads(source.read_text(encoding='utf-8'));check_hashes(proposal['sources'])
    if proposal.get('kind')!='consensus_recommendation':raise ValueError('Expected a saved consensus recommendation')
    survey_path=Path(proposal['survey'])
    if str(survey_path.resolve()) not in proposal['sources']:raise ValueError('Unsealed consensus report')
    survey=json.loads(survey_path.read_text(encoding='utf-8'));check_hashes(survey['sources'])
    design=proposal['recommendation'];validate(design,survey)
    if not np.isfinite(duplicate_distance) or duplicate_distance<=0:raise ValueError('Invalid duplicate distance')
    out=Path(output);out.mkdir(parents=True,exist_ok=False)
    anchors={a['anchor_id']:a for a in survey['anchors']};order=selected_anchors(design)
    _,pocket=_load_query(Path(survey['consensus_npz']))
    indices=np.array([anchors[a]['feature_index'] for a in order])
    observations={}
    for a in survey['anchors']:
        for o in a['observations']:
            observations.setdefault(o['query_id'],[]).append(dict(copy.deepcopy(o),anchor_id=a['anchor_id']))
    rows=[];records=[];unknown=[];timing=[];comparisons=[];weighted=[];extra_sources={}
    for meta in survey['prepared_complexes']:
        qid=meta['query_id'];path=Path(meta['query_npz']);_,q=_load_query(path)
        own=observations.get(qid,[])
        native=path.parent/'native.npz';manifest_path=native.with_suffix('.manifest.json')
        manifest=None
        if str(manifest_path.resolve()) in survey['sources'] and native.exists():
            manifest,native_query=_load_query(native)
            extra_sources.update(fingerprint([native,manifest_path]))
            # Alignment preserves feature ordering and chemistry.
            if not np.array_equal(native_query['feature_types'],q['feature_types']):
                raise ValueError('Native/aligned feature identity mismatch')
        identity,result,perf=diag.match_details(diag.query_tuple(q),q,
            [dict(query_id=qid,layer='feature_identity',feature_index=i) for i in range(len(q['feature_points']))])
        rows.extend(identity);timing.append(dict(query_id=qid,layer='feature_identity',**perf))
        for o in own:
            i=diag.source_index(o,manifest) if manifest else None
            members=manifest['source'].get('feature_atom_indices',[]) if manifest else []
            valid=i is not None and i<len(q['feature_points']) and i<len(members) and bool(members[i])
            if not valid:
                o['ligand_atom_indices']=[]
                unknown.append(dict(query_id=qid,source_anchor=o['source_anchor'],reason='missing_verified_source_mapping'))
            else:
                o.update(source_feature_index=i,ligand_atom_indices=members[i],
                    extraction=diag.extraction_metadata(o['key'][2]),mapping_status='verified_native_manifest')
                if 'source_feature_index' in next((v for v in anchors[o['anchor_id']]['observations'] if v['source_anchor']==o['source_anchor']),{}):
                    original=next(v for v in anchors[o['anchor_id']]['observations'] if v['source_anchor']==o['source_anchor'])
                    if original['source_feature_index']!=i or original.get('ligand_atom_indices')!=members[i]:
                        raise ValueError('Recorded atom mapping disagrees with sealed manifest')
                source_query=[np.array([o['point']]),np.array([o['key'][3]]),np.array([o['direction']]),np.array([o['kind']])]
                # Isolate the original feature: no competition with other hypotheses.
                candidate={k:q[k][[i]] for k in ('feature_points','feature_types','feature_directions','feature_direction_kinds')}
                detail,_,_=diag.match_details(source_query,candidate,[dict(query_id=qid,layer='source_observation',anchor_id=o['anchor_id'],source_anchor=o['source_anchor'],source_feature_index=i)])
                for item in detail:item['interaction_type']=o['key'][2]
                rows.extend(detail)
            records.append(o)
        consensus,result,perf=diag.match_details(diag.query_tuple(pocket,indices),q,
            [dict(query_id=qid,layer='selected_consensus',anchor_id=a,interaction_type=anchors[a]['feature_class']) for a in order])
        rows.extend(consensus);timing.append(dict(query_id=qid,layer='selected_consensus',**perf))
        terms=diag.contributions(result['anchor_scores'],result['assignments'],order,design,anchors)
        weighted.extend(dict(t,query_id=qid) for t in terms)
        alternative=diag.optional_group_assignment(result['pair_scores'],order,design)
        current=sum(t['numerator'] for t in terms)
        comparisons.append(dict(query_id=qid,production_optional_numerator=current,alternative=alternative,
            diagnostic_gain=alternative['optional_numerator']-current if 'optional_numerator' in alternative else None))
    relations=diag.duplicate_relations(records,duplicate_distance)
    identity=[r for r in rows if r['layer']=='feature_identity']
    recovery=[r for r in rows if r['layer']=='source_observation']
    passed=lambda values:bool(values) and all(abs(r['score']-1)<1e-6 for r in values)
    by_type={}
    for row in weighted:
        item=by_type.setdefault(row['interaction_type'],dict(contact_contribution_sum=0.,composite_contribution_sum=0.))
        item['contact_contribution_sum']+=row['contact_contribution']
        item['composite_contribution_sum']+=row['composite_contribution']
    result=dict(kind='contact_evidence_audit',status='complete',design=design,
        scope='Identity and geometric diagnostics only; no biological discrimination claim, library scan or production assignment change',
        checks=dict(feature_identity_passed=passed(identity),source_recovery_passed=passed(recovery) if not unknown else None,
            source_mapping_unresolved=len(unknown)),unresolved=unknown,
        distributions=diag.distributions(rows),match_rows=rows,source_observations=records,
        suspected_shared_evidence=relations,weighted_contributions=weighted,
        contributions_by_interaction_type=by_type,contribution_summary_scope='Sum over crystal poses; group ties split credit; unmatched contributions included as zero',
        assignment_diagnostics=comparisons,timings=timing,
        parameters=dict(sigma_angstrom=1.,angular_power=2.,cutoff_angstrom=4.5,duplicate_distance_angstrom=duplicate_distance),
        sources={**proposal['sources'],**survey['sources'],**extra_sources,**fingerprint([source,survey_path])},
        implementation_sources=fingerprint(sorted(Path(__file__).parent.glob('*.py'))))
    _atomic_json(out/'report.json',result)
    lines=['# Contact evidence diagnostics','',
        'Crystal geometry diagnostics only; no biological recall or affinity claim.','',
        '## Checks','',json.dumps(result['checks']), '',
        '## Matching distributions','',
        '| Layer | Feature | Direction | Count | Unmatched | Min | Median | Max |',
        '| --- | --- | --- | ---: | ---: | ---: | ---: | ---: |']
    for d in result['distributions']:
        quantiles=d['score']['quantiles']
        lines.append(f"| {d['layer']} | {d['feature_name']} | {d['direction_name']} | {d['count']} | {d['unmatched']} | {quantiles[0]:.6f} | {quantiles[2]:.6f} | {quantiles[4]:.6f} |")
    gains=[r['diagnostic_gain'] for r in comparisons if r['diagnostic_gain'] is not None]
    lines+=['','## Review boundaries','',
        f"Suspected shared-evidence pairs: {len(relations['pairs'])}. No pairs are automatically merged.",
        f"Optional-only assignment improves {sum(g>1e-9 for g in gains)} of {len(gains)} tested poses. Production matching is unchanged.",
        'See report.json for each source mapping, geometry factor, weighted contribution and alternative assignment.',
        'Unmatched scores are zero; missing geometry factors are null, not successful geometry observations.',
        'Passing self-tests is an implementation check, not proof that candidate angular penalties discriminate activity.']
    (out/'review.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--recommendation',required=True)
    parser.add_argument('--output',required=True)
    parser.add_argument('--duplicate-distance',type=float,default=.5)
    args=parser.parse_args();r=run(args.recommendation,args.output,args.duplicate_distance)
    print(json.dumps(dict(report=str((Path(args.output)/'report.json').resolve()),checks=r['checks'],
        suspected_pairs=len(r['suspected_shared_evidence']['pairs'])),indent=2))
    if not r['checks']['feature_identity_passed'] or r['checks']['source_recovery_passed'] is False:raise SystemExit(1)


if __name__=='__main__':main()
