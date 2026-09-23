"""Read-only crystal-pose score diagnostics; no library search or design adoption."""
import argparse
import copy
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from .contact_groups import contact_semantics, optional_weight_total, selected_anchors, resolve_regions
from .consensus_design import validate
from .expanded_wee1 import fingerprint
from .gaussian_batch import _atomic_json, _load_query
from .guided_filters import pose_rank, same_pose_gaussian
from .interaction_fast import match_batch
from .screening_selection import check_hashes


def run(recommendation, output):
    source=Path(recommendation)
    proposal=json.loads(source.read_text(encoding='utf-8'))
    if proposal.get('kind')!='consensus_recommendation':raise ValueError('Supply a consensus recommendation report')
    check_hashes(proposal['sources'])
    survey_path=Path(proposal['survey'])
    survey=json.loads(survey_path.read_text(encoding='utf-8'));check_hashes(survey['sources'])
    if str(survey_path.resolve()) not in proposal['sources']:raise ValueError('Recommendation must seal its consensus report')
    design=copy.deepcopy(proposal['recommendation']);validate(design,survey)
    ledger={}
    if survey.get('contact_evidence'):
        ledger_path=Path(survey['contact_evidence']['path'])
        if str(ledger_path.resolve()) not in survey['sources']:raise ValueError('Unsealed contact ledger')
        evidence=json.loads(ledger_path.read_text(encoding='utf-8'));check_hashes(evidence['sources'])
        ledger={c['query_id']:c for c in evidence['complexes']}
    if design.get('spatial_groups'):design['resolved_spatial_groups']=resolve_regions(design,survey)
    out=Path(output)
    out.mkdir(parents=True,exist_ok=False)
    order=selected_anchors(design);anchors={a['anchor_id']:a for a in survey['anchors']}
    _,expanded=_load_query(Path(survey['consensus_npz']))
    indices=np.array([anchors[a]['feature_index'] for a in order])
    match_query=tuple(expanded[k][indices] for k in ('feature_points','feature_types',
        'feature_directions','feature_direction_kinds','anchored_weights'))
    templates=[(t,_load_query(Path(t['query_npz']))[1]) for t in survey['templates'] if t['query_id'] in design['template_ids']]
    baseline=optional_weight_total(design) or 1.
    fixed=copy.deepcopy(design)
    fixed.update(optional_normalization='fixed_budget',optional_budget=design.get('optional_budget',baseline))
    rows=[];invariance=[];winner_counts={t['query_id']:0 for t,_ in templates}
    for meta in survey['prepared_complexes']:
        _,q=_load_query(Path(meta['query_npz']))
        _,assign,values=match_batch(match_query,q['feature_points'][None],q['feature_types'][None],
            q['feature_directions'][None],q['feature_direction_kinds'][None],np.eye(4)[None])
        candidate=SimpleNamespace(**q);seed=SimpleNamespace(transform_matrix=np.eye(4).reshape(-1))
        local=[]
        for t,query in templates:
            gaussian=float(same_pose_gaussian(query,candidate,[seed])[0])
            terms=pose_rank(values[0],assign[0],order,gaussian,q['shape_points'],np.eye(4),design)
            row=dict(query_id=meta['query_id'],template_id=t['query_id'],self_template=meta['query_id']==t['query_id'],**terms)
            rows.append(row);local.append(row)
        other=[r for r in local if not r['self_template']]
        if other:
            best=max(r['gaussian_same_pose'] for r in other)
            # Split exact ties, so a deterministic template ordering cannot create dominance.
            winners=[r for r in other if abs(r['gaussian_same_pose']-best)<=1e-12]
            for r in winners:winner_counts[r['template_id']]+=1/len(winners)
        before=pose_rank(values[0],assign[0],order,0.,q['shape_points'],np.eye(4),fixed)
        edited=copy.deepcopy(fixed)
        removed={a for a in edited['optional_weights'] if anchors[a]['feature_class']=='water_bridge'}
        for a in removed:del edited['optional_weights'][a]
        families=[g for g in edited.get('optional_groups',[]) if all(anchors[a]['feature_class']=='water_bridge' for a in g['anchor_ids'])]
        edited['optional_groups']=[g for g in edited.get('optional_groups',[]) if g not in families]
        after=pose_rank(values[0],assign[0],order,0.,q['shape_points'],np.eye(4),edited)
        scores={a:float(v) if i>=0 else 0. for a,v,i in zip(order,values[0],assign[0])}
        removed_numerator=sum(fixed['optional_weights'][a]*scores[a] for a in removed)+sum(
            g['weight']*max(scores[a] for a in g['anchor_ids']) for g in families)
        error=abs(before['optional_score']-after['optional_score']-removed_numerator/fixed['optional_budget'])
        invariance.append(dict(query_id=meta['query_id'],before=before['optional_score'],
            after_water_removal=after['optional_score'],unaffected_contribution_error=error))
    distributions=[]
    for t,q in templates:
        scores=[r['gaussian_same_pose'] for r in rows if r['template_id']==t['query_id'] and not r['self_template']]
        distributions.append(dict(template_id=t['query_id'],heavy_atoms=len(q['shape_points']),
            nonself_count=len(scores),nonself_quantiles=dict(zip(['min','q25','median','q75','max'],
                np.quantile(scores,[0,.25,.5,.75,1]).tolist())) if scores else {},
            nonself_best_template_fractional_count=winner_counts[t['query_id']]))
    selected=[]
    for aid in order:
        a=anchors[aid]
        observations=[]
        for observation in a['observations']:
            record={k:observation[k] for k in ('query_id','source_url','evidence') if k in observation}
            complex_record=ledger.get(observation['query_id'])
            if complex_record:
                atoms={atom['atom']:atom for atom in complex_record['ligand_atoms']}
                partner=observation.get('evidence',{}).get('protein_partner',{})
                record['nearby_atom_pairs']=[dict(ligand_atom=p['ligand_atom'],protein_atom=p['protein_atom'],
                    distance_angstrom=p['distance'],ligand_chemistry=atoms[p['ligand_atom']].get('chemistry',{}))
                    for p in complex_record['pairs'] if p['target_residue']==a['target_residue'] and
                    p['protein_atom'] in a['protein_atom'].split('/') and
                    (not partner.get('chain') or p['protein_chain']==partner['chain'])]
                record['chemistry_status']=complex_record['chemistry_status']
                record['solvent_exposure']='not_measured'
            observations.append(record)
        selected.append({k:a[k] for k in ('anchor_id','target_residue','protein_atom','feature_class',
            'distinct_structures','eligible_structures','distinct_chemotypes','pdb_ids')} |
            contact_semantics(a['protein_atom'],a['feature_class']) |
            dict(observations=observations))
    result=dict(kind='score_contract_audit',status='complete',scope='Crystal identity poses only; no pose search, library scan, gate acceptance or independent activity validation',
        design=design,selected_contact_evidence=selected,template_distributions=distributions,
        cross_template_calibration='not_performed; distributions are diagnostics on dependent crystal controls',
        fixed_budget_water_removal=dict(budget=fixed['optional_budget'],rows=invariance,
            passed=all(r['unaffected_contribution_error']<1e-10 for r in invariance)),
        pose_scores=rows,implementation_sources=fingerprint(sorted(Path(__file__).parent.glob('*.py'))),
        sources={**proposal['sources'],**survey['sources'],**fingerprint([source,survey_path])})
    _atomic_json(out/'report.json',result)
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--recommendation',required=True)
    parser.add_argument('--output',required=True)
    args=parser.parse_args();result=run(args.recommendation,args.output)
    print(json.dumps(dict(status=result['status'],report=str((Path(args.output)/'report.json').resolve()),
        fixed_budget_invariance_passed=result['fixed_budget_water_removal']['passed']),indent=2))
    if not result['fixed_budget_water_removal']['passed']:raise SystemExit(1)


if __name__=='__main__':main()
