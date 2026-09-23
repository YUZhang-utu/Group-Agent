"""Protein-contact policy; keep unmodeled mediator evidence out of screening."""
import argparse
import copy
import json
from pathlib import Path

from .contact_groups import optional_weight_total, selected_anchors
from .expanded_wee1 import fingerprint
from .gaussian_batch import _atomic_json
from .screening_selection import check_hashes

EXCLUDED_CLASSES={'water_bridge','metal_coordination'}


def require_protein_contacts(design, anchors):
    lookup={a['anchor_id']:a for a in anchors}
    excluded=[a for a in selected_anchors(design) if lookup[a].get('feature_class') in EXCLUDED_CLASSES]
    if excluded:
        raise ValueError('Unmodeled water/metal-mediated anchors are evidence-only; remove them from the design: '+', '.join(excluded))


def without_mediators(design, anchors):
    """Explicit migration; freeze original budget, retain all other coefficients."""
    result=copy.deepcopy(design)
    excluded={a['anchor_id'] for a in anchors if a.get('feature_class') in EXCLUDED_CLASSES}
    removed=sorted(set(selected_anchors(design)) & excluded)
    result['mandatory_anchors']=[a for a in design['mandatory_anchors'] if a not in excluded]
    result['alternative_groups']=[[a for a in g if a not in excluded] for g in design['alternative_groups']]
    result['alternative_groups']=[g for g in result['alternative_groups'] if g]
    result['optional_weights']={a:w for a,w in design['optional_weights'].items() if a not in excluded}
    if 'optional_groups' in design:
        for g in result['optional_groups']:g['anchor_ids']=[a for a in g['anchor_ids'] if a not in excluded]
        result['optional_groups']=[g for g in result['optional_groups'] if g['anchor_ids']]
    evidence={a['evidence_id'] for a in anchors if a['anchor_id'] in excluded}
    result['evidence_ids']=[e for e in design['evidence_ids'] if e not in evidence]
    if result.get('optional_normalization')!='fixed_budget':
        result.update(optional_normalization='fixed_budget',optional_budget=optional_weight_total(design) or 1.)
        total=sum(result.get(k,0.) for k in ('gaussian_weight','optional_weight','occupancy_weight'))
        for key in ('gaussian_weight','optional_weight','occupancy_weight'):
            if key in result:result[key]/=total
    result['rationale']='Protein-contact-only revision. Unmodeled water/metal-mediated anchors removed from scoring and eligibility; original contact budget retained. Prior scientific rationale requires review.'
    require_protein_contacts(result,anchors)
    return result,removed


def migrate(source, output):
    from .consensus_design import validate, scoring_contract
    source=Path(source);proposal=json.loads(source.read_text(encoding='utf-8'))
    if proposal.get('kind')!='consensus_recommendation':raise ValueError('Expected saved recommendation report')
    check_hashes(proposal['sources']);survey_path=Path(proposal['survey'])
    if str(survey_path.resolve()) not in proposal['sources']:raise ValueError('Unsealed survey')
    survey=json.loads(survey_path.read_text(encoding='utf-8'));check_hashes(survey['sources'])
    validate(proposal['recommendation'],survey)
    design,removed=without_mediators(proposal['recommendation'],survey['anchors'])
    validate(design,survey)
    out=Path(output);out.mkdir(parents=True,exist_ok=False)
    result=dict(kind='consensus_recommendation',status='complete',readiness='proposal_ready',
        recommendation=design,survey=str(survey_path.resolve()),scoring_contract=scoring_contract(design),
        migration=dict(source=str(source.resolve()),removed_anchor_ids=removed,
            original_rationale=proposal['recommendation']['rationale'],
            scope='New proposal only; no adoption or library scan. Existing explicit pose threshold retained; null is recalculated on adoption.'),
        sources={**proposal['sources'],**survey['sources'],**fingerprint([source,survey_path])})
    _atomic_json(out/'report.json',result)
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--recommendation',required=True);parser.add_argument('--output',required=True)
    args=parser.parse_args();result=migrate(args.recommendation,args.output)
    print(json.dumps(dict(report=str((Path(args.output)/'report.json').resolve()),
        removed_anchor_ids=result['migration']['removed_anchor_ids'],
        optional_budget=result['recommendation']['optional_budget']),indent=2))


if __name__=='__main__':main()
