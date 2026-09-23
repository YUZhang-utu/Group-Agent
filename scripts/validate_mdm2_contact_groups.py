"""E054 real-structure regression case; no MDM2 rules enter the generic engine."""
import argparse
from collections import Counter
import csv
import json
from pathlib import Path

import numpy as np

from aidd_agent.contact_groups import resolve_regions, grouped_terms
from aidd_agent.consensus_design import validate, adopt
from aidd_agent.expanded_wee1 import fingerprint


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--survey',type=Path,required=True)
    parser.add_argument('--audit',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--fresh-survey',type=Path)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    survey=json.loads(args.survey.read_text())
    if survey['target']['accession']!='Q00987':raise ValueError('This validation case is MDM2-specific')
    ledger_path=Path(survey['contact_evidence']['path']);ledger=json.loads(ledger_path.read_text())
    raw=list(csv.DictReader((args.audit/'atom_contacts.csv').open()))
    expected={(r['query_id'],r['ligand_atom'],int(r['canonical_residue']),r['protein_atom']):float(r['distance']) for r in raw}
    actual={(c['query_id'],p['ligand_atom'],p['target_residue'],p['protein_atom']):p['distance'] for c in ledger['complexes'] for p in c['pairs']}
    if set(actual)!=set(expected) or any(abs(actual[k]-expected[k])>1e-6 for k in expected):
        raise ValueError('General ledger differs from independent coordinate audit')
    positions=list(csv.DictReader((args.audit/'spatial_atoms.csv').open()))
    reference='5C5A:NUT:A:201'
    groups=[dict(id=site,reference_query=reference,ligand_atoms=[r['ligand_atom'] for r in positions
        if r['query_id']==reference and r['site_at_2A']==site],radius=2.,minimum_atoms=1)
        for site in ('A_Trp23','B_Leu26','C_Phe19')]
    # Deterministic regression fixture, not a live recommendation or user-adopted design.
    weights={'70725df96c9e77ca':1.,'6204d67aa4613d5f':.8,'db6f8400e16ed844':.75,
        'bf17c6a4b6ed6967':.7,'632e01bc9c9f9911':.7,'85bacca7ef8f747b':.65,
        '93872f9b7eb7cb3b':.6,'a6ad0b8bb793e330':.6,'ff41e62948bdd497':.55,
        'baf993010192c2fa':.2,'7974cd0635c375f2':.2,'f4b886c1f1653aa8':.2,
        'a6190c97bee1efc7':.5,'af3e105fd2419ceb':.45,'7cded8b9c11d5b02':.4,
        'a5645004d8ee0341':.4,'ba545d5498c51d6f':.4,'da95c0f6270aefa7':.35,
        '5176049acb85ae4d':.35,'a3be914565446540':.2}
    families=[dict(id='His_pi',anchor_ids=['pocket:'+a for a in ['632e01bc9c9f9911','85bacca7ef8f747b']],weight=.6),
        dict(id='His_polar',anchor_ids=['pocket:'+a for a in ['93872f9b7eb7cb3b','af3e105fd2419ceb']],weight=.55),
        dict(id='Ile_modes',anchor_ids=['pocket:'+a for a in ['70725df96c9e77ca','7cded8b9c11d5b02','a5645004d8ee0341','ba545d5498c51d6f']],weight=.4)]
    family_members={a for g in families for a in g['anchor_ids']}
    design=dict(mandatory_anchors=[],alternative_groups=[],
        optional_weights={'pocket:'+a:w for a,w in weights.items() if 'pocket:'+a not in family_members},
        optional_groups=families,spatial_groups=groups,occupancy_rewards=[0.,.1,.3,1.],occupancy_weight=.3,
        spatial_ambiguity=.5,template_ids=['5C5A:NUT:A:201','4JV9:1MO:A:201'],
        evidence_ids=['consensus:pocket:'+a for a in weights],
        rationale='Deterministic E054 geometry regression. Coefficients are exploratory fixture values, not adopted user thresholds.',
        exclusions=[],permissiveness=1.,gaussian_weight=.7,optional_weight=.3,minimum_pose_score=None)
    validate(design,survey)
    resolved=resolve_regions(design,survey)
    geometries=[]
    for c in ledger['complexes']:
        moved=np.array([a['point'] for a in c['ligand_atoms'] if not a['quality_issues']])
        terms=grouped_terms([],[],[],moved,dict(design,optional_groups=[],resolved_spatial_groups=resolved))
        geometries.append(dict(query_id=c['query_id'],**terms))
    proposal=dict(kind='consensus_recommendation',status='complete',readiness='proposal_ready',
        recommendation=design,survey=str(args.survey.resolve()),
        model=dict(mode='deterministic_regression_fixture_not_live'),
        sources={**survey['sources'],**fingerprint([args.survey])})
    proposal_path=args.output/'fixture-proposal.json';proposal_path.write_text(json.dumps(proposal,indent=2))
    adopted=adopt(proposal_path,args.output/'fixture-adopted')
    if adopted['readiness']!='needs_template_state_review' or '4JV9:1MO:A:201' not in adopted['failed_template_self_controls']:
        raise ValueError('Unsafe template did not trigger the state review gate')
    safe=adopt(proposal_path,args.output/'fixture-compatible-template',overrides={'template_ids':[reference]})
    if safe['readiness']!='ready_for_consensus_funnel':raise ValueError('Compatible reference template failed adoption')
    from aidd_agent.consensus_admission import read_structure,admit,POLICY
    cached={}
    def structure(row):
        path=row['structure_path']
        if path not in cached:cached[path]=read_structure(Path(path),survey['target'])
        return cached[path]
    refrow=next(r for r in survey['cohort']['admitted'] if r['query_id']==reference)
    reference_structure=structure(refrow)
    rechecked=[admit(reference_structure,refrow,refrow['admission']['target_chain'],structure(r),r,r['admission']['target_chain'])
               for r in survey['cohort']['admitted'] if r['query_id'] in {c['query_id'] for c in ledger['complexes']}]
    if any(r['status']=='admitted' for r in rechecked if r['query_id'] in {'4JV7:1MN:A:201','4JV9:1MO:A:201'}):
        raise ValueError('Known severe reference-state clashes were admitted')
    result=dict(status='complete',scope='Exploratory local regression, no production adoption or library scan',
        independent_pair_comparison=dict(expected=len(expected),actual=len(actual),identical=True),
        ligand_chemistry= dict(Counter(c['chemistry_status'] for c in ledger['complexes'])),
        contact_classes=dict(Counter(k for c in ledger['complexes'] for p in c['pairs'] for k in p['classes'])),
        regions=resolved,spatial_reference='Reviewed 5C5A atom selections based on the earlier p53 probe audit',
        reference_dependence='These are 5C5A-based regions, not the original 1YCR probe coordinates',
        geometry=geometries,spatial_counts=dict(Counter(g['occupied_group_count'] for g in geometries)),
        crystal_self_controls=dict(total=len(adopted['crystal_self_controls']),
            passed=sum(c['rule_passed'] and c['pocket_passed'] for c in adopted['crystal_self_controls']),
            derived_pose_threshold=adopted['design']['minimum_pose_score']),
        template_review=dict(incompatible_design_readiness=adopted['readiness'],
            failed_templates=adopted['failed_template_self_controls'],compatible_design_readiness=safe['readiness']),
        admission_recheck=dict(policy=POLICY,counts=dict(Counter(r['status'] for r in rechecked)),decisions=rechecked),
        sources=fingerprint([args.survey,ledger_path,args.audit/'atom_contacts.csv',args.audit/'spatial_atoms.csv',Path(__file__)]))
    if args.fresh_survey:
        fresh=json.loads(args.fresh_survey.read_text())
        chosen=sorted(fresh['anchors'],key=lambda a:(-a['frequency'],a['anchor_id']))[:20]
        ids=[a['anchor_id'] for a in chosen]
        fresh_design=dict(design,template_ids=[reference],optional_weights={a:1. for a in ids[2:]},
            optional_groups=[dict(id='fixture_alternatives',anchor_ids=ids[:2],weight=.6)],
            evidence_ids=[a['evidence_id'] for a in chosen],
            rationale='Fresh-cohort executable fixture; frequency-selected optional anchors only; not an adopted scientific recommendation.')
        fresh_proposal=dict(proposal,recommendation=fresh_design,survey=str(args.fresh_survey.resolve()),
            sources={**fresh['sources'],**fingerprint([args.fresh_survey])})
        fresh_path=args.output/'fresh-fixture-proposal.json';fresh_path.write_text(json.dumps(fresh_proposal,indent=2))
        fresh_adopted=adopt(fresh_path,args.output/'fresh-fixture-adopted')
        if fresh_adopted['readiness']!='ready_for_consensus_funnel':raise ValueError('Fresh compatible-cohort fixture not ready')
        result['fresh_cohort_fixture']=dict(prepared=len(fresh['prepared_complexes']),anchors=len(fresh['anchors']),
            templates=len(fresh['templates']),readiness=fresh_adopted['readiness'],
            crystal_controls_passed=sum(c['rule_passed'] and c['pocket_passed'] for c in fresh_adopted['crystal_self_controls']),
            cutoff=fresh_adopted['design']['minimum_pose_score'])
        result['sources'].update(fingerprint([args.fresh_survey]))
    (args.output/'report.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps({k:result[k] for k in ('status','independent_pair_comparison','ligand_chemistry','contact_classes','spatial_counts','crystal_self_controls')},indent=2))


if __name__=='__main__':main()
