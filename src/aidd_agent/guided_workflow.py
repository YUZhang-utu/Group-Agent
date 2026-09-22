"""Trusted coordinator adapters for the structure-guided chat sequence."""
from pathlib import Path
import copy
import json
import sqlite3

from . import library_acceptance as ev
from .expanded_wee1 import fingerprint
from .gaussian_batch import _atomic_json
from .screening_selection import check_hashes

SOURCE_ACTIONS={'anchor_recommend':'structure_survey','anchor_design':'anchor_recommend',
                'guided_funnel':'anchor_design','guided_select':('guided_funnel','consensus_funnel'),
                'structure_consensus':'structure_diversity','consensus_recommend':'structure_consensus',
                'consensus_design':'consensus_recommend','consensus_funnel':'consensus_design'}


def select_candidates(source, output, params):
    """Select actual retained poses; never intersect anchors across poses."""
    report=ev.read(source);check_hashes(report['output_hashes'])
    order=report['anchor_order'];required=params['required_anchors']
    if not required or not set(required)<=set(order):raise ValueError('Unknown selection anchors')
    if params['minimum_score'] != report['policy']['minimum_score']:
        raise ValueError('Changing the search threshold requires a new full-library design')
    if params.get('coarse_constraints') is not None:
        raise ValueError('Post-search selection cannot change the sealed coarse eligibility rules')
    if params.get('max_molecules') is not None:
        raise ValueError('This selection preserves all matching members; no implicit candidate cap')
    cols=[order.index(a) for a in required]
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    selected=output/'selected-poses.jsonl';ids=output/'selected-molecules.txt'
    count=0;poses=0
    with sqlite3.connect(Path(report['outputs']['candidates.sqlite']).resolve().as_uri()+'?mode=ro',uri=True) as db, \
         selected.open('w',encoding='utf-8') as stream, ids.open('w',encoding='utf-8') as members:
        last=None
        for mid,payload,cluster in db.execute('SELECT p.mid,p.payload,m.cluster FROM poses p JOIN members m ON p.mid=m.mid ORDER BY p.mid,p.mask'):
            row=json.loads(payload)
            hits=[row['assignments'][i]>=0 and row['scores'][i]>=params['minimum_score'] for i in cols]
            if not (all(hits) if params['match_mode']=='all' else any(hits)):continue
            row['cluster_id']=cluster;stream.write(json.dumps(row)+'\n');poses+=1
            if mid!=last:members.write(mid+'\n');count+=1;last=mid
    result=dict(kind='guided_selection',status='complete',matching_molecules=count,pose_records=poses,
        selection={k:params[k] for k in ('required_anchors','match_mode','minimum_score')},
        outputs=dict(poses=str(selected.resolve()),molecule_ids=str(ids.resolve())),
        sources=fingerprint([source,*report['output_hashes']]),
        scope='Same-threshold selection among adopted-design survivors; original mandatory rules remain in force. No new search, SDF export or docking.')
    _atomic_json(output/'report.json',result);return result


def execute(action,source,output,params,cfg,allow_compute):
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    if action=='structure_consensus':
        from .consensus_model import build
        return build(source,output,params.get('reference_query'),params.get('target_chain'),params.get('maximum_templates',8))
    if action in {'consensus_recommend','consensus_design'}:
        from .consensus_design import recommend,adopt
        return recommend(source,output,params['provider']) if action=='consensus_recommend' else adopt(source,output,params.get('design'))
    if action=='consensus_funnel':
        from .prompt_workflow import Blocked
        from .consensus_funnel import run
        if not allow_compute:raise Blocked('Enable compute for the full-library consensus scan')
        if not cfg.get('search'):raise Blocked('Configure the trusted library runtime')
        return run(source,output,cfg['search'])
    if action=='anchor_recommend':
        from .anchor_advice import recommend
        return recommend(source,output,params['provider'])
    if action=='anchor_design':
        from .anchor_advice import adopt
        return adopt(source,output,params.get('design'))
    if action=='guided_select':
        return select_candidates(source,output,params)
    if action!='guided_funnel':raise ValueError('Unknown guided action')
    from .prompt_workflow import Blocked
    if not allow_compute:raise Blocked('Enable compute for the requested full-library guided scan')
    search=cfg.get('search')
    if not search:raise Blocked('Configure the existing trusted library search runtime before scanning')
    adopted=ev.read(source);check_hashes(adopted['sources'])
    if adopted.get('readiness')!='ready_for_guided_funnel':raise ValueError('Adopt a valid search design first')
    batch=Path(search['batch'])
    _,acceptance=ev.accept_library(batch,full=False)
    q=copy.deepcopy(adopted['query'])
    q.update(artifact_catalog=str((batch/'artifacts/catalog.json').resolve()),
             chemical_companion=str((batch/'chemical/catalog.json').resolve()))
    inputs=output/'inputs';inputs.mkdir(exist_ok=True)
    sources={**adopted['sources'],**fingerprint([source,Path(q['artifact_catalog']),Path(q['chemical_companion'])])}
    classification=dict(kind='screening_evidence',status='complete',queries=[q],sources=sources,
        library={k:acceptance[k] for k in ('library_conformers','library_molecules')},
        classification='structure-guided-crystal-hypotheses-v1')
    classified=inputs/'classification.json';_atomic_json(classified,classification)
    design=adopted['design']
    policy=dict(required_anchors=adopted['anchor_order'],match_mode='any',minimum_score=design['minimum_score'],
                coarse_constraints=design['coarse_constraints'])
    selection=inputs/'selection.json'
    _atomic_json(selection,dict(kind='selection_preview',policy=policy,query=q,evidence_report=str(classified.resolve()),
                                guided_design=design,sources={**sources,**fingerprint([classified])}))
    from .preselection_full import run
    report=run(selection,output/'full-library',workers=search['workers'],chunk_size=search['refine_chunk'])
    report['kind']='guided_funnel';report['design']=design
    report['gaussian_scope']='Ranking over geometry/pocket-surviving seeds; differs from unrestricted Gaussian winner ranking'
    _atomic_json(output/'report.json',report)
    return report
