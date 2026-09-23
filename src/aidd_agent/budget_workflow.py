"""Trusted project-scoped Chat adapters for contact-first budget search."""
from pathlib import Path
from types import SimpleNamespace
import os

from . import library_acceptance as ev
from .gaussian_batch import _atomic_json
from .expanded_wee1 import fingerprint
from .screening_selection import check_hashes


def validate_budget(value, page=False):
    allowed={'export_molecules','start_rank'} if page else {'retrieval_molecules','export_molecules'}
    if not isinstance(value,dict) or set(value)-allowed:raise ValueError('Unknown budget fields')
    if any(type(v) is not int or not 1<=v<=1000000000 for v in value.values()):raise ValueError('Budgets must be positive integer molecule counts')
    if not page and value.get('export_molecules',100000)>value.get('retrieval_molecules',1000000):
        raise ValueError('Export budget exceeds retrieval molecule budget')
    return value


def execute(action,source,output,params,cfg,allow_compute):
    from .prompt_workflow import Blocked
    from .budget_screen import execute as screen,validate_regions
    from .budget_export import export
    if not allow_compute:raise Blocked('Enable --allow-compute for budget search and MOL2 handoff')
    search=cfg.get('search',{})
    if not search.get('batch'):raise Blocked('Configure search.batch for the current library')
    output=Path(output);output.mkdir(parents=True,exist_ok=True);source=Path(source)
    batch=Path(search['batch']);previous=ev.read(source);page=action=='budget_page'
    settings={k:v for k,v in params.items() if k!='source_run'};validate_budget(settings,page)
    if page:
        if previous.get('kind') not in {'consensus_budget','budget_page'}:raise ValueError('Expected completed budget search or page')
        check_hashes(previous['sources'])
        run=Path(previous['run_directory']);start=settings.get('start_rank',previous['end_rank']+1)
        count=settings.get('export_molecules',100000)
        if start<=previous['end_rank']:raise ValueError('Requested page overlaps the source delivery; choose a later start rank')
        target=previous['target'];migration=previous.get('migration')
    else:
        if previous.get('kind') not in {'consensus_recommendation','consensus_design'}:raise ValueError('Expected consensus recommendation or design')
        check_hashes(previous['sources']);prepared=source;migration=None
        if previous['kind']=='consensus_recommendation':
            from .contact_policy import migrate
            revised=output/'protein-only-proposal/report.json'
            if not revised.exists():migrate(source,revised.parent)
            prepared=revised;proposal=ev.read(revised);check_hashes(proposal['sources']);migration=proposal['migration']
            survey=ev.read(proposal['survey']);target=survey['target'];reference=survey['cohort']['reference']
        else:
            from .contact_policy import without_mediators
            from .consensus_design import scoring_contract
            target=previous['target'];reference=previous['reference']
            revised,removed=without_mediators(previous['design'],previous['anchors'])
            if removed:
                import copy
                migrated=copy.deepcopy(previous);migrated['design']=revised
                migrated['anchors']=[a for a in previous['anchors'] if a['anchor_id'] not in removed]
                migrated['anchor_order']=[a for a in previous['anchor_order'] if a not in removed]
                migrated['scoring_contract']=scoring_contract(revised)
                migrated['sources'].update(fingerprint([source]))
                migration=dict(source=str(source),removed_anchor_ids=removed,scope='Protein-only copy for budget ranking; original design preserved')
                migrated['migration']=migration;prepared=output/'protein-only-design.json'
                _atomic_json(prepared,migrated)
        definitions=Path(search['budget_regions']) if search.get('budget_regions') else None
        if definitions is None:
            fixture=Path(__file__).resolve().parents[2]/'to_human/E058_MDM2_DUAL_REGIONS.json'
            if fixture.is_file():
                try:validate_regions(dict(target=target,reference=reference),ev.read(fixture));definitions=fixture
                except ValueError:pass
        run=output/'search';start=1;count=settings.get('export_molecules',100000)
        os.environ.setdefault('AIDD_ASSIGNMENT_BACKEND','numba')
        args=SimpleNamespace(batch=batch,recommendation=prepared,definitions=definitions,output=run,
            retrieval_molecules=settings.get('retrieval_molecules',1000000),workers=search.get('workers',24),
            chunk_conformers=search.get('refine_chunk',64),nprobe=128,template_quota=count,rrf_k=60,retrieve_only=False)
        screen(args)
    result=export(batch,run,output/f'page-{start:09d}',start,count)
    ranked=ev.read(run/'report.json')
    summary=dict(kind='budget_page' if page else 'consensus_budget',status='complete',target=target,
        run_directory=str(run.resolve()),ranked_molecules=ranked['ranked_molecules'],ranking=ranked['ranking'],
        exported_molecules=result['exported_molecules'],requested_molecules=count,start_rank=start,end_rank=result['end_rank'],
        shortfall=result['shortfall'],review_required=result['review_required'],migration=migration,
        outputs=dict(ranking=str(run/'ranking.sqlite'),page=str(output/f'page-{start:09d}'),
            names=str(output/f'page-{start:09d}'/'molecules.csv'),evidence=str(output/f'page-{start:09d}'/'molecules.jsonl')),
        sources=fingerprint([source,run/'report.json',output/f'page-{start:09d}'/'report.json']),
        limitations=['ANN recall and activity enrichment are unvalidated for this target.',
            'Contact-first ranking is docking priority, not evidence of activity.'])
    _atomic_json(output/'report.json',summary);return summary
