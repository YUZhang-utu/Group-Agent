"""Read-only navigation across the structure-to-library workflow."""
from pathlib import Path

STAGES = {
    'pdb_fetch': [], 'structure_diversity': ['consensus'],
    'structure_consensus': ['recommend'],
    'consensus_recommend': ['adopt', 'design', 'budget'],
    'consensus_design': ['budget', 'guided'], 'consensus_funnel': [],
    'consensus_budget': ['budget_page'], 'budget_page': ['budget_page'],
}


def inspect(tools, args):
    from .project_context import ensure_within
    requested=args.get('task_id')
    jobs=[tools.app.task(tools.sid,requested)] if requested else tools.app.jobs(tools.sid)[-100:]
    rows=[]
    for job in jobs:
        row=dict(task_id=job['id'],task_status=job['status'],request=job['request'],stages=[])
        try:
            if job.get('plan'):
                path=ensure_within(Path(job['plan']),tools.project);plan=tools.json_file(path)
                row['run_id']=path.parent.name
                row['source_run_ids']=list(dict.fromkeys(s['params']['source_run']
                    for s in plan.get('steps',[]) if s.get('params',{}).get('source_run')))
            if job.get('report'):
                report=tools.json_file(ensure_within(Path(job['report']),tools.project))
                for step_id,step in report.get('steps',{}).items():
                    action=step.get('action')
                    if action not in STAGES:continue
                    stage=dict(step_id=step_id,action=action,status=step.get('status'),next_intents=[])
                    row['stages'].append(stage)
                    if step.get('status')!='complete':continue
                    child_path=step.get('result',{}).get('report')
                    child=tools.json_file(ensure_within(Path(child_path),tools.project)) if child_path else {}
                    stage['evidence']={k:child[k] for k in (
                        'target','reference_site_pdb','readiness','selection','quality_site_unique_ligands',
                        'coverage_curve','template_selection','proposed_template_ids','ranked_molecules',
                        'exported_molecules','start_rank','end_rank','shortfall','review_required',
                        'limitations','outputs') if k in child}
                    stage['reference_options']=child.get('cohort',{}).get('reference_options',[])
                    stage['proposed_references']=[r.get('query_id') for r in child.get('proposed_references',[])]
                    stage['artifacts']=tools.discover(child,tools.project)
                    if child_path:
                        stage['artifacts'].append(tools.register(child_path,tools.project))
                        if action=='structure_diversity':
                            matrix=tools.register(Path(child_path).with_name('ligand-similarity.csv'),tools.project)
                            if matrix:stage['artifacts'].append(matrix)
                    ready=child.get('readiness')
                    stage['next_intents']=STAGES[action] if job['status']=='complete' and child.get('status')=='complete' else []
                    if ready and ready.startswith('needs_'):
                        stage['next_intents']=[];stage['required_input']=ready
                    if action=='structure_consensus' and ready!='proposal_ready':stage['next_intents']=[]
        except (OSError,ValueError,KeyError,TypeError) as exc:
            row['error']=f'{type(exc).__name__}: {exc}'
            for stage in row['stages']:stage['next_intents']=[]
        if row['stages'] or requested or job['status'] in {'queued','planning','running'}:rows.append(row)
    by_run={r['run_id']:r['task_id'] for r in rows if r.get('run_id')}
    for row in rows:row['source_task_ids']=[by_run[r] for r in row.get('source_run_ids',[]) if r in by_run]
    return dict(tasks=rows,scope='Current-session tasks; no computation dispatched; source_run IDs preserve branches',
        instructions=[
            'Use the exact source task ID. Never join unrelated targets or branches. Ask if several branches fit.',
            'Diversity uses Morgan/Tanimoto chemical differences; consensus performs protein/pocket alignment and admission.',
            'needs_reference_instance requires reference_query and target_chain. After selection, rerun consensus from the original diversity source.',
            'Budget queries the full-library ANN index but scores retrieved molecules only, not every library pose. Use budget for budgeted delivery, guided only for a threshold funnel.',
            'Budget pages reuse rankings. Select the latest page of the intended branch to avoid repeated delivery.',
            'Supported next intents are not authorization. Preserve user requirements and all scientific prerequisites.'])
