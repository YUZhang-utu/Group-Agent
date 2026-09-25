import json
from pathlib import Path

from aidd_agent.chat_agent import ChatAgent
from aidd_agent.domain_tools import DomainTools
from aidd_agent.prompt_workflow import project_root


def add(app,sid,jid,action,child,source=None,run_id=None):
    ctx=app.context;root=project_root(Path(ctx['db']),ctx['user_id'],ctx['project_id'])/(run_id or jid)
    root.mkdir();plan=root/'plan.json';report=root/'report.json';detail=root/'detail.json'
    plan.write_text(json.dumps(dict(plan=dict(steps=[dict(action=action,params=dict(source_run=source) if source else {})]))))
    detail.write_text(json.dumps(child))
    report.write_text(json.dumps(dict(steps=dict(stage=dict(action=action,status='complete',result=dict(report=str(detail)))))))
    with app.connect() as db:
        db.execute('INSERT INTO jobs(id,session,status,plan,report,request,created) VALUES(?,?,?,?,?,?,?)',
            (jid,sid,'complete',str(plan),str(report),'MDM2 '+action,1))
    return root


def test_chain_preserves_sources_and_actual_diversity_evidence(tmp_path):
    app=ChatAgent(tmp_path,start=False);sid=app.new_session()
    root=add(app,sid,'div','structure_diversity',dict(status='complete',target='P49137',
        proposed_references=[dict(query_id='6Q9L:HTZ:A:201')],selection=dict(fingerprint='Morgan')))
    (root/'ligand-similarity.csv').write_text('query_id,6Q9L:HTZ:A:201\n6Q9L:HTZ:A:201,1\n')
    add(app,sid,'con','structure_consensus',dict(status='complete',readiness='proposal_ready'),'div')
    add(app,sid,'rec','consensus_recommend',dict(status='complete',readiness='proposal_ready'),'con')
    add(app,sid,'search','consensus_budget',dict(status='complete',exported_molecules=100000),'rec')
    tools=DomainTools(app,sid,'Continue my workflow','deepseek')
    rows={r['task_id']:r for r in tools.call('structure_workflow',{})['tasks']}
    assert rows['con']['source_task_ids']==['div']
    assert rows['search']['source_task_ids']==['rec']
    assert rows['rec']['stages'][0]['next_intents']==['adopt','design','budget']
    assert rows['search']['stages'][0]['next_intents']==['budget_page']
    stage=rows['div']['stages'][0]
    assert stage['proposed_references']==['6Q9L:HTZ:A:201']
    aid=next(a['artifact_id'] for a in stage['artifacts'] if a['name']=='ligand-similarity.csv')
    assert tools.call('read_artifact',dict(artifact_id=aid))['content']['items'][0]['query_id']=='6Q9L:HTZ:A:201'
    assert len(app.jobs(sid))==4


def test_unresolved_reference_or_failed_child_has_no_next_execution(tmp_path):
    app=ChatAgent(tmp_path,start=False);sid=app.new_session()
    add(app,sid,'con','structure_consensus',dict(status='complete',readiness='needs_reference_instance',
        cohort=dict(reference_options=[dict(query_id='actual-reference')])),'div')
    add(app,sid,'bad','consensus_recommend',dict(status='failed'))
    rows=DomainTools(app,sid,'inspect','deepseek').call('structure_workflow',{})['tasks']
    assert all(not r['stages'][0]['next_intents'] for r in rows)
    con=next(r for r in rows if r['task_id']=='con')
    assert con['stages'][0]['reference_options'][0]['query_id']=='actual-reference'
    other=app.new_session()
    assert not DomainTools(app,other,'inspect','deepseek').call('structure_workflow',{})['tasks']


def test_report_path_escape_is_reported_without_offering_transition(tmp_path):
    app=ChatAgent(tmp_path,start=False);sid=app.new_session()
    root=add(app,sid,'bad','structure_consensus',dict(status='complete',readiness='proposal_ready'))
    outside=tmp_path/'outside.json';outside.write_text('{}')
    (root/'report.json').write_text(json.dumps(dict(steps=dict(stage=dict(action='structure_consensus',
        status='complete',result=dict(report=str(outside)))))))
    row=DomainTools(app,sid,'inspect','deepseek').call('structure_workflow',{})['tasks'][0]
    assert 'error' in row and row['stages'][0]['next_intents']==[]
