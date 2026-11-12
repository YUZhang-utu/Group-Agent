"""Native saved-review summaries bypass JSON planning, with bounded local reads."""
import io
import json

import pytest

from aidd_agent.chat_agent import ChatAgent
from aidd_agent.domain_agent import run
from aidd_agent import project_knowledge as knowledge
from aidd_agent.saved_review_flow import SavedReviewFlow,matches
from aidd_agent.literature_recovery import progress,incomplete_message


@pytest.fixture
def project(tmp_path):
    profiles=tmp_path/'profiles';profiles.mkdir()
    (profiles/'deepseek.json').write_text(json.dumps(dict(base_url='http://localhost:1',model='fixture',json_mode=True)))
    app=ChatAgent(tmp_path/'work',config_dir=profiles,start=False);sid=app.new_session()
    knowledge.ingest(app,sid,'private.txt',b'PRIVATE_FLOW_SENTINEL')
    yield app,sid
    app.close()


def seed_review(app,sid):
    doc=knowledge.ingest(app,sid,'CIP2A.txt',b'Public CIP2A observation.','public','https://example.org/CIP2A',True)
    calls=[]
    def planner(*a,**k):
        calls.append(1)
        if len(calls)==1:return dict(kind='tool',tool='project_knowledge',arguments=dict(operation='read',document_id=doc['document_id'],page=1),purpose='Read paper'),{}
        return dict(kind='answer',answer='CIP2A unresolved mechanism and assay limitations. [E1]',evidence_ids=['E1']),{}
    run(app,sid,'Read public paper','deepseek',planner=planner)
    return doc['document_id']


def install_native_transport(monkeypatch,fail=False):
    from aidd_agent import domain_agent,prompt_plan
    monkeypatch.setattr(domain_agent,'chat_plan',lambda *a,**k:pytest.fail('Saved-review flow must not invoke the JSON planner'))
    captured=[]
    class Transport:
        def open(self,request,timeout):
            payload=json.loads(request.data);captured.append(payload)
            assert 'response_format' not in payload and 'tools' not in payload
            assert 'PRIVATE_FLOW_SENTINEL' not in request.data.decode()
            if fail:return io.BytesIO(b'not a valid provider envelope')
            context=json.loads(payload['messages'][1]['content'])
            citation=next(r['id'] for r in context['evidence'] if r['result'].get('summary'))
            assert 'Valid source citations for this answer: ['+citation+']' in payload['messages'][0]['content']
            answer='Saved CIP2A interpretation with unresolved limitations. ['+citation+']'
            return io.BytesIO(json.dumps(dict(choices=[dict(finish_reason='stop',message=dict(content=answer))])).encode())
    monkeypatch.setattr(prompt_plan,'build_opener',lambda *a:Transport())
    return captured


def test_native_summary_reads_then_synthesizes_without_json_decisions(project,monkeypatch):
    app,sid=project;seed_review(app,sid);other=app.new_session();wire=install_native_transport(monkeypatch)
    reply=run(app,other,"Use this project's saved literature reviews to summarize what we have learned about CIP2A. Do not launch calculations.",'deepseek')
    assert len(wire)==1 and 'Saved CIP2A interpretation' in reply and '[E' not in reply
    trace=json.loads(next((app.root/'agents'/other).glob('*.json')).read_text())
    assert trace['answer_mode']=='saved_review_text_synthesis'
    assert trace['saved_review_coordination']['summary_characters_read']>0
    assert [r['step']['arguments']['operation'] for r in trace['steps'] if r.get('step',{}).get('kind')=='tool']==['list','read']
    assert not app.jobs(other)


def test_native_summary_adopts_older_reviews_without_rereading_papers(project,monkeypatch):
    app,sid=project;seed_review(app,sid)
    with knowledge.connect(knowledge.folder(app,sid,'public')) as db:db.execute('DELETE FROM literature_reviews')
    other=app.new_session();wire=install_native_transport(monkeypatch)
    reply=run(app,other,'Summarize saved literature reviews for CIP2A','deepseek')
    trace=json.loads(next((app.root/'agents'/other).glob('*.json')).read_text())
    operations=[r['step']['arguments']['operation'] for r in trace['steps'] if r.get('step',{}).get('kind')=='tool']
    assert 'index' in operations and 'read' in operations and len(wire)==1
    assert 'Saved CIP2A interpretation' in reply and not app.jobs(other)


@pytest.mark.parametrize('archived',[False,True])
def test_empty_or_archived_reviews_never_trigger_a_summary_call(project,monkeypatch,archived):
    app,sid=project
    if archived:
        doc=seed_review(app,sid);knowledge.control(app,sid,dict(operation='archive',collection='public',document_id=doc))
    other=app.new_session();wire=install_native_transport(monkeypatch)
    reply=run(app,other,'Summarize saved literature reviews for CIP2A','deepseek')
    assert not wire and 'No saved summary text was read' in reply
    assert 'Located 0 public document(s)' not in reply and not app.jobs(other)
    trace=json.loads(next((app.root/'agents'/other).glob('*.json')).read_text())
    assert trace['failure_stage']=='literature_retrieval'


def test_final_provider_failure_reports_real_saved_text_progress(project,monkeypatch):
    app,sid=project;seed_review(app,sid);other=app.new_session();wire=install_native_transport(monkeypatch,fail=True)
    reply=run(app,other,'Summarize saved literature reviews for CIP2A','deepseek')
    assert len(wire)==1 and 'Read 1 saved public review(s)' in reply
    assert 'Located 0 public document(s)' not in reply
    trace=json.loads(next((app.root/'agents'/other).glob('*.json')).read_text())
    assert trace['failure_stage']=='answer_synthesis'


def test_flow_caps_characters_and_records_and_preserves_partial_receipts():
    flow=SavedReviewFlow('Summarize saved reviews')
    first=flow.next_call();flow.accept(first,dict(reviews=[dict(review_id=str(i)*32) for i in range(1,4)],has_more=True),'complete')
    count=0
    while (args:=flow.next_call()) is not None:
        count+=1
        flow.accept(args,dict(summary='x'*5000,has_more=True,next_offset=args['offset']+5000),'complete')
    report=flow.report()
    assert count==4 and report['summary_characters_read']==20000 and len(report['selected_review_ids'])==3
    assert report['partial_review_ids'] and report['more_review_records_available']


def test_indexing_is_bounded_and_failed_reads_are_not_retried():
    flow=SavedReviewFlow('Summarize saved reviews');calls=[]
    while (args:=flow.next_call()) is not None:
        calls.append(args)
        flow.accept(args,dict(reviews=[],has_more=True,next_offset=args.get('offset',0)+10),'complete')
    assert sum(a['operation']=='index' for a in calls)==2 and flow.report()['older_traces_unchecked']
    flow=SavedReviewFlow('Summarize saved reviews');args=flow.next_call()
    flow.accept(args,dict(reviews=[dict(review_id='a'*32)]),'complete')
    args=flow.next_call();flow.accept(args,{},'failed')
    assert flow.next_call() is None and flow.report()['failed_review_ids']==['a'*32]


def test_explicit_review_id_and_indexing_opt_out_are_respected():
    flow=SavedReviewFlow('Summarize saved review '+'a'*32)
    assert flow.next_call()==dict(operation='read',review_id='a'*32,offset=0)
    flow=SavedReviewFlow('Summarize saved literature reviews. Do not index old traces.')
    calls=[]
    while (args:=flow.next_call()) is not None:
        calls.append(args);flow.accept(args,dict(reviews=[]),'complete')
    assert all(a['operation']=='list' for a in calls)


@pytest.mark.parametrize('prompt',[
    'List saved literature reviews.', 'Run AF3 and summarize saved literature reviews.',
    'Summarize this newly imported paper.', 'Archive previous reviews and summarize the rest.',
    'Summarize saved literature reviews and find latest papers.',
])
def test_general_and_mixed_requests_keep_normal_route(prompt):assert not matches(prompt)


def test_explicit_readonly_requests_match_in_both_languages():
    assert matches('Summarize saved literature reviews. Do not run calculations or change project goals.')
    assert matches('\u603b\u7ed3\u5df2\u4fdd\u5b58\u7684\u6587\u732e\u5ba1\u9605\uff0c\u4e0d\u8981\u8fd0\u884c\u8ba1\u7b97\u3002')


def test_list_only_progress_is_not_mislabeled_as_reading_failure():
    ledger=progress([dict(step=dict(tool='literature_reviews',arguments=dict(operation='list')),status='complete',
        result=dict(collection='public',operation='list',reviews=[dict(review_id='a'*32)]))])
    text=incomplete_message(ledger,'Malformed JSON')
    assert '1 available review record(s)' in text and 'scientific synthesis has not started' in text
    assert 'Located 0 public document(s)' not in text
