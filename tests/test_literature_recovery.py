"""Public literature retains paging, evidence and safe final-answer recovery."""
import io
import json
from pathlib import Path

import pytest

from aidd_agent.chat_agent import ChatAgent
from aidd_agent import project_knowledge as knowledge
from aidd_agent.domain_agent import run,validate_step
from aidd_agent.prompt_plan import chat_plan,LLMResponseError,LLMRequestError
from aidd_agent.literature_recovery import synthesize,progress


@pytest.fixture
def library(tmp_path):
    app=ChatAgent(tmp_path,start=False);sid=app.new_session()
    text=('Published CIP2A binding evidence. '+('x'*200)+'\n')*170
    public=knowledge.ingest(app,sid,'paper.txt',text.encode(),'public','https://example.org/paper',True)
    private=knowledge.ingest(app,sid,'private.txt',b'PRIVATE_SENTINEL_873 unpublished data')
    yield app,sid,public,private
    app.close()


def test_model_paging_preserves_complete_passages_and_cursor(library):
    app,sid,public,_=library;offset=0;parts=[]
    while True:
        result=knowledge.model_tool(app,sid,dict(operation='read',document_id=public['document_id'],limit=20,offset=offset))
        assert len(json.dumps(result,ensure_ascii=False))<12000
        assert result['returned_count']==len(result['results'])<20
        for row in result['results']:
            assert 'text' in row and len(row['text'])<=1800
            parts.append(row['part'])
        if not result['has_more']:break
        assert result['next_offset']==offset+len(result['results'])
        offset=result['next_offset']
    assert parts==list(range(public['chunks']))


def failing_reader(public,seen=None,error=None):
    calls=[0]
    def planner(prompt,profile,**kwargs):
        if seen is not None:seen.append(json.loads(prompt))
        calls[0]+=1
        if calls[0]==1:return dict(kind='tool',tool='project_knowledge',arguments=dict(operation='read',document_id=public['document_id'],limit=20),purpose='Read public evidence'),{}
        raise error or LLMResponseError('Malformed final JSON',dict(content_characters=1000,finish_reason='stop'))
    return planner,calls


def test_malformed_json_recovers_without_tools_or_private_evidence(library):
    app,sid,public,private=library;planner,calls=failing_reader(public);captured=[]
    def finalizer(context,profile):
        captured.append(context)
        assert 'PRIVATE_SENTINEL_873' not in json.dumps(context) and private['document_id'] not in json.dumps(context)
        assert context['public_read_progress']['documents'][0]['passages_read']>0
        return dict(kind='answer',answer='Partial public evidence from page 1; remaining passages unread. [E1]',evidence_ids=['E1']),{}
    reply=run(app,sid,'Review imported public literature','deepseek',planner=planner,finalizer=finalizer)
    assert 'Partial public evidence' in reply and 'step/time limit' not in reply
    assert calls[0]==3 and len(captured)==1 and not app.jobs(sid)
    trace=json.loads(next((app.root/'agents'/sid).glob('*.json')).read_text())
    assert trace['answer_mode']=='read_only_text_recovery' and trace['status']=='complete'
    assert trace['steps'][1]['diagnostics']['finish_reason']=='stop'


def test_failed_synthesis_returns_progress_not_import_failure(library):
    app,sid,public,_=library;planner,_=failing_reader(public)
    def finalizer(*args):raise ValueError('Summary format still invalid')
    reply=run(app,sid,'Review literature','deepseek',planner=planner,finalizer=finalizer)
    assert 'Public knowledge access succeeded' in reply and 'scientific literature summary is incomplete' in reply
    assert public['document_id'] in reply and 'next_offset' in reply and 'do not re-import' in reply
    trace=json.loads(next((app.root/'agents'/sid).glob('*.json')).read_text())
    assert trace['status']=='incomplete' and not app.jobs(sid)


def test_provider_failure_keeps_public_progress_without_retry(library):
    app,sid,public,_=library;planner,calls=failing_reader(public,error=LLMRequestError('LLM HTTP 429'))
    reply=run(app,sid,'Review literature','deepseek',planner=planner,finalizer=lambda *a:pytest.fail('Do not retry unavailable provider'))
    assert calls[0]==2 and 'LLM HTTP 429' in reply and 'Public knowledge access succeeded' in reply


def test_finalizer_cannot_dispatch_or_cite_unknown_evidence(library):
    app,sid,public,_=library;planner,_=failing_reader(public)
    def finalizer(*args):return dict(kind='answer',answer='Invented',evidence_ids=['E99']),{}
    reply=run(app,sid,'Review literature','deepseek',planner=planner,finalizer=finalizer)
    assert 'Invalid final synthesis evidence contract' in reply and not app.jobs(sid)


def test_text_synthesis_uses_actual_wire_payload_without_json_format():
    context=dict(request='Review public paper',evidence=[dict(id='E1',tool='project_knowledge',status='complete',result=dict(collection='public',results=[dict(text='Public binding observation',page=1)]))])
    class Opener:
        def open(self,request,timeout):
            payload=json.loads(request.data)
            assert 'response_format' not in payload and 'tools' not in payload
            assert 'plain Markdown text' in payload['messages'][0]['content']
            text='The paper reports a "binding observation" on page 1. [E1]\nCoverage is partial.'
            return io.BytesIO(json.dumps(dict(choices=[dict(finish_reason='stop',message=dict(content=text))])).encode())
    result,_=synthesize(context,dict(base_url='http://localhost:1',model='fixture',json_mode=True),Opener())
    assert '"binding observation"' in result['answer'] and result['evidence_ids']==['E1']


def test_whole_json_fence_is_accepted_but_prose_fragments_are_not():
    class Opener:
        content='```json\n{"kind":"answer","answer":"No evidence read","evidence_ids":[]}\n```'
        def open(self,*a,**k):return io.BytesIO(json.dumps(dict(choices=[dict(finish_reason='stop',message=dict(content=self.content))])).encode())
    profile=dict(base_url='http://localhost:1',model='fixture')
    assert chat_plan('test',profile,Opener(),validator=validate_step)[0]['kind']=='answer'
    opener=Opener();opener.content='Some prose '+opener.content
    with pytest.raises(LLMResponseError) as exc:chat_plan('test',profile,opener,validator=validate_step)
    assert exc.value.diagnostics['content_characters']>0


def test_completed_read_removes_previous_continuation_cursor(library):
    app,sid,public,_=library;steps=[];offset=0
    while True:
        args=dict(operation='read',document_id=public['document_id'],limit=20,offset=offset)
        result=knowledge.model_tool(app,sid,args)
        steps.append(dict(step=dict(tool='project_knowledge',arguments=args),status='complete',result=result))
        if not result['has_more']:break
        offset=result['next_offset']
    ledger=progress(steps)
    assert ledger['continuation_cursors']==[] and ledger['documents'][0]['passages_read']==public['chunks']


def test_native_domain_path_recovers_through_text_transport(library,monkeypatch):
    from aidd_agent import domain_agent,prompt_plan
    app,sid,public,private=library;planner,calls=failing_reader(public)
    monkeypatch.setattr(domain_agent,'chat_plan',planner)
    captured=[]
    class Opener:
        def open(self,request,timeout):
            payload=json.loads(request.data);captured.append(payload)
            assert 'response_format' not in payload
            assert 'PRIVATE_SENTINEL_873' not in request.data.decode()
            return io.BytesIO(json.dumps(dict(choices=[dict(finish_reason='stop',message=dict(content='Public page 1 contains binding evidence. Coverage is partial. [E1]'))])).encode())
    monkeypatch.setattr(prompt_plan,'build_opener',lambda *a:Opener())
    from aidd_agent.llm_profiles import select_llm_profile
    profile=json.loads(select_llm_profile('deepseek',config_dir=app.config_dir).read_text())
    if profile.get('api_key_env'):monkeypatch.setenv(profile['api_key_env'],'fixture-secret-never-in-content')
    reply=run(app,sid,'Review the public paper','deepseek')
    assert 'Public page 1 contains' in reply and len(captured)==1 and calls[0]==3
    assert not app.jobs(sid)


def test_repeated_invalid_planning_stops_without_claiming_import_failure(library):
    app,sid,public,_=library;calls=[]
    def invalid(*a,**k):calls.append(1);raise ValueError('Malformed LLM JSON response')
    reply=run(app,sid,'Review public literature','deepseek',planner=invalid)
    assert len(calls)==2 and 'No tool ran in this turn' in reply
    assert knowledge.retrieve(app,sid,dict(operation='list',collection='public'))['results']


def test_step_budget_reserves_synthesis_without_replaying_read(library):
    app,sid,public,_=library;planner,calls=failing_reader(public);summaries=[]
    def finalizer(context,profile):
        summaries.append(context)
        return dict(kind='answer',answer='Partial page 1 evidence. [E1]',evidence_ids=['E1']),{}
    reply=run(app,sid,'Review public literature','deepseek',planner=planner,finalizer=finalizer,max_steps=2)
    assert 'Partial page 1 evidence' in reply and calls[0]==1 and len(summaries)==1
    assert not app.jobs(sid)


def test_listing_alone_never_triggers_scientific_synthesis(library):
    app,sid,public,_=library
    def planner(*a,**k):
        return dict(kind='tool',tool='project_knowledge',arguments=dict(operation='list'),purpose='Find imported documents'),{}
    reply=run(app,sid,'Review public literature','deepseek',planner=planner,max_steps=2,
        finalizer=lambda *a:pytest.fail('Listing is not paper evidence'))
    assert 'metadata only' in reply and 'scientific literature summary is incomplete' in reply
    assert public['document_id'] in reply
