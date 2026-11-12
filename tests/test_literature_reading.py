"""Explicit page completion and bounded same-conversation review continuity."""
import json
from types import SimpleNamespace

import pytest

from aidd_agent.chat_agent import ChatAgent
from aidd_agent import project_knowledge as knowledge
from aidd_agent.domain_agent import run
from aidd_agent.literature_reading import coverage,next_page_call,previous_summary
from aidd_agent.prompt_plan import LLMResponseError


@pytest.fixture
def paper(tmp_path,monkeypatch):
    import pypdf
    texts=['CIP2A abstract.']+['']*11+['Methods\n'+('CIP2A experimental construct details. '*220)]
    monkeypatch.setattr(pypdf,'PdfReader',lambda _:SimpleNamespace(is_encrypted=False,
        pages=[SimpleNamespace(extract_text=lambda t=t:t) for t in texts]))
    app=ChatAgent(tmp_path,start=False);sid=app.new_session()
    doc=knowledge.ingest(app,sid,'fixture.pdf',b'Synthetic PDF extraction fixture','public','https://example.org/paper',True)
    knowledge.ingest(app,sid,'private.txt',b'PRIVATE_REVIEW_SENTINEL')
    yield app,sid,doc['document_id']
    app.close()


def tool(**args):return dict(kind='tool',tool='project_knowledge',arguments=args,purpose='Read public literature'),{}
def answer():return dict(kind='answer',answer='Partial review grounded in retrieved passages. [E1]',evidence_ids=['E1']),{}


def test_search_enforces_document_and_page_filters(paper):
    app,sid,doc=paper
    other=knowledge.ingest(app,sid,'other.txt',b'CIP2A unrelated construct study','public','https://example.org/other',True)
    for operation in (knowledge.retrieve,knowledge.model_tool):
        result=operation(app,sid,dict(operation='search',query='CIP2A',document_id=doc,page=13,limit=20))
        assert result['results'] and all(r['id']==doc and r['page']==13 for r in result['results'])
        assert other['document_id'] not in json.dumps(result)
        assert not operation(app,sid,dict(operation='search',query='CIP2A',document_id=doc,page=12))['results']
        assert not operation(app,sid,dict(operation='search',query='CIP2A',document_id='f'*32))['results']


@pytest.mark.parametrize('extra',[dict(document_id='../other'),dict(page=0),dict(page=True)])
def test_invalid_search_scope_is_rejected(paper,extra):
    app,sid,_=paper
    with pytest.raises(ValueError):knowledge.model_tool(app,sid,dict(operation='search',query='CIP2A',**extra))


def test_selected_page_drains_without_extra_model_calls(paper):
    app,sid,doc=paper;calls=[]
    def planner(prompt,*a,**k):
        calls.append(json.loads(prompt))
        if len(calls)==1:return tool(operation='read',document_id=doc,page=13,limit=1)
        return answer()
    reply=run(app,sid,'Read page 13 of '+doc,'deepseek',planner=planner)
    trace=json.loads(next((app.root/'agents'/sid).glob('*.json')).read_text())
    reads=[r for r in trace['steps'] if r.get('step',{}).get('tool')=='project_knowledge']
    assert len(reads)>3 and len(calls)==2
    assert [r['step']['arguments'].get('offset',0) for r in reads]==list(range(len(reads)))
    assert all(r['model'].get('coordinator')=='explicit_page_pagination_no_model_call' for r in reads[1:])
    assert trace['indexed_reading_coverage']['pages'][0]['indexed_page_complete']
    assert 'complete page reads [13]' in reply and not app.jobs(sid)


def test_methods_search_only_answer_gets_one_reading_nudge(paper):
    app,sid,doc=paper;contexts=[]
    def planner(prompt,*a,**k):
        contexts.append(json.loads(prompt))
        if len(contexts)==1:return tool(operation='search',query='construct',document_id=doc,page=13,limit=1)
        if len(contexts)==3:
            assert 'only retrieved fragments' in contexts[-1]['previous_error']
            return tool(operation='read',document_id=doc,page=13,limit=1)
        return answer()
    reply=run(app,sid,'Read and review the Methods in '+doc,'deepseek',planner=planner)
    assert len(contexts)==4 and 'complete page reads [13]' in reply
    assert not app.jobs(sid)


def test_ignored_nudge_is_bounded_and_exposes_fragment_coverage(paper):
    app,sid,doc=paper;calls=[]
    def planner(*a,**k):
        calls.append(1)
        if len(calls)==1:return tool(operation='search',query='construct',document_id=doc,page=13,limit=1)
        return answer()
    reply=run(app,sid,'Review the Methods in '+doc,'deepseek',planner=planner)
    assert len(calls)==3 and 'fragment/partial pages [13]' in reply
    assert 'section/full-paper completion is not established' in reply


def record(doc,offset,count,more,operation='read',status='complete'):
    return dict(status=status,step=dict(tool='project_knowledge',arguments=dict(operation=operation,document_id=doc,page=13,offset=offset,limit=1)),
        result=dict(collection='public',has_more=more,next_offset=offset+count if more else None,
            results=[dict(id=doc,page=13,part=i,text='Fixture passage') for i in range(offset,offset+count)]))


def test_last_page_chunk_does_not_prove_complete_page_or_search_coverage():
    steps=[record('a'*32,0,1,True),record('a'*32,2,1,False)]
    assert not coverage(steps)['pages'][0]['indexed_page_complete']
    assert not coverage([record('a'*32,0,3,False,operation='search')])['pages'][0]['indexed_page_complete']
    steps.append(record('a'*32,1,1,True))
    assert coverage(steps)['pages'][0]['indexed_page_complete']


def test_failed_or_stalled_page_cursor_is_not_retried():
    first=record('a'*32,0,1,True)
    assert next_page_call([first])['offset']==1
    assert next_page_call([first,record('a'*32,1,1,True,status='failed')]) is None
    first['result']['next_offset']=0
    assert next_page_call([first]) is None
    assert next_page_call([record('a'*32,0,1,True,operation='search')]) is None


def test_page_budget_stops_with_partial_coverage(paper):
    app,sid,doc=paper;seen=[]
    def finalizer(context,profile):seen.append(context);return answer()
    reply=run(app,sid,'Read page 13','deepseek',max_steps=2,
        planner=lambda *a,**k:tool(operation='read',document_id=doc,page=13,limit=1),finalizer=finalizer)
    assert not seen[0]['indexed_reading_coverage']['pages'][0]['indexed_page_complete']
    assert seen[0]['public_read_progress']['continuation_cursors'][0]['next_offset']==1
    assert 'fragment/partial pages [13]' in reply


def test_previous_full_summary_reaches_planner_and_recovery(paper):
    app,sid,doc=paper
    previous=doc+'\n'+('Previous review wording. '*100)+'TAIL_CORRECTION_TARGET [E99]'
    app.message(sid,'assistant',previous+'\n\nEvidence sources:\nOld receipt\n\nAgent trace: old.json')
    calls=[];summaries=[]
    def planner(prompt,*a,**k):
        context=json.loads(prompt);calls.append(context)
        if len(calls)==1:return tool(operation='read',document_id=doc,page=1)
        assert context['previous_literature_summary']['text']==previous
        raise LLMResponseError('Malformed JSON')
    def finalizer(context,profile):
        summaries.append(context)
        assert context['previous_literature_summary']['text']==previous
        assert all(row['id']!='E99' for row in context['evidence'])
        assert 'PRIVATE_REVIEW_SENTINEL' not in json.dumps(context)
        return answer()
    run(app,sid,'Correct the prior summary for '+doc,'deepseek',planner=planner,finalizer=finalizer)
    assert len(summaries)==1


def test_previous_context_selects_same_document_and_never_silently_truncates():
    doc='a'*32;steps=[record(doc,0,1,False)]
    history=[dict(role='assistant',text=doc+' previous'),dict(role='assistant',text='b'*32+' unrelated')]
    assert previous_summary(history,steps)['text']==doc+' previous'
    assert previous_summary([dict(role='assistant',text=doc+'x'*19000)],steps)['status']=='omitted_over_budget'
    assert previous_summary([],steps)['status']=='not_found_in_recent_same_conversation'


def test_another_conversation_does_not_supply_prior_summary(paper):
    app,sid,doc=paper;other=app.new_session()
    app.message(other,'assistant',doc+' OTHER_CHAT_SUMMARY')
    contexts=[]
    def planner(prompt,*a,**k):
        contexts.append(json.loads(prompt))
        if len(contexts)==1:return tool(operation='read',document_id=doc,page=1)
        return answer()
    run(app,sid,'Review '+doc,'deepseek',planner=planner)
    assert contexts[-1]['previous_literature_summary']['status']=='not_found_in_recent_same_conversation'
    assert 'OTHER_CHAT_SUMMARY' not in json.dumps(contexts)


def test_native_recovery_wire_includes_complete_previous_summary(paper,monkeypatch):
    import io
    from aidd_agent import domain_agent,prompt_plan
    from aidd_agent.llm_profiles import select_llm_profile
    app,sid,doc=paper
    previous=doc+' '+('Earlier summary. '*150)+'END_OF_PREVIOUS_REVIEW'
    app.message(sid,'assistant',previous)
    calls=[];wire=[]
    def planner(*a,**k):
        calls.append(1)
        if len(calls)==1:return tool(operation='read',document_id=doc,page=1)
        raise LLMResponseError('Malformed JSON')
    class Transport:
        def open(self,request,timeout):
            payload=json.loads(request.data);wire.append(payload)
            context=json.loads(payload['messages'][1]['content'])
            assert context['previous_literature_summary']['text']==previous
            assert 'PRIVATE_REVIEW_SENTINEL' not in request.data.decode()
            assert 'tools' not in payload
            return io.BytesIO(json.dumps(dict(choices=[dict(finish_reason='stop',message=dict(content='Updated review from page 1. [E1]'))])).encode())
    monkeypatch.setattr(domain_agent,'chat_plan',planner)
    monkeypatch.setattr(prompt_plan,'build_opener',lambda *a:Transport())
    profile=json.loads(select_llm_profile('deepseek',config_dir=app.config_dir).read_text())
    if profile.get('api_key_env'):monkeypatch.setenv(profile['api_key_env'],'fixture-only')
    reply=run(app,sid,'Correct the prior summary for '+doc,'deepseek')
    assert 'Updated review' in reply and len(wire)==1 and not app.jobs(sid)
