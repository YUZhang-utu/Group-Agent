"""Paper-level display and source-checked project review reuse."""
import json
from pathlib import Path

import pytest

from aidd_agent.chat_agent import ChatAgent
from aidd_agent.domain_agent import run
from aidd_agent import project_knowledge as knowledge,literature_reviews as reviews
from aidd_agent.literature_citations import present,references,label
from aidd_agent.literature_reading import previous_summary
from aidd_agent.prompt_plan import LLMResponseError


@pytest.fixture
def project(tmp_path):
    app=ChatAgent(tmp_path,start=False);sid=app.new_session()
    doc=knowledge.ingest(app,sid,'target_paper.pdf.txt',b'CIP2A construct assay public result.','public','https://example.org/paper',True)
    private=knowledge.ingest(app,sid,'private.txt',b'PRIVATE_LIBRARY_MARKER')
    yield app,sid,doc['document_id'],private['document_id']
    app.close()


def make_review(app,sid,doc,text='CIP2A assay interpretation; unresolved mechanism. [E1]'):
    calls=[]
    def planner(*a,**k):
        calls.append(1)
        if len(calls)==1:return dict(kind='tool',tool='project_knowledge',arguments=dict(operation='read',document_id=doc,page=1),purpose='Read paper'),{}
        return dict(kind='answer',answer=text,evidence_ids=['E1']),{}
    reply=run(app,sid,'Review public paper '+doc,'deepseek',planner=planner)
    path=max((app.root/'agents'/sid).glob('*.json'),key=lambda p:p.stat().st_mtime)
    return reply,json.loads(path.read_text())


def test_display_metadata_and_audit_keep_separate_roles(project):
    app,sid,doc,_=project
    knowledge.control(app,sid,dict(operation='bibliography',collection='public',document_id=doc,
        display_title='CIP2A scaffolds',authors='Example et al.',year='2025'))
    reply,audit=make_review(app,sid,doc)
    assert '[Example et al., 2025 - CIP2A scaffolds](https://example.org/paper)' in reply
    assert '[E1]' not in reply and doc not in reply and 'Agent trace:' not in reply
    assert '```sources' in reply and 'complete page reads [1]' in reply
    assert audit['answer'].endswith('[E1]') and audit['steps'][0]['evidence_id']=='E1'
    assert audit['saved_literature_review']['status']=='saved'
    assert label(dict(title='file.pdf'))=='file'


def test_one_paper_repeated_receipts_are_deduplicated(project):
    app,sid,doc,_=project;_,audit=make_review(app,sid,doc)
    first=audit['steps'][0];second=dict(first,evidence_id='E2')
    reply,refs,_=present('Result [E1][E2].',[first,second],audit['indexed_reading_coverage'])
    assert len(refs)==1 and reply.split('Sources:',1)[0].strip()=='Result.'
    assert reply.split('```sources',1)[0].count('https://example.org/paper')==1


def test_multiple_papers_with_same_title_remain_distinct(project):
    app,sid,doc,_=project;_,audit=make_review(app,sid,doc)
    other=knowledge.ingest(app,sid,'target_paper.pdf.txt',b'Other conflicting public observation.','public','https://example.org/other',True)
    result=knowledge.model_tool(app,sid,dict(operation='read',document_id=other['document_id'],page=1))
    steps=[audit['steps'][0],dict(step=dict(tool='project_knowledge'),status='complete',result=result,evidence_id='E2')]
    display,refs,_=present('Compare [E1, E2].',steps,{})
    assert len(refs)==2 and label(refs[0])!=label(refs[1])
    assert 'https://example.org/other' in display and '[E1' not in display


def test_saved_review_survives_restart_and_is_used_from_new_chat(project):
    app,sid,doc,_=project;_,audit=make_review(app,sid,doc);identifier=audit['saved_literature_review']['review_id']
    app.close();reopened=ChatAgent(app.root.parent,start=False)
    try:
        other=reopened.new_session();contexts=[]
        def planner(prompt,*a,**k):
            context=json.loads(prompt);contexts.append(context)
            if len(contexts)==1:
                assert context['project_memory']['literature_reviews']['reviews'][0]['review_id']==identifier
                return dict(kind='tool',tool='literature_reviews',arguments=dict(operation='read',review_id=identifier),purpose='Reuse saved review'),{}
            assert context['evidence'][0]['result']['verification'].startswith('Model-generated')
            return dict(kind='answer',answer='The saved review records an unresolved mechanism. [E1]',evidence_ids=['E1']),{}
        reply=run(reopened,other,'What did our CIP2A paper review conclude?','deepseek',planner=planner)
        assert 'original paper text was not reread this turn' in reply and '[E1]' not in reply
        assert len(reviews.handle(reopened,other,dict(operation='list'))['reviews'])==1
        assert not reopened.jobs(other)
    finally:reopened.close()


@pytest.mark.parametrize('change',['archive','internal','source','review'])
def test_unavailable_or_modified_sources_never_surface_saved_review(project,change):
    app,sid,doc,_=project;_,audit=make_review(app,sid,doc);identifier=audit['id']
    directory=knowledge.folder(app,sid,'public')
    if change=='archive':knowledge.control(app,sid,dict(operation='archive',collection='public',document_id=doc))
    elif change=='internal':knowledge.ingest(app,sid,'group.txt',b'CIP2A construct assay public result.')
    elif change=='source':next((directory/'documents'/doc).glob('source.*')).write_bytes(b'Modified content')
    else:
        with knowledge.connect(directory) as db:db.execute('UPDATE literature_reviews SET payload=?',('{"private":"MUTATED"}',))
    assert reviews.handle(app,sid,dict(operation='list'))['reviews']==[]
    with pytest.raises(ValueError,match='unavailable'):reviews.handle(app,sid,dict(operation='read',review_id=identifier))


def test_review_isolation_and_private_metadata(project):
    from aidd_agent.research_projects import handle
    app,sid,doc,private=project;_,audit=make_review(app,sid,doc)
    knowledge.control(app,sid,dict(operation='bibliography',document_id=private,display_title='PRIVATE_TITLE'))
    payload=reviews.handle(app,sid,dict(operation='read',review_id=audit['id']))
    assert 'PRIVATE' not in json.dumps(payload) and private not in json.dumps(payload)
    other=handle(app,sid,dict(operation='create',name='Another target',target='Other',objective='Independent work'))['session']
    assert reviews.handle(app,other,dict(operation='list'))['reviews']==[]
    with pytest.raises(ValueError):reviews.handle(app,other,dict(operation='read',review_id=audit['id']))


def test_model_cannot_edit_citation_metadata(project):
    app,sid,doc,_=project
    with pytest.raises(ValueError):knowledge.model_tool(app,sid,dict(operation='bibliography',document_id=doc,display_title='invented'))


def test_long_review_pages_exactly_and_stays_inside_evidence_budget(project):
    app,sid,doc,_=project;text=('Scientific interpretation with limitations. '*220)+' [E1]'
    _,audit=make_review(app,sid,doc,text);offset=0;parts=[]
    while True:
        result=reviews.handle(app,sid,dict(operation='read',review_id=audit['id'],offset=offset));parts.append(result['summary'])
        assert len(json.dumps(result,ensure_ascii=False))<12000
        if not result['has_more']:break
        assert result['next_offset']==offset+len(result['summary']);offset=result['next_offset']
    with knowledge.connect(knowledge.folder(app,sid,'public')) as db:
        stored=json.loads(db.execute('SELECT payload FROM literature_reviews WHERE id=?',(audit['id'],)).fetchone()['payload'])
    assert ''.join(parts)==stored['summary']


def test_saved_review_supports_existing_synthesis_recovery(project):
    app,sid,doc,_=project;_,audit=make_review(app,sid,doc);other=app.new_session();calls=[];seen=[]
    def planner(*a,**k):
        calls.append(1)
        if len(calls)==1:return dict(kind='tool',tool='literature_reviews',arguments=dict(operation='read',review_id=audit['id']),purpose='Read saved review'),{}
        raise LLMResponseError('Malformed JSON')
    def finalizer(context,profile):
        seen.append(context)
        assert context['public_read_progress']['saved_review_passages_read']==1
        assert 'PRIVATE_LIBRARY_MARKER' not in json.dumps(context)
        citation=next(r['id'] for r in context['evidence'] if r['result'].get('summary'))
        return dict(kind='answer',answer='Saved interpretation, not fresh reading. ['+citation+']',evidence_ids=[citation]),{}
    reply=run(app,other,'Recall earlier review','deepseek',planner=planner,finalizer=finalizer)
    assert len(seen)==1 and 'Saved interpretation' in reply and '[E1]' not in reply


def test_index_old_reviews_is_project_scoped_and_idempotent(project):
    app,sid,doc,_=project;_,audit=make_review(app,sid,doc)
    with knowledge.connect(knowledge.folder(app,sid,'public')) as db:db.execute('DELETE FROM literature_reviews')
    other=app.new_session()
    assert reviews.handle(app,other,dict(operation='index',limit=10))['indexed_reviews']==1
    assert reviews.handle(app,other,dict(operation='index',limit=10))['indexed_reviews']==0
    assert reviews.handle(app,other,dict(operation='list'))['reviews'][0]['review_id']==audit['id']


def test_friendly_paper_links_preserve_same_chat_comparison(project):
    app,sid,doc,_=project;reply,audit=make_review(app,sid,doc)
    previous=previous_summary([dict(role='assistant',text=reply)],audit['steps'])
    assert previous['status']=='available' and 'CIP2A assay interpretation' in previous['text']
    assert '```sources' not in previous['text']


def test_error_cursor_json_is_not_mistaken_for_citation(project):
    app,sid,doc,_=project;_,audit=make_review(app,sid,doc)
    cursor=json.dumps([dict(evidence_id='E1',document_id=doc,next_offset=3)])
    shown,_,_=present(cursor,audit['steps'],{},replace_document_ids=False)
    assert cursor in shown


def test_saved_review_works_with_native_text_synthesis_contract(project):
    import io
    from aidd_agent.literature_recovery import synthesize
    app,sid,doc,private=project;_,audit=make_review(app,sid,doc)
    result=reviews.handle(app,sid,dict(operation='read',review_id=audit['id']))
    context=dict(request='Recall prior review',evidence=[dict(id='E1',tool='literature_reviews',status='complete',result=result)])
    class Transport:
        def open(self,request,timeout):
            payload=json.loads(request.data)
            assert 'PRIVATE_LIBRARY_MARKER' not in request.data.decode() and private not in request.data.decode()
            assert 'tools' not in payload and 'response_format' not in payload
            assert 'saved_model_interpretation_not_fresh_paper_reading' in request.data.decode()
            return io.BytesIO(json.dumps(dict(choices=[dict(finish_reason='stop',message=dict(content='The saved interpretation remains unverified. [E1]'))])).encode())
    step,_=synthesize(context,dict(base_url='http://localhost:1',model='fixture'),Transport())
    assert step['evidence_ids']==['E1']


@pytest.mark.parametrize('case',['uncited','mixed_tools'])
def test_ineligible_answers_are_not_promoted_to_public_review_memory(project,case):
    app,sid,doc,_=project;reply,audit=make_review(app,sid,doc)
    with knowledge.connect(knowledge.folder(app,sid,'public')) as db:db.execute('DELETE FROM literature_reviews')
    if case=='uncited':audit['evidence_ids']=[]
    else:audit['steps'].append(dict(step=dict(kind='tool',tool='read_artifact',arguments={}),status='complete'))
    assert reviews.save(app,sid,audit,reply,audit['paper_references']) is None
    assert reviews.handle(app,sid,dict(operation='list'))['reviews']==[]
