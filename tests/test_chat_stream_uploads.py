import io
import json
from pathlib import Path
import threading
from urllib.error import HTTPError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

import pytest

from aidd_agent.chat_agent import ChatAgent
from aidd_agent.chat_stream import answer_prefix, read_sse
from aidd_agent.chat_uploads import upload, attachment, attach_to_prompt
from aidd_agent.chat_web import make_server
from aidd_agent.prompt_plan import chat_plan


def sse(parts, finish='stop'):
    rows = [dict(id='fixture', choices=[dict(index=0, delta=dict(content=part))]) for part in parts]
    rows.append(dict(choices=[dict(index=0, delta={}, finish_reason=finish)]))
    return b''.join(('data: '+json.dumps(row)+'\n\n').encode() for row in rows)+b'data: [DONE]\n\n'


def test_partial_answer_does_not_expose_tools_or_nested_fields():
    value = dict(kind='answer', answer='A "quote", a newline\nand '+chr(0x1f9ec), evidence_ids=[])
    text = json.dumps(value)
    prefixes = [answer_prefix(text[:i]) for i in range(len(text)+1)]
    assert all(value['answer'].startswith(prefix) for prefix in prefixes)
    assert prefixes[-1] == value['answer']
    assert answer_prefix('{"kind":"tool","arguments":{"answer":"secret"}}') == ''
    assert answer_prefix('{"answer":"hello","kind":"answer"}') == 'hello'
    assert answer_prefix('{"answer":"hello') == ''


def test_provider_streaming_preserves_actual_deltas_and_validates_final_json():
    plan = dict(kind='answer', answer='Validated answer', evidence_ids=[])
    text = json.dumps(plan)
    fragments = [text[:19], text[19:34], text[34:]]
    seen = []
    class Opener:
        def open(self, request, timeout):
            assert json.loads(request.data)['stream'] is True
            return io.BytesIO(sse(fragments))
    result, meta = chat_plan('Inspect results', dict(base_url='http://localhost:1234/v1', model='fixture'),
                            Opener(), validator=lambda value: value, on_content=seen.append)
    assert result == plan and seen == fragments and meta['response_id'] == 'fixture'
    class Incomplete:
        def open(self, *a, **k):return io.BytesIO(sse([text], finish='length'))
    with pytest.raises(ValueError, match='incomplete'):
        chat_plan('Inspect', dict(base_url='http://localhost:1', model='fixture'), Incomplete(),
                  validator=lambda value: value, on_content=lambda _: None)


def test_stream_disabled_provider_uses_normal_response():
    class Opener:
        def open(self, request, timeout):
            assert 'stream' not in json.loads(request.data)
            return io.BytesIO(json.dumps(dict(choices=[dict(finish_reason='stop', message=dict(content='{"ok":true}'))])).encode())
    value, _ = chat_plan('Inspect', dict(base_url='http://localhost:1', model='fixture', stream=False),
                         Opener(), validator=lambda value:value, on_content=lambda _:pytest.fail('Not a stream'))
    assert value == {'ok':True}


def test_sse_provider_error_is_not_success():
    with pytest.raises(ValueError, match='streaming error'):
        read_sse(io.BytesIO(b'data: {"error":{"message":"private"}}\n\n'), lambda _:None)


def test_upload_bytes_ownership_and_tamper_checks(tmp_path):
    app = ChatAgent(tmp_path, start=False)
    from test_research_projects import create
    sid=app.new_session();other=create(app,sid,'Different target')['session']
    payload = b'@<TRIPOS>MOLECULE\r\nligand\r\n'
    row = upload(app, sid, 'ligand.mol2', io.BytesIO(payload), len(payload))
    assert Path(row['path']).read_bytes() == payload
    assert attachment(app,app.new_session(),row['id'])['path']==row['path']
    prompt = attach_to_prompt(app, sid, 'Inspect this ligand', [row['id']])
    assert json.loads(prompt.split('\n')[-1])[0]['path'] == row['path'] and '@<TRIPOS>' not in prompt
    with pytest.raises(ValueError, match='unavailable'):attachment(app, other, row['id'])
    Path(row['path']).write_bytes(b'changed')
    with pytest.raises(ValueError, match='changed'):attachment(app, sid, row['id'])
    app.close()


@pytest.mark.parametrize('name', ['../a.pdb', 'a\\b.pdb', 'a.py', 'a.zip', 'x\n.pdb'])
def test_bad_upload_names(tmp_path, name):
    app = ChatAgent(tmp_path, start=False)
    with pytest.raises(ValueError):upload(app, app.new_session(), name, io.BytesIO(b'x'), 1)
    app.close()


def test_truncated_and_oversized_uploads_are_not_registered(tmp_path):
    app = ChatAgent(tmp_path, start=False);sid = app.new_session()
    with pytest.raises(ValueError, match='interrupted'):upload(app,sid,'a.pdb',io.BytesIO(b'a'),3)
    with pytest.raises(ValueError, match='50 MiB'):upload(app,sid,'a.pdb',io.BytesIO(b'a'),51*1024*1024)
    assert not list(tmp_path.rglob('attachment.json'))
    assert not list(tmp_path.rglob('upload.part'))
    app.close()


def test_attachment_requests_validate_ids_and_preserve_command_syntax(tmp_path):
    app = ChatAgent(tmp_path, start=False); sid = app.new_session()
    with pytest.raises(ValueError, match='strings'):
        attach_to_prompt(app, sid, 'Inspect', [{}])
    with pytest.raises(ValueError, match='natural-language'):
        app.ask(sid, '/status', 'gpt', attachments=['a'*32])
    assert app.snapshot(sid)['messages'] == []
    app.close()


def test_basic_chat_snapshot_does_not_import_scientific_runtime(tmp_path):
    import os
    import subprocess
    import sys
    script = ('import sys; from aidd_agent.chat_agent import ChatAgent; '
              'a=ChatAgent(sys.argv[1],start=False); a.snapshot(a.new_session()); '
              'assert not any(n=="rdkit" or n.startswith("rdkit.") for n in sys.modules); a.close()')
    subprocess.run([sys.executable, '-c', script, str(tmp_path)], check=True,
                   env={**os.environ, 'PYTHONPATH': str(Path(__file__).resolve().parents[1]/'src')}, timeout=30)


def test_http_stream_arrives_before_completion_and_upload_is_authenticated(tmp_path, monkeypatch):
    app = ChatAgent(tmp_path,start=False);sid = app.new_session()
    release = threading.Event()
    def ask(session, text, provider, *, on_event=None, attachments=None):
        assert session == sid
        on_event(dict(type='answer_delta',text='First '))
        assert release.wait(5)
        on_event(dict(type='answer_delta',text='second'))
        return 'First second'
    monkeypatch.setattr(app,'ask',ask)
    server = make_server(app,0,token='fixture-token')
    thread = threading.Thread(target=server.serve_forever,daemon=True);thread.start()
    base = 'http://127.0.0.1:'+str(server.server_port)
    headers = {'Authorization':'Bearer fixture-token','Content-Type':'application/json'}
    try:
        request = Request(base+'/api/message/stream',data=json.dumps(dict(session=sid,text='Inspect')).encode(),headers=headers)
        with urlopen(request,timeout=10) as response:
            assert json.loads(response.readline())['type']=='status'
            assert json.loads(response.readline())==dict(type='answer_delta',text='First ')
            release.set()
            events=[json.loads(line) for line in response]
            assert events[-1]==dict(type='done',text='First second')
        url=base+'/api/uploads?'+urlencode(dict(session=sid,name='ligands.csv'))
        with pytest.raises(HTTPError) as error:urlopen(Request(url,data=b'id,smiles\na,CC\n'))
        assert error.value.code==401
        with urlopen(Request(url,data=b'id,smiles\na,CC\n',headers=headers)) as response:
            assert response.status==201
            row=json.load(response)
        assert attachment(app,sid,row['id'])['bytes']==15
        with pytest.raises(HTTPError) as error:
            urlopen(Request(url,data=b'a',headers={**headers,'Origin':'https://foreign.test'}))
        assert error.value.code==403
    finally:
        release.set();server.shutdown();server.server_close();app.close()


def test_invalid_streamed_answer_is_reset_before_retry(tmp_path, monkeypatch):
    from aidd_agent import domain_agent
    calls=[];events=[]
    def planner(prompt, profile, **kwargs):
        calls.append(1)
        if len(calls)==1:
            kwargs['on_content']('{"kind":"answer","answer":"Unvalidated claim","evidence_ids":["absent"]}')
            return dict(kind='answer',answer='Unvalidated claim',evidence_ids=['absent']),{}
        kwargs['on_content']('{"kind":"answer","answer":"No evidence yet","evidence_ids":[]}')
        return dict(kind='answer',answer='No evidence yet',evidence_ids=[]),{}
    monkeypatch.setattr(domain_agent,'chat_plan',planner)
    app=ChatAgent(tmp_path,start=False);sid=app.new_session()
    answer=app.ask(sid,'Inspect results','gpt',on_event=events.append)
    assert answer.startswith('No evidence yet')
    assert any(row['type']=='answer_reset' for row in events)
    saved = app.snapshot(sid)['messages'][-1]['text']
    assert saved.startswith('No evidence yet')
    assert 'Unvalidated claim' not in saved
    app.close()


def test_http_upload_prompt_tool_stream_and_reload(tmp_path, monkeypatch):
    """Exercise the same HTTP path as the composer, with a real read-only tool."""
    from aidd_agent import domain_agent
    contexts = []
    def planner(prompt, profile, **kwargs):
        context = json.loads(prompt); contexts.append(context)
        assert 'Attached local file references' in context['request']
        assert 'id,smiles' not in context['request']
        if not context['evidence']:
            return dict(kind='tool', tool='tasks', arguments={}, purpose='Read local tasks'), {}
        assert context['evidence'][0]['tool'] == 'tasks'
        assert context['evidence'][0]['status'] == 'complete'
        value = dict(kind='answer', answer='**Local tasks checked.**', evidence_ids=['E1'])
        raw = json.dumps(value)
        for part in [raw[:35], raw[35:]]: kwargs['on_content'](part)
        return value, {}
    monkeypatch.setattr(domain_agent, 'chat_plan', planner)
    app = ChatAgent(tmp_path, start=False)
    server = make_server(app, 0, token='fixture-token')
    threading.Thread(target=server.serve_forever, daemon=True).start()
    base = 'http://127.0.0.1:'+str(server.server_port)
    headers = {'Authorization': 'Bearer fixture-token', 'Content-Type': 'application/json'}
    def request(path, body):
        with urlopen(Request(base+path, data=json.dumps(body).encode(), headers=headers), timeout=10) as response:
            return json.load(response)
    try:
        for asset in ['/chat.js', '/markdown.js', '/chat.css']:
            with urlopen(base+asset) as response: assert response.status == 200 and response.read()
        sid = request('/api/session', {})['session']
        payload = b'id,smiles\nfixture,CCO\n'
        url = base+'/api/uploads?'+urlencode(dict(session=sid, name='test.csv'))
        with urlopen(Request(url, data=payload, headers=headers)) as response: row = json.load(response)
        assert Path(row['path']).read_bytes() == payload
        data = dict(session=sid, text='Inspect my local tasks and file reference', provider='gpt', attachments=[row['id']])
        with urlopen(Request(base+'/api/message/stream', data=json.dumps(data).encode(), headers=headers), timeout=10) as response:
            events = [json.loads(line) for line in response]
        assert len(contexts) == 2
        assert any(e['type']=='tool_start' and e['tool']=='tasks' for e in events)
        assert any(e['type']=='tool_end' and e['status']=='complete' for e in events)
        assert ''.join(e['text'] for e in events if e['type']=='answer_delta') == '**Local tasks checked.**'
        assert events[-1]['type']=='done'
        with urlopen(Request(base+'/api/state?session='+sid, headers=headers)) as response: state = json.load(response)
        assert state['messages'][-1]['text'] == events[-1]['text']
        assert state['tasks'] == []
        assert state['agent_runs'][0]['steps'][0]['tool'] == 'tasks'
        with urlopen(Request(base+'/api/agent/trace?'+urlencode(dict(session=sid,id=state['agent_runs'][0]['id'])), headers=headers)) as response:
            assert json.load(response)['status']=='complete'
    finally:
        server.shutdown(); server.server_close(); app.close()
