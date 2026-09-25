import json
from pathlib import Path
import subprocess
import sys
import threading
import time
from urllib.request import Request, urlopen
from urllib.error import HTTPError

import pytest

from aidd_agent.chat_agent import ChatAgent, validate_route
from aidd_agent.chat_web import make_server
from aidd_agent.prompt_workflow import create_plan


def route(intent, request="", task_id=None, message=""):
    return dict(intent=intent,request=request,task_id=task_id,message=message)


def wait_for(agent,sid,status):
    deadline=time.monotonic()+15
    while time.monotonic()<deadline:
        task=agent.task(sid)
        if task["status"]==status:return task
        time.sleep(.1)
    raise AssertionError(agent.snapshot(sid))


def test_dialogue_clarification_routing_and_session_ownership(tmp_path):
    seen=[]
    def router(prompt, profile):
        seen.append(json.loads(prompt))
        return route("clarify",message="Which organism?") if len(seen)==1 else route("run",request="Prepare human WEE1 input only.")
    app=ChatAgent(tmp_path,router=router,start=False)
    sid=app.new_session()
    app.ask(sid,"Prepare WEE1 input","deepseek")
    assert not app.jobs(sid)
    app.ask(sid,"Human","deepseek")
    assert any(m["text"]=="Which organism?" for m in seen[1]["conversation"])
    jid=app.task(sid)["id"]
    assert app.task(sid)["request"]=="Prepare human WEE1 input only."
    app.ask(sid,"/status","gpt")
    assert len(app.jobs(sid))==1
    other=app.new_session()
    with pytest.raises(ValueError,match="No matching"):
        app.ask(other,"/cancel "+jid,"deepseek")
    assert app.task(sid)["cancel"]==0
    app.close()
    reopened=ChatAgent(tmp_path,start=False)
    assert reopened.snapshot(sid)["messages"]
    assert reopened.task(sid)["id"]==jid
    reopened.update(jid,status="running")
    reopened.close()
    recovered=ChatAgent(tmp_path,start=False)
    assert recovered.task(sid)["status"]=="interrupted"
    recovered.close()


def test_results_are_read_from_actual_report_not_model(tmp_path):
    app=ChatAgent(tmp_path,router=lambda *a:route("run",request="Run WEE1 search"),start=False)
    sid=app.new_session();app.ask(sid,"Run WEE1 search","deepseek")
    job=app.task(sid)
    report=tmp_path/"report.json"
    report.write_text(json.dumps(dict(status="complete",steps={"search":dict(status="complete",result=dict(report="/actual/search/report.json",biological_quality="not_evaluated"))})))
    app.update(job["id"],status="complete",report=str(report))
    answer=app.ask(sid,"/results","deepseek")
    assert "/actual/search/report.json" in answer and "not_evaluated" in answer
    assert len(app.jobs(sid))==1
    plan=tmp_path/"plan.json";plan.write_text('{}')
    child=tmp_path/"search-report.json"
    child.write_text(json.dumps(dict(queries=[dict(query_id="fixture-query",execution_seconds=12.5,gaussian_seconds=10,annotation_seconds=1)])))
    report.write_text(json.dumps(dict(steps={"search":dict(status="complete",result=dict(report=str(child)))})))
    app.update(job["id"],plan=str(plan))
    answer=app.ask(sid,"/results","deepseek")
    assert "fixture-query: execution 12.5 s" in answer
    app.close()


def test_real_worker_handles_clarification_report_without_network(tmp_path):
    fixture=tmp_path/"clarification.json"
    fixture.write_text(json.dumps(dict(version=1,summary="Missing organism",clarifications=["Which organism?"],steps=[])))
    def planner(db,user,project,prompt,**kwargs):
        return create_plan(db,user,project,prompt,response_file=fixture)
    app=ChatAgent(tmp_path,planner=planner,router=lambda *a:route("run",request="Prepare protein input"))
    try:
        sid=app.new_session();app.ask(sid,"Prepare protein input","deepseek")
        job=wait_for(app,sid,"blocked")
        assert Path(job["report"]).is_file()
        assert any(m["text"]=="Which organism?" for m in app.snapshot(sid)["messages"])
    finally: app.close()


def test_background_cancel_preserves_plan_and_does_not_block_status(tmp_path,monkeypatch):
    fixture=tmp_path/"clarification.json"
    fixture.write_text(json.dumps(dict(version=1,summary="Fixture",clarifications=["Fixture?"],steps=[])))
    def planner(db,user,project,prompt,**kwargs):
        return create_plan(db,user,project,prompt,response_file=fixture)
    real_popen=subprocess.Popen
    def sleeping_child(cmd,**kwargs):
        if "aidd_agent.prompt_workflow" in cmd:
            return real_popen([sys.executable,"-u","-c","import time; print('fixture running',flush=True); time.sleep(30)"],**kwargs)
        return real_popen(cmd,**kwargs)
    monkeypatch.setattr(subprocess,"Popen",sleeping_child)
    app=ChatAgent(tmp_path,planner=planner,router=lambda *a:route("run",request="Fixture slow task"))
    try:
        sid=app.new_session();app.ask(sid,"Run fixture","deepseek")
        job=wait_for(app,sid,"running")
        assert "running" in app.ask(sid,"/status","deepseek")
        app.ask(sid,"/cancel","deepseek")
        cancelled=wait_for(app,sid,"cancelled")
        assert Path(cancelled["plan"]).is_file()
        assert app.process is None
    finally:app.close()


def test_http_token_origin_and_assets(tmp_path):
    app=ChatAgent(tmp_path,start=False)
    server=make_server(app,0,token="fixture-token")
    thread=threading.Thread(target=server.serve_forever,daemon=True);thread.start()
    base=f"http://127.0.0.1:{server.server_port}"
    try:
        with urlopen(base+"/") as r:
            assert b"AIDD Workbench" in r.read()
            assert "frame-ancestors 'none'" in r.headers["Content-Security-Policy"]
        with pytest.raises(HTTPError) as err:urlopen(base+"/api/state")
        assert err.value.code==401
        def post(path,body,origin=None):
            headers={"Authorization":"Bearer fixture-token","Content-Type":"application/json"}
            if origin:headers["Origin"]=origin
            with urlopen(Request(base+path,data=json.dumps(body).encode(),headers=headers)) as r:return json.load(r)
        with pytest.raises(HTTPError) as err:post('/api/session',{},'https://other.example')
        assert err.value.code==403
        sid=post('/api/session',{})['session']
        result=post('/api/message',dict(session=sid,text='/capabilities',provider='deepseek'))
        assert "WEE1" in result['message']
        assert "Glide preparation/execution adapter; workstation validation pending" in result['message']
        assert "PLANTS execution and cross-docking validation remain pending" in result['message']
        with urlopen(Request(base+'/api/state?session='+sid,headers={'Authorization':'Bearer fixture-token'})) as r:
            state=json.load(r)
            assert len(state['messages'])==2 and not state['tasks']
        viewer=app.root/'viewers'/sid;viewer.mkdir(parents=True)
        identifier='123-abcdefabcdef';artifact=viewer/'contacts.csv'
        artifact.write_text('distance\n3.2\n')
        receipt=viewer/(identifier+'.result.json')
        receipt.write_text(json.dumps(dict(artifacts=dict(contacts_csv=str(artifact)))))
        endpoint=base+'/api/viewer/artifact?session='+sid+'&id='+identifier+'&artifact=contacts_csv'
        with pytest.raises(HTTPError) as err:urlopen(endpoint)
        assert err.value.code==401
        with urlopen(Request(endpoint,headers={'Authorization':'Bearer fixture-token'})) as r:
            assert b'3.2' in r.read()
        receipt.write_text(json.dumps({'artifacts': {'contacts_csv': str(tmp_path/'outside.csv')}}))
        with pytest.raises(HTTPError) as err:urlopen(Request(endpoint,headers={'Authorization':'Bearer fixture-token'}))
        assert err.value.code==400
    finally:
        server.shutdown();server.server_close();app.close()


@pytest.mark.parametrize('bad',[
    dict(intent='shell',message='',request='',task_id=None),
    dict(intent='status',message='',request='Run AF3',task_id=None),
    dict(intent='run',message='',request='',task_id=None),
    dict(intent='clarify',message=chr(0x4e2d),request='',task_id=None)])
def test_invalid_router_output(bad):
    with pytest.raises(ValueError):validate_route(bad)


def test_attach_existing_owned_result_without_rerunning(tmp_path):
    app=ChatAgent(tmp_path,start=False)
    fixture=tmp_path/"clarification.json"
    fixture.write_text(json.dumps(dict(version=1,summary="Existing fixture",clarifications=["Fixture?"],steps=[])))
    ctx=app.context
    plan=create_plan(Path(ctx["db"]),ctx["user_id"],ctx["project_id"],"fixture",response_file=fixture)
    sid=app.new_session()
    jid=app.attach(sid,str(plan))
    assert app.task(sid,jid)["status"]=="interrupted"
    assert str(plan) in app.ask(sid,"/results","gpt")
    with pytest.raises(ValueError):app.attach(sid,str(fixture))
    plan.write_text('{}')
    with pytest.raises(ValueError):app.attach(sid,str(plan))
    app.close()
