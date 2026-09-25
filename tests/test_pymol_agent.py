import copy
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from aidd_agent.pymol_program import compile_program, run_program
from aidd_agent.pymol_agent import run_agent, skill_context, validate_plan
from aidd_agent.chat_agent import validate_route


@pytest.mark.parametrize('code',[
    'import os', 'cmd.do("quit")', 'cmd.load("x")', 'cmd.quit()',
    'cmd.color(__import__("os"), "all")', 'x = 1',
    'for x in []:\n cmd.hide("everything")', 'cmd.color("red", ") or all (")',
    'cmd.delete("v001")', 'cmd.set("surface_quality",10000)',
    'cmd.distance("ai_d","all","all",mode=3)',
])
def test_reject_unbounded_program(code):
    with pytest.raises((ValueError,SyntaxError)):compile_program(code)


class Cmd:
    def __init__(self):self.state={'color':'green'};self.calls=[];self.fail=False
    def get_names(self,kind):return ['v001','manual'] if kind=='objects' else []
    def get_model(self,selection,state=1):
        return SimpleNamespace(atom=[SimpleNamespace(chain='A',resi='23',resn='TRP',name='CA')])
    def count_atoms(self,selection):return 0 if 'none' in selection else 3
    def get_session(self):return copy.deepcopy(self.state)
    def set_session(self,state):self.state=copy.deepcopy(state)
    def save(self,path):Path(path).write_text('session')
    def png(self,path,**kwargs):Path(path).write_bytes(b'png')
    def color(self,**kwargs):
        self.calls.append(kwargs);self.state['color']=kwargs['color']
        if self.fail:raise ValueError('bad color')


CATALOG=[dict(id='v001',label='complex',kind='pdb')]


def test_scoping_artifacts_rollback_and_undo(tmp_path):
    cmd=Cmd()
    message=dict(id='123-abcdef012345',view=dict(operation='program',code='cmd.color("red", "all")'))
    result=run_program(cmd,message,CATALOG,tmp_path)
    assert result['status']=='complete'
    assert cmd.calls[0]['selection']=='(v001) and (all)'
    assert all(Path(p).exists() for p in result['artifacts'].values())
    cmd.fail=True
    result=run_program(cmd,dict(message,view=dict(operation='program',code='cmd.color("blue")')),CATALOG,tmp_path)
    assert result['rollback']=='complete' and cmd.state['color']=='red'
    run_program(cmd,dict(view=dict(operation='undo')),CATALOG,tmp_path)
    assert cmd.state['color']=='green'


def test_zero_count_is_valid(tmp_path):
    result=run_program(Cmd(),dict(id='123-abcdef012345',view=dict(operation='program',
        code='cmd.count_atoms("none")')),CATALOG,tmp_path)
    assert result['status']=='complete'
    assert result['calls'][0]['result']==0


def test_pinned_skill():
    content,source=skill_context()
    assert 'PyMOL' in content and len(source['commit'])==40


@pytest.mark.parametrize('pending',[False,True])
def test_agent_repair_and_pending(tmp_path,pending):
    submitted=[];plans=[]
    def dispatch(root,view,catalog,executable):
        submitted.append(view)
        return dict(status='queued',operation_id='123-abcdef012345',result='receipt')
    def waiter(queued):
        if len(submitted)==1:return dict(status='complete',scene=dict(objects=CATALOG))
        if pending:return dict(status='queued',pending=True)
        if len(submitted)==2:return dict(status='failed',error='unknown color',rollback='complete')
        return dict(status='complete',artifacts={})
    def planner(prompt,profile,**kwargs):
        plans.append(json.loads(prompt))
        assert 'Upstream skill' in kwargs['system_prompt']
        return dict(explanation='Color the protein.',code='cmd.color("red")'),{}
    answer=run_agent('color protein',{},tmp_path,CATALOG,planner=planner,dispatch=dispatch,waiter=waiter)
    assert len(submitted)==(2 if pending else 3)
    assert answer['status']==('queued' if pending else 'complete')
    if not pending:assert plans[1]['previous_attempt']['error']=='unknown color'


def test_agent_route():
    assert validate_route(dict(intent='pymol_agent',message='',task_id=None,request='show pocket'))['intent']=='pymol_agent'
    with pytest.raises(ValueError):validate_route(dict(intent='pymol_agent',message='',task_id=None,request=''))


def test_pocket_plan_executes_and_exports_residues(tmp_path):
    class PocketCmd(Cmd):
        def __getattr__(self,name):
            return lambda *args,**kwargs:self.calls.append((name,args,kwargs))
        def color(self,*args,**kwargs):self.calls.append(('color',args,kwargs))
    cmd=PocketCmd();submitted=[]
    catalog=[dict(CATALOG[0],ligand=dict(chain='A',resi='201',resn='HTZ'))]
    def dispatch(root,view,catalog,executable):
        submitted.append(view)
        return run_program(cmd,dict(id='pocket-test',view=view),catalog,root)|dict(operation_id='pocket-test')
    def planner(prompt,profile,**kwargs):
        assert json.loads(prompt)['request']=='show residues around ligand in 5 angstrom for v001'
        assert 'pocket_view' in kwargs['capabilities']['methods']
        return dict(explanation='Show the ligand neighborhood.',calls=[dict(method='pocket_view',
                    arguments=dict(objects='v001',radius=5))]),dict(source='test_planner')
    answer=run_agent('show residues around ligand in 5 angstrom for v001',{},tmp_path,catalog,
                     dispatch=dispatch,waiter=lambda value:value,planner=planner)
    assert answer['status']=='complete'
    pocket=answer['result']['calls'][0]['result']
    assert pocket['objects']==['v001'] and pocket['residues']['v001']
    assert Path(pocket['artifacts']['residues_json']).exists()
    selections=[str(call) for call in cmd.calls if isinstance(call,tuple) and call[0]=='show']
    assert any('byres' in text and 'within 5' in text and 'resn HTZ' in text for text in selections)
    assert ('enable',('v001',),{}) in cmd.calls


def test_pocket_adapter_validates_arguments():
    with pytest.raises(ValueError):compile_program('cmd.pocket_view(objects="v001",radius=50)')
    with pytest.raises(ValueError):compile_program('cmd.pocket_view(objects="all",radius=5)')


def test_history_reaches_planner_and_final_error_is_actionable(tmp_path):
    def dispatch(*args):return dict(operation_id='history-test',result='unused')
    def planner(prompt,profile,**kwargs):
        assert json.loads(prompt)['conversation']==[dict(role='user',text='Keep v001 visible')]
        raise ValueError('Unsupported label expression: example')
    answer=run_agent('show its pocket',{},tmp_path,CATALOG,planner=planner,dispatch=dispatch,
        waiter=lambda queued:dict(status='complete',scene=dict(objects=CATALOG)),
        conversation=[dict(role='user',text='Keep v001 visible')])
    assert 'Unsupported label expression' in answer['result']['message']


def test_structured_calls_quote_identifiers_and_large_scene():
    plan=validate_plan(dict(explanation='Style multiple complexes.',calls=[dict(method='color',
        arguments=dict(color='cyan',selection='v001 and resn 6ZT')) for _ in range(120)]))
    calls=compile_program(plan['code'])
    assert len(calls)==120 and calls[0][1]['selection']=='v001 and resn 6ZT'
    with pytest.raises(ValueError):
        validate_plan(dict(explanation='',calls=[dict(method='do',arguments=dict(command='quit'))]))
    with pytest.raises(ValueError,match='512'):
        validate_plan(dict(explanation='',calls=[{}]*513))


def test_invalid_source_preserved_for_repair(tmp_path):
    submitted=[];contexts=[]
    bad='cmd.color("red", 6ZT)'
    def dispatch(root,view,catalog,executable):
        submitted.append(view)
        return dict(status='queued',operation_id='123-abcdef012345',result='receipt')
    def waiter(queued):
        return dict(status='complete',scene=dict(objects=CATALOG))
    def planner(prompt,profile,**kwargs):
        contexts.append(json.loads(prompt))
        raw=(dict(explanation='Invalid source.',code=bad) if len(contexts)==1 else
            dict(explanation='Color ligand.',calls=[dict(method='color',arguments=dict(color='red',selection='v001'))]))
        return kwargs['validator'](raw),{}
    answer=run_agent('color ligand',{},tmp_path,CATALOG,planner=planner,dispatch=dispatch,waiter=waiter)
    assert answer['status']=='complete'
    previous=contexts[1]['previous_attempt']
    assert previous['plan']['code']==bad and previous['line']==1
    assert len(submitted)==2  # No program dispatched for the invalid plan.
    assert json.loads((tmp_path/'123-abcdef012345-agent.json').read_text())['attempts'][0]['plan']['code']==bad
