import copy
import json
from types import SimpleNamespace

import numpy as np
import pytest

from aidd_agent import structure_survey as survey
from aidd_agent.anchor_advice import validate_advice, recommend
from aidd_agent.guided_filters import GuidedFilter, rule_passes
from aidd_agent import preselection_full as pre
from test_preselection_full import geometry


def evidence():
    return dict(queries=[dict(query_id='Q', evidence_id='crystal:Q', pdb_id='ABCD', ccd_id='LIG',
        mapped_residues=[20], anchors=[dict(anchor_id=a,feature_index=i,feature_class='HBD',
        target_residue=20,evidence=dict(protein_partner=dict(atom_name='N'))) for i,a in enumerate(['A','B'])])],
        recurrence=[])


def advice():
    return dict(query_id='Q',mandatory_anchors=['A'],alternative_groups=[],optional_anchors=['B'],
                evidence_ids=['crystal:Q'],rationale='Exploratory structural hypothesis.')


def test_recurrence_counts_structures_not_copies_or_contact_rows():
    q=evidence()['queries'][0]
    other=copy.deepcopy(q);other.update(query_id='Q2',pdb_id='EFGH',ccd_id='OTHER')
    result=survey.recurrence([q,copy.deepcopy(q),other])
    assert len(result)==1
    assert result[0]['distinct_structures']==2
    assert result[0]['distinct_ligands']==2
    assert result[0]['eligible_structures']==2


def test_grounded_advice_and_same_pose_rule():
    s=evidence();a=advice()
    assert validate_advice(a,s)==a
    assert rule_passes(1,['A','B'],a)
    assert not rule_passes(2,['A','B'],a)
    a['alternative_groups']=[['B']];a['optional_anchors']=[]
    assert not rule_passes(1,['A','B'],a)
    assert rule_passes(3,['A','B'],a)
    for change in [dict(mandatory_anchors=['invented']),dict(evidence_ids=['invented']),
                   dict(optional_anchors=['A']),dict(mandatory_anchors=[],alternative_groups=[])]:
        with pytest.raises(ValueError):validate_advice(dict(a,**change),s)


def test_pocket_filter_uses_transformed_coordinates(tmp_path):
    q,_,_=geometry()
    path=tmp_path/'receptor.npz';np.savez(path,points=np.array([[10.,0,0]]))
    design=dict(advice(),minimum_score=.5,receptor_npz=str(path),
                pocket=dict(minimum_heavy_atom_distance=1.2,maximum_clashing_fraction=0.))
    gate=GuidedFilter(evidence()['queries'][0],q,design)
    matrices=np.array([np.eye(4),np.eye(4)]);matrices[1,0,3]=10
    assert gate.pocket_mask(np.array([[0.,0,0]]),matrices).tolist()==[True,False]


def test_target_mapping_requires_accession_and_sequence_not_author_number():
    pytest.importorskip('Bio')
    sequence='ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMN'
    protein=dict(accession='P1',sequence='WWW'+sequence+'WWW')
    data={'_struct_ref.id':['1'],'_struct_ref.pdbx_db_accession':['P1'],
          '_struct_ref.entity_id':['1'],'_struct_ref.db_name':['UNP'],
          '_entity_poly.entity_id':['1'],'_entity_poly.pdbx_seq_one_letter_code_can':[sequence],
          '_atom_site.group_PDB':['ATOM'],'_atom_site.label_entity_id':['1'],
          '_atom_site.label_seq_id':['1'],'_atom_site.auth_asym_id':['A'],
          '_atom_site.auth_seq_id':['999'],'_atom_site.pdbx_PDB_ins_code':['?']}
    assert survey.target_residues(data,protein)[0]=={('A','999'):4}
    assert survey.target_residues(data,dict(protein,accession='P2'))[0]=={}


@pytest.mark.parametrize('failed',[False,True])
def test_missing_structures_fall_back_without_inventing_coordinates(tmp_path,failed):
    def search(*args,**kwargs):
        if failed:raise TimeoutError()
        return {},[]
    services=SimpleNamespace(pdb_search=search,fetch_json=lambda url:dict(hitCount=1,
        resultList=dict(result=[dict(id='1',source='MED',title='Evidence',abstractText='Binding study')])) )
    protein=dict(accession='P1',organism='human',sequence_sha256='hash',source_url='fixture',genes=[])
    report=survey.survey(protein,tmp_path/'survey',services)
    assert report['readiness']=='needs_structure_input' and not report['queries']
    assert report['pdb_search_status']==('failed' if failed else 'complete')
    proposal=recommend(tmp_path/'survey/report.json',tmp_path/'advice','gpt',
        adviser=lambda _:dict(summary='Evidence needs review.',evidence_ids=['paper:MED:1'],needed_inputs=['Coordinates']))
    assert proposal['readiness']=='needs_structure_input'
    assert 'recommendation' not in proposal


@pytest.mark.parametrize('allowed',[[False,False],[True,False]])
def test_gaussian_only_receives_pocket_survivors(tmp_path,monkeypatch,allowed):
    query,features,seeds=geometry()
    candidate=SimpleNamespace(**vars(features),shape_points=query['shape_points'])
    policy=dict(required_anchors=['A','B'],minimum_score=.5,match_mode='any',coarse_constraints={})
    q=dict(condition_policy=policy,guided_design=advice(),anchors=[dict(anchor_id='A',score_column=0),dict(anchor_id='B',score_column=1)])
    monkeypatch.setattr(pre.full,'_STATE',(SimpleNamespace(get=lambda gid:candidate),query,query,q))
    monkeypatch.setattr(pre.full,'_FILTER_READER',SimpleNamespace(get=lambda gid:features))
    monkeypatch.setattr(pre.full,'_BOUND',SimpleNamespace(check=lambda f:(True,None)))
    monkeypatch.setattr(pre,'coarse_stages',lambda *a,**k:4)
    monkeypatch.setattr(pre,'prepare_seeds',lambda *a,**k:(seeds,2))
    monkeypatch.setattr(pre,'possible_seed_mask',lambda *a:np.ones(2,dtype=bool))
    gate=SimpleNamespace(geometry_possible=lambda f:True,geometry_seeds=lambda f,s,p:p,
                         pocket_mask=lambda p,t:np.array(allowed))
    monkeypatch.setattr(pre,'_GUIDED_CACHE',(json.dumps(q['guided_design'],sort_keys=True),gate))
    received=[]
    def score(*a,**k):
        received.extend(k['prepared'][0][1])
        return {pre.full.OBJECTIVE+'__objective':np.array([.9])}
    monkeypatch.setattr(pre,'_score_ids',score)
    report=pre.compute((0,1,str(tmp_path/'chunk.npz')))
    assert len(received)==sum(allowed)
    assert report['counts']['gaussian_evaluated_conformers']==int(any(allowed))
    assert report['counts']['pocket_rejected_seeds']==2-sum(allowed)


def test_adoption_checks_real_reference_and_applies_manual_edits(tmp_path,monkeypatch):
    from aidd_agent import anchor_advice as module
    from aidd_agent.expanded_wee1 import fingerprint
    query,_,_=geometry()
    receptor=tmp_path/'receptor.npz';np.savez(receptor,points=np.array([[20.,20,20]]))
    s=evidence();s['queries'][0].update(receptor_npz=str(receptor),query_npz='fixture')
    source=tmp_path/'survey.json';source.write_text(json.dumps(s))
    proposal=tmp_path/'proposal.json'
    proposal.write_text(json.dumps(dict(readiness='proposal_ready',survey=str(source),recommendation=advice(),
        defaults=module.DEFAULTS,sources=fingerprint([source,receptor]))))
    monkeypatch.setattr(module,'_load_query',lambda p:({},query))
    result=module.adopt(proposal,tmp_path/'accepted',dict(mandatory_anchors=['B'],optional_anchors=['A']))
    assert result['design']['mandatory_anchors']==['B']
    assert result['reference_control']['anchor_rule_passed']
    np.savez(receptor,points=np.array([[0.,0,0]]))
    raw=json.loads(proposal.read_text());raw['sources']=fingerprint([source,receptor]);proposal.write_text(json.dumps(raw))
    with pytest.raises(ValueError,match='reference fails pocket'):
        module.adopt(proposal,tmp_path/'blocked')


def test_post_search_selection_does_not_union_different_poses(tmp_path):
    import sqlite3
    from aidd_agent.guided_workflow import select_candidates
    from aidd_agent.expanded_wee1 import fingerprint
    database=tmp_path/'candidates.sqlite'
    with sqlite3.connect(database) as db:
        pre.schema(db)
        db.execute('INSERT INTO members VALUES (?,?,?,?,?,?)',('m0','cluster','C','',3,2))
        for mask,scores in [(1,[1.,0.]),(2,[0.,1.])]:
            pre.merge_pose(db,dict(molecule_id='m0',mask=mask,min_matched_score=1.,global_id=0,
                                  scores=scores,assignments=[0,1]))
    source=tmp_path/'report.json';source.write_text(json.dumps(dict(anchor_order=['A','B'],
        policy=dict(minimum_score=.5),outputs={'candidates.sqlite':str(database)},output_hashes=fingerprint([database]))))
    params=dict(required_anchors=['A','B'],match_mode='all',minimum_score=.5)
    result=select_candidates(source,tmp_path/'all',params)
    assert result['matching_molecules']==0
    result=select_candidates(source,tmp_path/'any',dict(params,match_mode='any'))
    assert result['matching_molecules']==1 and result['pose_records']==2
    with pytest.raises(ValueError,match='threshold'):
        select_candidates(source,tmp_path/'invalid',dict(params,minimum_score=.4))


def test_chat_guided_actions_are_session_owned_and_sequential(tmp_path):
    from aidd_agent.chat_agent import ChatAgent
    app=ChatAgent(tmp_path,start=False,router=lambda *args:dict(intent='run',request='Survey human WEE1',task_id=None,message=''))
    sid=app.new_session();app.ask(sid,'Survey human WEE1','deepseek')
    job=app.task(sid);folder=tmp_path/'PROMPT-0123456789abcdef';folder.mkdir()
    plan=folder/'plan.json';plan.write_text('{}')
    report=folder/'report.json';report.write_text(json.dumps(dict(steps=dict(survey=dict(action='structure_survey',status='complete')))))
    app.update(job['id'],status='complete',plan=str(plan),report=str(report))
    app.ask(sid,'/recommend '+job['id'],'deepseek')
    queued=app.task(sid)
    assert json.loads(__import__('pathlib').Path(queued['plan']).read_text())['plan']['steps'][0]['action']=='anchor_recommend'
    with pytest.raises(ValueError):app.ask(app.new_session(),'/recommend '+job['id'],'deepseek')
    with pytest.raises(ValueError):app.ask(sid,'/guided '+job['id'],'deepseek')
    app.close()
