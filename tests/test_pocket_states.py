import json
from pathlib import Path

import numpy as np
import pytest

from aidd_agent import pocket_states as ps


def atom(point, name='CB', residue='ALA', index=1, element='C'):
    return dict(xyz=np.array(point,float),auth_atom_id=name,auth_comp_id=residue,
                canonical_residue=index,type_symbol=element,label_alt_id='',occupancy='1',
                pdbx_PDB_ins_code='')


def test_alignment_and_small_noise_preserve_pocket_state():
    p=ps.policy(); grid=ps.pocket_grid([[0,0,0],[2,0,0]],p)
    atoms=[atom([x,y,z],index=i) for i,(x,y,z) in enumerate([
        [-3,0,0],[3,0,0],[0,-3,0],[0,3,0],[0,0,-3],[0,0,3]])]
    fixed=np.array([a['xyz'] for a in atoms]); rot=np.array([[0,-1,0],[1,0,0],[0,0,1.]])
    moved=fixed@rot.T+[10,20,30]
    transform,_,_=ps.fit_protein(moved,fixed)
    aligned=[dict(a,xyz=x@transform[:3,:3].T+transform[:3,3]) for a,x in zip(atoms,moved)]
    a,fa=ps.describe(grid,atoms,p); b,fb=ps.describe(grid,aligned,p)
    assert np.array_equal(a,b) and np.array_equal(fa,fb)
    noisy=[dict(a,xyz=a['xyz']+.04) for a in atoms]
    c,fc=ps.describe(grid,noisy,p)
    distance,*_=ps.compare(np.array([a,c]),np.array([fa,fc]),grid,[[0,0,0]],p)
    rows=[dict(id='A',pdb_id='1AAA'),dict(id='B',pdb_id='2AAA')]
    assert len(ps.clusters(distance,rows,p['cluster_distance']))==1


def test_material_occlusion_and_chemistry_are_measured():
    p=ps.policy(); grid=ps.pocket_grid([[0,0,0]],p)
    a=np.ones(len(grid),bool); b=grid[:,0]<0
    f=np.zeros((2,6,len(grid)),bool)
    distance,overlap,_,local=ps.compare(np.array([a,b]),f,grid,[[0,0,0]],p)
    assert overlap[0,1]<.5 and local[0,1]>.5
    rows=[dict(id='A',pdb_id='1AAA'),dict(id='B',pdb_id='2AAA')]
    groups=ps.clusters(distance,rows,p['cluster_distance'])
    assert len(groups)==2 and all(g['small_support'] for g in groups)
    f[0,0,:]=True; f[1,1,:]=True
    distance,overlap,chemical,_=ps.compare(np.array([a,a]),f,grid,[],p)
    assert overlap[0,1]==1 and chemical[0,1]==1 and distance[0,1]>0


def test_support_is_not_chain_multiplicity_and_rare_states_survive():
    rows=[dict(id=str(i),pdb_id='1AAA' if i<3 else '2AAA') for i in range(4)]
    distances=np.zeros((4,4));distances[3,:3]=distances[:3,3]=.8
    groups=ps.clusters(distances,rows,.3)
    assert len(groups)==2
    assert [g['distinct_pdb_support'] for g in groups]==[1,1]
    assert sum(len(g['members']) for g in groups)==4


def test_chain_duplication_does_not_overweight_linkage():
    basic=np.array([[0,.25,.8],[.25,0,.3],[.8,.3,0]])
    owners=[0]+[1]*20+[2]
    rows=[dict(id=str(i),pdb_id=str(owner)) for i,owner in enumerate(owners)]
    groups=ps.clusters(basic[np.ix_(owners,owners)],rows,.4)
    assert len(groups)==2
    assert sorted(g['distinct_pdb_support'] for g in groups)==[1,2]


def test_missing_atoms_are_not_interpreted_as_open_pocket():
    atoms=[atom([i,0,0],name=name) for i,name in enumerate(['N','CA','C','O'])]
    assert ps.quality(atoms,{1})==['missing_pocket_atoms:1:CB']
    atoms.append(atom([4,0,0]))
    assert not ps.quality(atoms,{1})
    # Missing remote atoms are not a reason to discard the observed pocket.
    assert not ps.quality(atoms[:-1],{1},{1:{'N','CA','C','O'}})


def test_adoption_provenance_and_state_restricted_consensus(tmp_path):
    raw=tmp_path/'structure.cif';raw.write_text('experimental fixture')
    ids=['pocket-000000000001','pocket-000000000002']
    rows=[dict(id=code+':A',pdb_id=code,target_chain='A',transform=np.eye(4).tolist(),
               queries=[dict(query_id=code+':LIG:B:1')]) for code in ('1AAA','2AAA')]
    report=dict(kind='pocket_states',readiness='needs_user_adoption',sources=ps.fingerprint([raw]),
                diversity_report='diversity.json',policy=ps.policy(),limitations=[],structures=rows,
                clusters=[dict(id=sid,members=[row['id']],representative=row['id']) for sid,row in zip(ids,rows)])
    path=tmp_path/'report.json';path.write_text(json.dumps(report))
    with pytest.raises((KeyError,ValueError)):ps.state_cohort(path)
    adoption=ps.adopt(path,tmp_path/'adopt')
    adopted_path=tmp_path/'adopt/report.json'
    assert adoption['readiness']=='ready_for_state_consensus'
    with pytest.raises(ValueError,match='Choose one'):ps.state_cohort(adopted_path)
    _,cohort,_=ps.state_cohort(adopted_path,ids[1])
    assert [row['query_id'] for row in cohort['admitted']]==['2AAA:LIG:B:1']
    raw.write_text('changed coordinates')
    with pytest.raises(ValueError):ps.state_cohort(adopted_path,ids[1])
    with pytest.raises(ValueError):ps.adopt(path,tmp_path/'changed')


def test_new_diversity_cannot_bypass_review(tmp_path):
    from aidd_agent.consensus_model import build
    (tmp_path/'protein.json').write_text('{}')
    source=tmp_path/'report.json';source.write_text(json.dumps(dict(pocket_state_review_required=True)))
    with pytest.raises(ValueError,match='explicit adoption'):build(source,tmp_path/'consensus')


def test_route_and_plan_validate_pocket_steps():
    from aidd_agent.chat_agent import validate_route
    from aidd_agent.prompt_plan import validate_plan
    assert validate_route(dict(intent='pockets',message='',task_id=None,request='',reference={}))
    for action,params in [('pocket_states',{'cluster_distance':.3}),('pocket_adopt',{}),
                          ('pocket_consensus',{'pocket_state_id':'pocket-000000000001'})]:
        assert validate_plan(dict(version=1,summary='Pocket workflow',clarifications=[],
            steps=[dict(id='pocket',action=action,params=dict(source_run='PROMPT-'+'a'*16,**params))]))
    with pytest.raises(ValueError):ps.policy({'cluster_distance':float('nan')})


def test_pymol_catalog_uses_actual_representative_files(tmp_path):
    from aidd_agent.structure_review import structures
    structure=tmp_path/'1AAA.cif';structure.write_text('experimental coordinates fixture')
    path=tmp_path/'pockets.json'
    path.write_text(json.dumps(dict(kind='pocket_states',sources=ps.fingerprint([structure]),
        clusters=[dict(representative='1AAA:A')],structures=[dict(id='1AAA:A',
            structure_path=str(structure),transform=np.eye(4).tolist())])))
    report=tmp_path/'execution.json';report.write_text(json.dumps(dict(steps={
        'pocket':dict(action='pocket_states',status='complete',result=dict(report=str(path)))})))
    job=dict(status='complete',report=str(report),plan=str(tmp_path/'plan.json'))
    catalog=structures(job,tmp_path)
    assert len(catalog)==1
    structure.write_text('modified')
    with pytest.raises(ValueError,match='provenance mismatch'):structures(job,tmp_path)
