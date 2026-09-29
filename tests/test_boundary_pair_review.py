import json
import sqlite3

import numpy as np
import pytest

from aidd_agent.boundary_pair_review import aligned_profile,compare,rmsd,select_panel


def test_balanced_selection_keeps_all_candidate_rotations():
    with sqlite3.connect(':memory:') as db:
        db.execute('CREATE TABLE candidate(cid,original_group,maximum_deviation,candidates)')
        matches=json.dumps([dict(hard_group='core',rotation=0),dict(hard_group='core',rotation=1)])
        db.executemany('INSERT INTO candidate VALUES(?,?,?,?)',[
            ('a','same',34,matches),('b','same',34,matches),('c','other',40,matches),
            ('d','third',60,matches),('excluded','none',34,'[]')])
        panel,populations=select_panel(db,100)
        assert len(panel)==3 and populations['up_to_35']==2
        assert all(len(r['matches'])==2 for r in panel)
        assert panel==select_panel(db,100)[0]


def test_proper_rmsd_and_profile_atom_mapping():
    points=np.random.default_rng(18).normal(size=(9,3))
    rot=np.array([[0,-1,0],[1,0,0],[0,0,1]])
    assert rmsd(points,points@rot+[2,4,8])<1e-12
    mirrored=points.copy();mirrored[:,0]*=-1
    assert rmsd(points,mirrored)>.01
    units=[[0,1,2,3],[4,5,6,7]]
    value=dict(units=units,omega_states=['cis','trans'],properties=[[1],[2]],
               typed_spatial_moments=[[3],[4]],sterics=[{'id':1},{'id':2}])
    result=aligned_profile(value,units[::-1])
    assert result['omega_states']==['trans','cis'] and result['properties']==[[2],[1]]
    with pytest.raises(ValueError):aligned_profile(value,[[0,1,2,3],[0,1,2,3]])


def test_candidate_rotation_applies_to_query_not_reference():
    xyz=np.random.default_rng(7).normal(size=(3,3,3))
    sterics=[dict(side_heavy_atoms=i+1,proximal_branch_excess=0,local_lower_vdw=[-1]*3,
                   local_upper_vdw=[i+1]*3,spatial_eigenvalues=[1]*3,maximum_reach=i+2) for i in range(3)]
    reference=dict(properties=[[i]*11 for i in range(3)],sterics=sterics,backbone_coordinates=xyz.tolist())
    query={k:np.roll(v,1,axis=0).tolist() if k!='sterics' else v[-1:]+v[:-1] for k,v in reference.items()}
    result=compare(reference,query,1)
    assert result['property_rms']==0 and result['steric_rms']==0
    assert result['backbone_rmsd_angstrom']<1e-12
