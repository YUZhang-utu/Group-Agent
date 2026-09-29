import math

import pytest

from aidd_agent.boundary_block_candidates import VERSION,candidates,digest,omega_angles,recover_chirality


def row(states,angles):
    vector=[]
    for angle in angles:
        x=math.radians(angle);vector.extend([0,0,math.sin(x),1,1,math.cos(x)])
    return dict(descriptor=vector,provenance=dict(omega_states=states))


def test_angles_and_chirality_recovery():
    states=['boundary','trans'];labels=['R','S']
    group=digest([VERSION,'backbone',2,labels,states])
    recovered,trials=recover_chirality({group:states},'backbone')
    assert recovered[group]==labels and trials==9
    assert omega_angles(row(states,[-40,-170]),'backbone')==pytest.approx([40,170])
    with pytest.raises(ValueError):omega_angles(row(states,[0,170]),'backbone')
    assert recover_chirality({group:states},'backbone',maximum_trials=1)[0]=={}


def test_candidates_preserve_definite_states_and_handedness():
    target=digest([VERSION,'backbone',2,['R','S'],['cis','trans']])
    matches,deviation=candidates(['boundary','trans'],[40,170],['R','S'],'backbone',{target})
    assert matches==[dict(hard_group=target,rotation=0)] and deviation==40
    assert candidates(['boundary','trans'],[40,170],['S','S'],'backbone',{target})[0]==[]
    assert candidates(['boundary','cis'],[40,20],['R','S'],'backbone',{target})[0]==[]
    assert candidates(['boundary','trans'],[90,170],['R','S'],'backbone',{target})[0]==[]
    assert candidates(['boundary','trans'],[40,170],None,'backbone',{target})[0]==[]
    assert candidates(['trans','boundary'],[170,40],['S','R'],'backbone',{target})[0]==[dict(hard_group=target,rotation=1)]
