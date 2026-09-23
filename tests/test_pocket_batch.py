import numpy as np
import pytest
from scipy.spatial import cKDTree
from aidd_agent.guided_filters import GuidedFilter


def scalar(gate,points,matrices):
    result=[]
    for matrix in matrices:
        moved=points@matrix[:3,:3].T+matrix[:3,3]
        distances,_=gate.tree.query(moved,k=1)
        physical=float(np.mean(distances<gate.cutoff-1e-8))<=gate.fraction
        excluded=any(r['mode']=='hard' and np.any(np.linalg.norm(moved-np.array(r['center']),axis=1)<r['radius']) for r in gate.design['exclusions'])
        result.append(physical and not excluded)
    return np.asarray(result,bool)


@pytest.mark.parametrize('fraction',[0.,.2,1.])
def test_batched_pocket_matches_scalar_across_blocks(fraction):
    rng=np.random.default_rng(62)
    gate=GuidedFilter.__new__(GuidedFilter)
    gate.tree=cKDTree(rng.normal(size=(500,3))*5)
    gate.cutoff=1.2;gate.fraction=fraction
    gate.design=dict(exclusions=[dict(mode='hard',center=[4,0,0],radius=.7),dict(mode='soft',center=[0,0,0],radius=10)])
    points=rng.normal(size=(61,3))
    matrices=np.repeat(np.eye(4)[None],517,axis=0)
    matrices[:,:3,:3]=np.linalg.qr(rng.normal(size=(517,3,3)))[0]
    matrices[:,:3,3]=rng.normal(size=(517,3))*5
    np.testing.assert_array_equal(gate.pocket_mask(points,matrices),scalar(gate,points,matrices))
    assert gate.pocket_mask(points,matrices[:0]).shape==(0,)


def test_pocket_strict_boundary_unchanged():
    gate=GuidedFilter.__new__(GuidedFilter);gate.tree=cKDTree([[0.,0,0]])
    gate.cutoff=1.2;gate.fraction=0.;gate.design=dict(exclusions=[])
    matrices=np.repeat(np.eye(4)[None],3,axis=0)
    matrices[:,0,3]=[1.2-1e-8-1e-10,1.2-1e-8,1.2-1e-8+1e-10]
    np.testing.assert_array_equal(gate.pocket_mask(np.zeros((1,3)),matrices),[False,True,True])
