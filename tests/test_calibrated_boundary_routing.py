import numpy as np
import pytest

from aidd_agent.calibrated_boundary_routing import distances,diverse_indices,fit_class,score_candidate,priority


def test_distinct_disjoint_molecule_sets_and_sparse_guard():
    records=[(str(i),'m'+str(i),np.zeros((3,29))) for i in range(64)]
    model=fit_class(records)
    assert model['status']=='routable' and model['check_retention']==1
    parts=[set(model[k]) for k in ('fit_ids','calibration_ids','check_ids')]
    assert len(set.union(*parts))==64 and not any(parts[i]&parts[j] for i in range(3) for j in range(i))
    assert set(model['prototype_ids'])<=parts[0]
    assert fit_class(records[:8])['status']=='insufficient_distinct_references'
    bad=records.copy();bad[1]=(bad[1][0],bad[0][1],bad[1][2])
    with pytest.raises(ValueError):fit_class(bad)


def test_joint_support_cannot_mix_different_prototypes():
    a=np.zeros((2,29));b=a.copy();a[:,:11]=10;b[:,11:23]=10
    model=dict(status='routable',prototype_vectors=[a.tolist(),b.tolist()],
               prototype_ids=['property_bad','steric_bad'],scales=[1,1,1],cutoff=1)
    result=score_candidate(np.zeros((2,29)),model,0)
    assert result['score']==10
    assert score_candidate(a,model,0)['score']==0
    with pytest.raises(ValueError):score_candidate(a,model,2)


def test_diversity_selection_and_directed_rotation():
    points=np.array([np.zeros((2,29)),np.ones((2,29)),np.full((2,29),10)])
    assert diverse_indices(points,2)==[0,2]
    source=np.zeros((2,29));source[0,:]=1
    model=dict(status='routable',prototype_vectors=[source.tolist()],prototype_ids=['ref'],scales=[1,1,1],cutoff=1)
    assert score_candidate(np.roll(source,1,axis=0),model,1)['score']==0
    np.testing.assert_allclose(distances(source,np.array([source])),0)


def test_failed_disjoint_check_prevents_routing():
    mids=sorted(['m'+str(i) for i in range(64)],key=priority)
    records=[(mid,mid,np.zeros((2,29)) if i<48 else np.ones((2,29))) for i,mid in enumerate(mids)]
    model=fit_class(records)
    assert model['status']=='independent_check_failed' and model['check_retention']==0
    assert score_candidate(np.zeros((2,29)),model,0) is None
