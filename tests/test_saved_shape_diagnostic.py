"""Kernel bounds, rigid refinement and transform composition."""
import sys
from pathlib import Path
import time
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import diagnose_saved_shapes as shape


class ShapeDiagnosticTests(unittest.TestCase):
    def fixture(self):
        return np.array([[0.,0,0],[1.4,0,0],[2,1.5,0],[0,2,1],[1,1,2],[3,2,2]])

    def test_bound_and_rigid_invariance(self):
        q=self.fixture();p=shape.primitives(q,q)
        self.assertAlmostEqual(p['shape'],1.)
        self.assertAlmostEqual(p['norm_upper_bound'],1.)
        m=shape.delta_transform(np.array([.2,.1,-.1,1,2,-1]),q.mean(axis=0))
        moved=shape.native.move(q,m);other=shape.primitives(q,moved)
        self.assertAlmostEqual(p['candidate_self'],other['candidate_self'])
        self.assertLessEqual(other['shape'],other['norm_upper_bound']+1e-10)

    def test_refinement_and_composition(self):
        q=self.fixture();original=shape.delta_transform(np.array([.1,-.1,.1,.5,-.3,.3]),q.mean(axis=0))
        c=shape.native.move(q,original)
        before=shape.primitives(q,c)['shape']
        result=shape.refine(q,c,time.monotonic()+10)
        self.assertGreater(result['shape'],before)
        self.assertGreater(result['shape'],.999)
        matrix=np.array(result['delta'])
        self.assertTrue(np.allclose(shape.native.move(c,matrix),shape.native.move(q,matrix@original)))
        self.assertLessEqual(result['calls'],180)

    def test_endpoint_recheck_uses_all_gates(self):
        import test_native_state_pose_audit as fixtures
        import check_refined_shape_poses as checker
        from scipy.spatial import cKDTree
        state,features,candidate=fixtures.NativeStateTests().fixture()
        good=checker.assess(candidate,features,state,np.eye(4))
        self.assertTrue(good['original_policy_pass'])
        state['tree']=cKDTree(candidate.shape_points)
        bad=checker.assess(candidate,features,state,np.eye(4))
        self.assertAlmostEqual(bad['shape'],1.)
        self.assertFalse(bad['original_policy_pass'])
        self.assertFalse(bad['contacts_shape_geometry'])

    def test_deadline_and_asymmetric_coverage(self):
        q=self.fixture();c=np.vstack([q,[[100,100,100]]])
        coverage=shape.point_coverage(q,c)
        self.assertEqual(coverage['reference_near_candidate'],1.)
        self.assertLess(coverage['candidate_near_reference'],1.)
        result=shape.refine(q,c,0)
        self.assertEqual(result['status'],'deadline')
        self.assertEqual(result['calls'],0)


if __name__=='__main__':unittest.main()
