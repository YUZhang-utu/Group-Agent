"""Reference-normalized coverage and same-pose diagnostic regression."""
import sys
from pathlib import Path
from types import SimpleNamespace
import time
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import diagnose_contact_features as diag
import test_native_state_pose_audit as fixtures


class CoverageTests(unittest.TestCase):
    def fixture(self):
        return fixtures.NativeStateTests().fixture()

    def test_remote_extra_features_do_not_penalize_coverage(self):
        state,features,candidate=self.fixture()
        before=diag.match_contact_features(state,features,np.eye(4))[2]
        features.feature_points=np.vstack([features.feature_points,[[100,100,100]]])
        features.feature_types=np.append(features.feature_types,1)
        features.feature_directions=np.vstack([features.feature_directions,[[0,0,0]]])
        features.feature_kinds=np.append(features.feature_kinds,0)
        after=diag.match_contact_features(state,features,np.eye(4))[2]
        self.assertAlmostEqual(before,1.0)
        self.assertAlmostEqual(after,before)
        from aidd_agent.gaussian_overlay import score_overlay
        data=state['aligned']
        score=score_overlay(data['shape_points'],candidate.shape_points,data['feature_points'],
            data['feature_types'],features.feature_points,features.feature_types,np.eye(4).ravel(),cutoff=None)
        self.assertLess(score['color']['tanimoto'],1.0)

    def test_missing_reference_and_one_to_one(self):
        state,features,_=self.fixture()
        features.feature_points=features.feature_points[:1]
        features.feature_types=features.feature_types[:1]
        features.feature_directions=features.feature_directions[:1]
        features.feature_kinds=features.feature_kinds[:1]
        state['aligned']['feature_points'][:]=features.feature_points[0]
        ok,_,coverage,detail=diag.match_contact_features(state,features,np.eye(4))
        self.assertFalse(ok)
        self.assertAlmostEqual(coverage,1/3)
        self.assertEqual(sum(d['candidate_feature']>=0 for d in detail),1)

    def test_shared_contact_labels_not_double_weighted(self):
        state,features,_=self.fixture()
        features.feature_points[2]+=100
        before=diag.match_contact_features(state,features,np.eye(4))[2]
        state['anchors'].append(dict(state['anchors'][0],anchor_id='duplicate',target_residue=67))
        after=diag.match_contact_features(state,features,np.eye(4))[2]
        self.assertAlmostEqual(before,after)

    def test_search_control_and_baseline(self):
        state,_,_=self.fixture()
        candidate,features=diag.controlled_candidate(state,'T')
        result=diag.search(candidate,features,{'T':state},time.monotonic()+30)
        self.assertTrue(result['complete'])
        self.assertTrue(result['original_pass'])
        self.assertTrue(result['replacement_0.9'])
        self.assertAlmostEqual(result['best']['shape']['shape'],1.0)
        self.assertTrue(result['best']['coverage']['feature_assignments'])
        self.assertIn('aligned_heavy_atom_coordinates',result['best']['shape'])
        from aidd_agent.gaussian_overlay import score_overlay
        pose=result['best']['shape'];data=state['aligned']
        score=score_overlay(data['shape_points'],candidate.shape_points,data['feature_points'],
            data['feature_types'],candidate.feature_points,candidate.feature_types,pose['transform'],cutoff=None)
        self.assertAlmostEqual(score['color']['tanimoto'],pose['color'])
        expired=diag.search(candidate,features,{'T':state},0)
        self.assertFalse(expired['complete'])
        self.assertIsNone(diag.summarize([dict(expired,kind='library')])['original_pass_fraction'])

    def test_direction_and_contact_compatibility(self):
        state,features,_=self.fixture()
        state['aligned']['feature_direction_kinds'][:]=1
        state['aligned']['feature_directions'][:]=[1,0,0]
        features.feature_kinds[:]=1
        features.feature_directions[:]=[-1,0,0]
        self.assertEqual(diag.match_contact_features(state,features,np.eye(4))[2],0.0)
        features.feature_directions[:]=[1,0,0]
        a=diag.match_contact_features(state,features,np.eye(4))
        b=diag.native.native_contacts(state,features,np.eye(4))
        self.assertEqual(a[:2],b[:2])


if __name__=='__main__':unittest.main()
