"""Native-state pairing, contact identity and fresh-pose gate regression."""
import importlib.util
from pathlib import Path
import sys
from types import SimpleNamespace
import time
import unittest
from unittest.mock import patch
import numpy as np
from scipy.spatial import cKDTree

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
import audit_native_state_poses as audit


class NativeStateTests(unittest.TestCase):
    def fixture(self):
        xyz=np.array([[x,y,0.] for x in [-5.,0.,5.] for y in [-.2,0.,.2]])
        fp=np.array([[-5.,0.,0.],[0.,0.,0.],[5.,0.,0.]])
        data=dict(shape_points=xyz,feature_points=fp,feature_types=np.full(3,5),
                  feature_directions=np.zeros((3,3)),feature_direction_kinds=np.zeros(3,dtype=np.uint8),
                  anchored_weights=np.ones(3),anchor_feature_indices=np.arange(3))
        anchors=[dict(anchor_id='a'+str(i),feature_index=i,target_residue=r,protein_atom=atom,
                      feature_class='hydrophobic') for i,(r,atom) in enumerate([(61,'CD1'),(54,'CD1'),(93,'CG1')])]
        tree=cKDTree(xyz+[0,0,3]);contact,regions=audit.base.make_regions(xyz,tree)
        state=dict(aligned=data,tree=tree,contact=contact,regions=regions,anchors=anchors,
                   record=dict(supported_for_requested_contact_rule=True))
        feature=SimpleNamespace(feature_points=fp,feature_types=np.full(3,5),feature_directions=np.zeros((3,3)),
                                feature_kinds=np.zeros(3,dtype=np.uint8),conformer_id='c1',molecule_id='m1')
        candidate=SimpleNamespace(shape_points=xyz,feature_points=fp,feature_types=np.full(3,5),
                                  conformer_id='c1',molecule_id='m1')
        return state,feature,candidate

    def test_pairing_preserves_native_distances(self):
        receptor=np.array([[0,0,3],[2,0,3],[0,2,3]])
        ligand=np.array([[0,0,0],[2,0,0],[0,2,0]])
        m=np.eye(4);m[:3,:3]=[[0,-1,0],[1,0,0],[0,0,1]];m[:3,3]=[8,3,-4]
        moved,d=audit.paired_receptor(receptor,ligand,audit.move(ligand,m),m)
        self.assertAlmostEqual(d,3)
        self.assertTrue(np.allclose(moved,audit.move(receptor,m)))
        with self.assertRaises(ValueError):audit.paired_receptor(receptor,ligand,ligand,m)
        m[0,0]=2
        with self.assertRaises(ValueError):audit.rigid(m)

    def test_native_contacts_and_missing_mode(self):
        state,feature,_=self.fixture()
        self.assertTrue(audit.native_contacts(state,feature,np.eye(4))[0])
        state['anchors'][2]['target_residue']=54
        self.assertFalse(audit.native_contacts(state,feature,np.eye(4))[0])
        state['anchors']=[]
        self.assertEqual(audit.native_contacts(state,feature,np.eye(4)),(False,[],[]))

    def test_feature_assignment_not_contact_list_presence(self):
        state,feature,_=self.fixture()
        feature.feature_points=feature.feature_points+100
        self.assertFalse(audit.native_contacts(state,feature,np.eye(4))[0])

    def test_contact_extraction_chain_atom_and_type_filters(self):
        state,_,_=self.fixture()
        atoms=[dict(auth_asym_id='A',auth_seq_id='61',canonical_residue=61,auth_atom_id='CD1')]
        row=dict(anchor_id='x',feature_index=0,feature_class='hydrophobic',
                 evidence=dict(protein_partner=dict(chain='A',residue_number='61',atom_name='CD1')))
        result=audit.filtered_contacts([row,row],state['anchors'],atoms,'A')
        self.assertEqual(len(result),1)
        row['evidence']['protein_partner']['chain']='B'
        self.assertEqual(audit.filtered_contacts([row],state['anchors'],atoms,'A'),[])

    def test_fresh_seed_search_and_timeout(self):
        state,feature,candidate=self.fixture()
        reader=SimpleNamespace(get=lambda gid:candidate)
        fr=SimpleNamespace(get=lambda gid:feature)
        coarse=dict(heavy_atoms=[1,100],feature_minimum={},extent_lower=[0,0,0],extent_upper=[100,100,100])
        audit._FRESH=(reader,fr,{'T':state},coarse,time.monotonic()+30)
        result=audit.fresh_one(dict(gid='1',cid='c1',mid='m1'))
        self.assertTrue(result['complete'])
        self.assertTrue(result['regions'])
        self.assertGreater(result['tested_seeds'],0)
        self.assertTrue(result['passing_poses'])
        audit._FRESH=(reader,fr,{'T':state},coarse,0)
        expired=audit.fresh_one(dict(gid='1',cid='c1',mid='m1'))
        self.assertFalse(expired['complete'])
        self.assertEqual(expired['tested_seeds'],0)

    def test_correct_receptor_required_for_fresh_acceptance(self):
        state,feature,candidate=self.fixture()
        state['tree']=cKDTree(candidate.shape_points)
        coarse=dict(heavy_atoms=[1,100],feature_minimum={},extent_lower=[0,0,0],extent_upper=[100,100,100])
        audit._FRESH=(SimpleNamespace(get=lambda gid:candidate),SimpleNamespace(get=lambda gid:feature),
                      {'T':state},coarse,time.monotonic()+30)
        from aidd_agent import gaussian_batch
        seed=SimpleNamespace(transform_matrix=np.eye(4).ravel())
        with patch.object(gaussian_batch,'prepare_seeds',return_value=([seed],0)):
            result=audit.fresh_one(dict(gid='1',cid='c1',mid='m1'))
        self.assertTrue(result['contacts'])
        self.assertFalse(result['steric'])
        self.assertFalse(result['regions'])


if __name__=='__main__':unittest.main()
