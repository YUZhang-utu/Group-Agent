"""Checks for rigid refinement, analytic derivatives and unchanged hard contacts."""
from pathlib import Path
import sys
from types import SimpleNamespace
import unittest
from unittest.mock import patch
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
import kras_contact_refine as refine
import run_kras_contact_refine as runner
import export_kras_soft_rank as exporter


class ContactRefinementTests(unittest.TestCase):
    def test_exact_source_pose_check(self):
        from rdkit import Chem
        state=refine.soft.run.kw.load_state(ROOT/'configs/targets/kras_sip_e116','KRAS_G12D')
        q=state['query'];chem=SimpleNamespace(**state['record']['chemistry'],molecule_id='M',conformer_id='C')
        candidate=SimpleNamespace(molecule_id='M',conformer_id='C',shape_points=q['shape_points'],feature_types=q['feature_types'])
        graph=Chem.RWMol()
        for z in chem.atomic_numbers:graph.AddAtom(Chem.Atom(int(z)))
        molecule=graph.GetMol();conf=Chem.Conformer(len(chem.atomic_numbers))
        for i,point in enumerate(q['shape_points']):conf.SetAtomPosition(i,point)
        molecule.AddConformer(conf)
        with patch('aidd_agent.gaussian_batch.ArtifactCatalogReader',return_value=SimpleNamespace(get=lambda gid:candidate)), \
             patch('aidd_agent.chemical_companion.ChemicalCompanionReader',return_value=SimpleNamespace(get=lambda gid:chem)):
            check=exporter.source_pose_checker(ROOT,ROOT/'configs/targets/kras_sip_e116')
        pose=dict(gid=0,mid='M',cid='C',transform=np.eye(4).ravel().tolist())
        valid,distance=check(molecule,'KRAS_G12D',pose)
        self.assertTrue(valid);self.assertGreater(distance,2.)
        matrix=np.eye(4);matrix[:3,3]=state['geometry']['receptor'][0]-q['shape_points'][0]
        pose['transform']=matrix.ravel().tolist()
        valid,distance=check(molecule,'KRAS_G12D',pose)
        self.assertFalse(valid);self.assertLess(distance,1e-8)

    def test_native_perturbation_recovery(self):
        results=runner.native_checks(ROOT/'configs/targets/kras_sip_e116')
        self.assertEqual(len(results),2)
        for result in results.values():
            self.assertTrue(result['all_contacts'])
            self.assertGreaterEqual(result['minimum_distance'],2.)

    def test_jacobian_and_rigid_geometry(self):
        state=refine.soft.run.kw.load_state(ROOT/'configs/targets/kras_sip_e116','KRAS_G12D')
        q=state['query'];chem=SimpleNamespace(**state['record']['chemistry'])
        sets=refine.soft.run.kw.chemical_sets(chem,q['feature_types'])
        initial=np.eye(4);initial[:3,3]=[.7,-.4,.3]
        objective=refine.Objective(q['shape_points'],initial,state,sets)
        x=np.array([.11,-.08,.06,.22,-.12,.31]);values,jac,matrix,_,_=objective.calculate(x,False)
        numerical=np.empty_like(jac);eps=1e-6
        for i in range(6):
            delta=np.zeros(6);delta[i]=eps
            numerical[:,i]=(objective.calculate(x+delta,False)[0]-objective.calculate(x-delta,False)[0])/(2*eps)
        np.testing.assert_allclose(jac,numerical,atol=2e-5,rtol=2e-5)
        before=q['shape_points'];after=before@matrix[:3,:3].T+matrix[:3,3]
        np.testing.assert_allclose(np.linalg.norm(before[:,None]-before[None,:],axis=2),
                                   np.linalg.norm(after[:,None]-after[None,:],axis=2),atol=1e-10)

    def test_contact_loss_is_not_accepted(self):
        state=refine.soft.run.kw.load_state(ROOT/'configs/targets/kras_sip_e116','KRAS_G12V')
        q=state['query'];chem=SimpleNamespace(**state['record']['chemistry'])
        sets=refine.soft.run.kw.chemical_sets(chem,q['feature_types'])
        initial=np.eye(4);initial[:3,3]=[100,100,100]
        result,evidence=refine.optimize_pose(q['shape_points'],initial,state,sets,max_nfev=2)
        self.assertIsNone(result)
        self.assertGreater(evidence['endpoint_minimum_distance'],2.)

    def test_spread_order_covers_every_conformer_once(self):
        tasks=list(runner.spread_tasks(103,4));ids=[gid for _,values in tasks for gid in values]
        self.assertEqual(sorted(ids),list(range(103)))
        self.assertEqual(len({index for index,_ in tasks}),len(tasks))


if __name__=='__main__':unittest.main()
