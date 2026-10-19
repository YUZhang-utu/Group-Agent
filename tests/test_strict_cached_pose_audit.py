"""Regression checks for strict cached-pose audit semantics and I/O."""
import importlib.util
import json
from pathlib import Path
import sqlite3
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch
import numpy as np
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('strict_audit', ROOT/'scripts/audit_strict_cached_poses.py')
a = importlib.util.module_from_spec(spec)
spec.loader.exec_module(a)


class StrictAuditTests(unittest.TestCase):
    def test_distinct_residues_not_atom_count(self):
        order = [a.ILE61, 'leu1', 'leu2']
        anchors = [dict(anchor_id=k, target_residue=r) for k,r in zip(order,[61,54,54])]
        self.assertFalse(a.contact_gate([1,1,1],[0,1,2],order,anchors)[0])
        anchors[-1]['target_residue'] = 93
        self.assertTrue(a.contact_gate([0.5,1,1],[0,1,2],order,anchors)[0])
        self.assertFalse(a.contact_gate([0.49,1,1],[0,1,2],order,anchors)[0])
        self.assertFalse(a.contact_gate([1,1,1],[-1,1,2],order,anchors)[0])

    def test_no_cross_pose_compensation(self):
        metric = dict(minimum_distance=3, envelope_coverage=1, regions_supported=True,
                      region_coverages=[1,1,1])
        contact_only = a.gates(True, 0.1, 0.9, metric)
        shape_only = a.gates(False, 0.9, 0.9, metric)
        self.assertFalse(a.accepted(contact_only))
        self.assertFalse(a.accepted(shape_only))
        self.assertFalse(a.accepted(a.gates(True,0.9,0.29,metric)))
        self.assertTrue(a.accepted(a.gates(True,0.5,0.3,metric)))

    def test_regions_and_sterics(self):
        reference = np.array([[x,y,0] for x in [-5,0,5] for y in [-0.2,0,0.2]])
        tree = cKDTree(reference + [0,0,3])
        contact, regions = a.make_regions(reference,tree)
        self.assertEqual([len(r) for r in regions], [3,3,3])
        good = a.geometry_metrics(reference,tree,contact,regions)
        self.assertTrue(a.accepted(a.gates(True,1,1,good)))
        partial = a.geometry_metrics(reference[:6],tree,contact,regions)
        self.assertFalse(a.gates(True,1,1,partial)['regions'])
        collision = a.geometry_metrics(reference+[0,0,2],tree,contact,regions)
        self.assertFalse(a.gates(True,1,1,collision)['steric'])
        self.assertIsNone(a.make_regions(reference[:3],tree)[1])

    def test_invalid_coordinates(self):
        for value in [[], [[0,0,float('nan')]], [[0,0]]]:
            with self.assertRaises(ValueError):
                a.points(value)

    def test_end_to_end_includes_original_failures_and_seals_outputs(self):
        from aidd_agent import gaussian_batch, block_sampling
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            prep, ref, batch = root/'prepared',root/'reference',root/'batch'
            for directory in [prep,ref,batch]: directory.mkdir()
            xyz = np.array([[x,y,0] for x in [-5,0,5] for y in [-0.2,0,0.2]])
            data = dict(shape_points=xyz, feature_points=xyz, feature_types=np.full(9,5))
            template, receptor = root/'template.npz',root/'receptor.npz'
            np.savez(template,**data)
            np.savez(receptor,points=xyz+[0,0,3])
            order = sorted([a.ILE61,'other54','other93'])
            anchors = [dict(anchor_id=k,target_residue=61 if k==a.ILE61 else int(k[-2:])) for k in order]
            design = dict(mandatory_anchors=[],alternative_groups=[],optional_weights={k:1 for k in order},
                          receptor_npz=str(receptor))
            query = dict(design=design,anchor_order=order,anchors=anchors,
                         templates=[dict(query_id='T',query_npz=str(template))],
                         sources={str(template):a.sha(template),str(receptor):a.sha(receptor)},
                         crystal_self_controls=[dict(query_id='T',scores=[1,1,1])])
            a.save(prep/'query.json',query)
            code = {n:a.sha(ROOT/'src/aidd_agent'/n) for n in ['gaussian_overlay.py','gaussian_batch.py','contact_groups.py']}
            a.save(prep/'report.json',dict(status='complete',batch=str(batch),query_sha256=a.sha(prep/'query.json'),
                                         registry_snapshot='fixture',code=code))
            with sqlite3.connect(ref/'schedule.sqlite') as db:
                db.execute('CREATE TABLE selected(gid INTEGER PRIMARY KEY,cid TEXT,mid TEXT)')
                db.executemany('INSERT INTO selected VALUES(?,?,?)',[(1,'c1','m1'),(2,'c2','m2')])
            db.close()
            pose = dict(global_id=1,conformer_id='c1',molecule_id='m1',transform=np.eye(4).ravel().tolist(),
                        scores=[1,1,1],assignments=[0,1,2],gaussian_same_pose=1)
            with sqlite3.connect(ref/'scores.sqlite') as db:
                db.execute('CREATE TABLE poses(gid INTEGER,template TEXT,payload TEXT,PRIMARY KEY(gid,template))')
                db.execute('INSERT INTO poses VALUES(?,?,?)',(1,'T',json.dumps(pose)))
            db.close()
            a.save(ref/'report.json',dict(status='complete',mode='full',include_training=False,
                                         prepared_report_sha256=a.sha(prep/'report.json'),query_sha256=a.sha(prep/'query.json'),
                                         output_hashes={n:a.sha(ref/n) for n in ['schedule.sqlite','scores.sqlite']}))
            candidate=SimpleNamespace(**data,conformer_id='c1',molecule_id='m1')
            class Reader:
                def __init__(self,path): pass
                def get(self,gid): return candidate
            with patch.object(gaussian_batch,'ArtifactCatalogReader',Reader), patch.object(gaussian_batch,'_load_query',return_value=({},data)), patch.object(block_sampling,'snapshot',return_value='fixture'):
                result=a.audit(SimpleNamespace(prepared=prep,reference=ref,output=root/'out',sample=2,seed=1,seconds=60))
            self.assertEqual(result['strict_pass_fraction'],0.5)
            self.assertEqual(result['old_pass_conformers'],1)
            self.assertEqual(result['template_controls_passed'],1)
            for name,digest in result['output_hashes'].items(): self.assertEqual(a.sha(root/'out'/name),digest)
            with patch.object(gaussian_batch,'ArtifactCatalogReader',Reader), patch.object(gaussian_batch,'_load_query',return_value=({},data)), patch.object(block_sampling,'snapshot',return_value='fixture'), patch.object(a.time,'monotonic',side_effect=[0,0,0,0,100,100,100]):
                partial=a.audit(SimpleNamespace(prepared=prep,reference=ref,output=root/'partial',sample=2,seed=1,seconds=60))
            self.assertEqual(partial['status'],'time_budget_exhausted')
            self.assertEqual(partial['processed_conformers'],0)
            self.assertIsNone(partial['strict_pass_fraction'])
            with self.assertRaises(FileExistsError):
                a.audit(SimpleNamespace(prepared=prep,reference=ref,output=root/'out',sample=2,seed=1,seconds=60))


if __name__ == '__main__':
    unittest.main()
