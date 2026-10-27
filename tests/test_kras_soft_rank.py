"""Run with unittest or pytest; includes native coordinates and source integrity."""
import json
from pathlib import Path
import sqlite3
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import kras_soft_rank as score
import export_kras_soft_rank as exporter


class SoftRankingTests(unittest.TestCase):
    def test_cli_pool_and_full_resume(self):
        from test_kras_weekend import tiny_library
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); batch = tiny_library(root/'library')
            prepared=ROOT/'configs/targets/kras_sip_e116'; direct=root/'direct'; direct.mkdir()
            original=score.run.seal_inputs(batch,prepared,128,0)
            score.run.initialize(batch,prepared)
            reader,features,states,chemistry=score.run.WORK
            rows=[]; totals={t:dict(direct_poses=0,distinct_molecules=0) for t in exporter.TARGETS}
            for gid in range(4):
                candidate=reader.get(gid)
                for target,state in states.items():
                    result=score.run.kw.evaluate(candidate,features.get(gid),state,chemistry.get(gid))
                    n=result['direct_region_pose_counts'].get('all',0)
                    totals[target]['direct_poses']+=n
                    if n: rows.append((target,candidate.molecule_id))
            connection=sqlite3.connect(direct/'counts.sqlite')
            connection.execute('CREATE TABLE molecules(target TEXT,mid TEXT,PRIMARY KEY(target,mid))')
            connection.executemany('INSERT OR IGNORE INTO molecules VALUES (?,?)',rows)
            for target in totals:
                totals[target]['distinct_molecules']=connection.execute('SELECT COUNT(*) FROM molecules WHERE target=?',(target,)).fetchone()[0]
            connection.commit();connection.close()
            exporter.save(direct/'protocol.json',dict(original=original))
            exporter.save(direct/'report.json',dict(status='complete',full_catalog=True,
                scope='direct_contacts_only_not_final_screening_hits',targets=totals))
            command=[sys.executable,str(ROOT/'scripts/kras_soft_rank.py'),'--direct',str(direct),
                     '--workers','2','--chunk-size','2']
            subprocess.run(command+['--output',str(root/'pilot'),'--sample','4'],check=True,capture_output=True,text=True)
            full=command+['--output',str(root/'full'),'--sample','0','--pilot',str(root/'pilot/report.json')]
            subprocess.run(full,check=True,capture_output=True,text=True)
            report=exporter.read(root/'full/report.json')
            self.assertEqual(report['status'],'complete')
            self.assertTrue(all(report['direct_pose_reconciliation'].values()))
            subprocess.run(full,check=True,capture_output=True,text=True)
            self.assertEqual(exporter.read(root/'full/report.json')['invocation_new_conformers'],0)
            # Release Windows fixture memory maps before deleting the directory.
            score.run.WORK=None
            del reader,features,states,chemistry,candidate
            import gc
            gc.collect()

    def test_native_controls_and_collision(self):
        controls = score.controls(ROOT/'configs/targets/kras_sip_e116')
        self.assertEqual(set(controls), set(exporter.TARGETS))
        for target in exporter.TARGETS:
            self.assertGreater(controls[target]['pose']['shape'], .99)
            state = score.run.kw.load_state(ROOT/'configs/targets/kras_sip_e116', target)
            q = state['query']
            candidate = SimpleNamespace(global_id=0, molecule_id='M', conformer_id='C',
                shape_points=q['shape_points'], feature_points=q['feature_points'], feature_types=q['feature_types'])
            features = SimpleNamespace(feature_points=q['feature_points'], feature_types=q['feature_types'],
                feature_directions=q['feature_directions'], feature_kinds=q['feature_direction_kinds'])
            # Add all ligand atoms as receptor atoms: every pose passing through
            # the reference position must now fail the unchanged distance floor.
            state['tree'] = score.run.kw.cKDTree(q['shape_points'])
            identity = SimpleNamespace(transform_matrix=score.np.eye(4))
            with patch.object(score.run.kw, 'prepare_seeds', return_value=([identity], None)):
                pose, counts = score.score_candidate(candidate, features, SimpleNamespace(**state['record']['chemistry']), state)
            self.assertIsNone(pose)
            self.assertGreater(counts['direct_poses'], 0)
            self.assertEqual(counts['geometry_poses'], 0)

    def test_atomic_resume_and_deterministic_ties(self):
        db = sqlite3.connect(':memory:'); score.create_tables(db)
        pose = dict(mid='M', cid='C1', gid=1, score=.5, bucket=1)
        output = {t: dict(rows=[pose], seeds=1) for t in exporter.TARGETS}
        score.commit_chunk(db, (0, 1, output))
        with self.assertRaises(ValueError): score.commit_chunk(db, (0, 1, output))
        for i in (4, 3, 2, 0):
            new = dict(pose, cid='C'+str(i), gid=i)
            score.commit_chunk(db, (i+1, 1, {t:dict(rows=[new],seeds=1) for t in exporter.TARGETS}))
        poses = json.loads(db.execute('SELECT payload FROM molecules LIMIT 1').fetchone()[0])
        self.assertEqual([p['gid'] for p in poses], [0,1,2])
        self.assertEqual(db.execute('SELECT COUNT(*) FROM molecules').fetchone()[0], 2)
        db.close()

    def test_diversity_selection(self):
        db = sqlite3.connect(':memory:'); score.create_tables(db)
        for i in range(12):
            db.execute('INSERT INTO molecules VALUES (?,?,?,?,?)', ('KRAS_G12D',str(i),12-i,0 if i<10 else i,'[]'))
        rows = list(exporter.selected_rows(db,'KRAS_G12D',10))
        self.assertEqual(len(rows), 10)
        self.assertEqual(len({r[0] for r in rows}), 10)
        self.assertEqual([r[0] for r in rows[:8]], list(map(str,range(8))))
        self.assertIn('10', [r[0] for r in rows])
        self.assertNotIn('9', [r[0] for r in rows])
        db.close()

    def test_export_batches_and_source_hash(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); batch=root/'batch'; batch.mkdir(); ranking=root/'ranking'; ranking.mkdir()
            source=root/'source.mol2'
            raw='@<TRIPOS>MOLECULE\nwater\n3 2 1 0 0\nSMALL\nUSER_CHARGES\n\n@<TRIPOS>ATOM\n1 O 0 0 0 O.3 1 WAT -0.8\n2 H1 0.96 0 0 H 1 WAT 0.4\n3 H2 -0.24 0.93 0 H 1 WAT 0.4\n@<TRIPOS>BOND\n1 1 2 1\n2 1 3 1\n'
            source.write_text(raw*5,encoding='utf-8')
            registry=sqlite3.connect(batch/'registry.sqlite3')
            registry.execute('CREATE TABLE conformer(id TEXT,molecule_id TEXT,source_path TEXT,source_record_index INTEGER,content_sha256 TEXT)')
            db=sqlite3.connect(ranking/'ranking.sqlite'); score.create_tables(db)
            digest=exporter.hashlib.sha256(raw.encode()).hexdigest()
            for i in range(5):
                mid,cid=f'M{i}',f'C{i}'
                registry.execute('INSERT INTO conformer VALUES(?,?,?,?,?)',(cid,mid,str(source),i,digest))
                pose=dict(mid=mid,cid=cid,gid=i,score=.9-i*.1,shape=.2,color=.1,feature_match=.3,
                          envelope=.4,minimum_distance=2.1,transform=score.np.eye(4).ravel().tolist())
                for target in exporter.TARGETS:
                    db.execute('INSERT INTO molecules VALUES(?,?,?,?,?)',(target,mid,pose['score'],i,json.dumps([pose])))
            registry.commit();registry.close();db.commit();db.close()
            exporter.save(ranking/'protocol.json',dict(common=dict(original=dict(batch=str(batch),prepared=str(root/'prepared')))))
            exporter.save(ranking/'report.json',dict(status='complete',full_catalog=True,
                ranking_sha256=exporter.sha(ranking/'ranking.sqlite'),
                targets={t:dict(eligible_molecules=5) for t in exporter.TARGETS}))
            result=exporter.export(ranking,root/'export',count=5,batch_size=2)
            self.assertTrue(result['quota_met'])
            for target in exporter.TARGETS:
                self.assertEqual([b['molecules'] for b in result['targets'][target]['batches']],[2,2,1])
                records=list(exporter.iter_mol2_blocks(root/'export'/target/'batch_001.mol2'))
                self.assertEqual(records[0][1].splitlines()[2:],raw.splitlines()[2:])
            source.write_text(raw.replace('-0.8','-0.7')*5,encoding='utf-8')
            with self.assertRaisesRegex(ValueError,'hash mismatch'):
                exporter.export(ranking,root/'bad-export',count=5,batch_size=2)


if __name__ == '__main__':
    unittest.main()
