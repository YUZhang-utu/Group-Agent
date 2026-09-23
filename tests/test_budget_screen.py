import json
import sqlite3
from types import SimpleNamespace

import numpy as np
import pytest

from aidd_agent import budget_screen as screen,preselection_full as pre
from aidd_agent.budget_export import hydrogen_roundtrip,transformed_mol2
from test_library_acceptance import library
from test_preselection_full import geometry


def test_unique_molecule_retrieval_and_deterministic_fusion():
    rows=screen.unique_neighbors([4,3,2,1],[b'a',b'a',b'b',b'c'],[.3,.2,.2,.9],3)
    assert rows==[('a',.2,3),('b',.2,2),('c',.9,1)]
    assert screen.fuse([['a','b'],['b','c']],3)==['b','a','c']
    with pytest.raises(ValueError):screen.fuse([['a','a']])


def test_expand_molecules_across_shards(library):
    batch,_,mids=library
    keys=np.unique(mids)[[0,2]]
    assert screen.collect_conformers(batch,keys).tolist()==np.flatnonzero(np.isin(mids,keys)).tolist()


def test_budget_retains_zero_contacts_and_below_old_threshold():
    query,features,seeds=geometry();features.feature_points+=100
    design=dict(mandatory_anchors=['A'],alternative_groups=[],optional_weights={'A':1,'B':1},
                gaussian_weight=.7,optional_weight=.3,minimum_pose_score=.9)
    args=(features,seeds,np.ones(2,bool),query,[0,1],.5,design,['A','B'],np.array([.1,.2]),query['shape_points'])
    assert pre.pose_representatives(*args)[0]==[]
    rows,count=pre.pose_representatives(*args,budget_mode=True)
    assert len(rows)==1 and count==2
    assert rows[0]['mask']==0 and rows[0]['seed_index']==1
    assert rows[0]['legacy_score_threshold_passed'] is False


def test_budget_compute_does_not_call_empirical_filters(tmp_path,monkeypatch):
    query,features,seeds=geometry();candidate=SimpleNamespace(**vars(features),shape_points=query['shape_points'])
    q=dict(selection_mode='budget',condition_policy=dict(required_anchors=['A','B'],minimum_score=.5),
        anchors=[dict(anchor_id='A',score_column=0),dict(anchor_id='B',score_column=1)])
    monkeypatch.setattr(pre.full,'_STATE',(SimpleNamespace(get=lambda _:candidate),query,query,q))
    monkeypatch.setattr(pre.full,'_FILTER_READER',SimpleNamespace(get=lambda _:features))
    def fail(*args,**kwargs):raise AssertionError('Empirical filter called')
    monkeypatch.setattr(pre,'coarse_stages',fail);monkeypatch.setattr(pre,'possible_seed_mask',fail)
    monkeypatch.setattr(pre.full,'_BOUND',SimpleNamespace(check=fail))
    monkeypatch.setattr(pre,'prepare_seeds',lambda *a,**kw:(seeds,2))
    monkeypatch.setattr(pre,'_score_ids',lambda *a,**kw:{pre.full.OBJECTIVE+'__objective':np.array([.1])})
    path=tmp_path/'chunk.npz';pre.compute((0,1,str(path),[7]))
    assert screen.verified_chunk(path,np.array([7]))
    with pytest.raises(ValueError,match='schedule'):screen.verified_chunk(path,np.array([8]))
    path.with_suffix('.poses.jsonl').write_text('corrupted')
    with pytest.raises(ValueError):screen.verified_chunk(path,np.array([7]))


def test_hydrogen_reconstruction_and_original_pose_transform():
    from rdkit import Chem
    for smiles in ('c1cc[nH]c1','[NH4+]','c1ccccc1'):
        audit=hydrogen_roundtrip(Chem.AddHs(Chem.MolFromSmiles(smiles)))
        assert audit['status']=='ok'
        assert any(a['explicit_h_neighbors'] for a in audit['atoms'])
    raw='@<TRIPOS>MOLECULE\noriginal name\n1 0\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n1 N 0 1 2 N.ar 1 MOL 0\n@<TRIPOS>BOND\n'
    matrix=np.eye(4);matrix[:3,3]=[2,3,4]
    moved=transformed_mol2(raw,matrix)
    assert 'original name' in moved and '2.000000 4.000000 6.000000 N.ar' in moved


def test_ann_grows_depth_until_unique_molecule_budget(library,tmp_path,monkeypatch):
    import sys
    batch,raw,mids=library;calls=[]
    class Index:
        ntotal=12;d=60;nlist=4
        def search(self,q,k):
            calls.append(self.nprobe)
            ids=np.arange(12) if self.nprobe==4 else np.array([0,1])
            return np.zeros((1,len(ids))),ids[None]
    monkeypatch.setitem(sys.modules,'faiss',SimpleNamespace(omp_set_num_threads=lambda _:None,read_index=lambda _:Index()))
    monkeypatch.setattr(screen,'query_vector',lambda _:raw[0])
    out=tmp_path/'retrieval';out.mkdir()
    r=screen.retrieve(batch,out,dict(templates=[dict(query_id='T')]),5,1,1)
    assert r['selected_molecules']==5 and calls==[1,2,4]
    assert len(np.unique(np.load(out/'selected-molecules.npy')))==5


def test_rank_merge_deduplicates_conformers_and_keeps_reserve(tmp_path,monkeypatch):
    from aidd_agent import gaussian_batch
    from aidd_agent.expanded_wee1 import fingerprint
    monkeypatch.setattr(gaussian_batch,'ArtifactCatalogReader',lambda _:SimpleNamespace(get=lambda gid:SimpleNamespace(molecule_id={1:'a',2:'a',3:'b',4:'c'}[gid],shape_points=np.array([[0.,0,0]]))))
    design=dict(templates=[dict(query_id='T1'),dict(query_id='T2')])
    for ti,values in enumerate([[('a',.1,1),('a',.8,2),('b',.7,3),('c',.2,4)],[('b',.9,3),('c',.5,4)]]):
        folder=tmp_path/'chunks'/f'{ti:02d}';folder.mkdir(parents=True)
        path=folder/'000.poses.jsonl'
        path.write_text(''.join(json.dumps(dict(molecule_id=mid,composite_score=score,global_id=gid,transform=np.eye(4).reshape(-1).tolist()))+'\n' for mid,score,gid in values))
        path.with_name('000.receipt.json').write_text(json.dumps(dict(files=fingerprint([path]))))
    region=[dict(id='A',points=[[0,0,0]],radius=1.,minimum_atoms=1)]
    definitions=dict(definitions=dict(one=region,two=region),ambiguity_margin=.5)
    assert screen.merge_and_rank(tmp_path,design,definitions,tmp_path,1,60)==3
    with sqlite3.connect(tmp_path/'ranking.sqlite') as db:
        assert db.execute('SELECT mid FROM ranking ORDER BY rank').fetchall()==[('b',),('a',),('c',)]
        assert db.execute("SELECT gid FROM poses WHERE mid='a'").fetchone()==(2,)
        assert json.loads(db.execute('SELECT payload FROM poses LIMIT 1').fetchone()[0])['region_occupancy']['agreed_regions']==1


def test_export_pages_original_names_hashes_and_resume(tmp_path,monkeypatch):
    import hashlib
    from aidd_agent import chemical_companion as companion
    from aidd_agent.budget_export import export
    from aidd_agent.expanded_wee1 import fingerprint
    from aidd_agent.mol2 import load_rdkit_mol2,iter_mol2_blocks
    from test_mol2_rdkit import _block
    batch=tmp_path/'batch';batch.mkdir();run=tmp_path/'run';run.mkdir()
    source=tmp_path/'source.mol2'
    raw=_block('1 C 0 0 0 C.3 1 M 0\n2 N 1 0 0 N.3 1 M 0\n','1 1 2 1\n',2,1)
    source.write_text(raw+raw.replace('\ntest\n','\nsecond\n'),newline='')
    blocks=list(iter_mol2_blocks(source));mol,_=load_rdkit_mol2(raw,'fixture')
    def get(gid):
        return SimpleNamespace(molecule_id=f'M{gid}',conformer_id=f'C{gid}',
            atomic_numbers=np.array([a.GetAtomicNum() for a in mol.GetAtoms()]),
            formal_charges=np.array([a.GetFormalCharge() for a in mol.GetAtoms()]),atom_flags=np.zeros(2,dtype=np.uint8),
            bonds=np.array([(0,1,1,0,0)],dtype=companion.BOND_DTYPE))
    monkeypatch.setattr(companion,'ChemicalCompanionReader',lambda _:SimpleNamespace(get=get))
    with sqlite3.connect(batch/'registry.sqlite3') as db:
        db.execute('CREATE TABLE molecule(id,source_name)')
        db.execute('CREATE TABLE conformer(id,molecule_id,source_path,source_record_index,content_sha256,source_record_name)')
        for i,(index,text) in enumerate(blocks):
            db.execute('INSERT INTO molecule VALUES(?,?)',(f'M{i}','original name '+str(i)))
            db.execute('INSERT INTO conformer VALUES(?,?,?,?,?,?)',(f'C{i}',f'M{i}',str(source),index,hashlib.sha256(text.encode()).hexdigest(),'record '+str(i)))
    with sqlite3.connect(run/'ranking.sqlite') as db:
        db.execute('CREATE TABLE ranking(rank,mid,rrf,template_support)')
        db.execute('CREATE TABLE poses(mid,template,payload)');db.execute('CREATE TABLE template_ranks(mid,template,rank)')
        for i in range(2):
            db.execute('INSERT INTO ranking VALUES(?,?,?,?)',(i+1,f'M{i}',.1,1))
            db.execute('INSERT INTO template_ranks VALUES(?,?,?)',(f'M{i}','T',i+1))
            db.execute('INSERT INTO poses VALUES(?,?,?)',(f'M{i}','T',json.dumps(dict(global_id=i,transform=np.eye(4).reshape(-1).tolist()))))
    (run/'report.json').write_text(json.dumps(dict(status='complete',outputs=fingerprint([run/'ranking.sqlite']))))
    first=export(batch,run,tmp_path/'page1',1,1);second=export(batch,run,tmp_path/'page2',2,1)
    assert first['exported_molecules']==second['exported_molecules']==1
    assert (tmp_path/'page1/original/M0.mol2').read_text()==blocks[0][1]
    assert not (tmp_path/'page2/original/M0.mol2').exists()
    assert export(batch,run,tmp_path/'page1',1,1)==first
    assert 'original name 0' in (tmp_path/'page1/molecules.csv').read_text(encoding='utf-8-sig')
    assert export(batch,run,tmp_path/'page3',3,1)['shortfall']==1
