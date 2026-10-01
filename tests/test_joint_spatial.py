"""Spatial invariance, real MOL2 provenance, resumability and joint ownership tests."""
import json
from pathlib import Path
import sqlite3

import numpy as np
import pytest
from rdkit import Chem, RDConfig, rdBase

from aidd_agent.joint_spatial_descriptor import geometry,compute,WIDTH
from aidd_agent.macrocycle_descriptors import describe_peptide
from aidd_agent.mol2 import load_rdkit_mol2,parse_mol2_block
from aidd_agent.final_work_blocks import sha
from aidd_agent.joint_spatial_profiles import prepare
from aidd_agent.joint_spatial_blocks import run,scaling,GROUPS,partition


def molecule():
    mol=Chem.MolFromSmiles('N1[C@@H](C)C(=O)N[C@@H](C)C(=O)N[C@@H](C)C1=O')
    conf=Chem.Conformer(mol.GetNumAtoms())
    for i,p in enumerate(np.random.default_rng(81).normal(size=(mol.GetNumAtoms(),3))):conf.SetAtomPosition(i,p)
    mol.AddConformer(conf);return mol


def vector(mol):
    d=describe_peptide(mol,None,'backbone')
    return geometry(mol,d['units'],np.array(d['descriptor']).reshape(-1,6))


def test_proper_rigid_and_cyclic_origin_invariance():
    mol=molecule();a=vector(mol);moved=Chem.Mol(mol)
    rotation=np.array([[0,-1,0],[1,0,0],[0,0,1]])
    for i,p in enumerate(np.array(mol.GetConformer().GetPositions())@rotation+[7,-9,2]):moved.GetConformer().SetAtomPosition(i,p)
    np.testing.assert_allclose(a,vector(moved),rtol=1e-5,atol=1e-6)
    d=describe_peptide(mol,None,'backbone')
    units=np.roll(d['units'],1,axis=0);torsions=np.roll(np.array(d['descriptor']).reshape(-1,6),1,axis=0)
    np.testing.assert_allclose(a,geometry(mol,units,torsions),atol=1e-6)


def test_sidechain_motion_and_reflection_are_visible():
    mol=molecule();a=vector(mol);changed=Chem.Mol(mol)
    atom=next(x.GetIdx() for x in mol.GetAtoms() if x.GetAtomicNum()==6 and x.GetDegree()==1)
    p=np.array(changed.GetConformer().GetAtomPosition(atom));changed.GetConformer().SetAtomPosition(atom,p+[3,2,1])
    b=vector(changed)
    np.testing.assert_allclose(a[:24],b[:24],atol=1e-6)
    assert np.linalg.norm(a[24:]-b[24:])>.1
    reflected=Chem.Mol(mol)
    for i,p in enumerate(np.array(mol.GetConformer().GetPositions())*[-1,1,1]):reflected.GetConformer().SetAtomPosition(i,p)
    assert np.linalg.norm(a-vector(reflected))>.1


def test_atom_renumbering_is_not_a_new_geometry():
    mol=molecule();order=np.random.default_rng(7).permutation(mol.GetNumAtoms()).tolist()
    np.testing.assert_allclose(vector(mol),vector(Chem.RenumberAtoms(mol,order)),rtol=1e-5,atol=1e-6)


def test_same_molecule_cannot_supply_both_fit_and_check():
    x=np.zeros((10000,WIDTH),dtype=np.float32);x[5000:]=1
    leaves,evidence=partition(x,['same-molecule']*len(x))
    assert len(leaves)==1 and evidence[0]['reason']=='insufficient_internal_check'


def test_small_admitted_cycles_do_not_fail_on_repeated_lags():
    mol=Chem.MolFromSmiles('N1[C@@H](C)C(=O)N[C@@H](C)C1=O')
    conf=Chem.Conformer(mol.GetNumAtoms())
    for i,p in enumerate(np.random.default_rng(71).normal(size=(mol.GetNumAtoms(),3))):conf.SetAtomPosition(i,p)
    mol.AddConformer(conf)
    assert np.isfinite(vector(mol)).all()


def mol2(mol,name):
    mol=Chem.AddHs(mol,addCoords=True)
    lines=['@<TRIPOS>MOLECULE',name,f'{mol.GetNumAtoms()} {mol.GetNumBonds()} 0 0 0','SMALL','USER_CHARGES','','@<TRIPOS>ATOM']
    for atom in mol.GetAtoms():
        i=atom.GetIdx();p=mol.GetConformer().GetAtomPosition(i)
        kind={'C':'C.2' if atom.GetHybridization()==Chem.HybridizationType.SP2 else 'C.3','N':'N.am','O':'O.2','H':'H'}[atom.GetSymbol()]
        lines.append(f'{i+1} {atom.GetSymbol()}{i+1} {p.x:.8f} {p.y:.8f} {p.z:.8f} {kind} 1 MOL 0.0')
    lines.append('@<TRIPOS>BOND')
    for b in mol.GetBonds():lines.append(f'{b.GetIdx()+1} {b.GetBeginAtomIdx()+1} {b.GetEndAtomIdx()+1} {int(b.GetBondTypeAsDouble())}')
    return '\n'.join(lines)+'\n'


def fixture(tmp_path,count=4):
    source=tmp_path/'source.mol2';profiles=tmp_path/'old';build=tmp_path/'build';blocks=tmp_path/'backbone'
    for p in (profiles,build,blocks):p.mkdir()
    texts=[mol2(molecule(),f'sample_conf{i+1}') for i in range(count)]
    source.write_text(''.join(texts))
    rows=[]
    for i,text in enumerate(texts):
        record=parse_mol2_block(source,i,text);mol,mode=load_rdkit_mol2(text,record.name);assert mode=='strict'
        d=describe_peptide(mol,None,'backbone')
        row=dict(conformer_id=f'cid-{i}',molecule_id='same-molecule',hard_group=d['hard_group'],descriptor=d['descriptor'],
            provenance=dict(content_sha256=record.content_sha256,source_record_name=record.name,source_name=record.molecule_name,
                source_path=str(source),source_record_index=i,flags=['cross_conformer_alignment_unverified'],units=d['units']))
        rows.append(row)
    with sqlite3.connect(build/'descriptors.sqlite') as db:
        db.execute('CREATE TABLE descriptor(cid TEXT PRIMARY KEY,payload TEXT)')
        db.executemany('INSERT INTO descriptor VALUES(?,?)',[(r['conformer_id'],json.dumps(r)) for r in rows])
    with sqlite3.connect(blocks/'work_blocks.sqlite') as db:
        db.execute('CREATE TABLE block(block_id TEXT,n INTEGER,kind TEXT)')
        db.executemany('INSERT INTO block VALUES(?,?,?)',[('parent',count,'single_class'),('special-exhaustive',1,'special_exhaustive')])
    with sqlite3.connect(profiles/'profiles.sqlite') as db:
        db.execute('CREATE TABLE files(id INTEGER PRIMARY KEY,path TEXT)');db.execute('INSERT INTO files VALUES(1,?)',(str(source),))
        db.execute('CREATE TABLE item(cid TEXT PRIMARY KEY,parent TEXT,file_id INTEGER,idx INTEGER,vector BLOB)')
        db.executemany('INSERT INTO item VALUES(?,?,?,?,?)',[(r['conformer_id'],'parent',1,i,np.zeros(69,dtype='<f4').tobytes()) for i,r in enumerate(rows)])
    br=dict(status='complete',source_conformers=count+1,regular_conformers=count,special_conformers=1,
            output_hashes={'work_blocks.sqlite':sha(blocks/'work_blocks.sqlite')})
    (blocks/'report.json').write_text(json.dumps(br))
    sr=dict(status='complete',schema=dict(rdkit=rdBase.rdkitVersion,variant='backbone',feature_definition_sha256=sha(Path(RDConfig.RDDataDir)/'BaseFeatures.fdef')),
            output_hashes={'descriptors.sqlite':sha(build/'descriptors.sqlite')})
    (build/'report.json').write_text(json.dumps(sr))
    pr=dict(status='complete',width=69,coverage=1,total_regular=count,blocks=str(blocks),
        signature=dict(blocks_report_sha256=sha(blocks/'report.json'),build_report_sha256=sha(build/'report.json')),
        output_hashes={'profiles.sqlite':sha(profiles/'profiles.sqlite')})
    (profiles/'report.json').write_text(json.dumps(pr))
    receipt=tmp_path/'e094-validation.json';receipt.write_text(json.dumps(dict(structural_gate='passed',property_gate='passed',conformers=count+1,regular_conformers=count)))
    return profiles,build,receipt,rows,texts


def test_real_mol2_extraction_joint_build_and_resume(tmp_path):
    profiles,build,receipt,rows,texts=fixture(tmp_path)
    output=tmp_path/'spatial';report=prepare(profiles,build,output,workers=1)
    assert report['width']==WIDTH and report['total_regular']==4
    assert prepare(profiles,build,output,workers=1,resume=True)==report
    with sqlite3.connect(output/'profiles.sqlite') as db:
        assert db.execute('SELECT count(DISTINCT mid) FROM item').fetchone()[0]==1
        assert db.execute('SELECT min(length(vector)) FROM item').fetchone()[0]==WIDTH*4
    final=run(output,tmp_path/'joint',receipt)
    assert final['regular_conformers']==4 and final['special_conformers']==1
    assert final['regular_work_blocks']==1 and final['check_grouping']=='source_molecule_id_crc32_mod5'
    assert final['output_hashes']['property_blocks.sqlite']==run(output,tmp_path/'joint',receipt,resume=True)['output_hashes']['property_blocks.sqlite']


def test_source_identity_failure_is_not_accepted(tmp_path):
    profiles,build,receipt,rows,texts=fixture(tmp_path)
    row=rows[0];row['provenance']['content_sha256']='0'*64
    with pytest.raises(ValueError,match='identity/hash'):
        compute((str(tmp_path/'source.mol2'),0,texts[0],json.dumps(row),'backbone',np.zeros(69,dtype='<f4').tobytes()))


def test_partial_extraction_resume(tmp_path,monkeypatch):
    import aidd_agent.joint_spatial_profiles as module
    profiles,build,receipt,rows,texts=fixture(tmp_path)
    monkeypatch.setattr(module,'BATCH_SIZE',2);original=module.compute;calls=0
    def fail_after_commit(task):
        nonlocal calls
        calls+=1
        if calls==3:raise RuntimeError('Injected interruption')
        return original(task)
    monkeypatch.setattr(module,'compute',fail_after_commit)
    with pytest.raises(RuntimeError):prepare(profiles,build,tmp_path/'spatial',workers=1)
    assert not (tmp_path/'spatial/report.json').exists()
    with sqlite3.connect(tmp_path/'spatial/profiles.sqlite') as db:assert db.execute('SELECT count(*) FROM item WHERE vector IS NOT NULL').fetchone()[0]==2
    monkeypatch.setattr(module,'compute',original)
    assert prepare(profiles,build,tmp_path/'spatial',workers=1,resume=True)['total_regular']==4


def test_global_group_weights_have_fixed_budgets():
    with sqlite3.connect(':memory:') as db:
        db.execute('CREATE TABLE item(vector BLOB)')
        db.executemany('INSERT INTO item VALUES(?)',[(np.full(WIDTH,v,dtype='<f4').tobytes(),) for v in (-1,1)])
        transform=scaling(db,2)
    weights=np.array(transform['weights'])
    for group,budget in GROUPS:assert np.sum(weights[group]**2)==pytest.approx(budget)


def test_parallel_workers_match_serial_vectors(tmp_path):
    profiles,build,receipt,rows,texts=fixture(tmp_path)
    prepare(profiles,build,tmp_path/'serial',workers=1)
    prepare(profiles,build,tmp_path/'parallel',workers=2)
    with sqlite3.connect(tmp_path/'serial/profiles.sqlite') as a,sqlite3.connect(tmp_path/'parallel/profiles.sqlite') as b:
        assert a.execute('SELECT cid,mid,vector FROM item ORDER BY cid').fetchall()==b.execute('SELECT cid,mid,vector FROM item ORDER BY cid').fetchall()
