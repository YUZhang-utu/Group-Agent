from types import SimpleNamespace
import numpy as np
import pytest

from aidd_agent.chemical_companion import BOND_DTYPE
from aidd_agent.macrocycle_blocks import describe,partition,amide_state


def example():
    n=18;t=np.arange(n)*2*np.pi/n
    xyz=np.column_stack((3*np.cos(t),3*np.sin(t),.6*np.sin(3*t)+.1*np.cos(2*t)))
    chem=SimpleNamespace(atomic_numbers=np.array([7,6,6]*6),formal_charges=np.zeros(n,int),
        bonds=np.array([(i,(i+1)%n,1,0,0) for i in range(n)],dtype=BOND_DTYPE))
    return chem,xyz


def test_rigid_transform_and_atom_reindexing():
    chem,xyz=example();a=describe(chem,xyz)
    rotation=np.array([[0,-1,0],[1,0,0],[0,0,1]])
    b=describe(chem,xyz@rotation+np.array([9,3,-2]))
    assert a['hard_group']==b['hard_group']
    np.testing.assert_allclose(a['descriptor'],b['descriptor'],atol=1e-10)
    perm=np.random.default_rng(4).permutation(18);inverse=np.argsort(perm)
    bonds=chem.bonds.copy();bonds['begin']=inverse[bonds['begin']];bonds['end']=inverse[bonds['end']]
    reordered=SimpleNamespace(atomic_numbers=chem.atomic_numbers[perm],formal_charges=chem.formal_charges[perm],bonds=bonds)
    c=describe(reordered,xyz[perm])
    assert a['hard_group']==c['hard_group']
    np.testing.assert_allclose(a['descriptor'],c['descriptor'],atol=1e-10)


def test_ambiguous_ring_and_explicit_map():
    chem,xyz=example();chem.bonds=np.concatenate((chem.bonds,np.array([(0,9,1,0,0)],dtype=BOND_DTYPE)))
    with pytest.raises(ValueError):describe(chem,xyz)
    assert describe(chem,xyz,list(range(18)))['ring_mapping']=='supplied'
    with pytest.raises(ValueError):describe(chem,xyz,[0]*18)


def test_hard_groups_capacity_and_determinism():
    rows=[dict(global_id=i,molecule_id='M'+str(i//3),conformer_id='C'+str(i),
        hard_group='cis' if i%2 else 'trans',descriptor=[np.sin(i),np.cos(i)]) for i in range(30)]
    members,blocks=partition(rows,4)
    assert len(members)==30 and all(b['conformers']<=4 for b in blocks)
    assert members==partition(list(reversed(rows)),4)[0]
    groups={r['conformer_id']:r['hard_group'] for r in rows}
    for block in blocks:
        assert {groups[r['conformer_id']] for r in members if r['block_id']==block['block_id']}=={block['hard_group']}
    assert len({r['molecule_id'] for r in members})==10


def test_cis_trans_boundary_and_invalid_geometry():
    assert [amide_state(a) for a in [0,30,31,149,150,-179]]==['cis','cis','boundary','boundary','trans','trans']
    chem,xyz=example();xyz[0]=np.nan
    with pytest.raises(ValueError):describe(chem,xyz)


def test_catalog_pilot_preserves_registry_identity(tmp_path,monkeypatch):
    import json
    import sqlite3
    import hashlib
    import aidd_agent.gaussian_batch as geometry
    import aidd_agent.chemical_companion as chemistry
    from aidd_agent.macrocycle_blocks import run
    batch=tmp_path/'batch';(batch/'artifacts').mkdir(parents=True);(batch/'chemical').mkdir()
    catalog=batch/'artifacts/catalog.json';catalog.write_text(json.dumps(dict(library_id='L',conformers=3)))
    (batch/'chemical/catalog.json').write_text(json.dumps(dict(library_id='L',artifact_v1_catalog_sha256=hashlib.sha256(catalog.read_bytes()).hexdigest())))
    chem,xyz=example()
    class Shape:
        def __init__(self,path):self.catalog=json.loads(path.read_text())
        def get(self,gid):return SimpleNamespace(molecule_id='M',conformer_id='C'+str(gid),shape_points=xyz)
    class Chemical:
        def __init__(self,path):pass
        def get(self,gid):return SimpleNamespace(**vars(chem),molecule_id='M',conformer_id='C'+str(gid))
    monkeypatch.setattr(geometry,'ArtifactCatalogReader',Shape);monkeypatch.setattr(chemistry,'ChemicalCompanionReader',Chemical)
    with sqlite3.connect(batch/'registry.sqlite3') as db:
        db.executescript('CREATE TABLE molecule(id,source_name,library_id); CREATE TABLE conformer(id,molecule_id,source_record_name,source_path,source_record_index,content_sha256,topology_sha256,conformer_index);')
        db.execute('INSERT INTO molecule VALUES(?,?,?)',('M','original-cis','L'))
        for i in range(3):db.execute('INSERT INTO conformer VALUES(?,?,?,?,?,?,?,?)',('C'+str(i),'M','original-cis_conf'+str(i+1),'source.mol2',i,'sha'+str(i),'topology',i+1))
    output=tmp_path/'output';report=run(batch,output,3,2)
    assert report['assigned_conformers']==3
    rows=[json.loads(line) for line in (output/'descriptors.jsonl').read_text().splitlines()]
    assert {r['conformer_id'] for r in rows}=={'C0','C1','C2'}
    assert {r['source_name'] for r in rows}=={'original-cis'}
    with pytest.raises(ValueError):run(batch,output,3,2)
