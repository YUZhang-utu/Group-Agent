import json
import sqlite3
from types import SimpleNamespace

import numpy as np
import pytest

from aidd_agent.interaction_matching import _maximum_weight_assignment, _assignment_reference
from aidd_agent.pilot_review import summarize, compare_topology
from aidd_agent.spatial_consistency import compare


def test_compiled_assignment_preserves_ties_zeros_and_rectangles(monkeypatch):
    pytest.importorskip('numba');monkeypatch.setenv('AIDD_ASSIGNMENT_BACKEND','numba')
    rng=np.random.default_rng(58)
    for rows,columns in [(0,0),(3,0),(1,9),(9,1),(14,30),(30,14)]:
        for _ in range(12):
            values=rng.integers(0,3,(rows,columns)).astype(float)
            weights=rng.integers(0,3,rows).astype(float)
            assert np.array_equal(_maximum_weight_assignment(values,weights),_assignment_reference(values,weights))
            values=rng.random((rows,columns));weights=rng.random(rows)
            assert np.array_equal(_maximum_weight_assignment(values,weights),_assignment_reference(values,weights))


def test_review_counts_exclusive_templates_and_failed_scaffolds():
    db=sqlite3.connect(':memory:')
    db.execute('CREATE TABLE poses(mid,template,payload)');db.execute('CREATE TABLE members(mid,cluster,error)')
    for mid,t in [('a','x'),('b','x'),('b','y')]:
        db.execute('INSERT INTO poses VALUES(?,?,?)',(mid,t,json.dumps(dict(matched_anchor_count=1,composite_score=.4))))
    db.executemany('INSERT INTO members VALUES(?,?,?)',[('a','s',''),('b','unresolved-b','KekulizeException')])
    result=summarize(db,4)
    assert result['matching_molecules']==2 and result['valid_scaffold_groups']==1
    assert result['templates']==[dict(template='x',hits=2,exclusive_hits=1),dict(template='y',hits=1,exclusive_hits=0)]


def region(name,point):return dict(id=name,points=[point],radius=1.,minimum_atoms=1)


def test_dual_definition_states_and_no_double_counting():
    a=[region('A',[0,0,0]),region('B',[5,0,0])]
    b=[region('A',[0,0,0]),region('B',[8,0,0])]
    r=compare([[0,0,0],[5,0,0]],{'one':a,'two':b})
    assert r['states']==dict(A='agreement',B='boundary')
    assert compare([],{'one':a,'two':b})['states']==dict(A='unoccupied',B='unoccupied')
    assert compare([],{'one':a,'two':None})['states']==dict(A='unknown',B='unknown')
    close=[region('A',[0,0,0]),region('B',[.2,0,0])]
    assert compare([[.1,0,0]],{'one':close,'two':close})['agreed_regions']==0
    shifted=[[10,2,3],[15,2,3]]
    def shift(groups):return [dict(g,points=(np.array(g['points'])+[10,2,3]).tolist()) for g in groups]
    assert compare(shifted,{'one':shift(a),'two':shift(b)})==r


def test_heavy_graph_can_match_while_hydrogen_metadata_is_missing():
    from rdkit import Chem
    from aidd_agent.chemical_companion import BOND_DTYPE,_bond_order
    mol=Chem.MolFromSmiles('c1cc[nH]c1')
    chem=SimpleNamespace(atomic_numbers=np.array([a.GetAtomicNum() for a in mol.GetAtoms()]),
        formal_charges=np.array([a.GetFormalCharge() for a in mol.GetAtoms()]),atom_flags=np.ones(5,dtype=np.uint8),
        bonds=np.array([(b.GetBeginAtomIdx(),b.GetEndAtomIdx(),_bond_order(b),1,0) for b in mol.GetBonds()],dtype=BOND_DTYPE))
    r=compare_topology(mol,chem)
    assert all(r[k] for k in ('atomic_numbers_equal','charges_equal','aromatic_flags_equal','bonds_equal'))
    assert r['aromatic_heteroatom_hydrogens'][0]['hydrogens']==1
    assert r['hydrogen_metadata_in_companion'] is False
