import copy

import numpy as np
from rdkit import Chem

from aidd_agent.macrocycle_property_profiles import profile,property_difference


def molecule():
    mol=Chem.MolFromSmiles('N1[C@@H](C)C(=O)N(C)[C@@H](C)C(=O)N2CCC[C@H]2C1=O')
    conf=Chem.Conformer(mol.GetNumAtoms())
    for i,point in enumerate(np.random.default_rng(81).normal(size=(mol.GetNumAtoms(),3))):conf.SetAtomPosition(i,point)
    mol.AddConformer(conf);return mol


def test_sterics_rigid_invariance_and_sidechain_response():
    mol=molecule();a=profile(mol,'c--A-Anme-P-c');moved=Chem.Mol(mol)
    rot=np.array([[0,-1,0],[1,0,0],[0,0,1]])
    for i,p in enumerate(np.array(moved.GetConformer().GetPositions())@rot+[7,-9,2]):
        moved.GetConformer().SetAtomPosition(i,p)
    b=profile(moved,'c--A-Anme-P-c');distance=property_difference(a,b)
    # Square roots of nearly zero covariance eigenvalues amplify roundoff.
    assert distance['property_rms']==0 and distance['steric_rms']<1e-8
    atom=next(x.GetIdx() for x in moved.GetAtoms() if x.GetAtomicNum()==6 and x.GetDegree()==1)
    point=np.array(moved.GetConformer().GetAtomPosition(atom));moved.GetConformer().SetAtomPosition(atom,point+[5,2,0])
    c=profile(moved,'c--A-Anme-P-c');distance=property_difference(b,c)
    assert distance['property_rms']==0 and distance['steric_rms']>.01


def test_properties_are_continuous_and_rotation_is_directed():
    a=profile(molecule(),'c--A-Anme-P-c');b=copy.deepcopy(a)
    b['properties'][0][0]+=.1
    assert 0<property_difference(a,b)['property_rms']<.1
    b=copy.deepcopy(a)
    for key in ('properties','sterics'):b[key]=b[key][-1:]+b[key][:-1]
    assert property_difference(a,b,rotation=1)['steric_rms']==0
