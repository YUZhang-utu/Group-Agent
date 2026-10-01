"""Backbone shape and chemically typed spatial moments anchored to directed residues."""
import json
from pathlib import Path

import numpy as np
from rdkit import Chem

from .macrocycle_descriptors import FAMILIES, describe_peptide, feature_factory
from .mol2 import parse_mol2_block, load_rdkit_mol2

VERSION = 'directed-residue-joint-spatial-v1'
BACKBONE_WIDTH = 24
SPATIAL_WIDTH = 4 * 8 * 7
WIDTH = 69 + BACKBONE_WIDTH + SPATIAL_WIDTH


def geometry(mol, units, torsions):
    """Finite spatial summary, not an injective encoding or receptor interaction score."""
    mol=Chem.RemoveHs(Chem.Mol(mol))
    xyz=np.asarray(mol.GetConformer().GetPositions(),dtype=float)
    units=np.asarray(units,dtype=int);n=len(units)
    if n<2 or units.shape!=(n,4) or not np.isfinite(xyz).all():raise ValueError('Invalid cyclic backbone')
    ca=xyz[units[:,1]]
    frames=[]
    for nitrogen,alpha,carbon,_ in units:
        e1=xyz[carbon]-xyz[alpha];norm=np.linalg.norm(e1)
        if norm<1e-8:raise ValueError('Degenerate backbone frame')
        e1/=norm;e2=xyz[nitrogen]-xyz[alpha];e2-=e1*np.dot(e1,e2);norm=np.linalg.norm(e2)
        if norm<1e-8:raise ValueError('Degenerate backbone frame')
        e2/=norm;frames.append(np.stack((e1,e2,np.cross(e1,e2)),axis=1))
    cuts=Chem.RWMol(mol)
    for i,(_,_,carbon,_) in enumerate(units):cuts.RemoveBond(int(carbon),int(units[(i+1)%n,0]))
    fragments=Chem.GetMolFrags(cuts.GetMol(),sanitizeFrags=False)
    features=feature_factory().GetFeaturesForMol(mol)
    clouds=[]
    for unit in units:
        fragment=set(next(f for f in fragments if int(unit[1]) in f))
        side=fragment-set(map(int,unit))
        families=[]
        for family in FAMILIES:
            # Include residue backbone donor/acceptor sites as well as side-chain sites.
            points=[xyz[list(f.GetAtomIds())].mean(0) for f in features
                    if f.GetFamily()==family and set(f.GetAtomIds())<=fragment]
            families.append(np.asarray(points,dtype=float).reshape(-1,3))
        families.append(xyz[sorted(side)].reshape(-1,3))
        clouds.append(families)
    torsions=np.asarray(torsions,dtype=float)
    if torsions.shape!=(n,6):raise ValueError('Expected backbone sine/cosine channels')
    backbone=[*torsions.mean(0),*torsions.std(0)]
    for lag in (1,2,3,n//2):
        distances=np.linalg.norm(ca-np.roll(ca,-lag,axis=0),axis=1)/5
        backbone.extend([distances.mean(),distances.std()])
    backbone_atoms=xyz[units[:,:3].ravel()];centered=(backbone_atoms-backbone_atoms.mean(0))/5
    backbone.extend(np.sqrt(np.maximum(np.linalg.eigvalsh(centered.T@centered/len(centered)),0)))
    edges=(np.roll(ca,-1,axis=0)-ca)/5
    backbone.append(np.einsum('ij,ij->i',np.cross(edges,np.roll(edges,-1,axis=0)),np.roll(edges,-2,axis=0)).mean())
    spatial=[]
    for lag in (0,1,2,n//2):
        for family in range(8):
            points=[];counts=[]
            for i in range(n):
                cloud=clouds[(i+lag)%n][family];counts.append(len(cloud))
                if len(cloud):points.append((cloud-ca[i])@frames[i]/5)
            if points:
                values=np.concatenate(points)
                spatial.extend([np.mean(counts)/(10 if family==7 else 4),*values.mean(0),*values.std(0)])
            else:spatial.extend([0.]*7)
    result=np.array([*backbone,*spatial],dtype='<f4')
    if result.shape!=(BACKBONE_WIDTH+SPATIAL_WIDTH,) or not np.isfinite(result).all():raise ValueError('Invalid joint spatial vector')
    return result


def compute(task):
    path,index,text,payload,variant,old_vector=task
    row=json.loads(payload);p=row['provenance']
    record=parse_mol2_block(Path(path),index,text)
    if (record.content_sha256,record.name,record.molecule_name)!=(p['content_sha256'],p['source_record_name'],p['source_name']):raise ValueError('Raw spatial source identity/hash mismatch')
    if str(Path(path).resolve())!=str(Path(p['source_path']).resolve()) or index!=p['source_record_index']:raise ValueError('Source locator mismatch')
    mol,mode=load_rdkit_mol2(text,record.name)
    if mode!='strict':raise ValueError('Strict spatial source chemistry required')
    name=None if 'cross_conformer_alignment_unverified' in p['flags'] else record.molecule_name
    description=describe_peptide(mol,name,variant)
    if description['hard_group']!=row['hard_group'] or description['units']!=p['units'] or not np.allclose(description['descriptor'],row['descriptor'],rtol=1e-9,atol=1e-9):raise ValueError('Stored backbone descriptor not reproduced')
    old=np.frombuffer(old_vector,dtype='<f4')
    if old.shape!=(69,) or not np.isfinite(old).all():raise ValueError('Invalid reusable property vector')
    torsions=np.asarray(row['descriptor']).reshape(len(p['units']),-1)[:,:6]
    new=geometry(mol,p['units'],torsions)
    mid=row.get('molecule_id')
    if not isinstance(mid,str) or not mid:raise ValueError('Molecule identity required for disjoint checks')
    return row['conformer_id'],mid,np.concatenate((old,new)).astype('<f4').tobytes()
