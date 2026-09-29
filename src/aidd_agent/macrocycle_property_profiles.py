"""Continuous side-chain properties and coordinate-based steric proxies for review."""
import argparse
import json
from pathlib import Path

import numpy as np
from rdkit import Chem

from .macrocycle_descriptors import describe_peptide
from .mol2 import iter_mol2_records, load_rdkit_mol2


def steric_profile(mol, units):
    """VDW-expanded local-axis extents, not union volume or a clash energy."""
    mol=Chem.RemoveHs(Chem.Mol(mol));xyz=np.asarray(mol.GetConformer().GetPositions())
    cuts=Chem.RWMol(mol)
    for i,(_,_,carbon,_) in enumerate(units):cuts.RemoveBond(carbon,units[(i+1)%len(units)][0])
    fragments=Chem.GetMolFrags(cuts.GetMol(),sanitizeFrags=False);periodic=Chem.GetPeriodicTable()
    result=[]
    for n,ca,carbon,oxygen in units:
        side=sorted(set(next(f for f in fragments if ca in f))-{n,ca,carbon,oxygen})
        e1=xyz[carbon]-xyz[ca];norm=np.linalg.norm(e1)
        if norm<1e-8:raise ValueError('Degenerate backbone frame')
        e1/=norm;e2=xyz[n]-xyz[ca];e2-=e1*np.dot(e1,e2);norm=np.linalg.norm(e2)
        if norm<1e-8:raise ValueError('Degenerate backbone frame')
        e2/=norm;frame=np.stack((e1,e2,np.cross(e1,e2)),axis=1)
        if side:
            local=(xyz[side]-xyz[ca])@frame
            radii=np.array([periodic.GetRvdw(mol.GetAtomWithIdx(i).GetAtomicNum()) for i in side])
            lower=(local-radii[:,None]).min(0);upper=(local+radii[:,None]).max(0)
            centered=local-local.mean(0);eigen=np.linalg.eigvalsh(centered.T@centered/len(side))
            reach=float((np.linalg.norm(local,axis=1)+radii).max())
            # Side-chain branch degree near the alpha carbon, excluding backbone atoms.
            branches=sum(max(0,sum(a.GetAtomicNum()>1 for a in mol.GetAtomWithIdx(i).GetNeighbors())-2)
                         for i in side if mol.GetBondBetweenAtoms(ca,i) is not None)
        else:lower=upper=eigen=np.zeros(3);reach=0.;branches=0
        result.append(dict(side_heavy_atoms=len(side),proximal_branch_excess=branches,
            local_lower_vdw=lower.tolist(),local_upper_vdw=upper.tolist(),
            spatial_eigenvalues=eigen.tolist(),maximum_reach=reach))
    return result


def profile(mol,name=None):
    base=describe_peptide(mol,name,'typed');count=len(base['units'])
    matrix=np.asarray(base['descriptor']).reshape(count,59)
    return dict(version='sidechain-property-steric-review-v1',units=base['units'],
        omega_states=base['omega_states'],hard_group=base['hard_group'],
        properties=matrix[:,6:17].tolist(),typed_spatial_moments=matrix[:,17:].tolist(),
        sterics=steric_profile(mol,base['units']),
        limitations=['Continuous similarity inputs, not exact sequence identities or validated activity thresholds',
                     'Local VDW extents and spread are steric proxies, not molecular union volumes or docking clashes'])


def property_difference(first,second,rotation=0):
    """Compare corresponding directed positions; identity admissibility is external."""
    a=np.asarray(first['properties'],float);b=np.asarray(second['properties'],float)
    if a.shape!=b.shape or a.ndim!=2 or a.shape[1]!=11:raise ValueError('Incompatible property profiles')
    if not isinstance(rotation,int) or not 0<=rotation<len(a):raise ValueError('Invalid directed rotation')
    b=np.roll(b,-rotation,axis=0)
    def steric_vector(p):
        return np.asarray([[r['side_heavy_atoms']/10,r['proximal_branch_excess']/3,
                            *np.asarray(r['local_lower_vdw'])/5,*np.asarray(r['local_upper_vdw'])/5,
                            *np.sqrt(np.maximum(r['spatial_eigenvalues'],0))/5,r['maximum_reach']/5]
                           for r in p['sterics']])
    sa=steric_vector(first);sb=np.roll(steric_vector(second),-rotation,axis=0)
    if not all(np.isfinite(v).all() for v in (a,b,sa,sb)):raise ValueError('Nonfinite profile')
    return dict(property_rms=float(np.sqrt(np.mean((a-b)**2))),
                steric_rms=float(np.sqrt(np.mean((sa-sb)**2))),
                maximum_local_steric_difference=float(np.linalg.norm(sa-sb,axis=1).max()),
                note='Separate diagnostics; no calibrated combined cutoff or automatic merge')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--mol2',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--limit',type=int,default=200)
    a=p.parse_args()
    if a.limit<1:raise ValueError('Positive review limit required')
    with a.output.open('x',encoding='utf-8') as stream:
        for i,record in enumerate(iter_mol2_records(a.mol2)):
            if i>=a.limit:break
            mol,mode=load_rdkit_mol2(record.raw_text,record.name)
            if mode!='strict':raise ValueError('Strict molecular identity required')
            result=profile(mol)
            stream.write(json.dumps(dict(source_record=record.name,source_sha256=record.content_sha256,**result))+'\n')


if __name__=='__main__':main()
