"""Within-job AF3 ligand pose consistency after protein C-alpha alignment."""
from pathlib import Path
import numpy as np


def coordinates(path,protein_chains,ligand_chains):
    from Bio.PDB import MMCIFParser
    models=list(MMCIFParser(QUIET=True).get_structure('sample',str(path)))
    if len(models)!=1:raise ValueError('Expected one AF3 model')
    protein={};ligands={}
    for chain in models[0]:
        for residue in chain:
            for atom in residue:
                key=(chain.id,residue.id,atom.name)
                if chain.id in protein_chains and atom.name=='CA':protein[key]=np.asarray(atom.coord,dtype=float)
                if chain.id in ligand_chains and atom.element not in {'H','D'}:ligands[key]=np.asarray(atom.coord,dtype=float)
    if len(protein)<3 or not ligands:raise ValueError('Insufficient protein/ligand coordinates for pose comparison')
    if not np.isfinite(list(protein.values())).all() or not np.isfinite(list(ligands.values())).all():raise ValueError('Nonfinite model coordinates')
    return protein,ligands


def summarize(samples,input_data,cutoff=2.):
    protein_chains=[e['protein']['id'] for e in input_data['sequences'] if 'protein' in e]
    ligand_chains=[e['ligand']['id'] for e in input_data['sequences'] if 'ligand' in e]
    rows=[];clusters=[];reference=None;keys=None
    for sample in samples:
        row=dict(sample=sample['sample'],status='not_assessed')
        try:
            if not sample.get('model'):raise ValueError('Sample model missing')
            protein,ligand=coordinates(Path(sample['model']),protein_chains,ligand_chains)
            if reference is None:
                reference=(protein,ligand);keys=sorted(ligand)
            common=sorted(set(reference[0]) & set(protein))
            if len(common)<3 or sorted(ligand)!=keys:raise ValueError('Atom correspondence differs between samples')
            x=np.array([protein[k] for k in common]);y=np.array([reference[0][k] for k in common])
            center=x.mean(axis=0);target=y.mean(axis=0)
            u,_,vt=np.linalg.svd((x-center).T@(y-target))
            correction=np.eye(3);correction[2,2]=np.linalg.det(u@vt)
            rotation=u@correction@vt
            aligned=(np.array([ligand[k] for k in keys])-center)@rotation+target
            rmsds=[float(np.sqrt(np.mean(np.sum((aligned-c)**2,axis=1)))) for c in clusters]
            eligible=[i for i,d in enumerate(rmsds) if d<=cutoff]
            cluster=eligible[0] if eligible else len(clusters)
            if not eligible:clusters.append(aligned)
            row.update(status='compared',cluster=cluster+1,aligned_ca_atoms=len(common),
                protein_alignment_rmsd=float(np.sqrt(np.mean(np.sum(((x-center)@rotation+target-y)**2,axis=1)))),
                ligand_rmsd_to_first=float(np.sqrt(np.mean(np.sum((aligned-np.array([reference[1][k] for k in keys]))**2,axis=1)))))
        except (ValueError,KeyError,TypeError,OSError) as exc:row['error']=str(exc)
        rows.append(row)
    return dict(status='compared' if rows and all(r['status']=='compared' for r in rows) else 'incomplete',
        samples=rows,clusters=len(clusters),cutoff_angstrom=cutoff,
        limitations=['Exploratory greedy clusters in sample order; not validated binding modes',
                    'Atom-name correspondence, no symmetry correction; whole-protein C-alpha alignment',
                    'Pose consistency is not affinity or binding evidence'])
