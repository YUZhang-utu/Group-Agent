"""Chain-resolved, protein-aligned admission into a reference pocket cohort."""
from collections import defaultdict
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from .chemistry_prep import _mmcif_dict, _column
from .structure_survey import target_residues

POLICY=dict(minimum_pocket_coverage=.8,minimum_aligned_residues=12,
            maximum_pocket_ca_rmsd=1.5,maximum_centroid_distance=5.,
            minimum_ligand_envelope_overlap=.4,ligand_envelope_distance=3.,
            minimum_shared_contact_residues=3,contact_distance=4.5)


def fit_protein(moving, fixed):
    moving=np.asarray(moving,float);fixed=np.asarray(fixed,float)
    if moving.shape!=fixed.shape or moving.ndim!=2 or len(moving)<3:raise ValueError('Insufficient protein correspondence')
    a=moving.mean(0);b=fixed.mean(0)
    if np.linalg.matrix_rank(moving-a)<2:raise ValueError('Degenerate protein alignment')
    u,_,vt=np.linalg.svd((moving-a).T@(fixed-b));r=vt.T@u.T
    if np.linalg.det(r)<0:vt[-1]*=-1;r=vt.T@u.T
    t=b-r@a;m=np.eye(4);m[:3,:3]=r;m[:3,3]=t
    distances=np.linalg.norm(moving@r.T+t-fixed,axis=1)
    return m,float(np.sqrt(np.mean(distances**2))),distances.tolist()


def read_structure(path, protein):
    d=_mmcif_dict(path);mapping,notes=target_residues(d,protein)
    n=len(d['_atom_site.group_PDB'])
    names=('group_PDB','label_entity_id','label_asym_id','auth_asym_id','auth_seq_id','label_seq_id',
           'auth_comp_id','auth_atom_id','type_symbol','pdbx_PDB_model_num','label_alt_id','occupancy','pdbx_PDB_ins_code')
    columns={k:_column(d,'_atom_site.'+k,n,'1' if k in {'pdbx_PDB_model_num','occupancy'} else '') for k in names}
    xyz=np.column_stack([np.asarray(d['_atom_site.'+k],float) for k in ('Cartn_x','Cartn_y','Cartn_z')])
    chains={};atoms=[]
    for i in range(n):
        if columns['pdbx_PDB_model_num'][i]!='1' or columns['type_symbol'][i] in {'H','D'}:continue
        row={k:v[i] for k,v in columns.items()};row.update(xyz=xyz[i],index=i)
        row['canonical_residue']=mapping.get((row['auth_asym_id'],row['auth_seq_id']))
        atoms.append(row)
        if row['group_PDB']=='ATOM' and row['canonical_residue'] is not None:
            chains.setdefault(row['auth_asym_id'],[]).append(row)
    # Organism assertions refer to the biological source, never expression host.
    taxonomy={}
    for category,field in [('entity_src_gen','pdbx_gene_src_ncbi_taxonomy_id'),('entity_src_nat','pdbx_ncbi_taxonomy_id')]:
        ids=d.get('_'+category+'.entity_id',[])
        if isinstance(ids,str):ids=[ids]
        for eid,tax in zip(ids,_column(d,'_'+category+'.'+field,len(ids))):
            if tax:taxonomy.setdefault(eid,set()).add(tax)
    expected=str(protein.get('organism',{}).get('taxonId',''))
    for chain in list(chains):
        observed=taxonomy.get(chains[chain][0]['label_entity_id'],set())
        if not observed or not expected or observed!={expected}:del chains[chain]
    return dict(data=d,atoms=atoms,chains=chains,mapping_notes=notes,path=str(path))


def ligand_atoms(structure, row):
    return [a for a in structure['atoms'] if a['auth_asym_id']==row['chain'] and
            a['auth_seq_id']==row['residue'] and a['auth_comp_id']==row['ccd_id'] and a['group_PDB']=='HETATM']


def contact_residues(atoms, ligand, distance=4.5):
    if not atoms or not ligand:return set()
    tree=cKDTree(np.array([a['xyz'] for a in atoms]));out=set()
    for neighbors in tree.query_ball_point(np.array([a['xyz'] for a in ligand]),distance):
        out.update(atoms[i]['canonical_residue'] for i in neighbors)
    return out-{None}


def assembly_compatible(structure, chain, ligand, assembly='1'):
    """Accept explicit identity assembly membership; transformed interfaces need review."""
    d=structure['data'];n=len(d.get('_pdbx_struct_assembly_gen.assembly_id',[]))
    groups=list(zip(_column(d,'_pdbx_struct_assembly_gen.assembly_id',n),
        _column(d,'_pdbx_struct_assembly_gen.oper_expression',n),_column(d,'_pdbx_struct_assembly_gen.asym_id_list',n)))
    required={a['label_asym_id'] for a in structure['chains'][chain]}|{a['label_asym_id'] for a in ligand}
    # Operator 1 must actually be identity; its numeric label is not proof.
    ops=d.get('_pdbx_struct_oper_list.id',[])
    if isinstance(ops,str):ops=[ops]
    identity=set()
    for i,op in enumerate(ops):
        try:
            r=np.array([[float(_column(d,f'_pdbx_struct_oper_list.matrix[{a}][{b}]',len(ops))[i]) for b in range(1,4)] for a in range(1,4)])
            t=np.array([float(_column(d,f'_pdbx_struct_oper_list.vector[{a}]',len(ops))[i]) for a in range(1,4)])
            if np.allclose(r,np.eye(3),atol=1e-6) and np.allclose(t,0,atol=1e-6):identity.add(op)
        except ValueError:continue
    return any(aid==assembly and op.strip('()') in identity and required<=set(chains.split(',')) for aid,op,chains in groups)


def ca_map(atoms):
    out={}
    for a in atoms:
        if a['auth_atom_id']=='CA' and not a['label_alt_id']:
            out[a['canonical_residue']]=a['xyz']
    return out


def admit(reference, ref_row, ref_chain, candidate, row, chain, policy=None, assembly='1'):
    p=dict(POLICY,**(policy or {}));result=dict(query_id=row['query_id'],target_chain=chain,status='pending',reasons=[])
    if chain not in candidate['chains']:result['reasons']=['target_identity_or_organism_unverified'];return result
    lig=ligand_atoms(candidate,row);rl=ligand_atoms(reference,ref_row)
    if not lig:result['reasons']=['ligand_instance_missing'];return result
    if not assembly_compatible(candidate,chain,lig,assembly):result['reasons']=['assembly_or_packing_context_unresolved'];return result
    ra=reference['chains'][ref_chain];ca=candidate['chains'][chain]
    rc=contact_residues(ra,rl,p['contact_distance']);cc=contact_residues(ca,lig,p['contact_distance'])
    rm=ca_map(ra);cm=ca_map(ca);common=sorted(set(rm)&set(cm))
    result.update(reference_contacts=sorted(rc),candidate_contacts=sorted(cc),shared_contacts=sorted(rc&cc),
                  pocket_coverage=len(rc&set(cm))/max(len(rc),1),alignment_residues=common)
    if len(common)<p['minimum_aligned_residues'] or result['pocket_coverage']<p['minimum_pocket_coverage']:
        result['reasons']=['insufficient_observed_pocket_or_domain'];return result
    # Fit the pocket and its local backbone environment, not unrelated domains.
    rc_xyz=np.array([rm[r] for r in rc if r in rm]);tree=cKDTree(rc_xyz)
    local=[r for r in common if tree.query(rm[r])[0]<=8.]
    if len(local)<p['minimum_aligned_residues']:result['reasons']=['insufficient_local_alignment'];return result
    m,rmsd,residuals=fit_protein([cm[r] for r in local],[rm[r] for r in local])
    xyz=np.array([a['xyz'] for a in lig])@m[:3,:3].T+m[:3,3];rxyz=np.array([a['xyz'] for a in rl])
    center=float(np.linalg.norm(xyz.mean(0)-rxyz.mean(0)))
    overlap=float(np.mean(cKDTree(rxyz).query(xyz)[0]<=p['ligand_envelope_distance']))
    local_atoms=[a for a in ca if a['canonical_residue'] in rc]
    quality_uncertain=any(a['label_alt_id'] or float(a['occupancy'])<.9 or a['pdbx_PDB_ins_code'] for a in local_atoms)
    result.update(transform=m.tolist(),local_ca_rmsd=rmsd,alignment_residues=local,residuals=residuals,
                  centroid_distance=center,envelope_overlap=overlap)
    if quality_uncertain:result['reasons'].append('pocket_alternate_insertion_or_low_occupancy')
    if rmsd>p['maximum_pocket_ca_rmsd']:result['reasons'].append('distinct_or_unresolved_receptor_state')
    if len(rc&cc)<p['minimum_shared_contact_residues'] or center>p['maximum_centroid_distance'] or overlap<p['minimum_ligand_envelope_overlap']:
        result['status']='rejected';result['reasons'].append('different_or_unresolved_site')
    if not result['reasons']:result['status']='admitted'
    return result


def cohort(diversity_root, output, protein, reference_pdb, reference_query=None, target_chain=None, assembly='1'):
    root=Path(diversity_root);output=Path(output);output.mkdir(parents=True,exist_ok=True)
    rows=json.loads((root/'ligand-instances.json').read_text());structures={}
    def get(code):
        if code not in structures:structures[code]=read_structure(root/'structures'/f'{code}.cif',protein)
        return structures[code]
    options=[]
    for row in rows:
        if row['pdb_id']!=reference_pdb or not row['quality_passed']:continue
        if reference_query and row['query_id']!=reference_query:continue
        s=get(row['pdb_id']);lig=ligand_atoms(s,row)
        for chain,atoms in s['chains'].items():
            if target_chain and chain!=target_chain:continue
            if len(contact_residues(atoms,lig))>=3 and assembly_compatible(s,chain,lig,assembly):options.append((row,chain))
    if len(options)!=1:
        return dict(status='complete',readiness='needs_reference_instance',reference_options=[dict(query_id=r['query_id'],target_chain=c) for r,c in options],
                    reason='Choose one ligand and target chain; polymer/interface references require explicit supported preparation')
    rr,rc=options[0];reference=get(rr['pdb_id']);decisions=[];accepted=[]
    for row in rows:
        if not row['quality_passed']:
            decisions.append(dict(query_id=row['query_id'],status='rejected',reasons=row['quality_rejections']));continue
        s=get(row['pdb_id']);checks=[admit(reference,rr,rc,s,row,c,assembly=assembly) for c in s['chains']]
        passed=[c for c in checks if c['status']=='admitted']
        if len(passed)==1:
            decision=passed[0];accepted.append(dict(row,admission=decision))
        elif len(passed)>1:decision=dict(query_id=row['query_id'],status='pending',reasons=['multiple_target_chain_contacts'],chain_checks=checks)
        else:decision=dict(query_id=row['query_id'],status='pending',reasons=['no_unique_admitted_target_site'],chain_checks=checks)
        decisions.append(decision)
    result=dict(status='complete',readiness='cohort_ready' if accepted else 'no_admitted_complexes',
        reference=dict(query_id=rr['query_id'],target_chain=rc,assembly_id=assembly,coordinate_frame=rr['query_id']),
        policy=POLICY,policy_status='Versioned exploratory same-state admission, not universal calibration',
        admitted=accepted,decisions=decisions,limitations=['Transformed biological assemblies and ambiguous interfaces are pending, not silently pooled',
        'Pocket mutations/insertions and missing correspondence can exclude legitimate alternative states',
        'Density metrics remain unassessed; cohort is geometrically admitted, not experimentally validated binding evidence'])
    return result
