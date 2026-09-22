"""Separate polymer-ligand contact and modified-residue sequence diagnostics."""
import csv
import json
from pathlib import Path

import numpy as np

from .structure_diversity import save, diverse_indices
from .structure_survey import target_residues


def token_similarity(left, right):
    """Normalized edit similarity of deposited CCD residue tokens, not chemistry."""
    previous=list(range(len(right)+1))
    for i,a in enumerate(left,1):
        current=[i]
        for j,b in enumerate(right,1):
            current.append(min(current[-1]+1,previous[j]+1,previous[j-1]+(a!=b)))
        previous=current
    return 1-previous[-1]/max(len(left),len(right),1)


def analyze(root):
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict
    from scipy.spatial import cKDTree
    root=Path(root);protein=json.loads((root/'protein.json').read_text())
    entries=json.loads((root/'entries.json').read_text())['entries']
    report=json.loads((root/'report.json').read_text());site=set(report['reference_site_residues'])
    results=[];failures=[]
    for entry in entries:
        path=root/'structures'/(entry['pdb_id']+'.cif')
        if not path.exists() or not entry['other_polymers']:continue
        try:
            d=MMCIF2Dict(str(path));mapping,notes=target_residues(d,protein)
            n=len(d['_atom_site.group_PDB'])
            def col(k,default=''):return d.get('_atom_site.'+k,[default]*n)
            entity,label_chain,auth_chain,auth_res,model,element,occ,alt=[col(k,v) for k,v in [
                ('label_entity_id',''),('label_asym_id',''),('auth_asym_id',''),('auth_seq_id',''),
                ('pdbx_PDB_model_num','1'),('type_symbol',''),('occupancy','1'),('label_alt_id','.')]]
            xyz=np.column_stack([np.asarray(col(k),float) for k in ('Cartn_x','Cartn_y','Cartn_z')])
            target=[i for i in range(n) if model[i]=='1' and element[i] not in {'H','D'} and (auth_chain[i],auth_res[i]) in mapping]
            if not target:raise ValueError('No verified target mapping')
            tree=cKDTree(xyz[target])
            seqrows=list(zip(d.get('_entity_poly_seq.entity_id',[]),d.get('_entity_poly_seq.num',[]),d.get('_entity_poly_seq.mon_id',[])))
            for partner in entry['other_polymers']:
                if partner['polymer_type']!='Protein':continue
                eid=partner['entity_id'];tokens=[r[2] for r in seqrows if r[0]==eid]
                for chain in sorted({label_chain[i] for i in range(n) if entity[i]==eid and model[i]=='1'}):
                    indices=[i for i in range(n) if label_chain[i]==chain and model[i]=='1' and element[i] not in {'H','D'}]
                    if not indices:continue
                    contacts=set()
                    for neighbors in tree.query_ball_point(xyz[indices],4.5):
                        contacts.update(mapping[auth_chain[target[j]],auth_res[target[j]]] for j in neighbors)
                    links=[]
                    for j,kind in enumerate(d.get('_struct_conn.conn_type_id',[])):
                        if kind not in {'covale','disulf'}:continue
                        if any(d.get(f'_struct_conn.ptnr{p}_label_asym_id',[])[j]!=chain for p in (1,2)):continue
                        links.append(dict(kind=kind,partners=[dict(residue=d[f'_struct_conn.ptnr{p}_label_seq_id'][j],
                            component=d[f'_struct_conn.ptnr{p}_label_comp_id'][j],atom=d[f'_struct_conn.ptnr{p}_label_atom_id'][j]) for p in (1,2)]))
                    results.append(dict(pdb_id=entry['pdb_id'],query_id=f"{entry['pdb_id']}:polymer:{chain}",
                        entity_id=eid,label_chain=chain,auth_chains=sorted({auth_chain[i] for i in indices}),
                        description=partner['description'],resolution=entry['resolution'],r_free=entry['r_free'],
                        ccd_sequence=tokens,length=len(tokens),observed_heavy_atoms=len(indices),
                        minimum_occupancy=min(float(occ[i]) for i in indices),
                        alternate_locations=sorted({alt[i] for i in indices if alt[i] not in {'.','?',''}}),
                        target_contacts=sorted(contacts),same_reference_site=len(site&contacts)>=3,
                        deposited_covalent_links=links,reference_readiness='requires_validated_polymer_chemical_graph_and_query_preparation'))
        except Exception as exc:failures.append(dict(pdb_id=entry['pdb_id'],error=str(exc)))
    unique={}
    for r in sorted(results,key=lambda r:(r['resolution'] or 999,r['query_id'])):
        if r['same_reference_site'] and r['length']<=50 and r['resolution'] is not None and r['resolution']<=2.5:
            signature=(tuple(r['ccd_sequence']),json.dumps(r['deposited_covalent_links'],sort_keys=True))
            unique.setdefault(signature,[]).append(r)
    ligands=[rows[0] for rows in unique.values()]
    matrix=np.array([[token_similarity(a['ccd_sequence'],b['ccd_sequence']) for b in ligands] for a in ligands])
    selected=diverse_indices(matrix,list(range(len(ligands))),8,.6)
    with (root/'peptide-sequence-similarity.csv').open('w',newline='',encoding='utf-8') as f:
        writer=csv.writer(f);writer.writerow(['query_id',*[r['query_id'] for r in ligands]])
        for row,scores in zip(ligands,matrix):writer.writerow([row['query_id'],*scores])
    result=dict(scope='Polymer ligand contacts and CCD-token sequence diagnostics; not chemical or 3D similarity',
        instances=results,unique_sequence_link_variants=len(unique),sequence_diversity_examples=[ligands[i] for i in selected],
        failures=failures,limitations=['Residue sequence similarity does not capture macrocyclization or stereochemistry fully',
        'Atom completeness, covalent graph and query preparation remain unvalidated; not search-ready references'])
    save(root/'polymer-ligands.json',result)
    return result
