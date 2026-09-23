"""Independent pocket evidence and diverse ligand shape templates in one frame."""
import copy
import hashlib
import json
from pathlib import Path

import numpy as np

from .gaussian_batch import _atomic_json, _load_query, write_gaussian_query
from .expanded_wee1 import fingerprint
from .contact_groups import contact_semantics


def transform_query(query, matrix, output, source):
    values={k:query[k].copy() for k in ('shape_points','feature_points','feature_types',
        'feature_directions','feature_direction_kinds','anchored_weights','anchor_feature_indices')}
    for key in ('shape_points','feature_points'):
        values[key]=values[key]@matrix[:3,:3].T+matrix[:3,3]
    values['feature_directions']=values['feature_directions']@matrix[:3,:3].T
    # Shape ranking is independent of which interactions will later be adopted.
    values['anchored_weights']=np.ones(len(values['feature_points']))
    values['anchor_feature_indices']=np.arange(len(values['feature_points']))
    write_gaussian_query(output,**values,source=source)
    return values


def consensus(observations, prepared, distance=1.0):
    """Complete-link spatial/directional modes; no pooling of distant contacts."""
    groups=[]
    for item in sorted(observations,key=lambda x:(x['key'],x['query_id'],x['source_anchor'])):
        compatible=[]
        for group in groups:
            if group[0]['key']!=item['key']:continue
            if all(np.linalg.norm(np.array(o['point'])-item['point'])<=distance and
                   (item['kind']==0 or (abs(np.dot(o['direction'],item['direction'])) if item['kind']==2
                    else np.dot(o['direction'],item['direction']))>=.5) for o in group):compatible.append(group)
        if compatible:compatible[0].append(item)
        else:groups.append([item])
    anchors=[];features=[]
    for group in groups:
        first=group[0];residue,partner,interaction,ftype,kind=first['key']
        points=np.array([o['point'] for o in group]);medoid=int(np.argmin(np.linalg.norm(points[:,None]-points[None,:],axis=2).sum(1)))
        representative=group[medoid]
        pdbs=sorted({o['pdb_id'] for o in group});chemotypes=sorted({o['chemotype'] for o in group})
        eligible=sorted({q['pdb_id'] for q in prepared if all(f'{residue}:{atom}' in q['observed_partner_atoms'] for atom in partner.split('/'))})
        signature=json.dumps([first['key'],sorted({o['query_id']+':'+o['source_anchor'] for o in group})])
        aid='pocket:'+hashlib.sha256(signature.encode()).hexdigest()[:16]
        fraction=len(pdbs)/max(len(eligible),1)
        anchors.append(dict(anchor_id=aid,feature_index=len(features),score_column=len(features),
            target_residue=residue,protein_atom=partner,feature_class=interaction,
            **contact_semantics(partner,interaction),
            evidence_id='consensus:'+aid,distinct_structures=len(pdbs),eligible_structures=len(eligible),
            frequency=fraction,distinct_chemotypes=len(chemotypes),pdb_ids=pdbs,eligible_pdb_ids=eligible,
            mandatory_proposal_eligible=len(pdbs)>=3 and len(chemotypes)>=2 and fraction>=.5,
            observations=group,representative_query=representative['query_id'],
            interpretation='Observed mode-specific geometric recurrence, not proof of energetic necessity'))
        features.append(representative)
    return anchors,features


def build(source, output, reference_query=None, target_chain=None, maximum_templates=8):
    from rdkit import Chem, DataStructs
    from rdkit.Chem import rdFingerprintGenerator
    from rdkit.Chem.Scaffolds import MurckoScaffold
    from .consensus_admission import cohort, read_structure
    from .structure_diversity import PublicCache, diverse_indices
    from .anchor_extraction import extract_and_write_query_manifest
    from .screening_selection import anchor_evidence
    from .classified_features import crystal_environment, extra_hypotheses
    from .gaussian_batch import prepare_gaussian_query
    root=Path(source).parent;out=Path(output);out.mkdir(parents=True,exist_ok=True)
    survey=json.loads(Path(source).read_text());protein=json.loads((root/'protein.json').read_text())
    admitted=cohort(root,out,protein,survey['reference_site_pdb'],reference_query,target_chain)
    _atomic_json(out/'cohort.json',admitted)
    report=dict(kind='structure_consensus',status='complete',target=protein,cohort=admitted,
        readiness=admitted['readiness'],sources=fingerprint([source,root/'protein.json',root/'ligand-instances.json']))
    if admitted['readiness']!='cohort_ready':
        _atomic_json(out/'report.json',report);return report
    prepared=[];observations=[];failures=[];cache=PublicCache(out/'raw')
    for row in admitted['admitted']:
        qid=row['query_id'];folder=out/'queries'/qid.replace(':','_');folder.mkdir(parents=True,exist_ok=True)
        cif=root/'structures'/f'{row["pdb_id"]}.cif';ccd=folder/(row['ccd_id']+'.cif')
        try:
            ccd.write_bytes(cache.get(f'https://files.rcsb.org/ligands/download/{row["ccd_id"]}.cif',binary=True))
            manifest=folder/'query_manifest.json'
            extract_and_write_query_manifest(cif,ccd,row['ccd_id'],qid,manifest)
            native=folder/'native.npz';prepare_gaussian_query(cif,ccd,manifest,native)
            anchors,paths,query=anchor_evidence(native)
            atoms,waters,metals=crystal_environment(cif,qid)
            chain=row['admission']['target_chain'];atoms=[a for a in atoms if a['chain']==chain]
            anchors+=extra_hypotheses(qid,query,atoms,waters,metals)
            structure=read_structure(cif,protein)
            mapping={(a['auth_asym_id'],a['auth_seq_id']):a['canonical_residue'] for a in structure['chains'][chain]}
            matrix=np.array(row['admission']['transform']);aligned=folder/'aligned.npz'
            transformed=transform_query(query,matrix,aligned,dict(query_id=qid,protein_alignment=row['admission']))
            molecule=Chem.MolFromSmiles(row['smiles'])
            scaffold=MurckoScaffold.MurckoScaffoldSmiles(mol=molecule,includeChirality=False) or Chem.MolToSmiles(molecule)
            observed=sorted({f'{a["canonical_residue"]}:{a["auth_atom_id"]}' for a in structure['chains'][chain]
                             if not a['label_alt_id'] and float(a['occupancy'])>=.9})
            prepared.append(dict(query_id=qid,pdb_id=row['pdb_id'],smiles=row['smiles'],chemotype=scaffold,
                resolution=row['resolution'],query_npz=str(aligned.resolve()),observed_partner_atoms=observed))
            for anchor in anchors:
                partner=anchor['evidence'].get('protein_partner',{});res=mapping.get((partner.get('chain'),partner.get('residue_number')))
                if res is None:continue
                i=anchor['feature_index'];kind=int(transformed['feature_direction_kinds'][i]);ft=int(transformed['feature_types'][i])
                interaction=anchor.get('interaction_class',anchor['feature_class'])
                observations.append(dict(key=[res,partner.get('atom_name',''),interaction,ft,kind],
                    query_id=qid,pdb_id=row['pdb_id'],chemotype=scaffold,source_anchor=anchor['anchor_id'],
                    point=transformed['feature_points'][i].tolist(),direction=transformed['feature_directions'][i].tolist(),kind=kind,
                    evidence=anchor['evidence'],source_url=f'https://www.rcsb.org/structure/{row["pdb_id"]}'))
            report['sources'].update(fingerprint([cif,ccd,manifest,*paths,aligned,aligned.with_suffix('.manifest.json')]))
        except Exception as exc:
            failures.append(dict(query_id=qid,error=type(exc).__name__,detail=str(exc)[:400]))
    anchors,features=consensus(observations,prepared)
    if not anchors:
        report.update(readiness='needs_prepared_complexes',failures=failures)
        _atomic_json(out/'report.json',report);return report
    ref=admitted['reference'];refrow=next(r for r in admitted['admitted'] if r['query_id']==ref['query_id'])
    receptor=read_structure(root/'structures'/f'{refrow["pdb_id"]}.cif',protein)
    receptor_path=out/'receptor.npz'
    np.savez_compressed(receptor_path,points=np.array([a['xyz'] for a in receptor['chains'][ref['target_chain']]]))
    pocket=out/'pocket.npz'
    write_gaussian_query(pocket,shape_points=[f['point'] for f in features],feature_points=[f['point'] for f in features],
        feature_types=[f['key'][3] for f in features],feature_directions=[f['direction'] for f in features],
        feature_direction_kinds=[f['kind'] for f in features],anchored_weights=np.ones(len(features)),anchor_feature_indices=range(len(features)),
        source=dict(reference_frame=ref,scope='Consensus features only; never used as a ligand shape template'))
    unique={}
    for q in sorted(prepared,key=lambda q:(q['resolution'],q['query_id'])):unique.setdefault(q['smiles'],q)
    templates=list(unique.values());generator=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048,includeChirality=True)
    fps=[generator.GetFingerprint(Chem.MolFromSmiles(q['smiles'])) for q in templates]
    similarity=np.array([DataStructs.BulkTanimotoSimilarity(fp,fps) for fp in fps])
    chosen=diverse_indices(similarity,list(range(len(templates))),maximum_templates,.6)
    uncovered=[templates[i]['query_id'] for i in range(len(templates)) if similarity[i,chosen].max()<.6]
    report.update(readiness='proposal_ready',anchors=anchors,templates=templates,prepared_complexes=prepared,failures=failures,
        proposed_template_ids=[templates[i]['query_id'] for i in chosen],consensus_npz=str(pocket.resolve()),receptor_npz=str(receptor_path.resolve()),
        template_selection=dict(maximum=maximum_templates,similarity_threshold=.6,unique_ligands=len(templates),
            covered=len(templates)-len(uncovered),uncovered=uncovered,scope='Morgan chemical similarity coverage, not 3D retrieval recall'),
        mandatory_support_policy=dict(minimum_frequency=.5,minimum_distinct_pdbs=3,minimum_chemotypes=2),
        limitations=admitted['limitations']+['Nonpolymer CCD-backed references only; polymer chemistry requires separate preparation',
            'Only the admitted receptor-state cohort is pooled; alternative states remain pending',
            'Support thresholds and spatial mode radius are exploratory', 'Known active retention is not established by crystal self-controls'])
    report['sources'].update(fingerprint([out/'cohort.json',pocket,pocket.with_suffix('.manifest.json'),receptor_path]))
    from .contact_evidence import attach
    attach(report,out/'contacts.json')
    _atomic_json(out/'report.json',report)
    from .consensus_report import render
    render(out)
    return report
