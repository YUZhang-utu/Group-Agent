"""Independent ChEMBL binding controls, kept outside catalog coverage counts."""
from collections import Counter
from pathlib import Path
from types import SimpleNamespace
from urllib.parse import urlencode
import numpy as np

from .gaussian_batch import _atomic_json, _load_query, prepare_seeds
from .structure_diversity import PublicCache


def binding_records(protein, fetch, maximum=24):
    """Bounded diagnostic panel; never merge incompatible assay values into a potency."""
    base='https://www.ebi.ac.uk/chembl/api/data/'
    query=urlencode({'target_components__accession':protein['accession'],'format':'json','limit':100})
    targets=fetch(base+'target.json?'+query).get('targets',[])
    eligible=[t for t in targets if t.get('target_type')=='SINGLE PROTEIN' and
        str(t.get('tax_id'))==str(protein['organism']['taxonId']) and
        {c.get('accession') for c in t.get('target_components',[])}=={protein['accession']}]
    rows=[];seen=set();assays={};rejected=Counter()
    for target in eligible:
        query=urlencode(dict(target_chembl_id=target['target_chembl_id'],assay_type='B',standard_relation='=',
            standard_units='nM',standard_value__lte=1000,standard_type__in='Ki,Kd,IC50',limit=200,format='json'))
        page=fetch(base+'activity.json?'+query)
        assay_ids=sorted({r['assay_chembl_id'] for r in page.get('activities',[])})
        if assay_ids:
            assay_page=fetch(base+'assay.json?'+urlencode(dict(assay_chembl_id__in=','.join(assay_ids),limit=200,format='json')))
            assays.update({a['assay_chembl_id']:a for a in assay_page.get('assays',[])})
        # Fixed first page, explicitly reported diagnostic sampling, not exhaustive activity evidence.
        for row in page.get('activities',[]):
            if row.get('data_validity_comment') or row.get('potential_duplicate'):rejected['flagged_record']+=1;continue
            if row.get('standard_relation')!='=' or row.get('standard_units')!='nM' or row.get('standard_type') not in {'Ki','Kd','IC50'}:continue
            try:value=float(row['standard_value'])
            except (TypeError,ValueError):continue
            if not 0<value<=1000:continue
            aid=row['assay_chembl_id']
            assay=assays.get(aid,{})
            if assay.get('confidence_score')!=9 or assay.get('target_chembl_id')!=target['target_chembl_id']:
                rejected['assay_identity_or_confidence']+=1;continue
            mid=row['molecule_chembl_id']
            if mid in seen:continue
            smiles=row.get('canonical_smiles')
            if not smiles:rejected['missing_structure']+=1;continue
            seen.add(mid);rows.append(dict(molecule_id=mid,smiles=smiles,activity_id=row['activity_id'],assay_id=aid,
                target_id=target['target_chembl_id'],type=row['standard_type'],value=value,units='nM',relation='=',
                document_id=row.get('document_chembl_id'),source_url=base+f'activity/{row["activity_id"]}.json'))
            if len(rows)>=maximum:break
        if len(rows)>=maximum:break
    return rows,dict(targets=[t['target_chembl_id'] for t in eligible],rejected=dict(rejected),
        sampling=f'First 200 matching records per target, at most {maximum} distinct molecules; no potency aggregation across assays')


def candidates(smiles, molecule_id):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from .conformer_artifacts import _feature_factory,_feature_template,_features_from_template,_heavy_coordinates
    from .chemical_companion import _feature_direction
    mol=Chem.MolFromSmiles(smiles)
    if mol is None or len(Chem.GetMolFrags(mol))!=1:raise ValueError('Unresolved or multicomponent chemical structure')
    mol=Chem.AddHs(mol);parameters=AllChem.ETKDGv3();parameters.randomSeed=20260922;parameters.numThreads=1;parameters.pruneRmsThresh=.5
    ids=list(AllChem.EmbedMultipleConfs(mol,numConfs=8,params=parameters))
    if not ids:raise ValueError('Independent conformer preparation failed')
    for index in ids:
        one=Chem.Mol(mol);one.RemoveAllConformers();one.AddConformer(mol.GetConformer(index),assignId=True)
        one=Chem.RemoveHs(one)  # Match the library heavy-neighbor direction convention.
        template=_feature_template(one,_feature_factory());points,types,_=_features_from_template(one,template)
        directions,kinds=zip(*[_feature_direction(one,family,atoms) for family,atoms in template])
        yield SimpleNamespace(molecule_id=molecule_id,conformer_id=str(index),shape_points=_heavy_coordinates(one)[0],
            feature_points=points,feature_types=types,feature_directions=np.array(directions),feature_kinds=np.array(kinds))


def run(adopted, output, fetch=None):
    from .guided_filters import GuidedFilter,adaptive_stage,same_pose_gaussian
    from .necessary_conditions import NecessaryConditions
    from .pose_feasibility import possible_seed_mask
    from .preselection_full import pose_representatives
    from .full_library_screen import POSE_PARAMETERS
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    protocol=dict(classification='exploratory',hypothesis='The adopted funnel retains independently sourced binding controls',
        preparation='ChEMBL canonical graph, eight deterministic ETKDGv3 conformers, heavy-neighbor feature directions; not the library MOL2 ingestion route',
        metrics='Molecule-level cumulative stage retention, union over templates and conformers; assays remain individually cited',
        acceptance='Diagnostic only; zero controls or preparation failures cannot establish retrieval recall')
    _atomic_json(output/'protocol.json',protocol)
    report=dict(status='running',protocol=protocol,records=[],failures=[],catalog_count_contribution=0)
    try:
        records,sampling=binding_records(adopted['target'],fetch or PublicCache(output/'raw').get)
        from rdkit import Chem
        def graph(smiles):
            mol=Chem.MolFromSmiles(smiles)
            return Chem.MolToSmiles(mol,isomericSmiles=False) if mol else None
        reference_graphs={graph(s) for s in adopted.get('reference_smiles',[])}
        report['reference_overlap_excluded']=[r['molecule_id'] for r in records if graph(r['smiles']) in reference_graphs]
        records=[r for r in records if graph(r['smiles']) not in reference_graphs]
        report['sampling']=sampling
        _,pocket=_load_query(Path(adopted['consensus_npz']));design=adopted['design'];order=adopted['anchor_order']
        lookup={a['anchor_id']:a['feature_index'] for a in adopted['anchors']};feature_indices=[lookup[a] for a in order]
        pocket=dict(pocket,anchor_feature_indices=np.array(sorted(feature_indices)))
        columns=[list(pocket['anchor_feature_indices']).index(i) for i in feature_indices]
        gate=GuidedFilter(dict(anchors=adopted['anchors']),pocket,design)
        bound=NecessaryConditions(pocket,feature_indices,'any',design['minimum_score'])
        templates=[_load_query(Path(q['query_npz']))[1] for q in adopted['templates']]
        for record in records:
            best=0
            try:
                for candidate in candidates(record['smiles'],record['molecule_id']):
                    stage=adaptive_stage(candidate,design['adaptive_coarse']);best=max(best,stage+1)
                    if stage<3 or not bound.check(candidate)[0] or not gate.geometry_possible(candidate):continue
                    best=max(best,5)
                    for template in templates:
                        seeds,_=prepare_seeds(candidate,template,**POSE_PARAMETERS)
                        possible=gate.geometry_seeds(candidate,seeds,possible_seed_mask(bound,candidate,seeds))
                        if not possible.any():continue
                        best=max(best,6);positions=np.flatnonzero(possible)
                        matrices=np.array([seeds[i].transform_matrix for i in positions]).reshape(-1,4,4)
                        possible[positions]=gate.pocket_mask(candidate.shape_points,matrices)
                        if not possible.any():continue
                        best=max(best,7);gs=np.full(len(seeds),np.nan);positions=np.flatnonzero(possible)
                        gs[positions]=same_pose_gaussian(template,candidate,[seeds[i] for i in positions])
                        rows,_=pose_representatives(candidate,seeds,possible,pocket,columns,.5,design,order,gs,candidate.shape_points)
                        if rows:best=8;break
                    if best==8:break
                report['records'].append(dict(record,last_stage=best))
            except Exception as exc:report['failures'].append(dict(molecule_id=record['molecule_id'],reason=type(exc).__name__,detail=str(exc)[:200]))
        labels=['prepared','size','feature_count','extent','anchor_geometry','seed_geometry','pocket','final_pose']
        denominator=len(records)
        report.update(status='complete',identified_controls=denominator,prepared_controls=len(report['records']),
            retention=[dict(stage=name,retained=sum(r['last_stage']>=i+1 for r in report['records']),
                fraction_over_identified=(sum(r['last_stage']>=i+1 for r in report['records'])/denominator if denominator else None)) for i,name in enumerate(labels)],
            acceptance='exploratory_measured' if denominator and not report['failures'] else 'not_validated_missing_or_failed_controls')
    except Exception as exc:report.update(status='unavailable',error=type(exc).__name__,acceptance='not_validated')
    _atomic_json(output/'report.json',report);return report
