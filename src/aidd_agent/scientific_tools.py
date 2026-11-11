"""Project-scoped coordinators for scientific quality, evaluation and AF3 batches."""
import json
from pathlib import Path

from .expanded_wee1 import fingerprint

FIELDS={
    'pose_qc': {'ligand_attachment_id','set_id','count','receptor_attachment_id','receptor_chains','exclude_residues','expected_smiles','required_contacts'},
    'screening_benchmark': {'reference_set_id','stage_set_ids','labels_attachment_id','top_k','minimum_recovery','protocol_note'},
    'af3_batch': {'name','protein_sequences','protein_task_id','protein_copies','seeds','seed_count','seed_start','ligands',
                  'ligand_attachment_id','set_id','count','shared_msa'},
}
DESCRIPTIONS={
    'pose_qc': 'Queue independent pose geometry/chemistry checks. ligand_attachment_id is one posed SDF/MOL2; optional receptor_attachment_id PDB/mmCIF, receptor_chains (default all supplied), exclude_residues [[chain,residue,name]], expected_smiles and required_contacts [{chain,residue,distance}]. Includes other protein chains, waters and cofactors present in the supplied assembly. Does not infer biological assembly. Pass means available checks passed, not affinity or experimental validation.',
    'screening_benchmark': 'Evaluate saved rankings without a new library scan. reference_set_id defines the complete reference universe; stage_set_ids is an ordered list of saved candidate sets. Optional labels_attachment_id CSV molecule_id,label (0/1), top_k (default 100), minimum_recovery (explicit 0..1), protocol_note. Reports positive losses, exact-reference recovery and enrichment only with complete labels. Same library and comparable target/query/assay protocols required; unlabelled molecules are never negative.',
    'af3_batch': 'Queue independent protein-ligand predictions with shared compatible protein MSA/templates. Required name and literal protein_sequences; optional protein_copies, seeds OR seed_count/seed_start, shared_msa (default true). Exactly one ligand source: ligands [{name,smiles}] copied from user; ligand_attachment_id SMILES table/SDF/MOL2; OR set_id and count (1..1000) from local candidate set. Each ligand gets a separate complex with all protein chains. Sequential GPU allocation, child receipts, failed-child resume, all-sample confidence and pose QC. Never interpret scores as affinity. Uncertain interrupted child dispatch requires operator reconciliation; do not duplicate tasks.',
}


def handle(app,sid,kind,args,request=''):
    from .chat_uploads import attachment
    from .candidate_sets import get
    from .research_tasks import queue
    from .af3_inputs import validate,literal_check
    if set(args)-FIELDS[kind]:raise ValueError('Unknown '+kind+' fields')
    root=Path(app.context_for(sid)['project_root']);values=dict(args);inputs={}
    def attach(identifier):
        item=attachment(app,sid,identifier);inputs[item['path']]=item['sha256'];return item['path']
    def candidate(identifier):
        item=get(root,identifier);inputs.update(item['outputs'])
        inputs.update(fingerprint([root/'candidate-sets'/identifier/'manifest.json']))
        return item
    if kind=='pose_qc':
        if ('ligand_attachment_id' in args)==('set_id' in args):raise ValueError('QC requires an attached ligand OR candidate set')
        if args.get('set_id'):
            candidate(args['set_id'])
            if type(args.get('count')) is not int or not 1<=args['count']<=1000:raise ValueError('Specify 1..1000 candidate poses for QC')
        else:values['ligand_path']=attach(args['ligand_attachment_id'])
        if args.get('receptor_attachment_id'):values['receptor_path']=attach(args['receptor_attachment_id'])
        if args.get('expected_smiles') and args['expected_smiles'] not in request:raise ValueError('Expected SMILES must come from the user')
    elif kind=='screening_benchmark':
        reference=candidate(args['reference_set_id']);values['reference_path']=reference['members']
        ids=args['stage_set_ids']
        if not isinstance(ids,list) or not 1<=len(ids)<=20:raise ValueError('Use 1..20 stage candidate sets')
        values['stages']=[]
        for identifier in ids:
            item=candidate(identifier)
            if (item['batch'],item['catalog_sha256'])!=(reference['batch'],reference['catalog_sha256']):raise ValueError('Benchmark library versions differ')
            values['stages'].append(dict(name=item['name'],path=item['members']))
        if args.get('labels_attachment_id'):values['labels_path']=attach(args['labels_attachment_id'])
    else:
        if args.get('protein_task_id'):
            from .library_acceptance import read,sha
            from .project_context import ensure_within
            if 'protein_sequences' in args:raise ValueError('Choose literal protein sequences OR an existing protein task')
            job=app.task(sid,args['protein_task_id'])
            if job['status']!='complete' or not job.get('report'):raise ValueError('Protein source task must be complete')
            report_path=ensure_within(Path(job['report']),root);report=read(report_path);inputs.update(fingerprint([report_path]))
            paths=[]
            for step in report.get('steps',{}).values():
                result=step.get('result',{})
                if step.get('status')=='complete' and result.get('input_json'):
                    path=ensure_within(Path(result['input_json']),root)
                    if sha(path)!=result.get('input_sha256'):raise ValueError('Source protein input changed since preparation')
                    paths.append(path)
            if len(paths)!=1:raise ValueError('Protein source must have one unambiguous prepared AF3 input')
            inputs.update(fingerprint(paths));data=read(paths[0])
            sequences=[]
            for entity in data['sequences']:
                if 'protein' in entity:
                    p=entity['protein'];copies=len(p['id']) if isinstance(p['id'],list) else 1
                    if p.get('modifications'):raise ValueError('Modified proteins are not supported by this batch adapter')
                    sequences.extend([p['sequence']]*copies)
            values['protein_sequences']=sequences
            if args.get('protein_copies',1)!=1:raise ValueError('Existing task already specifies assembly copies; do not multiply it again')
        if not isinstance(args.get('name'),str) or not 1<=len(args['name'])<=160:raise ValueError('Name the AF3 batch')
        if type(args.get('shared_msa',True)) is not bool:raise ValueError('shared_msa must be boolean')
        validate(dict(values,ligand_smiles=['placeholder']),{})
        literal_check({k:v for k,v in values.items() if k!='protein_sequences'} if args.get('protein_task_id') else args,request)
        if sum(key in args for key in ('ligands','ligand_attachment_id','set_id'))!=1:raise ValueError('Choose exactly one ligand source')
        if 'ligands' in args:
            rows=args['ligands']
            if not isinstance(rows,list) or not 1<=len(rows)<=1000:raise ValueError('Use 1..1000 separate ligands')
            from .af3_analysis import ligand_evidence
            for item in rows:
                if not isinstance(item,dict) or set(item)-{'name','smiles'} or not isinstance(item.get('smiles'),str) or item['smiles'] not in request:
                    raise ValueError('Ligand SMILES must be copied literally from the user')
                ligand_evidence(item['smiles'])
        elif 'ligand_attachment_id' in args:
            values['input_path']=attach(args['ligand_attachment_id'])
        else:
            candidate(args['set_id']);count=args.get('count')
            if type(count) is not int or not 1<=count<=1000:raise ValueError('Specify 1..1000 candidate molecules for AF3')
    return queue(app,sid,kind,values,inputs)


def resolve_batch(args,directory,root):
    """Read local chemistry only in the executor; freeze exact per-ligand input."""
    from rdkit import Chem
    args=dict(args)
    if args.get('set_id'):
        import csv
        from .candidate_sets import export
        from .library_acceptance import read
        from .screening_selection import check_hashes
        export_dir=Path(directory)/'candidate-inputs';export_dir.mkdir(parents=True,exist_ok=True)
        result=export(dict(set_id=args['set_id'],count=args['count']),export_dir,root)
        manifest=read(Path(result['export_manifest']));check_hashes(manifest['outputs'])
        if result['shortfall']:raise ValueError('Candidate set has fewer molecules than requested')
        rows=[]
        with Path(result['names']).open(encoding='utf-8-sig',newline='') as stream:
            for row in csv.DictReader(stream):
                if row['qc_status']!='ok':raise ValueError('Candidate chemistry requires review before AF3: '+row['molecule_id'])
                path=Path(result['export_manifest']).parent/row['original_mol2']
                mol=Chem.MolFromMol2File(str(path),removeHs=True)
                if mol is None:raise ValueError('Candidate MOL2 failed strict chemistry parsing')
                rows.append(dict(name=row['name'],molecule_id=row['molecule_id'],smiles=Chem.MolToSmiles(mol,isomericSmiles=True)))
        args['ligands']=rows
    elif args.get('input_path'):
        from .molecule_preparation import records
        rows=[]
        for name,original,mol in records(args):
            if len(rows)>=1000:raise ValueError('Split batches larger than 1000 ligands')
            if mol is None:raise ValueError('Unparseable ligand record: '+name)
            rows.append(dict(name=name,smiles=original or Chem.MolToSmiles(Chem.RemoveHs(mol),isomericSmiles=True)))
        if not rows:raise ValueError('No ligands in uploaded batch')
        args['ligands']=rows
    return args
