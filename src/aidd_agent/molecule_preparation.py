"""Bounded, reproducible chemical input preparation with per-record failures."""
import csv
import json
from pathlib import Path

from .expanded_wee1 import fingerprint


def handle(app,sid,args):
    from .chat_uploads import attachment
    from .research_tasks import queue
    if set(args)-{'attachment_id','smiles','name','conformers','seed','smiles_column','name_column','allow_unspecified_stereo'}:
        raise ValueError('Unknown molecule preparation fields')
    if ('attachment_id' in args)==('smiles' in args):raise ValueError('Provide an uploaded structure/table OR literal SMILES')
    values=dict(args);inputs={}
    if 'attachment_id' in args:
        row=attachment(app,sid,args['attachment_id']); values['input_path']=row['path'];inputs={row['path']:row['sha256']}
    elif not isinstance(args['smiles'],str) or not 1<=len(args['smiles'])<=6000:raise ValueError('Literal SMILES must be 1..6000 characters')
    for key,low,high in [('conformers',1,100),('seed',0,2147483647)]:
        val=args.get(key,1 if key=='conformers' else 47)
        if type(val) is not int or not low<=val<=high:raise ValueError('Invalid '+key)
    if type(args.get('allow_unspecified_stereo',False)) is not bool:raise ValueError('allow_unspecified_stereo must be boolean')
    return queue(app,sid,'molecule_prepare',values,inputs)


def records(args):
    from rdkit import Chem
    if 'smiles' in args:
        yield args.get('name','molecule'),args['smiles'],Chem.MolFromSmiles(args['smiles']);return
    path=Path(args['input_path']);suffix=path.suffix.lower()
    if suffix in {'.sdf','.mol'}:
        values=Chem.SDMolSupplier(str(path),removeHs=False) if suffix=='.sdf' else [Chem.MolFromMolFile(str(path),removeHs=False)]
        for i,mol in enumerate(values):yield mol.GetProp('_Name') if mol and mol.HasProp('_Name') else f'record-{i+1}',None,mol
    elif suffix=='.mol2':
        from .mol2 import iter_mol2_blocks
        for i,block in iter_mol2_blocks(path):
            mol=Chem.MolFromMol2Block(block,removeHs=False,sanitize=True)
            yield mol.GetProp('_Name') if mol and mol.HasProp('_Name') else f'record-{i+1}',None,mol
    elif suffix in {'.csv','.tsv'}:
        with path.open(encoding='utf-8-sig',newline='') as stream:
            rows=csv.DictReader(stream,delimiter='\t' if suffix=='.tsv' else ',')
            column=args.get('smiles_column','smiles')
            if column not in (rows.fieldnames or []):raise ValueError('Missing SMILES column: '+column)
            for i,row in enumerate(rows):
                value=row[column] or '';yield row.get(args.get('name_column','name')) or f'record-{i+1}',value,Chem.MolFromSmiles(value)
    elif suffix in {'.smi','.smiles','.txt'}:
        with path.open(encoding='utf-8-sig') as stream:
            for i,line in enumerate(stream):
                if not line.strip():continue
                parts=line.strip().split(maxsplit=1);yield parts[1] if len(parts)>1 else f'record-{i+1}',parts[0],Chem.MolFromSmiles(parts[0])
    else:raise ValueError('Preparation supports SMILES, CSV/TSV, SDF/MOL and MOL2')


def execute(args,directory):
    from rdkit import Chem,rdBase
    from rdkit.Chem import AllChem
    from .ligand_library import prepare_ligand
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=True)
    prepared=[];failures=[];queries=[];seen={};total=0
    count=args.get('conformers',1);seed=args.get('seed',47)
    for index,(name,original,mol) in enumerate(records(args)):
        if index>=5000 or (index+1)*count>10000:raise ValueError('Preparation limit: 5000 records and 10000 requested conformers; split this input')
        total+=1;record_id=f'mol-{index+1:06d}'
        try:
            if mol is None or mol.GetNumHeavyAtoms()<1:raise ValueError('Chemical graph could not be parsed')
            if len(Chem.GetMolFrags(mol))!=1:raise ValueError('Disconnected fragments require an explicit salt/component decision')
            if mol.GetNumHeavyAtoms()>500:raise ValueError('More than 500 heavy atoms; use a dedicated preparation workflow')
            Chem.AssignStereochemistry(mol,cleanIt=True,force=True)
            unassigned=[str(info.centeredOn) for info in Chem.FindPotentialStereo(mol) if str(info.specified)=='Unspecified']
            if unassigned and not args.get('allow_unspecified_stereo',False):raise ValueError('Unspecified stereochemistry; explicitly allow it or provide the intended isomer')
            canonical=Chem.MolToSmiles(Chem.RemoveHs(mol),isomericSmiles=True)
            duplicate=seen.get(canonical);seen.setdefault(canonical,record_id)
            preserved=mol.GetNumConformers()==1 and mol.GetConformer().Is3D()
            if preserved:
                ids=[mol.GetConformer().GetId()];energies=[None]
            else:
                mol=Chem.AddHs(mol);mol.RemoveAllConformers()
                parameters=AllChem.ETKDGv3();parameters.randomSeed=seed;parameters.numThreads=1
                parameters.pruneRmsThresh=.5;parameters.enforceChirality=True
                ids=list(AllChem.EmbedMultipleConfs(mol,numConfs=count,params=parameters))
                if not ids:raise ValueError('3D embedding failed')
                # Optimize only generated coordinates. Original 3D poses stay unchanged.
                if AllChem.MMFFHasAllMoleculeParams(mol):
                    optimized=AllChem.MMFFOptimizeMoleculeConfs(mol,numThreads=1,maxIters=200)
                    energies=[dict(converged=status==0,energy_kcal_mol=energy) for status,energy in optimized]
                else:energies=[dict(converged=False,reason='MMFF parameters unavailable') for _ in ids]
            folder=directory/record_id;folder.mkdir(exist_ok=True);outputs=[]
            for ordinal,(cid,energy) in enumerate(zip(ids,energies),1):
                sdf=folder/f'conformer-{ordinal:03d}.sdf'
                mol.SetProp('_Name',name)
                with Chem.SDWriter(str(sdf)) as writer:writer.write(mol,confId=cid)
                item=dict(conformer=ordinal,sdf=str(sdf),optimization=energy)
                # Only conformers meeting the search engine's geometric contract become queries.
                query=sdf.with_suffix('.npz');qid=record_id+f':conf-{ordinal:03d}'
                try:
                    prepare_ligand(sdf,query,qid)
                    queries.append(dict(query_id=qid,query_npz=str(query),sdf=str(sdf),name=name,record_id=record_id))
                except ValueError as exc:item['search_query_error']=str(exc)
                outputs.append(item)
            prepared.append(dict(record_id=record_id,name=name,input_smiles=original,canonical_isomeric_smiles=canonical,
                formal_charge=Chem.GetFormalCharge(mol),unspecified_stereo=unassigned,duplicate_of=duplicate,
                coordinates='original_3d_preserved' if preserved else 'generated_ETKDGv3',conformers=outputs))
        except (ValueError,RuntimeError) as exc:
            failures.append(dict(record_id=record_id,name=name,input_smiles=original,error=str(exc)))
    if total==0:raise ValueError('No molecular records found')
    details=directory/'molecules.json';details.write_text(json.dumps(prepared,indent=2),encoding='utf-8')
    failure_path=directory/'failures.csv'
    with failure_path.open('w',encoding='utf-8',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=['record_id','name','input_smiles','error']);writer.writeheader();writer.writerows(failures)
    return dict(kind='molecule_prepare',status='complete',input_records=total,prepared_molecules=len(prepared),
        failed_molecules=len(failures),prepared_queries=queries,details=str(details),failures=str(failure_path),
        policy=dict(seed=seed,requested_conformers=count,rdkit=rdBase.rdkitVersion,stereochemistry='preserve specified; no enumeration',
                    protonation='unchanged',tautomers='unchanged',salts='reject disconnected; never silently strip'),
        sources=fingerprint([Path(args['input_path'])]) if args.get('input_path') else {},
        outputs=fingerprint([details,failure_path,*directory.glob('mol-*/*.sdf'),*directory.glob('mol-*/*.npz'),*directory.glob('mol-*/*.manifest.json')]))
