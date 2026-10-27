"""Create separate AF3 protein-ligand inputs from a base JSON and a SMILES table."""
import argparse
import copy
import csv
import hashlib
import json
from pathlib import Path
import re


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path,value):
    Path(path).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')


def validate_smiles(smiles):
    from rdkit import Chem
    if not smiles or any(c.isspace() for c in smiles):
        raise ValueError('SMILES must be nonempty and contain no whitespace')
    mol=Chem.MolFromSmiles(smiles)
    if mol is None or mol.GetNumHeavyAtoms()==0:raise ValueError('Invalid SMILES: '+smiles)
    if len(Chem.GetMolFrags(mol))!=1:
        raise ValueError('Disconnected SMILES: supply the intended single ligand; salts/ions need explicit entities')
    # Preserve the literal input, including stereochemistry and formal charge.
    return mol.GetNumHeavyAtoms()


def validate_base(base,keep_ligands):
    if not isinstance(base,dict) or base.get('dialect')!='alphafold3':raise ValueError('AF3 dialect JSON required')
    if type(base.get('version')) is not int or base['version'] not in (1,2,3,4):raise ValueError('Unsupported AF3 version')
    def check_paths(value):
        if isinstance(value,dict):
            for k,v in value.items():
                if k.endswith('Path'):raise ValueError('External '+k+' must be inlined before packaging this batch')
                check_paths(v)
        elif isinstance(value,list):
            for item in value:check_paths(item)
    check_paths(base)
    used=set();proteins=0
    for entity in base.get('sequences',[]):
        if len(entity)!=1 or next(iter(entity)) not in ('protein','ligand'):
            raise ValueError('This helper supports protein chains and explicitly retained cofactors only')
        kind=next(iter(entity));value=entity[kind]
        if kind=='ligand' and not keep_ligands:
            raise ValueError('Base contains a ligand; remove it or explicitly use --keep-base-ligands for cofactors')
        ids=value.get('id');ids=ids if isinstance(ids,list) else [ids]
        if not ids:raise ValueError('Empty chain IDs')
        for chain in ids:
            if not isinstance(chain,str) or not re.fullmatch('[A-Z]+',chain) or chain in used:
                raise ValueError('Duplicate or invalid chain ID')
            used.add(chain)
        if kind=='protein':
            proteins+=len(ids)
            if not re.fullmatch('[ACDEFGHIKLMNPQRSTVWYX]+',value.get('sequence','')):
                raise ValueError('Supply a verified continuous uppercase amino-acid sequence; no placeholders')
        elif 'smiles' in value:validate_smiles(value['smiles'])
    if not proteins:raise ValueError('At least one protein chain required')
    return used


def next_ids(used,count):
    ids=[]
    for n in range(1,10000):
        number=n;chain=''
        while number:
            number,r=divmod(number-1,26);chain=chr(65+r)+chain
        if chain not in used:ids.append(chain)
        if len(ids)==count:return ids
    raise ValueError('Too many requested entities')


def prepare(base_path,table,out,seeds=None,copies=1,keep_ligands=False):
    base_path=Path(base_path);table=Path(table);out=Path(out)
    if out.exists():raise ValueError('Use a new output directory')
    if copies<1:raise ValueError('Positive ligand copy count required')
    base=json.loads(base_path.read_text(encoding='utf-8-sig'));used=validate_base(base,keep_ligands)
    seeds=seeds if seeds is not None else base.get('modelSeeds')
    if not isinstance(seeds,list) or not seeds or any(type(s) is not int or s<0 for s in seeds) or len(set(seeds))!=len(seeds):
        raise ValueError('Provide distinct nonnegative integer model seeds')
    ligand_ids=next_ids(used,copies);rows=[];seen=set()
    with table.open(encoding='utf-8-sig',newline='') as stream:
        reader=csv.DictReader(stream,delimiter='\t' if table.suffix.lower() in ('.tsv','.smi') else ',')
        if not reader.fieldnames or not {'id','smiles'}<=set(reader.fieldnames):
            raise ValueError('Table must have id and smiles columns (CSV or TSV with header)')
        for row in reader:
            name=row.get('id','');smiles=row.get('smiles','')
            if not re.fullmatch('[A-Za-z0-9][A-Za-z0-9_-]*',name):raise ValueError('Use safe nonempty ligand IDs')
            if name.lower() in seen:raise ValueError('Duplicate case-insensitive ligand ID: '+name)
            seen.add(name.lower());heavy=validate_smiles(smiles);rows.append((name,smiles,heavy))
    if not rows:raise ValueError('No ligand rows')
    (out/'inputs').mkdir(parents=True)
    jobs=[]
    for name,smiles,heavy in rows:
        job=copy.deepcopy(base);job['name']=name;job['modelSeeds']=seeds
        job['sequences'].append({'ligand':{'id':ligand_ids[0] if copies==1 else ligand_ids,'smiles':smiles}})
        path=out/'inputs'/(name+'.json');save(path,job)
        jobs.append(dict(id=name,input='inputs/'+path.name,sha256=sha(path),smiles=smiles,heavy_atoms=heavy))
    manifest=dict(status='prepared_not_predicted',base_sha256=sha(base_path),table_sha256=sha(table),
        modelSeeds=seeds,ligand_chain_ids=ligand_ids,keep_base_ligands=keep_ligands,jobs=jobs,
        limitations=['Each row is a separate complex prediction, not a ligand competition mixture.',
                    'No protonation or stereochemistry enumeration; literal SMILES retained.',
                    'RDKit parse success does not guarantee AF3 conformer generation or binding.',
                    'AF3 model is not an experimental co-crystal and confidence is not affinity.'])
    save(out/'manifest.json',manifest)
    print(f'Prepared {len(jobs)} jobs, {len(seeds)} seeds/job in {out}')
    return manifest


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--base',type=Path,required=True);p.add_argument('--ligands',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--seeds',type=int,nargs='+')
    p.add_argument('--ligand-copies',type=int,default=1);p.add_argument('--keep-base-ligands',action='store_true')
    a=p.parse_args();prepare(a.base,a.ligands,a.output,a.seeds,a.ligand_copies,a.keep_base_ligands)


if __name__=='__main__':main()
