"""Trusted Chat coordinator for explicitly selected receptor/reference searches."""
import hashlib
import json
from pathlib import Path
import re
import shutil
import time
import uuid

from .gaussian_batch import _atomic_json
from .library_acceptance import read, sha
from .project_context import ensure_within

OPERATIONS = {
    'prepare': {'name','receptor','pdb_id','reference_ligand','ligand_instance','receptor_chains','keep_residues'},
    'start': {'name','receptor','pdb_id','reference_ligand','ligand_instance','receptor_chains','keep_residues','required_contacts','policy'},
    'run': {'task_id','required_contacts','policy'},
}


def validate(args):
    if not isinstance(args,dict) or not isinstance(args.get('operation'),str) or args.get('operation') not in OPERATIONS:
        raise ValueError('target_search operation must be prepare, start or run')
    operation=args['operation']
    if set(args)-OPERATIONS[operation]-{'operation'}: raise ValueError('Unknown target search arguments')
    if operation!='run':
        if not isinstance(args.get('name'),str) or not re.fullmatch('[A-Za-z0-9][A-Za-z0-9_-]{0,63}',args['name']):
            raise ValueError('Use a short target/query name')
        if ('receptor' in args)+('pdb_id' in args)!=1: raise ValueError('Supply receptor file OR pdb_id')
        if 'pdb_id' in args and (not isinstance(args['pdb_id'],str) or not re.fullmatch('[0-9][A-Za-z0-9]{3}',args['pdb_id'])):
            raise ValueError('Invalid PDB ID')
        if ('reference_ligand' in args)+('ligand_instance' in args)!=1:
            raise ValueError('Supply a receptor-frame reference_ligand file OR explicit ligand_instance')
        chains=args.get('receptor_chains')
        if not isinstance(chains,list) or not chains or len(chains)>26 or any(not isinstance(c,str) or not re.fullmatch('[A-Za-z0-9_]{1,8}',c) for c in chains) or len(set(chains))!=len(chains):
            raise ValueError('Explicit unique receptor_chains are required')
        if 'ligand_instance' in args:
            row=args['ligand_instance']
            if not isinstance(row,dict) or set(row)!={'ccd','chain','residue'} or any(not isinstance(v,str) or not v for v in row.values()):
                raise ValueError('ligand_instance requires ccd, chain and author residue strings')
            if not re.fullmatch('[A-Z0-9]{1,5}',row['ccd']):raise ValueError('Invalid ligand CCD')
    elif not isinstance(args.get('task_id'),str):raise ValueError('run requires an existing target query task_id')
    for key in ('keep_residues','required_contacts'):
        rows=args.get(key,[])
        if not isinstance(rows,list) or len(rows)>20:raise ValueError(key+' must be a list of at most 20 residue records')
        for row in rows:
            allowed={'chain','residue','name'} if key=='keep_residues' else {'chain','residue','residue_name','feature_class','distance'}
            required={'chain','residue','name'} if key=='keep_residues' else {'chain','residue'}
            if not isinstance(row,dict) or not required<=set(row) or set(row)-allowed:
                raise ValueError('Invalid '+key+' record')
            if any(not isinstance(row[k],str) or not row[k] for k in set(row)-{'distance'}):raise ValueError('Residue identities must be strings')
            if 'distance' in row and (type(row['distance']) not in (int,float) or not 0<row['distance']<=8):raise ValueError('Contact distance must be in (0,8] A')
    if 'policy' in args and not isinstance(args['policy'],dict):raise ValueError('policy must be an object')


def import_file(app, project, supplied, destination, suffixes):
    if not isinstance(supplied,str) or not Path(supplied).is_absolute():raise ValueError('Use exact absolute input paths on the Chat host')
    path=Path(supplied).resolve()
    if not any(path.is_relative_to(root) for root in [project.resolve(),*app.result_import_roots]):
        raise ValueError('Input is outside the Project and configured --input-root directories')
    if path.suffix.lower() not in suffixes or not path.is_file():raise ValueError('Unsupported or missing molecular input file')
    if not 0<path.stat().st_size<=128*1024*1024:raise ValueError('Molecular input must be 1 byte..128 MiB')
    target=destination.with_suffix(path.suffix.lower()); before=sha(path)
    shutil.copyfile(path,target)
    if sha(target)!=before or sha(path)!=before:raise ValueError('Molecular input changed during import')
    return str(target.resolve())


def handle(app, sid, args):
    from .prompt_workflow import project_root, create_plan
    validate(args);app.session(sid)
    ctx=app.context;project=project_root(Path(ctx['db']),ctx['user_id'],ctx['project_id'])
    with app.lock:
        identifier=uuid.uuid4().hex[:24]
        folder=project/'target-search'/'requests'/identifier;folder.mkdir(parents=True)
        request=dict(args)
        if args['operation']=='run':
            task=app.task(sid,args['task_id'])
            if task['status']!='complete':raise ValueError('Wait for the target query task to complete')
            report=read(ensure_within(Path(task['report']),project))
            sources=[row['result']['report'] for row in report['steps'].values() if row['action']=='target_query_prepare' and row['status']=='complete']
            if len(sources)!=1:raise ValueError('Choose a task containing one completed target query')
            query=ensure_within(Path(sources[0]),project)
            request.update(query_report=str(query),query_report_sha256=sha(query))
        else:
            for field,suffixes in (('receptor',{'.pdb','.cif','.mmcif'}),('reference_ligand',{'.sdf','.mol2'})):
                if field in args:request[field]=import_file(app,project,args[field],folder/field,suffixes)
        request['input_hashes']={str(p.resolve()):sha(p) for p in folder.iterdir() if p.is_file()}
        path=folder/'request.json';_atomic_json(path,request)
        params=dict(request_id=identifier,request_sha256=sha(path))
        steps=[]
        if args['operation']!='run':steps.append(dict(id='query',action='target_query_prepare',params=params))
        if args['operation']!='prepare':steps.append(dict(id='screen',action='target_library_search',params=dict(params,**({'query_step':'query'} if steps else {}))))
        plan=dict(version=1,summary='General receptor/reference '+args['operation'],clarifications=[],steps=steps)
        plan_path=create_plan(Path(ctx['db']),ctx['user_id'],ctx['project_id'],json.dumps(args),local_plan=plan)
        jid=uuid.uuid4().hex[:16]
        with app.connect() as db:
            db.execute('INSERT INTO jobs(id,session,request,provider,profile,status,plan,report,log,created) VALUES(?,?,?,?,?,?,?,?,?,?)',
                (jid,sid,plan['summary'],'local-coordinator','','queued',str(plan_path),str(plan_path.parent/'execution/report.json'),str(app.root/(jid+'.log')),time.time()))
        return dict(status='queued',task_id=jid,operation=args['operation'],plan=str(plan_path),
            full_library_requested=args['operation']!='prepare',
            note='Target-independent query; full catalog scan uses finite poses, not exhaustive conformations. Compute must be enabled.')


def request_file(execution, params):
    project=execution.parent.parent.parent
    path=ensure_within(project/'target-search'/'requests'/params['request_id']/'request.json',project)
    if sha(path)!=params['request_sha256']:raise ValueError('Target search request changed')
    return path


def execute(action, params, directory, execution, cfg, allow_compute, services):
    from .screening_selection import check_hashes
    from .target_query import prepare, design
    request=read(request_file(execution,params))
    if request['input_hashes']:check_hashes(request['input_hashes'])
    if action=='target_query_prepare':
        values=dict(request)
        if 'pdb_id' in values:
            values['receptor']=services.pdb_fetch(values['pdb_id'],directory/'source')['path']
        if 'ligand_instance' in values:
            from .chemistry_prep import fetch_ccd, enumerate_ligand_instances, _ccd_molecule
            from rdkit import Chem
            source=Path(values['receptor'])
            if source.suffix.lower()=='.pdb':
                from Bio.PDB import PDBParser, MMCIFIO
                writer=MMCIFIO();writer.set_structure(PDBParser(QUIET=True).get_structure('query',str(source)))
                source=directory/'source.cif';writer.save(str(source))
            item=values['ligand_instance']
            matches=[r for r in enumerate_ligand_instances(source,[item['ccd']]) if r['chain_id']==item['chain'] and r['residue_number']==item['residue']]
            if len(matches)!=1 or matches[0]['altloc'] or matches[0]['insertion_code'] or matches[0]['model']!='1':
                raise ValueError('Ligand instance is absent or ambiguous; supply an explicitly prepared reference SDF')
            ccd=fetch_ccd(item['ccd'],directory/'ccd')
            mol=_ccd_molecule(ccd,matches[0]['atoms'])
            ligand=directory/'reference.sdf'
            with Chem.SDWriter(str(ligand)) as writer:writer.write(mol)
            values['reference_ligand']=str(ligand.resolve())
        result=prepare(values,directory/'query')
        return dict(status=result['kind'],report=str(directory/'query/report.json'),target=result['target'],anchors=len(result['query']['anchors']))
    if not allow_compute:
        from .workflow_errors import Blocked
        raise Blocked('Enable --allow-compute for the requested full-library receptor search')
    if not cfg.get('search'):
        from .workflow_errors import Blocked
        raise Blocked('Configure the target-independent search.batch/workers/refine_chunk runtime')
    if 'query_step' in params:source=execution/params['query_step']/'query/report.json'
    else:
        source=ensure_within(Path(request['query_report']),execution.parent.parent.parent)
        if sha(source)!=request['query_report_sha256']:raise ValueError('Source target query changed')
    design(source,directory/'design',request.get('required_contacts',[]),request.get('policy'))
    from .guided_workflow import execute as guided
    result=guided('guided_funnel',directory/'design/report.json',directory/'search',{},cfg,True)
    return dict(status='target_library_search',report=str(directory/'search/report.json'),full_coverage=result.get('full_coverage'),
                outputs=result.get('outputs'),note='All catalog conformers considered under the explicit finite-pose policy; no candidate Top-K')
