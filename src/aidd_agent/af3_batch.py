"""Sealed, sequential AF3 ligand batches with compatible shared protein features."""
import json
import tempfile
import time
from pathlib import Path

from .af3_inputs import resolve, seeds, validate
from .af3_shared import ready, transfer, load_data_output
from .expanded_wee1 import fingerprint
from .gaussian_batch import _atomic_json
from .library_acceptance import read, sha
from .prediction import (write_alphafold3_input,load_model_profile,validate_af3_installation,
                         prediction_command,inspect_alphafold3_output)
from .screening_selection import check_hashes
from .workflow_errors import Blocked


def sample_summary(output):
    """Retain every sample confidence and structure hash; never infer missing seeds."""
    import re
    root=Path(output).resolve();rows=[];seen=set()
    for path in sorted(root.rglob('*summary_confidences.json')):
        if not re.fullmatch(r'seed-\d+_sample-\d+',path.parent.name):continue
        if path.parent in seen:raise ValueError('Ambiguous summary confidence files in sample '+path.parent.name)
        seen.add(path.parent)
        from .project_context import ensure_within
        path=ensure_within(path,root)
        prefix=path.name.removesuffix('summary_confidences.json')
        model=ensure_within(path.with_name(prefix+'model.cif'),root)
        row=dict(sample=path.parent.name,confidence=read(path),confidence_path=str(path),confidence_sha256=sha(path),
                 model=str(model) if model.is_file() else None,model_sha256=sha(model) if model.is_file() else None)
        rows.append(row)
    ranking_files=list(Path(output).rglob('ranking_scores.csv'))+list(Path(output).rglob('*_ranking_scores.csv'))
    ranking=[]
    if len(ranking_files)==1:
        import csv
        with ranking_files[0].open(encoding='utf-8',newline='') as stream:ranking=list(csv.DictReader(stream))
    finite=[r for r in ranking if r.get('ranking_score') is not None]
    import math
    if any(not math.isfinite(float(r['ranking_score'])) for r in finite):raise ValueError('Nonfinite AF3 ranking score')
    best=max(finite,key=lambda r:float(r['ranking_score'])) if finite else None
    ranges={}
    for key in ('iptm','ptm','ranking_score'):
        values=[r['confidence'][key] for r in rows if type(r['confidence'].get(key)) in (float,int)]
        if values:ranges[key]=dict(min=min(values),max=max(values),count=len(values))
    return dict(samples=rows,sample_count=len(rows),ranking_rows=ranking,best_recorded_sample=best,
        sample_model_count=sum(bool(r['model']) for r in rows),
        missing_sample_models=[r['sample'] for r in rows if not r['model']],
        metric_ranges=ranges,ranking_source=str(ranking_files[0]) if len(ranking_files)==1 else None,
        sample_coverage='recorded_outputs_only',affinity_ranking='not_supported')


def validate_sample_coverage(summary,input_data):
    expected=set(input_data['modelSeeds'])
    rows=summary['ranking_rows']
    if not rows:
        summary['seed_coverage']=dict(status='not_verified',expected=sorted(expected),reason='No unambiguous ranking table')
        return
    if any('seed' not in row for row in rows):raise ValueError('AF3 ranking lacks seed identifiers')
    observed={int(row['seed']) for row in rows}
    if observed!=expected:raise ValueError('AF3 output seed coverage differs from the requested seeds')
    identities=[(row['seed'],row.get('sample')) for row in rows]
    if len(identities)!=len(set(identities)):raise ValueError('Duplicate AF3 ranking sample identifiers')
    import re
    actual=set()
    for sample in summary['samples']:
        match=re.fullmatch(r'seed-(\d+)_sample-(\d+)',sample['sample'])
        if match and sample.get('model'):actual.add((int(match[1]),int(match[2])))
    expected_samples={(int(row['seed']),int(row['sample'])) for row in rows if row.get('sample') is not None}
    if actual and actual!=expected_samples:raise ValueError('AF3 ranked sample structures/confidence are incomplete')
    summary['seed_coverage']=dict(status='verified_ranked_samples' if actual else 'ranking_only',
                                 expected=sorted(expected),observed=sorted(observed))


def prepare(args,directory):
    from .af3_analysis import ligand_evidence
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=True)
    validate(dict(args,ligand_smiles=['placeholder']),{})
    sequences,proteins=resolve(args,directory)
    jobs=[]
    for index,item in enumerate(args['ligands'],1):
        jid=f'ligand-{index:05d}';folder=directory/'inputs';folder.mkdir(exist_ok=True)
        path=folder/(jid+'.json')
        write_alphafold3_input(jid,sequences,path,directory,seeds=seeds(args))
        data=read(path);chemistry=ligand_evidence(item['smiles'])
        data['sequences'].append(dict(ligand=dict(id=chr(65+len(sequences)),smiles=item['smiles'])))
        _atomic_json(path,data)
        jobs.append(dict(id=jid,label=item.get('name',jid),molecule_id=item.get('molecule_id'),
                         input=str(path),input_sha256=sha(path),chemistry=chemistry,status='prepared'))
    return dict(kind='af3_batch',status='prepared_not_predicted',proteins=proteins,model_seeds=seeds(args),
                jobs=jobs,shared_msa=args.get('shared_msa',True),experimental_status='not_assessed')


def collect(output,input_path,profile):
    data=read(input_path)
    manifest=inspect_alphafold3_output(output,construct_name=data['name'],model_version=str(profile.get('version','local_unspecified')),
        chain_ids=[next(iter(e.values()))['id'] for e in data['sequences']],input_path=input_path)
    manifest['sample_summary']=sample_summary(output)
    validate_sample_coverage(manifest['sample_summary'],data)
    from .pose_consistency import summarize
    manifest['sample_summary']['pose_consistency']=summarize(manifest['sample_summary']['samples'],data)
    from .pose_quality import af3_ligands
    try:manifest['pose_quality']=af3_ligands(manifest['structure_path'],input_path)
    except Exception as exc:manifest['pose_quality']=dict(qc_status='not_assessed',error=str(exc))
    return manifest


def execute(args,directory,cfg,compute,services):
    from .prompt_workflow import file_lock
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=True)
    state_path=directory/'batch-state.json'
    if state_path.exists():
        state=read(state_path)
        for job in state['jobs']:
            if sha(Path(job['input']))!=job['input_sha256']:raise ValueError('Prepared AF3 input changed')
        if state.get('shared_data'):check_hashes(state['shared_data']['outputs'])
    else:
        state=prepare(args,directory);_atomic_json(state_path,state)
    if not compute:raise Blocked('Compute disabled; AF3 inputs retained at '+str(directory/'inputs'))
    if not cfg.get('af3_profile'):raise Blocked('Configure af3_profile; prepared batch inputs retained at '+str(directory/'inputs'))
    try:
        profile=load_model_profile(Path(cfg['af3_profile']))
        if profile['backend']!='alphafold3':raise ValueError('AF3 batch requires an alphafold3 profile')
        validate_af3_installation(profile)
    except (ValueError,OSError) as exc:raise Blocked(str(exc)) from None
    contract=dict(arguments=args,profile_sha256=sha(Path(cfg['af3_profile'])),service_mode=services.mode)
    if state.get('contract') and state['contract']!=contract:raise ValueError('Batch runtime/input contract changed')
    state['contract']=contract;state['status']='running';_atomic_json(state_path,state)
    # The Chat queue is sequential; the host lock also excludes another project/process.
    lock=file_lock(Path(tempfile.gettempdir())/'aidd-af3-gpu.lock')
    try:lock.__enter__()
    except OSError as exc:raise Blocked('AF3 resource is occupied; resume this task after the active job finishes') from exc
    try:
        prepared=None
        if state['shared_msa']:
            shared=state.get('shared_data')
            if shared:
                prepared=read(shared['path'])
            else:
                dispatch=state.get('shared_dispatch',{})
                uncertain=dispatch.get('status')=='running'
                attempt=dispatch.get('attempt',0)+(0 if uncertain else 1)
                folder=directory/f'shared-data-attempt-{attempt:03d}';folder.mkdir(exist_ok=True)
                marker=folder/'dispatch.json'
                if uncertain:
                    # A finished data output can reconcile an interrupted receipt.
                    try:prepared,data_path=load_data_output(folder)
                    except ValueError:raise Blocked('Uncertain shared-data dispatch; verify process exit and reconcile '+str(marker)) from None
                else:
                    command=prediction_command(profile,Path(state['jobs'][0]['input']),folder)
                    command+=['--run_inference=false','--run_data_pipeline=true']
                    _atomic_json(marker,dict(status='running',argv=command))
                    state['shared_dispatch']=dict(status='running',attempt=attempt,folder=str(folder))
                    _atomic_json(state_path,state)
                    try:services.run_command(command,folder/'execution.log')
                    except Exception:
                        _atomic_json(marker,dict(status='failed',argv=command))
                        state['shared_dispatch']['status']='failed';_atomic_json(state_path,state);raise
                    prepared,data_path=load_data_output(folder)
                if not ready(prepared):raise ValueError('AF3 data output lacks complete inline protein features')
                for job in state['jobs']:transfer(read(job['input']),prepared)
                state['shared_data']=dict(path=str(data_path),outputs=fingerprint([data_path]))
                state['shared_dispatch']['status']='complete'
                _atomic_json(state_path,state)
        for job in state['jobs']:
            if job['status']=='complete':check_hashes(job['outputs']);continue
            if job['status']=='running':
                raise Blocked('Uncertain AF3 child dispatch '+job['id']+'; reconcile its local batch-state receipt before retry')
            attempt=job.get('attempt',0)+1;folder=directory/job['id']/f'attempt-{attempt:03d}';folder.mkdir(parents=True)
            submitted=folder/'input.json';data=read(job['input'])
            _atomic_json(submitted,transfer(data,prepared) if prepared else data)
            output=folder/'prediction';output.mkdir()
            command=prediction_command(profile,submitted,output)
            if prepared:command+=['--run_data_pipeline=false','--run_inference=true']
            job.update(status='running',attempt=attempt,submitted_input=str(submitted),output=str(output),log=str(folder/'execution.log'))
            _atomic_json(folder/'command.json',dict(argv=command,source='trusted_local_profile'))
            _atomic_json(state_path,state);started=time.monotonic()
            try:
                services.run_command(command,folder/'execution.log')
                manifest=collect(output,submitted,profile)
                manifest_path=folder/'prediction-manifest.json';_atomic_json(manifest_path,manifest)
                files=[submitted,manifest_path,*output.rglob('*.cif'),*output.rglob('*.json'),*output.rglob('*.csv')]
                job.update(status='complete',manifest=str(manifest_path),model=manifest['structure_path'],
                           confidence=manifest['confidence'],pose_quality=manifest['pose_quality'],
                           sample_summary=manifest['sample_summary'],outputs=fingerprint(files))
                job.pop('error',None)
            except Exception as exc:
                job.update(status='failed',error=str(exc))
            job['seconds']=time.monotonic()-started;_atomic_json(state_path,state)
        state['status']='complete' if all(j['status']=='complete' for j in state['jobs']) else 'failed'
        _atomic_json(state_path,state)
        if state['status']!='complete':raise RuntimeError('AF3 batch has failed children; completed models retained. Resume retries failed children only: '+str(state_path))
        state['outputs']={p:digest for job in state['jobs'] for p,digest in job['outputs'].items()}
        state['batch_state']=str(state_path)
        return state
    finally:lock.__exit__(None,None,None)


def reconcile(path,job_id,confirm_stopped=False):
    """Operator-only recovery. Never exposed as an autonomous model tool."""
    if not confirm_stopped:raise ValueError('Confirm the old process is stopped before reconciliation')
    path=Path(path);state=read(path)
    if job_id=='shared-data':
        if state.get('shared_dispatch',{}).get('status')!='running':raise ValueError('No uncertain shared-data dispatch')
        state['shared_dispatch'].update(status='failed',reconciled_at=time.time())
        _atomic_json(path,state);return
    rows=[j for j in state['jobs'] if j['id']==job_id]
    if len(rows)!=1 or rows[0]['status']!='running':raise ValueError('Expected one uncertain running child')
    rows[0].update(status='failed',error='Operator confirmed previous process stopped; retry allowed',reconciled_at=time.time())
    _atomic_json(path,state)


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--state',type=Path,required=True);parser.add_argument('--job',required=True)
    parser.add_argument('--confirm-process-stopped',action='store_true')
    options=parser.parse_args();reconcile(options.state,options.job,options.confirm_process_stopped)
