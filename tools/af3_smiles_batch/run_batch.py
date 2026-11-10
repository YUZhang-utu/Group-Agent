"""Run prepared AF3 jobs sequentially through the existing school Apptainer runner."""
import argparse
from datetime import datetime
import json
import os
from pathlib import Path
import subprocess
import time
import uuid
from contextlib import contextmanager
from prepare_batch import sha,save
from msa_reuse import ready,transfer,load_data_output


def invoke(command,log_path):
    output_root=None
    with Path(log_path).open('w',encoding='utf-8') as log:
        proc=subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
        for line in proc.stdout:
            log.write(line);log.flush();print(line,end='',flush=True)
            if line.startswith('Output: '):output_root=Path(line[8:].strip())
        return proc.wait(),output_root


@contextmanager
def batch_lock(bundle):
    """Release the interprocess lock on exit, including process interruption."""
    with (Path(bundle)/'.batch.lock').open('a+b') as handle:
        handle.seek(0,2)
        if handle.tell()==0:handle.write(b'0');handle.flush()
        handle.seek(0)
        if os.name=='nt':
            import msvcrt
            msvcrt.locking(handle.fileno(),msvcrt.LK_NBLCK,1)
        else:
            import fcntl
            fcntl.flock(handle.fileno(),fcntl.LOCK_EX|fcntl.LOCK_NB)
        try:yield
        finally:
            if os.name=='nt':
                handle.seek(0);msvcrt.locking(handle.fileno(),msvcrt.LK_UNLCK,1)


def run(bundle,check=False,runner=None,optimize=True,resume=None):
    with batch_lock(Path(bundle).resolve()):
        return _run(bundle,check,runner,optimize,resume)


def _run(bundle,check=False,runner=None,optimize=True,resume=None):
    bundle=Path(bundle).resolve();runner=Path(runner or Path(__file__).with_name('run_af3.sh')).resolve()
    manifest=json.loads((bundle/'manifest.json').read_text(encoding='utf-8'))
    # Validate the entire batch before allocating work.
    for job in manifest['jobs']:
        path=(bundle/job['input']).resolve()
        if bundle not in path.parents or sha(path)!=job['sha256']:
            raise ValueError('Input path/checksum mismatch')
    contract=dict(manifest_sha256=sha(bundle/'manifest.json'),optimize=optimize,
        code={p.name:sha(p) for p in [Path(__file__),Path(__file__).with_name('prepare_batch.py'),Path(__file__).with_name('msa_reuse.py'),runner]},
        runtime={name:os.environ.get(name) for name in ('AF3_INSTALL','AF3_IMAGE','AF3_MODELS','AF3_DB_DIR')})
    prior=None;completed={}
    if resume is not None:
        if check:raise ValueError('Host checks cannot resume predictions')
        resume=Path(resume).resolve()
        if bundle not in resume.parents:raise ValueError('Resume report must be inside this batch')
        prior=json.loads(resume.read_text(encoding='utf-8'))
        if prior.get('mode')!='prediction' or prior.get('contract')!=contract:
            raise ValueError('Changed inputs/code/runtime or legacy report; use a fresh run')
        for row in prior['jobs']:
            if row.get('status')!='model_generated':continue
            if sha(Path(row['model']))!=row['model_sha256']:
                raise ValueError('Completed model changed; inspect the original run')
            if sha(Path(row['submitted_input']))!=row['submitted_input_sha256']:
                raise ValueError('Submitted AF3 input changed')
            if row.get('confidence_path') and sha(Path(row['confidence_path']))!=row['confidence_sha256']:
                raise ValueError('Completed confidence data changed')
            completed[row['id']]=row
    run_id=datetime.now().strftime('%Y%m%d_%H%M%S')+'_'+str(os.getpid())+'_'+uuid.uuid4().hex[:8]
    logs=bundle/('check_'+run_id if check else 'run_'+run_id);logs.mkdir()
    report=dict(status='running',mode='host_check' if check else 'prediction',total_jobs=len(manifest['jobs']),jobs=[],
                contract=contract,resumed_from=str(resume) if resume else None,reused_jobs=0)
    save(logs/'report.json',report)
    inference_inputs={}
    if optimize and not check and len(completed)<len(manifest['jobs']):
        try:
            first=bundle/manifest['jobs'][0]['input']
            source=json.loads(first.read_text(encoding='utf-8'))
            started=time.monotonic()
            shared=prior.get('shared_data') if prior else None
            if shared:
                data_path=Path(shared['source'])
                if sha(data_path)!=shared['sha256']:raise ValueError('Shared protein data changed')
                prepared=json.loads(data_path.read_text(encoding='utf-8'));reused_inline=shared['reused_inline']
            elif ready(source):
                prepared=source;data_path=first;reused_inline=True
            else:
                print('Preparing shared protein MSA/templates once for this batch.',flush=True)
                code,output=invoke(['bash',str(runner),str(first),'--data-only'],logs/'shared_data_pipeline.log')
                if code:raise RuntimeError('Shared data pipeline failed')
                prepared,data_path=load_data_output(output);reused_inline=False
            inputs=logs/'inference_inputs';inputs.mkdir()
            for job in manifest['jobs']:
                if job['id'] in completed:continue
                original=json.loads((bundle/job['input']).read_text(encoding='utf-8'))
                enriched=transfer(original,prepared)
                path=inputs/(job['id']+'.json');save(path,enriched);inference_inputs[job['id']]=path
            report['shared_data']=dict(source=str(data_path),sha256=sha(data_path),
                reused_inline=reused_inline,seconds=time.monotonic()-started,
                note='All ligand tasks skip MSA/template searches; ligand featurisation still runs.')
            save(logs/'report.json',report)
        except Exception as exc:
            report.update(status='failed',stage='shared_data',error=str(exc));save(logs/'report.json',report)
            raise
    for job in manifest['jobs']:
        if job['id'] in completed:
            row=dict(completed[job['id']],reused=True,original_seconds=completed[job['id']]['seconds'],seconds=0)
            report['jobs'].append(row);report['reused_jobs']+=1;save(logs/'report.json',report)
            print('Reused verified model '+job['id'],flush=True);continue
        input_path=inference_inputs.get(job['id'],bundle/job['input'])
        mode=['--check'] if check else ['--inference-only'] if optimize else []
        command=['bash',str(runner),str(input_path)]+mode
        started=time.monotonic();output_root=None
        print('Starting '+job['id'],flush=True)
        row=dict(id=job['id'],input_sha256=job['sha256'],submitted_input=str(input_path))
        try:
            row['submitted_input_sha256']=sha(input_path)
            code,output_root=invoke(command,logs/(job['id']+'.log'))
            row['process_returncode']=code
            row['output']=str(output_root) if output_root else None
            if code==0 and not check:
                models=[] if output_root is None else [p for p in output_root.rglob('*_model.cif')
                    if not any(part.startswith('seed-') for part in p.relative_to(output_root).parts)]
                if len(models)!=1:
                    raise ValueError('Expected exactly one final model')
                model=models[0];row.update(model=str(model),model_sha256=sha(model))
                summary=model.with_name(model.name[:-len('_model.cif')]+'_summary_confidences.json')
                if summary.exists():
                    row['confidence']=json.loads(summary.read_text())
                    row.update(confidence_path=str(summary),confidence_sha256=sha(summary))
        except (OSError, ValueError) as exc:
            code=1;row['error']=str(exc)
        row.update(returncode=code,seconds=time.monotonic()-started)
        row['status']='host_checks_passed' if check and code==0 else 'model_generated' if code==0 else 'failed'
        report['jobs'].append(row);save(logs/'report.json',report)
        if code:
            report['status']='failed';save(logs/'report.json',report)
            raise RuntimeError('Job failed; see '+str(logs/'report.json'))
    report['status']='host_checks_passed_not_inference' if check else 'complete'
    save(logs/'report.json',report);print('Batch report: '+str(logs/'report.json'),flush=True)
    return report


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('bundle',type=Path)
    mode=p.add_mutually_exclusive_group(required=True);mode.add_argument('--check',action='store_true');mode.add_argument('--run',action='store_true')
    p.add_argument('--runner',type=Path);p.add_argument('--legacy',action='store_true',help='Run each full pipeline separately')
    p.add_argument('--resume',type=Path,help='Resume a matching report; verify and reuse completed local models')
    a=p.parse_args();run(a.bundle,a.check,a.runner,not a.legacy,a.resume)
