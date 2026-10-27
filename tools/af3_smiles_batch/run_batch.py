"""Run prepared AF3 jobs sequentially through the existing school Apptainer runner."""
import argparse
from datetime import datetime
import json
import os
from pathlib import Path
import subprocess
import time
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


def run(bundle,check=False,runner=None,optimize=True):
    bundle=Path(bundle).resolve();runner=Path(runner or Path(__file__).with_name('run_af3.sh')).resolve()
    manifest=json.loads((bundle/'manifest.json').read_text(encoding='utf-8'))
    # Validate the entire batch before allocating work.
    for job in manifest['jobs']:
        path=(bundle/job['input']).resolve()
        if bundle not in path.parents or sha(path)!=job['sha256']:
            raise ValueError('Input path/checksum mismatch')
    run_id=datetime.now().strftime('%Y%m%d_%H%M%S')+'_'+str(os.getpid())
    logs=bundle/('check_'+run_id if check else 'run_'+run_id);logs.mkdir()
    report=dict(status='running',mode='host_check' if check else 'prediction',jobs=[])
    save(logs/'report.json',report)
    inference_inputs={}
    if optimize and not check:
        try:
            first=bundle/manifest['jobs'][0]['input']
            source=json.loads(first.read_text(encoding='utf-8'))
            started=time.monotonic()
            if ready(source):
                prepared=source;data_path=first;reused_inline=True
            else:
                print('Preparing shared protein MSA/templates once for this batch.',flush=True)
                code,output=invoke(['bash',str(runner),str(first),'--data-only'],logs/'shared_data_pipeline.log')
                if code:raise RuntimeError('Shared data pipeline failed')
                prepared,data_path=load_data_output(output);reused_inline=False
            inputs=logs/'inference_inputs';inputs.mkdir()
            for job in manifest['jobs']:
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
        input_path=inference_inputs.get(job['id'],bundle/job['input'])
        mode=['--check'] if check else ['--inference-only'] if optimize else []
        command=['bash',str(runner),str(input_path)]+mode
        started=time.monotonic();output_root=None
        print('Starting '+job['id'],flush=True)
        code,output_root=invoke(command,logs/(job['id']+'.log'))
        row=dict(id=job['id'],returncode=code,seconds=time.monotonic()-started,
                 input_sha256=job['sha256'],submitted_input_sha256=sha(input_path),
                 submitted_input=str(input_path),output=str(output_root) if output_root else None)
        if code==0 and not check:
            models=[] if output_root is None else [p for p in output_root.rglob('*_model.cif')
                if not any(part.startswith('seed-') for part in p.relative_to(output_root).parts)]
            if len(models)!=1:
                row['error']='Expected exactly one final model';code=1
            else:
                model=models[0];row.update(model=str(model),model_sha256=sha(model))
                summary=model.with_name(model.name[:-len('_model.cif')]+'_summary_confidences.json')
                if summary.exists():row['confidence']=json.loads(summary.read_text())
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
    a=p.parse_args();run(a.bundle,a.check,a.runner,not a.legacy)
