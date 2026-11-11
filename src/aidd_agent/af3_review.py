"""Read recorded AF3 outputs without changing tasks or launching calculations."""
from pathlib import Path

from .af3_batch import sample_summary,validate_sample_coverage
from .af3_analysis import assess_chain_mapping
from .library_acceptance import read,sha
from .project_context import ensure_within


def review(app,sid,args):
    root=Path(app.context_for(sid)['project_root']);job=app.task(sid,args['task_id'])
    if job['status']!='complete' or not job.get('report'):raise ValueError('Review requires a completed owned task')
    report_path=ensure_within(Path(job['report']),root);report=read(report_path)
    options=[(key,row['result']) for key,row in report.get('steps',{}).items()
             if row.get('status')=='complete' and row.get('action')=='af3_run' and isinstance(row.get('result'),dict)]
    if args.get('step_id'):options=[item for item in options if item[0]==args['step_id']]
    if len(options)!=1:raise ValueError('Select one af3_run step from the completed task')
    step,result=options[0];manifest={}
    if result.get('manifest'):
        manifest_path=ensure_within(Path(result['manifest']),root);manifest=read(manifest_path)
    folder=result.get('output_directory') or manifest.get('provenance',{}).get('output_dir')
    if not folder and result.get('model'):folder=str(Path(result['model']).parent)
    if not folder:raise ValueError('No recorded AF3 output directory')
    folder=ensure_within(Path(folder),root)
    if not folder.is_dir():raise ValueError('Recorded AF3 output directory is unavailable on this host')
    for pattern in ('*summary_confidences.json','*ranking_scores.csv','*_model.cif','model.cif'):
        for path in folder.rglob(pattern):ensure_within(path,root);ensure_within(path,folder)
    summary=sample_summary(folder)
    prepared=[]
    for row in report.get('steps',{}).values():
        value=row.get('result',{})
        if row.get('status')=='complete' and value.get('input_json'):
            path=ensure_within(Path(value['input_json']),root)
            if value.get('input_sha256') and sha(path)!=value['input_sha256']:raise ValueError('Prepared input hash changed')
            prepared.append(path)
    recorded_input=manifest.get('provenance',{}).get('input_path')
    if recorded_input:
        path=ensure_within(Path(recorded_input),root)
        if sha(path)!=manifest['provenance'].get('input_sha256'):raise ValueError('Manifest input hash changed')
        prepared=[path]
    if len(prepared)==1:validate_sample_coverage(summary,read(prepared[0]))
    else:summary['seed_coverage']=dict(status='not_verified',reason='No unique recorded input JSON')
    confidence=manifest.get('confidence',result.get('confidence',{}))
    confidence_path=manifest.get('provenance',{}).get('confidence_path')
    if confidence_path:
        confidence_path=ensure_within(Path(confidence_path),root)
        expected=manifest['provenance'].get('confidence_sha256')
        if expected and sha(confidence_path)!=expected:raise ValueError('Recorded confidence file changed')
        confidence=read(confidence_path)
    model_path=result.get('model') or manifest.get('structure_path');model_check='not_available'
    if model_path:
        path=ensure_within(Path(model_path),root)
        expected=result.get('structure_sha256') or manifest.get('structure_sha256')
        if expected and sha(path)!=expected:raise ValueError('Recorded final CIF changed')
        model_check='hash_verified' if expected else 'present_unsealed' if path.is_file() else 'missing'
    best=summary.get('best_recorded_sample');best_file=None
    if best:
        key=f"seed-{int(best['seed'])}_sample-{int(best['sample'])}"
        matching=[r for r in summary['samples'] if r['sample']==key]
        if len(matching)==1:best_file=matching[0]['model']
    return dict(task_id=job['id'],step_id=step,status='reviewed',operation='read_only_no_inference',
        source_report=str(report_path),source_report_sha256=sha(report_path),output_directory=str(folder),
        final_model_check=model_check,best_recorded_sample_model=best_file,
        confidence_chain_mapping=assess_chain_mapping(confidence),**summary)
