"""Run explicit golden prompts against the real local Chat queue and provider."""
import argparse
import json
import os
import time
from pathlib import Path

from .gaussian_batch import _atomic_json
from .library_acceptance import read,sha


def preflight(runtime=None,config_dir=None,provider='deepseek'):
    from .doctor import environment_report
    from .llm_profiles import select_llm_profile
    from .prediction import load_model_profile,validate_af3_installation
    from .prompt_workflow import load_runtime
    result=dict(environment=environment_report(),provider=dict(status='blocked'),af3=dict(status='blocked'),library=dict(status='blocked'))
    try:
        profile=select_llm_profile(provider,config_dir=config_dir);data=read(profile)
        variable=data.get('api_key_env')
        if variable and not os.environ.get(variable):raise ValueError('Configured provider credential environment variable is unset')
        result['provider']=dict(status='configured_not_contacted',profile_sha256=sha(profile))
    except (ValueError,OSError) as exc:result['provider']['error']=str(exc)
    try:cfg,_=load_runtime(runtime)
    except (ValueError,OSError) as exc:
        result['runtime_error']=str(exc);cfg={}
    try:
        if not cfg.get('af3_profile'):raise ValueError('No configured AF3 profile')
        profile=load_model_profile(Path(cfg['af3_profile']));validate_af3_installation(profile)
        result['af3']=dict(status='configured_not_executed',profile_sha256=sha(Path(cfg['af3_profile'])),
                          gpu_validation='not_run')
    except (ValueError,OSError) as exc:result['af3']['error']=str(exc)
    try:
        batch=Path(cfg['search']['batch']);catalog=batch/'artifacts/catalog.json'
        data=read(catalog)
        if not (batch/'registry.sqlite3').is_file():raise ValueError('Missing library registry')
        result['library']=dict(status='configured_not_scanned',catalog_sha256=sha(catalog),conformers=data.get('conformers'))
    except (ValueError,OSError,KeyError) as exc:result['library']['error']='Local catalog unavailable: '+str(exc)
    return result


def check_report(report,case):
    if report.get('status')!='complete':raise ValueError('Scientific task did not complete')
    if report.get('service_mode')!='live':raise ValueError('Fixture service cannot pass live acceptance')
    actions={s['action'] for s in report['steps'].values() if s['status']=='complete'}
    if set(case.get('expected_actions',[]))-actions:raise ValueError('Expected scientific actions were not executed')
    for action,expected in case.get('result_expect',{}).items():
        matches=[s['result'] for s in report['steps'].values() if s['action']==action and s['status']=='complete']
        if len(matches)!=1:raise ValueError('Expected exactly one result for '+action)
        if any(matches[0].get(k)!=v for k,v in expected.items()):raise ValueError('Result contract mismatch for '+action)
    for expectation in case.get('af3_inputs',[]):
        candidates=[]
        for step in report['steps'].values():
            result=step.get('result',{})
            if result.get('input_json'):candidates.append(Path(result['input_json']))
            candidates.extend(Path(j['input']) for j in result.get('jobs',[]) if j.get('input'))
        if not candidates:raise ValueError('No prepared AF3 inputs to validate')
        for path in candidates:
            data=read(path)
            proteins=[e['protein']['sequence'] for e in data['sequences'] if 'protein' in e]
            ligands=[e['ligand'].get('smiles') for e in data['sequences'] if 'ligand' in e]
            if 'protein_sequences' in expectation and proteins!=expectation['protein_sequences']:raise ValueError('AF3 input protein identity/stoichiometry mismatch')
            if 'seeds' in expectation and data['modelSeeds']!=expectation['seeds']:raise ValueError('AF3 input seed mismatch')
            if 'ligand_count' in expectation and len(ligands)!=expectation['ligand_count']:raise ValueError('AF3 ligand count mismatch')
    for pointer,expected in case.get('expect',{}).items():
        value=report
        for part in pointer.strip('/').split('/'):
            part=part.replace('~1','/').replace('~0','~')
            value=value[int(part)] if isinstance(value,list) else value[part]
        if value!=expected:raise ValueError('Acceptance mismatch at '+pointer)
    from .screening_selection import check_hashes
    for step in report['steps'].values():
        result=step.get('result',{})
        if isinstance(result.get('outputs'),dict):check_hashes(result['outputs'])
        if result.get('model') and sha(Path(result['model']))!=result.get('structure_sha256'):
            raise ValueError('Final structure checksum mismatch')
    return dict(status='passed',actions=sorted(actions),scientific_quality='Only explicitly asserted checks; no predictive superiority claim')


def run(suite,output,storage,runtime=None,config_dir=None,provider='deepseek',execute=False,timeout=3600,input_roots=()):
    from .chat_agent import ChatAgent
    from .chat_uploads import upload
    suite=Path(suite);output=Path(output);output.mkdir(parents=True,exist_ok=True)
    specification=read(suite);cases=specification['cases']
    ids=[c['id'] for c in cases]
    if not cases or len(ids)!=len(set(ids)):raise ValueError('Acceptance case IDs must be unique and nonempty')
    for case in cases:
        if not case.get('prompt') or not case.get('expected_actions') or not (case.get('expect') or case.get('result_expect')):
            raise ValueError('Every case needs prompt, expected_actions and exact report assertions')
    state_path=output/'acceptance.json'
    contract=dict(suite_sha256=sha(suite),storage=str(Path(storage).resolve()),provider=provider,
        input_roots=[str(Path(p).resolve()) for p in input_roots],
        runtime_sha256=sha(Path(runtime)) if runtime else None,
        code={p.name:sha(p) for p in sorted(Path(__file__).parent.glob('*.py'))})
    state=read(state_path) if state_path.exists() else dict(contract=contract,cases=[])
    if state['contract']!=contract:raise ValueError('Acceptance inputs/code changed; choose a fresh output directory')
    state['preflight']=preflight(runtime,config_dir,provider);state['mode']='live' if execute else 'preflight'
    state['status']='running' if execute else 'preflight_only';_atomic_json(state_path,state)
    if not execute:return state
    if 'REPLACE_' in json.dumps(cases):raise ValueError('Replace workstation-specific fields in the acceptance suite before execution')
    if state['preflight']['provider']['status']=='blocked':
        state['status']='blocked';_atomic_json(state_path,state);return state
    app=ChatAgent(storage,runtime=runtime,config_dir=config_dir,allow_compute=True,start=True,result_import_roots=input_roots)
    try:
        if 'session' not in state:
            state['session']=app.new_session();_atomic_json(state_path,state)
        sid=state['session']
        for case in cases:
            prior=next((r for r in state['cases'] if r['id']==case['id']),None)
            if prior and prior.get('status')=='passed':
                for job in prior['tasks']:check_report(read(Path(job['report'])),case)
                continue
            if any(state['preflight'][key]['status']=='blocked' for key in case.get('requires',[]) if key in state['preflight']):
                if prior is None:state['cases'].append(dict(id=case['id'],status='blocked',reason='Required runtime missing'))
                _atomic_json(state_path,state);continue
            if prior and prior.get('status')=='dispatching':
                state['status']='blocked';state['error']='Uncertain prompt submission; reconcile existing task IDs before resuming'
                _atomic_json(state_path,state);return state
            if prior is None or not prior.get('tasks'):
                if prior is not None:state['cases'].remove(prior)
                prior=dict(id=case['id'],status='dispatching',tasks=[]);state['cases'].append(prior)
                attachments=[]
                for source in case.get('attachments',[]):
                    path=(suite.parent/source).resolve()
                    with path.open('rb') as stream:row=upload(app,sid,path.name,stream,path.stat().st_size)
                    attachments.append(row['id'])
                before={j['id'] for j in app.jobs(sid)};_atomic_json(state_path,state)
                answer=app.ask(sid,case['prompt'],provider,attachments=attachments or None)
                jobs=[j for j in app.jobs(sid) if j['id'] not in before]
                prior.update(status='waiting',answer=answer,tasks=[dict(id=j['id'],report=j['report'],plan=j['plan']) for j in jobs])
                _atomic_json(state_path,state)
                if not jobs:
                    prior.update(status='failed',error='Prompt produced no scientific task');_atomic_json(state_path,state);continue
            else:
                # Reuse the sealed task. The worker still enforces code/config and
                # uncertain scientific-dispatch checks; never recreate the prompt.
                for saved in prior['tasks']:
                    job=app.task(sid,saved['id'])
                    if job['status'] in {'failed','blocked','interrupted','cancelled'}:
                        app.ask(sid,'/resume '+job['id'],provider)
            deadline=time.monotonic()+timeout
            while time.monotonic()<deadline:
                current=[app.task(sid,j['id']) for j in prior['tasks']]
                if all(j['status'] not in {'queued','running','planning','pending'} for j in current):break
                time.sleep(1)
            current=[app.task(sid,j['id']) for j in prior['tasks']]
            prior['tasks']=[dict(id=j['id'],report=j['report'],plan=j['plan']) for j in current]
            if any(j['status'] in {'queued','running','planning','pending'} for j in current):
                prior['status']='interrupted_by_acceptance_timeout';state['status']='pending'
                state['error']='Acceptance timeout stops this runner and its worker. Rerun this suite to resume existing sealed tasks; uncertain GPU dispatch requires reconciliation.'
                _atomic_json(state_path,state);return state
            try:
                for job in current:
                    if job['status']!='complete':raise ValueError('Task '+job['id']+': '+job['status']+'; inspect '+str(job.get('log')))
                    check_report(read(Path(job['report'])),case)
                prior['status']='passed'
            except (ValueError,KeyError,OSError) as exc:prior.update(status='failed',error=str(exc))
            _atomic_json(state_path,state)
        state['status']='passed' if all(c['status']=='passed' for c in state['cases']) else 'incomplete'
        _atomic_json(state_path,state);return state
    finally:app.close()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('suite','output','storage'):parser.add_argument('--'+name,type=Path,required=True)
    parser.add_argument('--runtime',type=Path);parser.add_argument('--config-dir',type=Path)
    parser.add_argument('--provider',choices=['gpt','deepseek'],default='deepseek')
    parser.add_argument('--input-root',type=Path,action='append',default=[])
    parser.add_argument('--run',action='store_true');parser.add_argument('--timeout',type=float,default=3600)
    args=parser.parse_args();result=run(args.suite,args.output,args.storage,args.runtime,args.config_dir,args.provider,args.run,args.timeout,args.input_root)
    print(json.dumps(dict(status=result['status'],report=str(args.output/'acceptance.json'))))
    return 0 if result['status']=='passed' else 2


if __name__=='__main__':raise SystemExit(main())
