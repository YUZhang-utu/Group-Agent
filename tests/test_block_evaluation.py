"""Synthetic panel/adapter contracts; not real PLANTS or full-library acceptance."""
import csv
import json
from pathlib import Path
import sqlite3
import subprocess

import pytest

from aidd_agent.block_sampling import sample, check_outputs
from aidd_agent.block_plants import prepare, run, analyze, parse_ranking, recover
from aidd_agent.final_work_blocks import sha
from aidd_agent.joint_spatial_profiles import save
from aidd_agent.mol2 import parse_mol2_block, iter_mol2_blocks
from aidd_agent.prompt_plan import validate_plan, CAPABILITIES
from aidd_agent.chat_agent import ChatAgent, validate_route


def raw(i):
    return f'@<TRIPOS>MOLECULE\nm{i//2}_conf{i%2+1}\n2 1 0 0 0\nSMALL\nUSER_CHARGES\n\n@<TRIPOS>ATOM\n1 C1 {i}.0 0.0 0.0 C.3 1 MOL 0.0\n2 O1 {i}.0 1.0 0.0 O.3 1 MOL 0.0\n@<TRIPOS>BOND\n1 1 2 1\n'


def fixture(root, n=12):
    paths={k:root/k for k in ('E094','E095','E096','profiles','descriptors')}
    for p in paths.values(): p.mkdir(parents=True)
    source=root/'source.mol2'; source.write_text(''.join(raw(i) for i in range(n)),encoding='utf-8')
    with sqlite3.connect(paths['descriptors']/'descriptors.sqlite') as db:
        db.execute('CREATE TABLE descriptor(cid TEXT PRIMARY KEY,payload TEXT)')
        for i in range(n):
            r=parse_mol2_block(source,i,raw(i))
            d=dict(conformer_id=f'id{i}',molecule_id=f'm{i//2}',provenance=dict(content_sha256=r.content_sha256,
                source_record_name=r.name,source_name=r.molecule_name,source_path=str(source),source_record_index=i))
            db.execute('INSERT INTO descriptor VALUES(?,?)',(f'id{i}',json.dumps(d)))
    save(paths['descriptors']/'report.json',dict(status='complete'))
    for s in ('E094','E095','E096'):
        groups={f'b{j}':[] for j in range(2)}
        for i in range(n): groups[f'b{i%2 if s=="E096" else int(i>=n//2)}'].append(f'id{i}')
        with (paths[s]/'blocks.csv').open('w',newline='') as f:
            w=csv.writer(f);w.writerow(['block_id','conformers']);w.writerows((b,len(v)) for b,v in groups.items())
        r=dict(status='complete',regular_conformers=n,output_hashes={'blocks.csv':sha(paths[s]/'blocks.csv')})
        if s!='E094':
            with sqlite3.connect(paths[s]/'property_blocks.sqlite') as db:
                db.execute('CREATE TABLE membership(cid TEXT PRIMARY KEY,parent_block TEXT,block_id TEXT)')
                db.executemany('INSERT INTO membership VALUES(?,?,?)',[(cid,f'b{int(int(cid[2:])>=n//2)}',bid) for bid,ids in groups.items() for cid in ids])
            r.update(parent_blocks=str(paths['E094']),membership_gate='passed_against_complete_cached_profiles',sources={str(paths['E094']/'report.json'):sha(paths['E094']/'report.json')})
        save(paths[s]/'report.json',r)
    with sqlite3.connect(paths['profiles']/'profiles.sqlite') as db:
        db.executescript('CREATE TABLE files(id INTEGER PRIMARY KEY,path TEXT); CREATE TABLE item(cid TEXT PRIMARY KEY,parent TEXT,file_id INTEGER,idx INTEGER,mid TEXT,vector BLOB);')
        db.execute('INSERT INTO files VALUES(1,?)',(str(source),))
        db.executemany('INSERT INTO item VALUES(?,?,?,?,?,?)',[(f'id{i}',f'b{int(i>=n//2)}',1,i,f'm{i//2}',b'') for i in range(n)])
    save(paths['profiles']/'report.json',dict(status='complete',coverage=1,total_regular=n,blocks=str(paths['E094']),signature=dict(sources={str(paths[k]/'report.json'):sha(paths[k]/'report.json') for k in ('E094','descriptors')})))
    receipt=root/'validation.json'; save(receipt,dict(status='complete',structural_gate='passed',property_gate='passed',regular_conformers=n))
    cfg=root/'sampling.json'; save(cfg,{k:str(p) for k,p in paths.items()}|dict(validation=str(receipt)))
    return cfg,source


def plants_profile(root,source):
    receptor=root/'receptor.mol2'; receptor.write_text(raw(0))
    executable=root/'plants'; executable.write_text('synthetic fixture executable; subprocess mocked')
    p=root/'plants.json';save(p,dict(executable=str(executable),ligand_chemistry_reviewed=True,workers=1,
                                   receptors=[dict(id='MDM2_fixture',mol2=str(receptor),center=[0,0,0],radius=10,reviewed=True,evidence='Synthetic fixture only')]))
    return p


def fake_plants(argv, *, cwd, stdout, stderr, timeout, shell):
    assert shell is False and argv[-3:]==['--mode','screen','plantsconfig']
    p=Path(cwd);cfg=(p/'plantsconfig').read_text();assert 'rigid_ligand 1' in cfg
    out=p/'poses';out.mkdir();rows=[];texts=[]
    for i,(_,text) in enumerate(iter_mol2_blocks(p/'ligands.mol2'),1):
        lines=text.splitlines(keepends=True);alias=lines[1].strip();name=f'{alias}_entry_{i:05d}_conf_01'
        lines[1]=name+'\n';texts.append(''.join(lines));rows.append([name,-float(int(alias[1:]))])
    with (out/'bestranking.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['LIGAND_ENTRY','TOTAL_SCORE']);w.writerows(rows)
    (out/'docked_ligands.mol2').write_text(''.join(texts))
    return subprocess.CompletedProcess(argv,0)


def test_sampling_reproducible_deduplicated_and_raw_preserved(tmp_path):
    cfg,source=fixture(tmp_path)
    a=sample(cfg,tmp_path/'a',count=3); b=sample(cfg,tmp_path/'b',count=3)
    assert a['slots']==18 and a['unique_conformers']<18
    assert a['output_hashes']['block_samples.csv']==b['output_hashes']['block_samples.csv']
    db=sqlite3.connect(tmp_path/'a/samples.sqlite')
    assert db.execute('SELECT min(n),max(n) FROM (SELECT count(*) n FROM sample GROUP BY scheme,block_id)').fetchone()==(3,3)
    exports=json.loads((tmp_path/'a/exports.json').read_text())['chunks'];seen=set()
    original=dict(iter_mol2_blocks(source))
    for chunk in exports:
        for (_,text),record in zip(iter_mol2_blocks(tmp_path/'a'/chunk['path']),chunk['records'],strict=True):
            assert text.splitlines()[2:]==original[record['source_index']].splitlines()[2:]
            assert record['cid'] not in seen;seen.add(record['cid'])
    assert len(seen)==a['unique_conformers']
    assert sample(cfg,tmp_path/'a',count=3,resume=True)==a
    with pytest.raises(ValueError,match='fresh output'): sample(cfg,tmp_path/'a',count=3,seed=1,resume=True)
    with pytest.raises(ValueError): sample(cfg,tmp_path/'bad',schemes=['E096','E096'])


def test_small_blocks_all_members_and_changed_source_rejected(tmp_path):
    cfg,source=fixture(tmp_path);a=sample(cfg,tmp_path/'all',count=100)
    assert a['slots']==36 and a['unique_conformers']==12
    source.write_text(source.read_text().replace('C1 0.0','C1 8.0'))
    with pytest.raises(ValueError,match='Raw identity'): sample(cfg,tmp_path/'broken',count=100)


def test_shared_docking_analysis_resume_and_pose_integrity(tmp_path,monkeypatch):
    cfg,source=fixture(tmp_path);s=tmp_path/'samples';sample(cfg,s,count=100)
    profile=plants_profile(tmp_path,source);p=tmp_path/'prepared';prepare(s/'report.json',profile,p)
    calls=[]
    def runner(*a,**kw):calls.append(kw['cwd']);return fake_plants(*a,**kw)
    monkeypatch.setattr('aidd_agent.block_plants.subprocess.run',runner)
    out=tmp_path/'docking';r=run(p/'report.json',out)
    assert r['status']=='complete' and r['scored']==12 and len(calls)==1
    assert run(p/'report.json',out)['scored']==12 and len(calls)==1
    result=analyze(out/'report.json',tmp_path/'analysis')
    assert len(result['comparisons'])==3
    assert all(x['sampled']==12 and x['scored']==12 for x in result['comparisons'])
    assert len({x['population_weighted_mean'] for x in result['comparisons']})==1
    poses=next(out.rglob('docked_ligands.mol2'));poses.write_text(poses.read_text()+'\n')
    with pytest.raises(ValueError,match='changed'): run(p/'report.json',out)


def test_failed_jobs_are_partial_and_retried(tmp_path,monkeypatch):
    cfg,source=fixture(tmp_path);s=tmp_path/'samples';sample(cfg,s,count=2)
    p=tmp_path/'prepared';prepare(s/'report.json',plants_profile(tmp_path,source),p)
    monkeypatch.setattr('aidd_agent.block_plants.subprocess.run',lambda *a,**kw:subprocess.CompletedProcess(a[0],1))
    out=tmp_path/'docking';r=run(p/'report.json',out)
    assert r['status']=='partial' and r['scored']==0
    result=analyze(out/'report.json',tmp_path/'partial-analysis')
    assert all(v['population_weighted_mean'] is None for v in result['comparisons'])
    monkeypatch.setattr('aidd_agent.block_plants.subprocess.run',fake_plants)
    assert run(p/'report.json',out)['status']=='complete'
    assert list(out.rglob('attempt-002'))


def test_worker_override_and_stopped_checkpoint_migration(tmp_path,monkeypatch):
    from concurrent.futures import ThreadPoolExecutor
    cfg,source=fixture(tmp_path);s=tmp_path/'samples';sample(cfg,s,count=2)
    profile=plants_profile(tmp_path,source);config=json.loads(profile.read_text())
    config['receptors'].append(dict(config['receptors'][0],id='second'))
    save(profile,config)
    p=tmp_path/'prepared';prepare(s/'report.json',profile,p)
    hashes={str(f.relative_to(p)):sha(f) for f in p.rglob('*') if f.is_file()}
    pools=[];calls=[]
    def pool(*,max_workers):
        pools.append(max_workers);return ThreadPoolExecutor(max_workers=max_workers)
    def engine(*a,**kw):calls.append(1);return fake_plants(*a,**kw)
    monkeypatch.setattr('aidd_agent.block_plants.ThreadPoolExecutor',pool)
    monkeypatch.setattr('aidd_agent.block_plants.subprocess.run',engine)
    out=tmp_path/'docking'
    assert run(p/'report.json',out,workers=20,max_jobs=1)['workers_this_invocation']==20
    old_report=(out/'report.json').read_bytes();old_scores=(out/'scores.csv').read_bytes()
    run(p/'report.json',out,workers=8)
    # Simulate a stopped run: later completed receipts exist but summary is still old.
    (out/'report.json').write_bytes(old_report);(out/'scores.csv').write_bytes(old_scores)
    with pytest.raises(ValueError,match='include-checkpoints'):
        recover(out/'report.json',tmp_path/'rejected')
    new=tmp_path/'recovered'
    assert recover(out/'report.json',new,include_checkpoints=True)['status']=='complete'
    assert run(p/'report.json',new,workers=24)['status']=='complete'
    assert pools==[20,8,24] and len(calls)==2
    assert hashes=={str(f.relative_to(p)):sha(f) for f in p.rglob('*') if f.is_file()}
    for value in (0,65,True,1.5):
        with pytest.raises(ValueError,match='workers must'):
            run(p/'report.json',new,workers=value)


def test_failed_first_job_stops_unscheduled_receptors(tmp_path,monkeypatch):
    cfg,source=fixture(tmp_path);s=tmp_path/'samples';sample(cfg,s,count=2)
    profile=plants_profile(tmp_path,source);c=json.loads(profile.read_text())
    c['receptors'].append(dict(c['receptors'][0],id='second_fixture'));save(profile,c)
    p=tmp_path/'prepared';prepare(s/'report.json',profile,p);calls=[]
    def fail(*a,**kw):calls.append(a);return subprocess.CompletedProcess(a[0],1)
    monkeypatch.setattr('aidd_agent.block_plants.subprocess.run',fail)
    report=run(p/'report.json',tmp_path/'docking')
    assert len(calls)==1 and report['jobs']==2 and report['attempted_jobs']==1 and report['first_job_gate']=='failed'


@pytest.mark.parametrize('name,score',[('other_entry_00001_conf_01','-1'),('c000000001_entry_00002_conf_01','-1'),('c000000001_entry_00001_conf_01','nan')])
def test_bad_output_identity_and_scores_rejected(tmp_path,name,score):
    p=tmp_path/'scores.csv';p.write_text('LIGAND_ENTRY,TOTAL_SCORE\n'+name+','+score+'\n')
    with pytest.raises(ValueError): parse_ranking(p,[dict(cid='a',alias='c000000001')])


def test_native_plants_omitted_name_header_and_reordered_entries(tmp_path):
    p=tmp_path/'scores.csv'
    p.write_text('TOTAL_SCORE,SCORE_RB_PEN,TIME\n'
                 'c000141837_entry_00002_conf_01,-82.2662,-82.2662,1.04165\n'
                 'c000132700_entry_00001_conf_01,-73.9855,-73.9855,1.14469\n')
    records=[dict(cid='first',alias='c000132700'),dict(cid='second',alias='c000141837')]
    rows=parse_ranking(p,records)
    assert rows['c000132700']['score']==-73.9855
    assert rows['c000141837']['cid']=='second'


@pytest.mark.parametrize('header,row',[
    ('TOTAL_SCORE,TIME','c000000001_entry_00001_conf_01,-1'),
    ('TOTAL_SCORE,TIME','-1,1'),
    ('LIGAND_ENTRY,TOTAL_SCORE','c000000001_entry_00001_conf_01,-1,extra'),
    ('TOTAL_SCORE,TOTAL_SCORE','c000000001_entry_00001_conf_01,-1,-2'),
    ('TOTAL_SCORE,TIME','c000000001_entry_00002_conf_01,-1,1'),
    ('TOTAL_SCORE,TIME','c000000001_entry_00001_conf_01,nan,1'),
])
def test_ranking_does_not_guess_malformed_column_or_identity(tmp_path,header,row):
    p=tmp_path/'scores.csv';p.write_text(header+'\n'+row+'\n')
    with pytest.raises(ValueError):parse_ranking(p,[dict(cid='a',alias='c000000001')])


def test_recover_native_outputs_without_engine_and_resume(tmp_path,monkeypatch):
    cfg,source=fixture(tmp_path);s=tmp_path/'samples';sample(cfg,s,count=2)
    p=tmp_path/'prepared';prepare(s/'report.json',plants_profile(tmp_path,source),p)
    def native(*a,**kw):
        result=fake_plants(*a,**kw)
        ranking=Path(kw['cwd'])/'poses/bestranking.csv'
        ranking.write_text(ranking.read_text().replace('LIGAND_ENTRY,TOTAL_SCORE','TOTAL_SCORE'))
        return result
    def old_parser(*args):raise ValueError('Unknown PLANTS ranking header')
    monkeypatch.setattr('aidd_agent.block_plants.subprocess.run',native)
    monkeypatch.setattr('aidd_agent.block_plants.parse_ranking',old_parser)
    out=tmp_path/'old';assert run(p/'report.json',out)['failed_jobs']==1
    before={str(f.relative_to(out)):sha(f) for f in out.rglob('*') if f.is_file()}
    monkeypatch.setattr('aidd_agent.block_plants.parse_ranking',parse_ranking)
    def forbidden(*a,**kw):pytest.fail('Recovery/resume must not launch the engine')
    monkeypatch.setattr('aidd_agent.block_plants.subprocess.run',forbidden)
    new=tmp_path/'recovered';r=recover(out/'report.json',new)
    assert r['failed_jobs']==0 and r['scored']>0 and r['first_job_gate']=='passed'
    assert before=={str(f.relative_to(out)):sha(f) for f in out.rglob('*') if f.is_file()}
    assert run(p/'report.json',new,max_jobs=1)['scored']==r['scored']
    with (new/'scores.csv').open() as f:
        assert all(Path(row['pose_file']).is_relative_to(new) for row in csv.DictReader(f))
    assert analyze(new/'report.json',tmp_path/'analysis')['status']=='complete'
    pose=next(out.rglob('docked_ligands.mol2'));pose.write_text(pose.read_text()+'\n')
    with pytest.raises(ValueError,match='hash/path mismatch'):recover(out/'report.json',tmp_path/'tampered')


def test_recovery_rejects_engine_failure(tmp_path,monkeypatch):
    cfg,source=fixture(tmp_path);s=tmp_path/'samples';sample(cfg,s,count=2)
    p=tmp_path/'prepared';prepare(s/'report.json',plants_profile(tmp_path,source),p)
    monkeypatch.setattr('aidd_agent.block_plants.subprocess.run',lambda *a,**kw:subprocess.CompletedProcess(a[0],1))
    out=tmp_path/'failed';run(p/'report.json',out)
    with pytest.raises(ValueError,match='known post-execution'):
        recover(out/'report.json',tmp_path/'recovered')


def test_model_schema_dependencies_and_no_model_paths():
    steps=[dict(id='sample',action='block_sample',params=dict(count=100)),
           dict(id='prepare',action='block_plants_prepare',params=dict(sample_step='sample')),
           dict(id='dock',action='block_plants_run',params=dict(prepared_step='prepare')),
           dict(id='analysis',action='block_analyze',params=dict(docking_step='dock'))]
    plan=dict(version=1,summary='Evaluate block panels',clarifications=[],steps=steps)
    assert validate_plan(plan)==plan and 'block_sample' in CAPABILITIES['actions']
    steps[1]['params']=dict(sample_step='dock')
    with pytest.raises(ValueError,match='dependency'): validate_plan(plan)
    steps[1]['params']=dict(sample_step='sample',profile='/invented/path')
    with pytest.raises(ValueError,match='parameters'): validate_plan(plan)


def test_confirmed_receptor_and_dock_chains_are_allowed_but_arbitrary_source_chains_are_not():
    for first,second,ref in [('pocket_adopt','plants_receptors','adoption_step'),('block_plants_run','block_analyze','docking_step')]:
        p=dict(version=1,summary='Continue explicitly authorized stages',clarifications=[],steps=[
            dict(id='first',action=first,params=dict(source_run='PROMPT-'+'0'*16)),
            dict(id='second',action=second,params={ref:'first'})])
        assert validate_plan(p)==p
        p['steps'][1]['params']={ref:'missing'}
        with pytest.raises(ValueError):validate_plan(p)
    p=dict(version=1,summary='Assess receptors',clarifications=[],steps=[dict(id='target',action='protein_from_pdb',params=dict(pdb_id='1ABC')),
        dict(id='assess',action='receptor_assess',params=dict(protein_step='target',reference_pdb='1ABC'))])
    assert validate_plan(p)==p
    p['steps'][1]['params']['reference_pdb']='../../bad'
    with pytest.raises(ValueError):validate_plan(p)


def test_natural_language_router_queues_owned_block_plan(tmp_path):
    decision=dict(intent='block_evaluation',message='Sample all three schemes and prepare PLANTS jobs.',task_id=None,request='',evaluation=dict(stage='sample_prepare',count=100))
    app=ChatAgent(tmp_path,router=lambda *args:decision,start=False)
    try:
        sid=app.new_session();app.ask(sid,'Take 100 random conformers from each of the three block schemes and prepare docking','deepseek')
        job=app.task(sid);plan=json.loads(Path(job['plan']).read_text())['plan']
        assert [s['action'] for s in plan['steps']]==['block_sample','block_plants_prepare']
        assert job['status']=='queued'
        with pytest.raises(ValueError): validate_route(dict(decision,evaluation=dict(stage='run',count=100)))
        with pytest.raises(ValueError): validate_route(dict(decision,evaluation=dict(stage='sample',path='/invented')))
    finally:app.close()


def test_real_project_executor_all_four_stages_and_sealed_followup(tmp_path,monkeypatch):
    from aidd_agent.prompt_workflow import initialize_context, create_plan, run_plan
    cfg,source=fixture(tmp_path/'inputs');profile=plants_profile(tmp_path,source)
    runtime=tmp_path/'runtime.json';save(runtime,dict(block_evaluation=dict(sampling_profile=str(cfg),plants_profile=str(profile))))
    context=initialize_context(tmp_path/'project','tester','Block evaluation')
    args=(Path(context['db']),context['user_id'],context['project_id'])
    steps=[dict(id='sample',action='block_sample',params=dict(count=2)),
           dict(id='prepare',action='block_plants_prepare',params=dict(sample_step='sample')),
           dict(id='dock',action='block_plants_run',params=dict(prepared_step='prepare')),
           dict(id='analysis',action='block_analyze',params=dict(docking_step='dock'))]
    plan=create_plan(*args,'Evaluate samples',local_plan=dict(version=1,summary='Evaluate samples',clarifications=[],steps=steps))
    monkeypatch.setattr('aidd_agent.prompt_workflow.platform.platform',lambda:'synthetic-test-platform')
    monkeypatch.setattr('aidd_agent.block_plants.subprocess.run',fake_plants)
    report,_=run_plan(*args,plan,runtime=runtime,allow_compute=True)
    assert report['status']=='complete', report
    assert report['steps']['analysis']['result']['comparisons']
    followup=create_plan(*args,'Analyze the completed docking',local_plan=dict(version=1,summary='Analyze docking',clarifications=[],
        steps=[dict(id='analysis',action='block_analyze',params=dict(source_run=plan.parent.name))]))
    report,_=run_plan(*args,followup,runtime=runtime)
    assert report['status']=='complete',report
