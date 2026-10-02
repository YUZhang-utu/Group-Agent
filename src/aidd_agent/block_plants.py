"""Prepared-receptor PLANTS jobs and block-aware score analysis. No shell execution."""
from __future__ import annotations

import argparse
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
import csv
import json
import math
from pathlib import Path
import re
import shutil
import subprocess
import time

import numpy as np

from .block_sampling import check_outputs, verified_report, readonly, exclusive_output
from .final_work_blocks import read, sha
from .joint_spatial_profiles import save
from .mol2 import iter_mol2_blocks

SOURCE_ACTIONS={'plants_receptors':'pocket_adopt','block_plants_prepare':'block_sample', 'block_plants_run':'block_plants_prepare',
                'block_analyze':'block_plants_run'}


def write_csv(path, rows, fields):
    with Path(path).open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields); w.writeheader(); w.writerows(rows)


def profile_read(path):
    path=Path(path).resolve(); c=read(path)
    if set(c)-{'executable','receptors','ligand_mode','search_speed','cluster_structures','workers','timeout_seconds','ligand_chemistry_reviewed'}:
        raise ValueError('Unknown PLANTS profile setting')
    c.setdefault('ligand_chemistry_reviewed',False)
    if type(c['ligand_chemistry_reviewed']) is not bool: raise ValueError('Invalid ligand review annotation')
    def resolve(p):
        p=Path(p); return str((p if p.is_absolute() else path.parent/p).resolve())
    c['executable']=resolve(c['executable'])
    c.setdefault('ligand_mode','rigid'); c.setdefault('search_speed','speed1'); c.setdefault('cluster_structures',1)
    c.setdefault('workers',4); c.setdefault('timeout_seconds',7200)
    if c['ligand_mode'] not in ('rigid','flexible') or c['search_speed'] not in ('speed1','speed2','speed4'): raise ValueError('Unsupported mode/speed')
    for k,lo,hi in [('cluster_structures',1,20),('workers',1,64),('timeout_seconds',1,604800)]:
        if type(c[k]) is not int or not lo<=c[k]<=hi: raise ValueError('Invalid '+k)
    receptors=c.get('receptors')
    if not isinstance(receptors,list) or not receptors: raise ValueError('Prepared receptors and reviewed binding sites required')
    seen=set()
    for r in receptors:
        if set(r)!={'id','mol2','center','radius','reviewed','evidence'} or r['reviewed'] is not True or not isinstance(r['evidence'],str) or not r['evidence'].strip():
            raise ValueError('Each receptor needs review and source/preparation evidence')
        if not isinstance(r['id'],str) or not re.fullmatch('[A-Za-z0-9_-]{1,64}',r['id']) or r['id'] in seen: raise ValueError('Unique safe receptor ID required')
        seen.add(r['id']); r['mol2']=resolve(r['mol2'])
        if not isinstance(r['center'],list) or len(r['center'])!=3 or any(type(v) not in (int,float) or not math.isfinite(v) for v in r['center']): raise ValueError('Invalid site center')
        if type(r['radius']) not in (int,float) or not 0<r['radius']<=100: raise ValueError('Invalid site radius')
        if sum(1 for _ in iter_mol2_blocks(Path(r['mol2'])))!=1: raise ValueError('Receptor must contain one prepared MOL2 record')
    return c


@exclusive_output
def prepare(sample_report, profile, output):
    sample_dir=Path(sample_report).resolve().parent; check_outputs(sample_dir)
    c=profile_read(profile); output=Path(output).resolve()
    signature=dict(sample_report=str(sample_dir/'report.json'),sample_sha256=sha(sample_dir/'report.json'),
                   profile=c,receptor_hashes={r['id']:sha(Path(r['mol2'])) for r in c['receptors']},
                   implementation_sha256=sha(Path(__file__)))
    if output.exists():
        if not (output/'signature.json').exists() or read(output/'signature.json')!=signature: raise ValueError('Different prepared inputs; use fresh output')
        if (output/'report.json').exists(): return check_outputs(output)
    else: output.mkdir(parents=True)
    save(output/'signature.json',signature); (output/'receptors').mkdir(exist_ok=True)
    for r in c['receptors']: shutil.copyfile(r['mol2'],output/'receptors'/(r['id']+'.mol2'))
    chunks=read(sample_dir/'exports.json')['chunks']
    jobs=[dict(id=f'{r["id"]}-{i:06d}',receptor=r['id'],chunk=chunk) for r in c['receptors'] for i,chunk in enumerate(chunks)]
    save(output/'jobs.json',dict(jobs=jobs)); save(output/'profile.json',c)
    files=['signature.json','jobs.json','profile.json',*[f'receptors/{r["id"]}.mol2' for r in c['receptors']]]
    r=dict(status='complete',kind='block_plants_prepare',readiness='prepared_not_docked',jobs=len(jobs),
           unique_conformers=verified_report(sample_dir/'report.json')['unique_conformers'],receptors=len(c['receptors']),
           sample_report=str(sample_dir/'report.json'),ligand_mode=c['ligand_mode'],ligand_chemistry_reviewed=c['ligand_chemistry_reviewed'],
           output_hashes={f:sha(output/f) for f in files},
           limitations=['Source-backed ligand atom types and states are preserved; protonation suitability is not inferred from successful export.',
                        'ChemPLP docking scores are not binding affinity or search recall.',
                        'Reference redocking/cross-docking validation is separate; this is an exploratory panel.'])
    save(output/'report.json',r); return r


def parse_ranking(path, records):
    """Original-title aliases plus 1-based entry indices prevent filename guessing."""
    by_alias={r['alias']:r for r in records}; found={}
    with Path(path).open(encoding='utf-8-sig',newline='') as f:
        reader=csv.reader(f)
        fields=next(reader,[])
        normalized=[re.sub('[^A-Z0-9]','',k.upper()) for k in fields]
        if len(set(normalized))!=len(normalized) or 'TOTALSCORE' not in normalized:
            raise ValueError('Unknown or duplicate PLANTS ranking header')
        names=[i for i,k in enumerate(normalized) if k in ('LIGANDENTRY','LIGAND','NAME','PLANTSLIGANDDESCRIPTION')]
        implicit=not names and bool(normalized) and normalized[0]=='TOTALSCORE'
        if len(names)>1 or (not names and not implicit): raise ValueError('Unknown PLANTS ranking header')
        name_index=names[0] if names else 0
        score_index=normalized.index('TOTALSCORE')+int(implicit)
        for row in reader:
            if not row: continue
            if len(row)!=len(fields)+int(implicit): raise ValueError('PLANTS ranking row width mismatch')
            name=row[name_index].strip(); m=re.fullmatch(r'(c\d{9})_entry_(\d+)_conf_(\d+)(?:\.mol2)?',name)
            if not m or m[1] not in by_alias: raise ValueError('Unrecognized PLANTS pose identity: '+name)
            entry=int(m[2]); record=by_alias[m[1]]
            if not 1<=entry<=len(records) or records[entry-1]['alias']!=m[1]: raise ValueError('PLANTS entry/title mismatch')
            value=float(row[score_index])
            if not math.isfinite(value) or m[1] in found: raise ValueError('Duplicate or nonfinite best-ranking score')
            found[m[1]]=dict(cid=record['cid'],alias=m[1],score=value,pose_name=name.removesuffix('.mol2'))
    return found


def collect_pose_rows(attempt, job):
    rows=[]
    found=parse_ranking(attempt/'poses/bestranking.csv',job['chunk']['records'])
    pose_locations={}
    for p in sorted((attempt/'poses').glob('*.mol2')):
        for idx,text in iter_mol2_blocks(p):
            title=text.splitlines()[1].strip().removesuffix('.mol2')
            if title in {r['pose_name'] for r in found.values()}:
                if title in pose_locations: raise ValueError('Ambiguous output pose title')
                pose_locations[title]=(str(p.resolve()),idx)
    for record in job['chunk']['records']:
        hit=found.get(record['alias']); pose=pose_locations.get(hit['pose_name']) if hit else None
        rows.append(dict(cid=record['cid'],alias=record['alias'],receptor=job['receptor'],job=job['id'],
                         status='ok' if pose else 'missing_pose' if hit else 'missing_score',
                         score=hit['score'] if pose else None,pose_name=hit['pose_name'] if hit else '',
                         pose_file=pose[0] if pose else '',pose_index=pose[1] if pose else None))
    return rows


def plants_config(c, receptor):
    lines=['scoring_function chemplp','search_speed '+c['search_speed'],'protein_file protein.mol2','ligand_file ligands.mol2',
           'output_dir poses','bindingsite_center '+' '.join(map(str,receptor['center'])),f'bindingsite_radius {receptor["radius"]}',
           f'rigid_ligand {int(c["ligand_mode"]=="rigid")}',f'cluster_structures {c["cluster_structures"]}',
           'cluster_rmsd 2.0','keep_original_mol2_description 1','write_multi_mol2 1','merge_multi_conf_output 0']
    return '\n'.join(lines)+'\n'


@exclusive_output
def run(prepared_report, output, *, max_jobs=None, workers=None):
    if workers is not None and (type(workers) is not int or not 1<=workers<=64):
        raise ValueError("workers must be an integer from 1 to 64")
    prepared=Path(prepared_report).resolve().parent; prep=check_outputs(prepared)
    sample_dir=Path(prep['sample_report']).parent; check_outputs(sample_dir)
    c=read(prepared/'profile.json'); output=Path(output).resolve()
    effective_workers=c['workers'] if workers is None else workers
    executable=Path(c['executable'])
    if not executable.is_file(): raise ValueError('Configured PLANTS executable is missing')
    signature=dict(prepared_report=str(prepared/'report.json'),prepared_sha256=sha(prepared/'report.json'),
                   executable=str(executable),executable_sha256=sha(executable),implementation_sha256=sha(Path(__file__)))
    if output.exists():
        if not (output/'signature.json').exists() or read(output/'signature.json')!=signature: raise ValueError('Run inputs changed; fresh output required')
    else: output.mkdir(parents=True)
    save(output/'signature.json',signature)
    jobs=read(prepared/'jobs.json')['jobs']
    if max_jobs is not None and (type(max_jobs) is not int or max_jobs<1): raise ValueError('max_jobs must be positive')
    selected=jobs if max_jobs is None else jobs[:max_jobs]
    def execute(job):
        directory=output/job['id']; directory.mkdir(exist_ok=True); receipt=directory/'receipt.json'
        if receipt.exists():
            saved=read(receipt)
            if saved['status']=='complete':
                for p,h in saved['output_hashes'].items():
                    if sha(directory/p)!=h: raise ValueError('Completed docking output changed')
                return saved
        # Fresh attempts avoid mixing PLANTS outputs from interrupted/failed runs.
        number=1
        while (directory/f'attempt-{number:03d}').exists(): number+=1
        attempt=directory/f'attempt-{number:03d}'; attempt.mkdir()
        receptor=next(r for r in c['receptors'] if r['id']==job['receptor'])
        shutil.copyfile(prepared/'receptors'/(job['receptor']+'.mol2'),attempt/'protein.mol2')
        ligand=sample_dir/job['chunk']['path']
        if sha(ligand)!=job['chunk']['sha256']: raise ValueError('Ligand chunk changed')
        shutil.copyfile(ligand,attempt/'ligands.mol2')
        (attempt/'plantsconfig').write_text(plants_config(c,receptor),encoding='ascii')
        start=time.monotonic(); rows=[]; error=None
        try:
            with (attempt/'plants.log').open('w',encoding='utf-8') as log:
                result=subprocess.run([str(executable),'--mode','screen','plantsconfig'],cwd=attempt,stdout=log,
                                      stderr=subprocess.STDOUT,timeout=c['timeout_seconds'],shell=False)
            if result.returncode: raise RuntimeError('PLANTS exit code '+str(result.returncode))
            rows=collect_pose_rows(attempt,job)
        except (OSError, ValueError, RuntimeError, subprocess.TimeoutExpired) as exc:
            error=str(exc)
            rows=[dict(cid=r['cid'],alias=r['alias'],receptor=job['receptor'],job=job['id'],status='failed',score=None,pose_name='',pose_file='',pose_index=None) for r in job['chunk']['records']]
        save(attempt/'scores.json',rows)
        files=[p for p in attempt.rglob('*') if p.is_file()]
        saved=dict(status='complete' if error is None and all(r['status']=='ok' for r in rows) else 'failed',job=job['id'],error=error,
                   rows=rows,seconds=time.monotonic()-start,output_hashes={str(p.relative_to(directory)):sha(p) for p in files})
        save(receipt,saved); return saved
    # Validate a real canary before submitting the remaining panel to the engine.
    first=execute(selected[0]) if selected else None
    gate=bool(first and first['status']=='complete')
    print('PLANTS first-job gate: '+('passed' if gate else 'failed'),flush=True)
    if gate:
        with ThreadPoolExecutor(max_workers=effective_workers) as pool:
            for i,result in enumerate(pool.map(execute,selected[1:]),2):
                print(f'PLANTS {i}/{len(selected)}: {result["job"]} {result["status"]}',flush=True)
                save(output/'progress.json',dict(completed_jobs=i,requested_jobs=len(selected),last_job=result['job']))
    return summarize_run(prepared,prep,c,output,jobs,gate,workers=effective_workers)


def summarize_run(prepared, prep, c, output, jobs, gate, *, workers=None):
    all_rows=[]; done=0; failures=0
    for job in jobs:
        path=output/job['id']/'receipt.json'
        if path.exists():
            receipt=read(path)
            for name,digest in receipt['output_hashes'].items():
                if sha(path.parent/name)!=digest: raise ValueError('Docking receipt hash mismatch')
            all_rows.extend(receipt['rows']); done+=1; failures+=receipt['status']!='complete'
        else:
            all_rows.extend(dict(cid=r['cid'],alias=r['alias'],receptor=job['receptor'],job=job['id'],status='not_run',score=None,pose_name='',pose_file='',pose_index=None) for r in job['chunk']['records'])
    fields=['cid','alias','receptor','job','status','score','pose_name','pose_file','pose_index']
    write_csv(output/'scores.csv',all_rows,fields)
    status='complete' if done==len(jobs) and not failures else 'partial'
    hashes={'scores.csv':sha(output/'scores.csv'),'signature.json':sha(output/'signature.json')}
    if (output/'recovery.json').exists(): hashes['recovery.json']=sha(output/'recovery.json')
    # Seal receipts and pose files too: reuse is bound to coordinates, not just scores.
    for job in jobs:
        receipt=output/job['id']/'receipt.json'
        if receipt.exists():
            hashes[str(receipt.relative_to(output))]=sha(receipt)
            hashes.update({str(Path(job['id'])/name):digest for name,digest in read(receipt)['output_hashes'].items()})
    r=dict(status=status,kind='block_plants_run',prepared_report=str(prepared/'report.json'),sample_report=prep['sample_report'],
           jobs=len(jobs),attempted_jobs=done,failed_jobs=failures,scored=sum(r['status']=='ok' for r in all_rows),
           first_job_gate='passed' if gate else 'failed',
           ligand_mode=c['ligand_mode'],score='ChemPLP TOTAL_SCORE (lower is better)',workers_this_invocation=workers,output_hashes=hashes,
           limitations=['Docking is stochastic; sample seed does not seed the PLANTS optimizer.',
                        'Missing/failed/not-run scores are excluded, never zero-imputed.'])
    save(output/'report.json',r); return r


@exclusive_output
def recover(run_report, output, *, include_checkpoints=False):
    from .prompt_workflow import file_lock
    root=Path(run_report).resolve().parent
    with file_lock(root.parent/('.'+root.name+'.e097.lock')):
        return _recover(run_report,output,include_checkpoints=include_checkpoints)


def _recover(run_report, output, *, include_checkpoints=False):
    """Reparse sealed successful-engine outputs into a fresh resumable run; never execute PLANTS."""
    source=Path(run_report).resolve(); old=read(source); root=source.parent
    if old.get('kind')!='block_plants_run' or old.get('status') not in ('complete','partial'):
        raise ValueError('Existing PLANTS run report required')
    def sealed(base, name, digest):
        path=(base/name).resolve()
        if not path.is_relative_to(base.resolve()) or not path.is_file() or sha(path)!=digest:
            raise ValueError('Recovery source hash/path mismatch: '+name)
        return path
    for name,digest in old['output_hashes'].items(): sealed(root,name,digest)
    signature=read(root/'signature.json'); prepared=Path(old['prepared_report']).resolve().parent
    if signature.get('prepared_report')!=str(prepared/'report.json') or signature.get('prepared_sha256')!=sha(prepared/'report.json'):
        raise ValueError('Recovery preparation signature mismatch')
    prep=check_outputs(prepared); check_outputs(Path(prep['sample_report']).parent)
    c=read(prepared/'profile.json'); executable=Path(c['executable'])
    if signature.get('executable')!=str(executable) or signature.get('executable_sha256')!=sha(executable):
        raise ValueError('Recovery executable changed')
    jobs=read(prepared/'jobs.json')['jobs']; plans=[]
    for job in jobs:
        receipt_path=root/job['id']/'receipt.json'
        if not receipt_path.exists(): continue
        rel=str(receipt_path.relative_to(root))
        if rel not in old['output_hashes'] and not include_checkpoints: raise ValueError('Unsealed job receipt; stopped-run migration requires --include-checkpoints')
        receipt=read(receipt_path)
        if rel not in old['output_hashes'] and receipt.get('status')!='complete':
            continue  # Uncommitted/failed work will be executed normally in the new run.
        if receipt.get('job')!=job['id']: raise ValueError('Recovery job identity mismatch')
        if receipt.get('status')!='complete' and receipt.get('error')!='Unknown PLANTS ranking header':
            raise ValueError('Recovery only accepts completed jobs or the known post-execution ranking-header failure')
        for name,digest in receipt['output_hashes'].items(): sealed(receipt_path.parent,name,digest)
        rankings=[Path(name) for name in receipt['output_hashes'] if Path(name).parts[-2:]==('poses','bestranking.csv')]
        if len(rankings)!=1: raise ValueError('Expected one sealed ranking attempt')
        attempt=rankings[0].parent.parent
        if len(attempt.parts)!=1 or not re.fullmatch(r'attempt-\d+',attempt.name): raise ValueError('Invalid attempt path')
        for name in (attempt/'protein.mol2',attempt/'ligands.mol2',attempt/'plantsconfig',attempt/'plants.log'):
            if str(name) not in receipt['output_hashes']: raise ValueError('Missing sealed engine input/log')
        if sha(receipt_path.parent/attempt/'protein.mol2')!=sha(prepared/'receptors'/(job['receptor']+'.mol2')):
            raise ValueError('Recovery receptor mismatch')
        if sha(receipt_path.parent/attempt/'ligands.mol2')!=job['chunk']['sha256']:
            raise ValueError('Recovery ligand mismatch')
        receptor=next(r for r in c['receptors'] if r['id']==job['receptor'])
        if (receipt_path.parent/attempt/'plantsconfig').read_text()!=plants_config(c,receptor):
            raise ValueError('Recovery docking configuration mismatch')
        # Require all MOL2 files read by the parser to be covered by the old receipt.
        for pose in (receipt_path.parent/attempt/'poses').glob('*.mol2'):
            if str(pose.relative_to(receipt_path.parent)) not in receipt['output_hashes']:
                raise ValueError('Unsealed pose output')
        rows=collect_pose_rows(receipt_path.parent/attempt,job)
        if not rows or any(row['status']!='ok' for row in rows): raise ValueError('Recovery has missing scores or poses')
        plans.append((job,receipt,attempt))
    if not plans: raise ValueError('No recoverable jobs')
    output=Path(output).resolve()
    if output.exists(): raise ValueError('Recovery requires a fresh output directory')
    output.mkdir(parents=True)
    signature['implementation_sha256']=sha(Path(__file__));save(output/'signature.json',signature)
    save(output/'recovery.json',dict(source_report=str(source),source_report_sha256=sha(source),
         recovered_jobs=[job['id'] for job,_,_ in plans],engine_executed=False,include_checkpoints=include_checkpoints))
    for job,receipt,attempt in plans:
        directory=output/job['id'];directory.mkdir()
        for name,digest in receipt['output_hashes'].items():
            src=sealed(root/job['id'],name,digest);dst=directory/name;dst.parent.mkdir(parents=True,exist_ok=True)
            shutil.copyfile(src,dst)
            if sha(dst)!=digest: raise ValueError('Recovery copy hash mismatch')
        rows=collect_pose_rows(directory/attempt,job)
        score_path=attempt/'scores.json';save(directory/score_path,rows)
        hashes=dict(receipt['output_hashes']);hashes[str(score_path)]=sha(directory/score_path)
        recovered=dict(receipt,status='complete',error=None,rows=rows,
                       output_hashes=hashes,
                       recovery_source_receipt=str(root/job['id']/'receipt.json'))
        save(directory/'receipt.json',recovered)
    gate=(output/jobs[0]['id']/'receipt.json').exists()
    return summarize_run(prepared,prep,c,output,jobs,gate)


@exclusive_output
def analyze(run_report, output):
    run_dir=Path(run_report).resolve().parent; r=read(run_report)
    if r.get('kind')!='block_plants_run' or r.get('status') not in ('complete','partial'): raise ValueError('PLANTS run report required')
    for name,h in r['output_hashes'].items():
        if sha(run_dir/name)!=h: raise ValueError('Docking result changed')
    samples=Path(r['sample_report']).parent; check_outputs(samples)
    output=Path(output).resolve(); signature=dict(run_report=str(Path(run_report).resolve()),sha256=sha(Path(run_report)),implementation_sha256=sha(Path(__file__)))
    if output.exists():
        if not (output/'signature.json').exists() or read(output/'signature.json')!=signature: raise ValueError('Use fresh analysis output')
    else: output.mkdir(parents=True)
    save(output/'signature.json',signature)
    with (run_dir/'scores.csv').open(encoding='utf-8',newline='') as f: scores=list(csv.DictReader(f))
    indexed={}
    for row in scores:
        key=(row['cid'],row['receptor'])
        if key in indexed: raise ValueError('Duplicate conformer/receptor result')
        if row['status']=='ok' and not math.isfinite(float(row['score'])): raise ValueError('Nonfinite successful score')
        indexed[key]=row
    receptors=sorted({r['receptor'] for r in scores}); distributions=defaultdict(list); totals=defaultdict(int); joined=[]
    with readonly(samples/'samples.sqlite') as db:
        known={r[0] for r in db.execute('SELECT cid FROM selected')}
        if any(cid not in known for cid,_ in indexed): raise ValueError('Unknown score conformer')
        for scheme,bid,cid in db.execute('SELECT scheme,block_id,cid FROM sample'):
            for receptor in receptors:
                key=(scheme,bid,receptor); totals[key]+=1
                row=indexed.get((cid,receptor),dict(cid=cid,receptor=receptor,status='missing',score='',pose_name='',pose_file='',pose_index=''))
                if row['status']=='ok': distributions[key].append(float(row['score']))
                joined.append({k:row.get(k,'') for k in ['cid','receptor','status','score','pose_name','pose_file','pose_index']}|dict(scheme=scheme,block_id=bid))
        populations={(s,b):n for s,b,n in db.execute('SELECT scheme,block_id,n FROM population')}
    summary=[]
    for (scheme,bid,receptor),n in sorted(totals.items()):
        values=distributions[(scheme,bid,receptor)]; pop=populations[(scheme,bid)]
        row=dict(scheme=scheme,block_id=bid,receptor=receptor,population=pop,sampled=n,scored=len(values),missing=n-len(values),sample_fraction=n/pop)
        row.update(dict(zip(['minimum','p10','median','p90','maximum'],map(float,np.quantile(values,[0,.1,.5,.9,1])))) if values else dict.fromkeys(['minimum','p10','median','p90','maximum'],''))
        row['mean']=float(np.mean(values)) if values else ''; summary.append(row)
    comparisons=[]
    for scheme in sorted({s for s,_,_ in totals}):
        for receptor in receptors:
            rows=[v for v in summary if v['scheme']==scheme and v['receptor']==receptor]; complete=all(v['missing']==0 for v in rows)
            comparisons.append(dict(scheme=scheme,receptor=receptor,blocks=len(rows),sampled=sum(v['sampled'] for v in rows),
                                    scored=sum(v['scored'] for v in rows),population_weighted_mean=sum(v['population']*v['mean'] for v in rows)/sum(v['population'] for v in rows) if complete else None,
                                    estimate_status='complete_panel' if complete else 'withheld_missing_scores'))
    write_csv(output/'block_summary.csv',summary,list(summary[0]) if summary else ['scheme','block_id'])
    write_csv(output/'block_scores.csv',joined,['scheme','block_id','cid','receptor','status','score','pose_name','pose_file','pose_index'])
    text=['# Block docking evaluation','',f'Engine score: {r["score"]}. Mode: {r["ligand_mode"]}.','',
          'Results are separated by receptor. Lower ChemPLP scores are more favorable under this model; they are not measured affinity.','']
    for x in comparisons: text.append(f'- {x["scheme"]} / {x["receptor"]}: {x["scored"]}/{x["sampled"]} scored across {x["blocks"]} blocks; population-weighted mean: {x["population_weighted_mean"]}.')
    text+=['','No block rejection or significance claim. Equal block sampling, reused conformers and shared molecules induce unequal weights/dependence.',
           'Use block_scores.csv pose_file/pose_index with cid/receptor/pose_name for subsequent N-E rescoring. No rescore model was executed.']
    (output/'review.md').write_text('\n'.join(text)+'\n',encoding='utf-8')
    result=dict(status='complete',kind='block_analyze',docking_status=r['status'],comparisons=comparisons,
                output_hashes={f:sha(output/f) for f in ['block_summary.csv','block_scores.csv','review.md','signature.json']},
                limitations=['Exploratory panel, not held-out retrieval recall or affinity validation.',
                             'No causal comparison or independent significance claim across reused samples.'])
    save(output/'report.json',result); return result


def main():
    p=argparse.ArgumentParser(description=__doc__); sub=p.add_subparsers(dest='action',required=True)
    for name in ('prepare','run','recover','analyze'):
        s=sub.add_parser(name); s.add_argument('--source',required=True); s.add_argument('--output',required=True)
        if name=='prepare': s.add_argument('--profile',required=True)
        if name=='recover': s.add_argument('--include-checkpoints',action='store_true',help='Import verified completed receipts after the last report; stop the old process first')
        if name=='run':
            s.add_argument('--max-jobs',type=int)
            s.add_argument('--workers',type=int,help='Concurrent PLANTS processes, 1..64; scheduling only')
    a=p.parse_args()
    result=prepare(a.source,a.profile,a.output) if a.action=='prepare' else run(a.source,a.output,max_jobs=a.max_jobs,workers=a.workers) if a.action=='run' else recover(a.source,a.output,include_checkpoints=a.include_checkpoints) if a.action=='recover' else analyze(a.source,a.output)
    print(json.dumps(result,indent=2))


if __name__=='__main__': main()
