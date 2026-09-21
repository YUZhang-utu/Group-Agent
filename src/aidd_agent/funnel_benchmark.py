"""Small spread-out library pilot: correctness before CPU/CUDA throughput claims."""
from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
import copy
import json
import multiprocessing
import platform
from pathlib import Path
import time

import numpy as np

from . import full_library_screen as full
from .pose_acceleration import hardware_options, array_module


def sample_ids(total, count=256):
    if total < 1 or count < 1: raise ValueError('Positive pilot size required')
    return np.unique(np.linspace(0,total-1,min(total,count),dtype=np.int64))


def work(ids):
    reader, original, expanded, query = full._STATE
    started = time.perf_counter()
    records = [full._FILTER_READER.get(int(gid)) for gid in ids]
    decisions = [full._BOUND.check(c) for c in records]
    prefilter_seconds = time.perf_counter()-started
    tick = time.perf_counter()
    # Score even bound rejects in this pilot to check false rejection on real rows.
    rigid = full._score_ids(reader,original,ids,backend=query['pose_backend'],**full.POSE_PARAMETERS)
    pose_seconds = time.perf_counter()-tick
    tick = time.perf_counter()
    _, assignments, scores = full.score_batched(Path(query['artifact_catalog']),Path(query['chemical_companion']),
        expanded,ids,rigid['molecule_ids'],rigid['conformer_ids'],[full.OBJECTIVE],
        {full.OBJECTIVE:rigid[full.OBJECTIVE+'__transform']},sigma=1.,cutoff=4.5,angular_power=2.)
    rigid[full.OBJECTIVE+'__anchor_scores'] = scores[full.OBJECTIVE]
    rigid[full.OBJECTIVE+'__anchor_assignments'] = assignments[full.OBJECTIVE]
    columns = sorted({a['score_column'] for a in query['anchors'] if a['anchor_id'] in query['condition_policy']['required_anchors']})
    from .necessary_conditions import pose_mask
    rigid['condition_passed'] = pose_mask(rigid,columns,query['condition_policy'],full.OBJECTIVE)
    keep = np.array([d[0] for d in decisions])
    full.ev.require(not np.any(rigid['condition_passed'] & ~keep),'Prefilter rejected a passing reference pose')
    return rigid, dict(prefilter_seconds=prefilter_seconds,pose_seconds=pose_seconds,
        annotation_seconds=time.perf_counter()-tick,prefilter_passed=int(keep.sum()),
        rejection_reasons=dict(Counter(d[1] for d in decisions if not d[0])))


def compare(reference, actual):
    differences = []
    max_error = 0.
    for key in reference:
        if key not in actual or reference[key].shape != actual[key].shape:
            differences.append(key);continue
        a,b = reference[key],actual[key]
        if a.dtype.kind in 'fc':
            max_error=max(max_error,float(np.max(np.abs(a-b),initial=0)))
        # Reference recomputation is expected to preserve all arrays exactly.
        if not np.array_equal(a,b): differences.append(key)
    return dict(passed=not differences,different_arrays=differences,max_absolute_error=max_error)


def _run(selection, output, *, count=256, repeats=2, include_gpu=True):
    output=Path(output).resolve(); selection=Path(selection).resolve()
    selected=full.ev.read(selection)
    full.ev.require(selected.get('kind')=='selection_preview','Benchmark needs a saved selection preview')
    full.check_hashes(selected['sources'])
    source=Path(selected['evidence_report']).resolve()
    full.ev.require(str(source) in selected['sources'],'Unsealed classification')
    evidence=full.ev.read(source)
    full.ev.require(evidence.get('kind')=='screening_evidence' and evidence.get('classification'),
                    'Use a preview from completed classification')
    policy=full.normalize_policy(selected['policy'])
    owners={a['anchor_id']:q for q in evidence['queries'] for a in q['anchors']}
    full.ev.require(set(policy['required_anchors']) <= owners.keys(),'Unknown rule anchors')
    qs=[owners[a] for a in policy['required_anchors']]
    full.ev.require(len({q['query_id'] for q in qs})==1,'Use one query per rule')
    q=copy.deepcopy(qs[0]);q['condition_policy']=policy
    for protected in (source.parent,selection.parent,Path(q['artifact_catalog']).resolve().parent,
                      Path(q['chemical_companion']).resolve().parent):
        full.ev.require(not output.is_relative_to(protected) and not protected.is_relative_to(output),
                        'Benchmark output overlaps inputs')
    full.ev.require(repeats>0 and count>0,'Positive pilot count/repeats required')
    output.mkdir(parents=True,exist_ok=True)
    full.ev.require(not (output/'report.json').exists(),'Use a new benchmark output directory')
    indices=sorted({a['feature_index'] for a in q['anchors']})
    for a in q['anchors']:a['score_column']=indices.index(a['feature_index'])
    ids=sample_ids(full.catalog_count(full.ev.read(q['artifact_catalog'])),count)
    _,cpu_workers=hardware_options('numpy',None)
    # Compare equal concurrency first, then the hardware-sized CPU pool.
    scenarios=[('reference',min(8,cpu_workers)),('numpy',min(8,cpu_workers))]
    if cpu_workers != min(8,cpu_workers):scenarios.append(('numpy',cpu_workers))
    report=dict(kind='hardware_benchmark',status='running',query_id=q['query_id'],condition_policy=policy,
        environment=dict(host=platform.node(),platform=platform.platform(),numpy=np.__version__),
        sample_count=len(ids),sample_ids=ids.tolist(),sample_scope='Evenly spread global IDs; exploratory, not a full-library count',
        scenarios=[],sources={**selected['sources'],**full.fingerprint([selection,source,Path(q['query_npz']),Path(q['artifact_catalog']),Path(q['chemical_companion'])])},
        code=full.fingerprint(sorted(Path(__file__).parent.glob('*.py'))),e031_changes_ranking=False,
        numerical_policy='Float64 seed batches; CPU reference recheck of contenders within 1e-10 relative/absolute objective guard; sampled equivalence only',
        biological_quality='not_evaluated',timing_scope='Worker startup, sample reads, bounds, all sample poses, E031 and IPC; no full integrity scan or output serialization')
    if include_gpu:
        try:
            cp=array_module('cupy')
            cp.zeros(1).sum().item()
            report['gpu']=dict(status='available',cupy=cp.__version__,device=str(cp.cuda.Device()))
            scenarios.append(('cupy',1))
        except Exception as exc:
            report['gpu']=dict(status='unavailable',error=str(exc))
    else:report['gpu']=dict(status='not_requested')
    reference=None
    for backend,workers in scenarios:
        row=dict(backend=backend,workers=workers,timings=[],checks=[],status='running')
        report['scenarios'].append(row)
        for repeat in range(repeats):
            q['pose_backend']=backend
            full.ev.log(f'Hardware pilot: {backend}, {workers} workers, repeat {repeat+1}/{repeats}, {len(ids)} spread IDs')
            started=time.perf_counter()
            batches=[ids[i:i+8] for i in range(0,len(ids),8)]
            if backend=='cupy':
                full.initialize(q)
                outputs=[work(batch) for batch in batches]
            else:
                with ProcessPoolExecutor(max_workers=workers,initializer=full.initialize,initargs=(q,),
                                         mp_context=multiprocessing.get_context('spawn')) as pool:
                    outputs=list(pool.map(work,batches))
            elapsed=time.perf_counter()-started
            arrays={key:(outputs[0][0][key] if key=='objective_names' else np.concatenate([r[0][key] for r in outputs])) for key in outputs[0][0]}
            if reference is None:reference=arrays
            row['timings'].append(elapsed);row['checks'].append(compare(reference,arrays))
            row['prefilter_passed']=sum(r[1]['prefilter_passed'] for r in outputs)
            row['matching_conformers']=int(arrays['condition_passed'].sum())
            reasons=Counter()
            for _,stats in outputs:reasons.update(stats['rejection_reasons'])
            row['rejection_reasons']=dict(reasons)
            row['worker_seconds']={k:sum(r[1][k] for r in outputs) for k in ('prefilter_seconds','pose_seconds','annotation_seconds')}
            full._atomic_json(output/'report.json',report)
        row['median_seconds']=float(np.median(row['timings']))
        row['conformers_per_second']=len(ids)/row['median_seconds']
        row['status']='passed' if all(c['passed'] for c in row['checks']) else 'failed'
    valid=[r for r in report['scenarios'] if r['status']=='passed']
    best=min(valid,key=lambda r:r['median_seconds']) if valid else None
    report['recommendation']=dict(backend=best['backend'],workers=best['workers'],scope='Fastest equivalent sampled configuration; validate larger pilot before full run') if best else None
    report['status']='complete' if all(r['status']=='passed' for r in report['scenarios']) else 'failed'
    full.check_hashes(report['sources'])
    full._atomic_json(output/'report.json',report)
    full._atomic_json(output/'RUN_STATUS.json',dict(status=report['status'],report_sha256=full.ev.sha(output/'report.json')))
    return report


def run(selection, output, **kwargs):
    report_path=Path(output)/'report.json'
    existed=report_path.exists()
    try:
        return _run(selection,output,**kwargs)
    except Exception as exc:
        # Preserve pre-existing runs and protected-input failures. Only update a
        # report created by this invocation, so a failed accelerator is visible.
        if not existed and report_path.exists():
            report=full.ev.read(report_path)
            report.update(status='failed',error=str(exc))
            full._atomic_json(report_path,report)
            full._atomic_json(report_path.parent/'RUN_STATUS.json',dict(status='failed',report_sha256=full.ev.sha(report_path)))
        raise


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--selection',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--count',type=int,default=256)
    p.add_argument('--repeats',type=int,default=2)
    p.add_argument('--no-gpu',action='store_true')
    a=p.parse_args();result=run(a.selection,a.output,count=a.count,repeats=a.repeats,include_gpu=not a.no_gpu)
    print(json.dumps(dict(status=result['status'],report=str(a.output/'report.json'))))
    return 0 if result['status']=='complete' else 1


if __name__=='__main__':raise SystemExit(main())
