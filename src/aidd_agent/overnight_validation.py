"""Uncapped whole-library condition run with before/after correctness audits."""
import argparse
import copy
import json
from pathlib import Path
import time

import numpy as np

from . import full_library_screen as full
from .funnel_benchmark import run as pilot
from .prompt_workflow import file_lock
from .survivor_validation import validate_survivors


def sample_representatives(path, limit=128):
    """Deterministic reservoir across the full uncapped molecule result list."""
    rng = np.random.default_rng(49)
    selected = []
    with Path(path).open(encoding='utf-8') as stream:
        for i,line in enumerate(stream):
            gid = int(json.loads(line)['global_id'])
            if i < limit: selected.append(gid)
            else:
                j = int(rng.integers(i+1))
                if j < limit: selected[j] = gid
    return sorted(selected)


def run(selection, output, workers=22, chunk_size=2048):
    selection, output = Path(selection).resolve(), Path(output).resolve()
    selected = full.ev.read(selection)
    full.ev.require(selected.get('kind')=='selection_preview','Use a selection_preview report, not an execution wrapper')
    policy = full.normalize_policy(selected['policy'])
    full.ev.require('coarse_constraints' in policy,'Full validation requires the explicit joint rule')
    full.check_hashes(selected['sources'])
    protected = [selection.parent, Path(selected['evidence_report']).resolve().parent,
                 Path(selected['query']['artifact_catalog']).resolve().parent,
                 Path(selected['query']['chemical_companion']).resolve().parent]
    full.ev.require(all(not output.is_relative_to(p) and not p.is_relative_to(output) for p in protected),
                    'Output overlaps source evidence or library')
    full.ev.require(workers>0 and chunk_size>0,'Positive workers and chunk size required')
    output.mkdir(parents=True,exist_ok=True)
    with file_lock(output/'suite.lock'):
        protocol=dict(selection=full.fingerprint([selection]),code=full.fingerprint(sorted(Path(__file__).parent.glob('*.py'))),
                      workers=workers,chunk_size=chunk_size,policy=policy)
        protocol_path=output/'suite-protocol.json'
        if protocol_path.exists():
            full.ev.require(full.ev.read(protocol_path)==protocol,'Changed inputs/code/settings: use a new output directory')
        else: full._atomic_json(protocol_path,protocol)
        report=dict(status='running',stage='preflight_validation',selection=str(selection),policy=policy,
            full_library=True,candidate_top_k=None,gaussian_top_n=None,stages={},
            scope='Current selected query and explicit rule; not exhaustive orientations, biology or docking validation')
        def save(): full._atomic_json(output/'suite-report.json',report)
        save()
        try:
            # Unique directories keep interrupted small audits inspectable; the heavy
            # full scan alone resumes verified chunks in its stable directory.
            audit=output/f'audit-{time.time_ns()}'
            before=pilot(selection,audit,count=10000,validate_survivors=True,include_gpu=False)
            report['stages']['before']=str(audit/'report.json')
            full.ev.require(before['status']=='complete','Pre-run correctness/control failure; inspect audit report')
            report.update(stage='full_library'); save()
            full.ev.log('Starting uncapped whole-library joint screening')
            result=full.run_funnel(selection,output/'full-library',workers=workers,chunk_size=chunk_size,
                                   backend='numpy',pose_feasibility=True)
            report['stages']['full_library']=str(output/'full-library'/'report.json')
            full.ev.require(result['status']=='complete' and all(q['full_coverage'] for q in result['queries']),
                            'Whole-library coverage incomplete')
            report['library']=result['library']
            report['queries']=[dict(query_id=q['query_id'],full_coverage=q['full_coverage'],counts=q['counts']) for q in result['queries']]
            report.update(stage='uncapped_molecule_list'); save()
            selection_result=full.preview(output/'full-library'/'report.json',output/'all-matches',policy)
            report['stages']['all_matches']=str(output/'all-matches'/'report.json')
            report['representatives_jsonl']=selection_result['representatives_jsonl']
            ids=sample_representatives(selection_result['representatives_jsonl'])
            report.update(stage='final_hit_validation'); save()
            if ids:
                q=copy.deepcopy(selected['query']); q['condition_policy']=policy
                indices=sorted({a['feature_index'] for a in q['anchors']})
                for a in q['anchors']: a['score_column']=indices.index(a['feature_index'])
                full.initialize(q)
                after_dir=output/f'final-hit-audit-{time.time_ns()}'; after_dir.mkdir()
                after=validate_survivors(q,ids,[],after_dir)
                full._atomic_json(after_dir/'report.json',after)
                report['stages']['final_hit_validation']=str(after_dir/'report.json')
                full.ev.require(after['status']=='passed' and after['reference_positive_count']==len(ids),
                                'Final hit sample failed reference reproduction')
                report['final_hit_validation']=dict(sampled_molecules=len(ids),status='passed',scope='Sampled final representatives only')
            else:
                report['final_hit_validation']=dict(sampled_molecules=0,status='not_covered_no_final_hits')
            report.update(status='complete',stage='finished'); save()
            return report
        except Exception as exc:
            report.update(status='failed',error=str(exc)); save()
            raise


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--selection',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--workers',type=int,default=22)
    parser.add_argument('--chunk-size',type=int,default=2048)
    args=parser.parse_args()
    run(args.selection,args.output,args.workers,args.chunk_size)


if __name__=='__main__': main()
