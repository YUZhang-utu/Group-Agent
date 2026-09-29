"""Matched-target, equal-query descriptive analysis; no membership acceptance."""
import argparse
from collections import defaultdict,Counter
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

METRICS=('property_rms','steric_rms','maximum_local_steric_difference','backbone_rmsd_angstrom')


def describe(values):
    a=np.array(values,float)
    return dict(n=len(a),median=float(np.median(a)),p10=float(np.quantile(a,.1)),p90=float(np.quantile(a,.9))) if len(a) else dict(n=0)


def table(path,rows):
    if not rows:return
    with path.open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)


def run(source,output,expected_hash=None):
    source=Path(source);output=Path(output);digest=hashlib.sha256(source.read_bytes()).hexdigest()
    if expected_hash and digest!=expected_hash:raise ValueError('CSV hash mismatch')
    with source.open(newline='',encoding='utf-8-sig') as f:rows=list(csv.DictReader(f))
    seen=set();controls=defaultdict(list);pairs=defaultdict(list);by_band=defaultdict(list)
    for r in rows:
        for m in METRICS:
            r[m]=float(r[m])
            if not np.isfinite(r[m]) or r[m]<0:raise ValueError('Invalid metric')
        identity=tuple(r[k] for k in ('kind','query','reference','target_group','rotation'))
        if identity in seen:raise ValueError('Duplicate pair row')
        seen.add(identity)
        if r['kind']=='within_definite_control':controls[r['target_group']].append(r)
        elif r['kind']=='boundary_vs_reference':
            pairs[(r['band'],r['query'],r['target_group'],r['rotation'])].append(r)
            by_band[r['band']].append(r)
        else:raise ValueError('Unknown pair type')
    control_stats={g:{m:dict(median=float(np.median([r[m] for r in rs])),maximum=max(r[m] for r in rs)) for m in METRICS}
                   for g,rs in controls.items()}
    matched=[];missing=[];queries=defaultdict(list)
    for (b,q,g,rotation),rs in sorted(pairs.items()):
        if g not in controls:
            missing.append(dict(band=b,query=q,target_group=g,rotation=rotation,references=len(rs)));continue
        value=dict(band=b,query=q,target_group=g,rotation=rotation,references=len(rs),control_pairs=len(controls[g]))
        within=True
        for m in METRICS:
            med=float(np.median([r[m] for r in rs]));stats=control_stats[g][m]
            value[m]=med;value[m+'_control_median']=stats['median'];value[m+'_delta']=med-stats['median']
            value[m+'_within_control_max']=med<=stats['maximum']
            within &= med<=stats['maximum']
        value['within_all_observed_control_maxima']=bool(within)
        matched.append(value);queries[b,q].append(value)
    query_rows=[]
    for (b,q),rs in sorted(queries.items()):
        query_rows.append(dict(band=b,query=q,matched_candidate_rotations=len(rs),
            matched_target_classes=len({r['target_group'] for r in rs}),
            any_within_observed_envelope=any(r['within_all_observed_control_maxima'] for r in rs),
            all_within_observed_envelope=all(r['within_all_observed_control_maxima'] for r in rs),
            **{m+'_median_delta':float(np.median([r[m+'_delta'] for r in rs])) for m in METRICS}))
    summary={}
    for b,rs in sorted(by_band.items()):
        qs=[r for r in query_rows if r['band']==b]
        summary[b]=dict(pairs=len(rs),queries=len({r['query'] for r in rs}),
            target_classes=len({r['target_group'] for r in rs}),matched_queries=len(qs),
            pooled={m:describe([r[m] for r in rs]) for m in METRICS},
            equal_query_matched_deltas={m:describe([r[m+'_median_delta'] for r in qs]) for m in METRICS},
            queries_with_any_observed_envelope_candidate=sum(r['any_within_observed_envelope'] for r in qs),
            queries_with_all_candidates_in_envelope=sum(r['all_within_observed_envelope'] for r in qs))
    target_sets=[{r['target_group'] for r in rs} for rs in by_band.values()]
    common=set.intersection(*target_sets) if target_sets else set()
    common_rows=[]
    for g in sorted(common):
        for b in sorted(by_band):
            values=[r for r in matched if r['target_group']==g and r['band']==b]
            if values:common_rows.append(dict(target_group=g,band=b,queries=len({r['query'] for r in values}),
                **{m+'_median_delta':float(np.median([r[m+'_delta'] for r in values])) for m in METRICS}))
    # Paired query-band comparisons where an identical target has observations in both bands.
    pairwise_bands=[]
    for a,b in __import__('itertools').combinations(sorted(by_band),2):
        ag={r['target_group'] for r in matched if r['band']==a};bg={r['target_group'] for r in matched if r['band']==b}
        for m in METRICS:
            differences=[]
            for g in sorted(ag&bg):
                av=np.median([r[m] for r in matched if r['band']==a and r['target_group']==g])
                bv=np.median([r[m] for r in matched if r['band']==b and r['target_group']==g])
                differences.append(float(av-bv))
            pairwise_bands.append(dict(first=a,second=b,metric=m,**describe(differences)))
    result=dict(status='complete',input_sha256=digest,rows=len(rows),
        unique_conformers=len({r[k] for r in rows for k in ('query','reference')}),
        same_molecule_pairs=sum(r['same_molecule']=='True' for r in rows),
        control_pairs=sum(map(len,controls.values())),control_classes=len(controls),
        control_pair_count_histogram=dict(Counter(map(len,controls.values()))),
        matched_candidate_rotations=len(matched),candidate_rotations_without_controls=len(missing),
        bands=summary,common_targets_all_bands=len(common),
        common_target_deltas={b:{m:describe([r[m+'_median_delta'] for r in common_rows if r['band']==b]) for m in METRICS} for b in sorted(by_band)},
        pairwise_band_target_matched=pairwise_bands,
        limitations=['Descriptive class-balanced panel, not population coverage or recall',
            'Control envelopes with at most three pairs are unstable diagnostics, not acceptance thresholds',
            'Queries reuse target controls; paired differences are not independent population replicates',
            'CSV hash verified; raw profiles and coordinates were not supplied for independent recomputation',
            'No reassignment or selection of a best rotation'])
    output.mkdir(parents=True,exist_ok=False)
    table(output/'matched_candidates.csv',matched);table(output/'query_balanced.csv',query_rows)
    table(output/'missing_controls.csv',missing);table(output/'common_targets.csv',common_rows)
    table(output/'band_target_comparisons.csv',pairwise_bands)
    result['output_hashes']={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in output.glob('*.csv')}
    (output/'report.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',required=True);p.add_argument('--output',required=True)
    p.add_argument('--expected-hash');a=p.parse_args();run(a.source,a.output,a.expected_hash)
