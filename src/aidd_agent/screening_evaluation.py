"""Evaluate saved candidate rankings without repeating scientific calculations."""
import csv
import math
from pathlib import Path

from .budget_screen import readonly
from .expanded_wee1 import fingerprint
from .library_acceptance import read


def ranking(path):
    with readonly(path) as db:
        tables={r[0] for r in db.execute("SELECT name FROM sqlite_master WHERE type='table'")}
        table='members' if 'members' in tables else 'ranking'
        result=[];seen=set()
        for expected,(rank,mid) in enumerate(db.execute('SELECT rank,mid FROM '+table+' ORDER BY rank'),1):
            if rank!=expected or not isinstance(mid,str) or mid in seen:
                raise ValueError('Ranking must contain contiguous ranks and unique string molecule IDs')
            seen.add(mid);result.append(mid)
        if not result:raise ValueError('Empty ranking cannot establish retrieval quality')
        return result


def evaluate(reference,stages,labels=None,top_k=100,minimum_recovery=None):
    """Full reference defines the universe. Unlabelled molecules are never negatives."""
    if type(top_k) is not int or top_k<1:raise ValueError('top_k must be a positive integer')
    if len(reference)!=len(set(reference)) or not reference:raise ValueError('Invalid reference universe')
    universe=set(reference);labels=labels or {}
    if set(labels)-universe:raise ValueError('Label IDs absent from the reference universe')
    if any(type(v) is not int or v not in (0,1) for v in labels.values()):raise ValueError('Labels must be explicit binary measurements')
    if minimum_recovery is not None and (type(minimum_recovery) not in (float,int) or not math.isfinite(minimum_recovery) or not 0<=minimum_recovery<=1):
        raise ValueError('minimum_recovery must be in [0,1]')
    positives={mid for mid,value in labels.items() if value==1};prior=universe;reports=[]
    truth=set(reference[:top_k])
    for name,ordered in stages:
        selected=set(ordered)
        if len(selected)!=len(ordered) or selected-universe:raise ValueError('Stage contains duplicate or unknown molecule IDs')
        top=ordered[:top_k];top_set=set(top);labeled_top=[mid for mid in top if mid in labels]
        hits=len(top_set & positives);retained=selected & positives
        recovery=len(top_set & truth)/len(truth)
        complete_labels=len(labels)==len(universe)
        ef=(hits/len(top))/(len(positives)/len(universe)) if complete_labels and positives and top else None
        reports.append(dict(name=name,molecules=len(ordered),fraction_of_reference=len(ordered)/len(reference),
            nested_in_previous=selected<=prior,added_since_previous=len(selected-prior),
            retained_positives=len(retained),positive_recall=len(retained)/len(positives) if positives else None,
            lost_positive_ids=sorted((prior & positives)-selected),top_k_requested=top_k,top_k_actual=len(top),
            top_positive_hits=hits,top_labeled_count=len(labeled_top),
            precision_among_labeled=hits/len(labeled_top) if labeled_top else None,
            enrichment_factor=ef,exact_reference_top_k_recovery=recovery,
            recovery_gate='not_configured' if minimum_recovery is None else 'passed' if recovery>=minimum_recovery else 'failed'))
        prior=selected
    return dict(status='complete',kind='screening_benchmark',reference_molecules=len(reference),
        labeled_molecules=len(labels),positive_molecules=len(positives),stages=reports,
        scientific_acceptance='not_established',
        limitations=['Exact score-ranking recovery is not biological activity enrichment',
          'Enrichment factor requires labels for the complete reference universe; unlabelled is not inactive',
          'User-provided labels require assay and held-out split review before predictive claims',
          'Rankings must describe the same target/reference protocol; catalog identity alone is insufficient',
          'No runtime speedup is inferred from fewer retained molecules'])


def execute(args,directory):
    labels={};sources=[Path(args['reference_path'])]
    if args.get('labels_path'):
        path=Path(args['labels_path']);sources.append(path)
        with path.open(encoding='utf-8-sig',newline='') as stream:
            for row in csv.DictReader(stream):
                mid=row['molecule_id']
                if not mid or mid in labels:raise ValueError('Duplicate or empty label ID')
                if row['label'] not in {'0','1'}:raise ValueError('Labels must be 0 or 1')
                labels[mid]=int(row['label'])
    stages=[]
    for item in args['stages']:
        path=Path(item['path']);sources.append(path);stages.append((item['name'],ranking(path)))
    result=evaluate(ranking(args['reference_path']),stages,labels,args.get('top_k',100),args.get('minimum_recovery'))
    result['sources']=fingerprint(sources)
    result['protocol_note']=args.get('protocol_note','No target/assay comparability statement supplied')
    return result
