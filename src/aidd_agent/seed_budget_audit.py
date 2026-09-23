"""Replay a frozen molecular panel from an existing retrieval in fresh outputs."""
import argparse
from collections import Counter
import os
from pathlib import Path
import sqlite3

import numpy as np

from . import library_acceptance as ev
from .budget_screen import collect_conformers, run_refinement, merge_and_rank
from .expanded_wee1 import fingerprint
from .gaussian_batch import _atomic_json
from .screening_selection import check_hashes


def variants():
    rows=[('reference-512',dict(max_pair_seeds=512,backend='reference'))]
    rows += [(f'generated-{cap}',dict(max_pair_seeds=cap,backend='batched')) for cap in (128,256,512,1024,2048)]
    rows += [(f'survivors-{target}',dict(max_pair_seeds=4096,survivor_target=target,backend='batched')) for target in (100,200)]
    return rows


def ranking(path):
    with sqlite3.connect(f'{path.resolve().as_uri()}?mode=ro',uri=True) as db:
        ordered=[row[0] for row in db.execute('SELECT mid FROM ranking ORDER BY rank')]
        scores={(mid,template):(contact,gaussian) for mid,template,contact,gaussian in
                db.execute('SELECT mid,template,contact,gaussian FROM poses')}
        payloads={(mid,template):payload for mid,template,payload in db.execute('SELECT mid,template,payload FROM poses')}
    return ordered,scores,payloads


def compare(reference, candidate, top):
    a,sa,pa=reference;b,sb,pb=candidate
    width=min(top,len(a));common=set(sa)&set(sb)
    return dict(reference_top_size=width,candidate_top_size=min(top,len(b)),
        top_overlap=len(set(a[:width])&set(b[:top])),
        reference_top_recall=len(set(a[:width])&set(b[:top]))/width if width else None,
        identical_molecule_order=a==b,exact_pose_payloads=pa==pb,
        missing_molecule_templates=len(set(sa)-set(sb)),new_molecule_templates=len(set(sb)-set(sa)),
        improved_contact_pairs=sum(sb[k][0]>sa[k][0] for k in common),
        decreased_contact_pairs=sum(sb[k][0]<sa[k][0] for k in common),
        maximum_absolute_contact_change=max((abs(sb[k][0]-sa[k][0]) for k in common),default=0.))


def summarize(folder):
    seconds=Counter();counts=Counter();stops=Counter();survivors=[];conformers=[]
    for receipt in sorted(folder.glob('chunks/*/*.receipt.json')):
        row=ev.read(receipt);check_hashes(row['files'])
        seconds.update(row['worker_seconds']);counts.update(row['counts'])
        for item in row['seed_search']['conformers']:
            survivors.append(item['surviving']);stops[item['stop_reason']]+=1
            conformers.append(dict(template_index=receipt.parent.name,**item))
    return dict(worker_seconds=dict(seconds),counts=dict(counts),stop_reasons=dict(stops),
        zero_survivor_conformer_templates=sum(x==0 for x in survivors),
        survivor_quantiles=dict(zip(('min','q25','median','q75','max'),
            np.quantile(survivors,[0,.25,.5,.75,1]).tolist())) if survivors else {},
        conformer_templates=conformers)


def run(args):
    source=args.run.resolve();out=args.output.resolve()
    protocol=ev.read(source/'protocol.json')
    check_hashes(protocol['sources'])
    artifact=[Path(p) for p in protocol['sources'] if Path(p).as_posix().endswith('/artifacts/catalog.json')]
    if len(artifact)!=1:raise ValueError('Expected one artifact catalog in source protocol')
    batch=artifact[0].parent.parent
    if any(out.is_relative_to(p) or p.is_relative_to(out) for p in (source,batch)):
        raise ValueError('Audit output must be separate from source run and library')
    if out.exists() and any(out.iterdir()):raise ValueError('Use an empty fresh audit output')
    retrieval=ev.read(source/'retrieval.json');check_hashes(retrieval['outputs'])
    design=ev.read(source/'adopted-design/report.json');check_hashes(design['sources'])
    keys=np.load(source/'selected-molecules.npy',allow_pickle=False)
    if not len(keys):raise ValueError('Source retrieval is empty')
    if args.molecules>len(keys):raise ValueError('Requested panel exceeds retrieved molecule count')
    rng=np.random.default_rng(args.seed)
    chosen=np.sort(rng.choice(keys,args.molecules,replace=False))
    ids=collect_conformers(batch,chosen)
    out.mkdir(parents=True,exist_ok=True)
    np.save(out/'sample-molecules.npy',chosen);np.save(out/'sample-conformers.npy',ids)
    inputs=[source/'protocol.json',source/'retrieval.json',source/'adopted-design/report.json',
            source/'selected-molecules.npy',out/'sample-molecules.npy',out/'sample-conformers.npy']
    frozen=dict(kind='seed_budget_audit',seed=args.seed,molecules=len(chosen),conformers=len(ids),
        workers=args.workers,chunk_conformers=args.chunk_conformers,top=args.top,
        assignment_backend=os.environ.get('AIDD_ASSIGNMENT_BACKEND','python'),
        variants=variants(),sources=fingerprint(inputs),code=fingerprint(sorted(Path(__file__).parent.glob('*.py'))),
        scope='Exploratory uniform sample of retrieved molecules, all stored conformers, all templates; not full-library recall',
        timing_scope='Includes worker initialization; fixed variant order and warm-cache effects prevent production speed claims')
    _atomic_json(out/'protocol.json',frozen)
    rows=[];rankings={}
    try:
        for name,settings in variants():
            _atomic_json(out/'status.json',dict(status='running',variant=name))
            folder=out/name;folder.mkdir()
            _atomic_json(folder/'protocol.json',dict(seed_search=settings,parent=fingerprint([out/'protocol.json'])))
            print(f'Seed audit: {name}',flush=True)
            timing=run_refinement(batch,folder,design,ids,args.workers,args.chunk_conformers,settings)
            count=merge_and_rank(folder,design,None,batch,protocol['template_quota'],protocol['rrf_k'])
            row=dict(name=name,settings=settings,timing=timing,ranked_molecules=count,**summarize(folder))
            _atomic_json(folder/'report.json',row)
            rankings[name]=ranking(folder/'ranking.sqlite')
            rows.append({k:v for k,v in row.items() if k!='conformer_templates'})
        for row in rows:
            row['vs_generated_2048']=compare(rankings['generated-2048'],rankings[row['name']],args.top)
        equivalent=compare(rankings['reference-512'],rankings['generated-512'],args.top)
        check_hashes(frozen['sources']);check_hashes(frozen['code'])
        check_hashes(protocol['sources']);check_hashes(design['sources'])
        report=dict(kind='seed_budget_audit',status='complete',variants=rows,backend_equivalence=equivalent,
            equivalence_passed=equivalent['identical_molecule_order'] and equivalent['exact_pose_payloads'],
            biological_validation='not_run',scope=frozen['scope'],
            outputs=fingerprint([out/'sample-molecules.npy',out/'sample-conformers.npy',
                *out.glob('*/report.json'),*out.glob('*/ranking.sqlite')]))
        _atomic_json(out/'report.json',report)
        _atomic_json(out/'status.json',dict(status='complete',equivalence_passed=report['equivalence_passed']))
        return report
    except BaseException as exc:
        _atomic_json(out/'status.json',dict(status='failed_or_interrupted',error=str(exc)));raise


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--molecules',type=int,default=64);p.add_argument('--workers',type=int,default=20)
    p.add_argument('--chunk-conformers',type=int,default=32);p.add_argument('--seed',type=int,default=20260924)
    p.add_argument('--top',type=int,default=20);args=p.parse_args()
    if min(args.molecules,args.workers,args.chunk_conformers,args.top)<=0:p.error('Positive budgets required')
    result=run(args)
    print(f"Audit complete; backend equivalence: {result['equivalence_passed']}",flush=True)
    if not result['equivalence_passed']:raise SystemExit(2)


if __name__=='__main__':main()
