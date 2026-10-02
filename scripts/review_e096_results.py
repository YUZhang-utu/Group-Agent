"""Small-file E096 review. Reads reports/CSV/transform only, never large SQLite or MOL2."""
import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
import math
from pathlib import Path


def read(path):return json.loads(Path(path).read_text(encoding='utf-8-sig'))
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def quantiles(values):
    values=sorted(values)
    def q(p):
        i=(len(values)-1)*p;lo=int(i);hi=math.ceil(i)
        return values[lo]*(hi-i)+values[hi]*(i-lo) if hi!=lo else values[lo]
    return {k:q(p) for k,p in [('p10',.1),('p50',.5),('p90',.9),('p95',.95),('p99',.99)]} if values else {}


def feature_self_test():
    """Check the installed RDKit factory AND residue-fragment mapping on synthetic controls."""
    import sys
    sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
    import numpy as np
    from rdkit import Chem, RDConfig, rdBase
    from aidd_agent.macrocycle_descriptors import describe_peptide
    from aidd_agent.joint_spatial_descriptor import geometry
    from aidd_agent import joint_spatial_descriptor
    cases=[]
    for name,side,index in [('neutral_amine','CCCCN',5),('protonated_amine','CCCC[NH3+]',5),('carboxylic_acid','CC(=O)O',6),('carboxylate','CC(=O)[O-]',6)]:
        mol=Chem.MolFromSmiles(f'N1[C@@H]({side})C(=O)N[C@@H](C)C(=O)N[C@@H](C)C1=O')
        conf=Chem.Conformer(mol.GetNumAtoms())
        for i,xyz in enumerate(np.random.default_rng(81).normal(size=(mol.GetNumAtoms(),3))):conf.SetAtomPosition(i,xyz)
        mol.AddConformer(conf);d=describe_peptide(mol,None,'backbone')
        v=geometry(mol,d['units'],np.asarray(d['descriptor']).reshape(-1,6))
        counts=[float(v[24+lag*56+index*7]*4) for lag in range(4)]
        cases.append(dict(name=name,counts=counts,passed=all(c>0 for c in counts)))
    return dict(rdkit=rdBase.rdkitVersion,feature_definition_sha256=sha(Path(RDConfig.RDDataDir)/'BaseFeatures.fdef'),
        descriptor_implementation_sha256=sha(joint_spatial_descriptor.__file__),cases=cases,passed=all(c['passed'] for c in cases),
        scope='Synthetic feature plumbing controls; not a census of the workstation library')


def review(root,baseline=None):
    root=Path(root);b=root/'blocks';r=read(b/'report.json');p=read(root/'profiles/report.json')
    t=read(b/'transform.json');decisions=read(b/'split-decisions.json')
    errors=[];warnings=[];hashes={}
    for name in ('blocks.csv','transform.json','split-decisions.json'):
        digest=sha(b/name);hashes[name]=digest
        if digest!=r.get('output_hashes',{}).get(name):errors.append('Hash mismatch: '+name)
    for source in (root/'profiles/report.json',):
        expected=r.get('sources',{}).get(str(source.resolve()))
        if expected is None:warnings.append('Profile report path not found in block source map')
        elif sha(source)!=expected:errors.append('Profile report source hash mismatch')
    with (b/'blocks.csv').open(encoding='utf-8-sig',newline='') as f:rows=list(csv.DictReader(f))
    groups=defaultdict(list);sizes=[];ids=set()
    for row in rows:
        n=int(row['conformers']);sizes.append(n);groups[row['parent_block']].append(n)
        if n<=0:errors.append('Nonpositive block population')
        if row['block_id'] in ids:errors.append('Duplicate block ID')
        ids.add(row['block_id'])
    checks={
        'profiles_complete':p.get('status')=='complete' and p.get('coverage')==1 and p.get('width')==317,
        'blocks_complete':r.get('status')=='complete',
        'reported_membership_gate':r.get('membership_gate')=='passed_against_complete_cached_profiles',
        'population_preserved':sum(sizes)==r.get('regular_conformers')==p.get('total_regular')==24530851,
        'special_pool_unchanged':r.get('special_conformers')==132885,
        'block_count_matches':len(rows)==r.get('regular_work_blocks'),
        'parent_count_matches':len(groups)==r.get('parent_regular_blocks')==384,
        'bounded_fragmentation':all(n>=5000 for n in sizes) and len(rows)<=1536,
        'child_policy':all(len(ns)<=4 and (len(ns)==1 or min(ns)>=max(5000,math.ceil(.2*sum(ns)))) for ns in groups.values()),
        'molecule_grouped_check':r.get('check_grouping')=='source_molecule_id_crc32_mod5',
        'reported_minimum_matches':bool(sizes) and min(sizes)==r.get('smallest_regular_block'),
        'reported_maximum_matches':bool(sizes) and max(sizes)==r.get('largest_regular_block'),
        'reported_no_small_blocks':r.get('blocks_below_minimum')==sum(n<5000 for n in sizes),
    }
    errors.extend(k for k,v in checks.items() if not v)
    chem=[i+j for i in (0,23,46) for j in range(11)]
    steric=[i+j for i in (0,23,46) for j in range(11,23)]
    groups_idx={'broad_chemistry':chem,'local_sterics':steric,'backbone_geometry':list(range(69,93)),'typed_spatial':list(range(93,317))}
    std=t['std'];weights=t['weights'];mean=t['mean']
    if any(len(a)!=317 or not all(math.isfinite(v) for v in a) for a in (std,weights,mean)):
        raise ValueError('Invalid 317-channel transform')
    if any(v<0 for v in std+weights):errors.append('Negative scale or weight')
    contributions={k:sum((std[i]*weights[i])**2 for i in indices) for k,indices in groups_idx.items()}
    denominator=sum(contributions.values())
    channels={k:dict(active=sum(weights[i]>0 for i in indices),total=len(indices),weighted_variance=contributions[k],
                    share_of_total_variance=contributions[k]/denominator if denominator else 0) for k,indices in groups_idx.items()}
    for k,v in channels.items():
        if v['active']!=r.get('active_channels_by_group',{}).get(k):errors.append('Active channel count mismatch: '+k)
        if not v['active']:warnings.append('Inactive entire feature group: '+k)
    families=['Donor','Acceptor','Aromatic','Hydrophobe','LumpedHydrophobe','PosIonizable','NegIonizable','SideHeavyAtoms']
    family_review={}
    for i,name in enumerate(families):
        indices=[93+lag*56+i*7+j for lag in range(4) for j in range(7)]
        count=mean[93+i*7]*(10 if i==7 else 4)
        family_review[name]=dict(mean_count=count,active_spatial_channels=sum(weights[j]>0 for j in indices),total_spatial_channels=28)
        if not math.isclose(count,r.get('mean_residue_feature_counts',{}).get(name,float('nan')),rel_tol=1e-7,abs_tol=1e-9):errors.append('Feature count report mismatch: '+name)
        if count==0:warnings.append(name+': zero factory-recognized residue-contained sites; raw chemistry and positive controls still required')
    reasons=Counter();shares=defaultdict(list);gain=[]
    for parent in decisions:
        for proposal in parent['proposals']:
            reasons[proposal['reason']]+=1
            if proposal.get('accepted'):
                axis=proposal['axis']
                if len(axis)!=317 or not all(math.isfinite(a) for a in axis):raise ValueError('Invalid split axis')
                norm=sum(a*a for a in axis)
                for k,indices in groups_idx.items():shares[k].append(sum(axis[i]**2 for i in indices)/norm if norm else 0)
                gain.append(proposal['check_gain'])
    if dict(reasons)!=r.get('decision_counts'):errors.append('Decision counts mismatch')
    if reasons['accepted_internal_check']!=len(rows)-len(groups):errors.append('Accepted split/tree leaf count mismatch')
    comparison=None
    if baseline:
        prior=read(Path(baseline)/'report.json')
        comparison={k:dict(E095=prior.get(k),E096=r.get(k)) for k in ['regular_work_blocks','smallest_regular_block','largest_regular_block','regular_conformers','special_conformers']}
    total=sum(sizes)
    return dict(status='review_failed' if errors else 'metadata_checks_passed_review_pending',scope='Small-file consistency diagnostic, not independent full membership validation',
        checks=checks,errors=errors,warnings=warnings,regular_blocks=len(rows),regular_conformers=total,
        size_quantiles=quantiles(sizes),minimum=min(sizes,default=0),maximum=max(sizes,default=0),
        largest_ten_fraction=sum(sorted(sizes,reverse=True)[:10])/total if total else None,
        blocks_below_5000=sum(n<5000 for n in sizes),children_per_parent=dict(Counter(len(ns) for ns in groups.values())),
        feature_groups=channels,feature_families=family_review,decision_counts=dict(reasons),check_gain_quantiles=quantiles(gain),
        accepted_axis_squared_loading_share={k:quantiles(v) for k,v in shares.items()},comparison=comparison,
        timings=dict(profiles_this_invocation=p.get('wall_seconds_this_invocation'),blocks_this_invocation=r.get('wall_seconds')),
        verified_small_file_hashes=hashes,limitations=['Large database hashes and membership rows were not reread; existing successful membership receipt is trusted.',
        'Weighted variance uses stored global statistics; axis loading shares are not causal feature importance or retrieval gain.',
        'E095/E096 summary comparison does not compare individual memberships.',
        'No raw-molecule ionizable-site census, docking, binding affinity, or held-out retrieval recall is established.'])


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--baseline',type=Path);p.add_argument('--output',type=Path,required=True);p.add_argument('--feature-self-test',action='store_true');a=p.parse_args()
    if a.output.exists():raise SystemExit('Choose a fresh diagnostic output file')
    if a.root.resolve() in a.output.resolve().parents:raise SystemExit('Place diagnostic output outside the frozen E096 directory')
    result=review(a.root,a.baseline)
    if a.feature_self_test:
        result['feature_self_test']=feature_self_test()
        profile=read(a.root/'profiles/report.json')
        result['feature_self_test']['matches_run_feature_definition']=result['feature_self_test']['feature_definition_sha256']==profile['signature']['feature_definition_sha256']
        result['feature_self_test']['matches_run_descriptor_implementation']=result['feature_self_test']['descriptor_implementation_sha256']==profile['signature']['code_hashes']['joint_spatial_descriptor.py']
        if not result['feature_self_test']['passed']:result['errors'].append('Synthetic ionizable feature self-test failed')
        if not result['feature_self_test']['matches_run_feature_definition'] or not result['feature_self_test']['matches_run_descriptor_implementation']:result['warnings'].append('Self-test environment differs from the original run; do not use it as a matching-environment check')
    if result['errors']:result['status']='review_failed'
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,indent=2),encoding='utf-8');print(json.dumps(result,indent=2))
    if result['errors']:raise SystemExit(1)
