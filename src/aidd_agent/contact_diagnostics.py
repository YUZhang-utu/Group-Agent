"""Source-scoped evidence relations and transparent matching diagnostics."""
from collections import defaultdict
import itertools
import re
import time

import numpy as np

from .interaction_matching import interaction_match, _maximum_weight_assignment


def extraction_metadata(interaction):
    direct=interaction in {'HBA','HBD'}
    return dict(extraction_method='direct_polar_geometry' if direct else 'extra_crystal_hypothesis',
        criteria_version='direct-hbond-atomcenter-v1' if direct else 'classified-extra-v1',
        criteria=dict(HBA='Complementary role; nearest protein partner <= 3.5 Angstrom',
            HBD='Complementary role; nearest protein partner <= 3.5 Angstrom',
            hydrophobic='Feature type 5; curated protein atom table; nearest per residue <= 4 Angstrom',
            salt_bridge='Opposite ionizable features; feature-center distance <= 5.5 Angstrom',
            pi_stacking='Ring centers <= 5.5 Angstrom; axial angle <= 30 or >= 60 degrees; offset <= 2 Angstrom',
            water_bridge='Two distances <= 3.5 Angstrom; unresolved water orientation').get(interaction,'See source evidence geometry'),
        unvalidated=['candidate_protein_contact','energetic_contribution','protonation']+
            (['explicit_hydrogen_angle'] if direct else []))


def source_index(observation, native_manifest):
    """Resolve only namespaced source IDs, never infer indices from proximity."""
    prefix=observation['query_id']+'/'
    name=observation['source_anchor']
    if not name.startswith(prefix):return None
    local=name[len(prefix):]
    found=re.match(r'^F(\d+):',local)
    index=int(found.group(1)) if found else native_manifest['source'].get('anchor_mapping',{}).get(local)
    return index if type(index) is int and index>=0 else None


def duplicate_relations(records, distance):
    if not np.isfinite(distance) or distance<=0:raise ValueError('Duplicate distance must be positive and finite')
    groups=defaultdict(list);unknown=[]
    for r in records:
        if not r.get('ligand_atom_indices'):
            unknown.append(dict(query_id=r['query_id'],source_anchor=r['source_anchor']))
        else:groups[r['query_id']].append(r)
    pairs=[]
    for qid,rows in groups.items():
        for a,b in itertools.combinations(rows,2):
            if a['source_anchor']==b['source_anchor']:continue
            common=set(a['ligand_atom_indices']) & set(b['ligand_atom_indices'])
            d=float(np.linalg.norm(np.array(a['point'])-b['point']))
            if not common or d>distance:continue
            pairs.append(dict(query_id=qid,left=a['anchor_id'],right=b['anchor_id'],
                left_source=a['source_anchor'],right_source=b['source_anchor'],shared_atom_indices=sorted(common),
                center_distance_angstrom=d,same_interaction_type=a['key'][2]==b['key'][2],
                same_protein_partner=a['key'][:2]==b['key'][:2],status='suspected_shared_evidence_requires_review'))
    return dict(distance_angstrom=distance,pairs=pairs,unresolved=unknown,
        policy='Same source instance only; pairwise relations, no transitive union or automatic scoring groups')


def optional_group_assignment(pair_scores, order, design):
    """Exact optional-only objective, one candidate per family. Diagnostic only.

    Hard rules, multiple hits within a family, and atom-level injectivity are not
    represented by this alternative model. No production assignment is replaced.
    """
    if design['mandatory_anchors'] or design['alternative_groups']:
        return dict(status='not_applicable_hard_rules')
    units=[dict(id=a,anchor_ids=[a],weight=w) for a,w in design['optional_weights'].items()]
    units+=design.get('optional_groups',[])
    if not units:return dict(status='no_optional_units')
    matrix=np.asarray(pair_scores);columns=matrix.shape[1]
    members=[[order.index(a) for a in g['anchor_ids']] for g in units]
    collapsed=np.array([matrix[indices].max(axis=0) for indices in members]).reshape(len(units),columns)
    assigned=_maximum_weight_assignment(collapsed,np.array([g['weight'] for g in units]))
    selections=[];total=0.
    for unit,indices,row,col in zip(units,members,collapsed,assigned):
        if col<0:continue
        winner=max(indices,key=lambda i:matrix[i,col])
        total+=unit['weight']*row[col]
        selections.append(dict(unit=unit['id'],anchor=order[winner],candidate_feature=int(col),score=float(row[col])))
    return dict(status='diagnostic_only',optional_numerator=float(total),selections=selections,
        scope='One candidate feature per optional unit; feature uniqueness is not atom uniqueness')


def query_tuple(q, indices=None):
    keys=('feature_points','feature_types','feature_directions','feature_direction_kinds')
    return [q[k] if indices is None else q[k][indices] for k in keys]


def match_details(query, candidate, labels):
    start=time.perf_counter()
    result=interaction_match(*query,np.ones(len(query[0])),*query_tuple(candidate))
    elapsed=time.perf_counter()-start
    rows=[]
    for i,label in enumerate(labels):
        j=int(result['assignments'][i]);kind=int(query[3][i]);spatial=angular=None
        if j>=0:
            d=float(np.linalg.norm(query[0][i]-candidate['feature_points'][j]))
            spatial=float(np.exp(-d*d/2))
            cosine=float(np.clip(np.dot(query[2][i],candidate['feature_directions'][j]),-1,1))
            angular=1. if kind==0 else (max(cosine,0.) if kind==1 else abs(cosine))**2
        rows.append(dict(label,feature_type=int(query[1][i]),direction_kind=kind,candidate_feature=j,
            matched=j>=0,spatial_factor=spatial,angular_factor=angular,score=float(result['anchor_scores'][i])))
    return rows,result,dict(seconds=elapsed,anchor_rows=len(query[0]),candidate_columns=len(candidate['feature_points']),
        assignment_columns=len(query[0])+len(candidate['feature_points']))


def distributions(rows):
    groups=defaultdict(list)
    for row in rows:groups[(row['layer'],row['feature_type'],row['direction_kind'])].append(row)
    result=[]
    for (layer,ft,kind),items in sorted(groups.items()):
        entry=dict(layer=layer,feature_type=ft,direction_kind=kind,count=len(items),unmatched=sum(not r['matched'] for r in items))
        from .conformer_artifacts import FEATURE_TYPES
        entry.update(feature_name=next((name for name,value in FEATURE_TYPES.items() if value==ft),'unknown'),
            direction_name={0:'none',1:'signed',2:'axial'}[kind])
        for name in ('score','spatial_factor','angular_factor'):
            values=[r[name] for r in items if r[name] is not None]
            entry[name]=dict(count=len(values),quantiles=np.quantile(values,[0,.25,.5,.75,1]).tolist()) if values else dict(count=0,quantiles=[])
        result.append(entry)
    return result


def contributions(values, assignments, order, design, anchors):
    scores={a:float(v) if j>=0 else 0. for a,v,j in zip(order,values,assignments)}
    denominator=design.get('optional_budget') if design.get('optional_normalization')=='fixed_budget' else (
        sum(design['optional_weights'].values())+sum(g['weight'] for g in design.get('optional_groups',[])))
    top=sum(design.get(k,0.) for k in ('gaussian_weight','optional_weight','occupancy_weight'))
    units=[dict(id=a,anchor_ids=[a],weight=w) for a,w in design['optional_weights'].items()]+design.get('optional_groups',[])
    rows=[]
    for unit in units:
        best=max(scores[a] for a in unit['anchor_ids'])
        winners=[a for a in unit['anchor_ids'] if scores[a]==best]
        for aid in winners:
            numerator=unit['weight']*best/len(winners)
            term=numerator/denominator if denominator else 0.
            rows.append(dict(unit=unit['id'],anchor_id=aid,interaction_type=anchors[aid]['feature_class'],
                numerator=numerator,contact_contribution=term,composite_contribution=term*design['optional_weight']/top,
                tie_credit=1/len(winners)))
    return rows
