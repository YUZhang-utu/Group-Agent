"""Evidence-bound spatial regions and optional alternatives, evaluated per pose."""
import json
from pathlib import Path
import re

import numpy as np
from scipy.spatial.distance import cdist

from .screening_selection import check_hashes

EXTENSION_FIELDS = {'optional_groups','spatial_groups','occupancy_rewards','occupancy_weight','spatial_ambiguity'}


def selected_anchors(design):
    return sorted(set(design['mandatory_anchors']) | set(design['optional_weights']) |
        {a for g in design['alternative_groups'] for a in g} |
        {a for g in design.get('optional_groups',[]) for a in g['anchor_ids']})


def number(value, low, high):
    return type(value) in (int,float) and np.isfinite(value) and low <= value <= high


def validate_extensions(design, survey):
    families=design.get('optional_groups',[]); spatial=design.get('spatial_groups',[])
    if not isinstance(families,list) or len(families)>20 or not isinstance(spatial,list) or len(spatial)>16:
        raise ValueError('Invalid grouped design lists')
    names=[];members=[];known={a['anchor_id'] for a in survey['anchors']}
    for group in families:
        if not isinstance(group,dict) or set(group)!={'id','anchor_ids','weight'}:
            raise ValueError('Optional family requires id, anchor_ids and weight')
        names.append(group['id']);ids=group['anchor_ids']
        if not isinstance(ids,list) or not ids or any(not isinstance(a,str) or a not in known for a in ids):
            raise ValueError('Optional family requires known anchor IDs')
        members.extend(ids)
        if not number(group['weight'],0,1):raise ValueError('Invalid optional family weight')
    occupied=set(design['mandatory_anchors']) | set(design['optional_weights']) | {a for g in design['alternative_groups'] for a in g}
    if len(set(members))!=len(members) or occupied.intersection(members):
        raise ValueError('Optional family members must not be double counted or required')
    for group in spatial:
        if not isinstance(group,dict) or set(group)!={'id','reference_query','ligand_atoms','radius','minimum_atoms'}:
            raise ValueError('Spatial group requires explicit reference atom selection and thresholds')
        names.append(group['id'])
        if not isinstance(group['reference_query'],str) or not isinstance(group['ligand_atoms'],list) or not group['ligand_atoms'] or any(not isinstance(a,str) or not a for a in group['ligand_atoms']):
            raise ValueError('Invalid spatial reference atom selection')
        if len(set(group['ligand_atoms']))!=len(group['ligand_atoms']):raise ValueError('Duplicate spatial reference atoms')
        if not number(group['radius'],.1,4) or type(group['minimum_atoms']) is not int or not 1<=group['minimum_atoms']<=100:
            raise ValueError('Invalid spatial geometry thresholds')
    if any(not isinstance(n,str) or not re.fullmatch(r'[A-Za-z][A-Za-z0-9_-]{0,63}',n) for n in names) or len(set(names))!=len(names):
        raise ValueError('Group identifiers must be unique short ASCII names')
    ambiguity=design.get('spatial_ambiguity',.5)
    if not number(ambiguity,0,2):raise ValueError('Invalid spatial ambiguity margin')
    rewards=design.get('occupancy_rewards',[]);weight=design.get('occupancy_weight',0.)
    if not number(weight,0,1):raise ValueError('Invalid occupancy weight')
    if not isinstance(rewards,list):raise ValueError('Occupancy rewards must be a list')
    if spatial:
        if len(rewards)!=len(spatial)+1 or any(not number(x,0,1) for x in rewards) or rewards[0]!=0 or any(a>b for a,b in zip(rewards,rewards[1:])):
            raise ValueError('Provide monotone occupancy rewards from zero for every group count')
        resolve_regions(design,survey)
    elif rewards or weight:
        raise ValueError('Occupancy rewards require explicit spatial groups')


def resolve_regions(design,survey):
    if not design.get('spatial_groups'):return []
    evidence=survey.get('contact_evidence')
    if not evidence:raise ValueError('Rebuild consensus contact evidence before selecting spatial regions')
    path=Path(evidence['path'])
    # Validate both the ledger itself and its upstream sources before using coordinates.
    if str(path.resolve()) not in survey['sources']:raise ValueError('Unsealed contact ledger')
    check_hashes({str(path.resolve()):survey['sources'][str(path.resolve())]})
    ledger=json.loads(path.read_text());check_hashes(ledger['sources'])
    if ledger['target']!=survey['target']['accession'] or ledger['coordinate_frame']!=survey['cohort']['reference']['coordinate_frame']:
        raise ValueError('Contact ledger target or coordinate frame mismatch')
    complexes={c['query_id']:c for c in ledger['complexes']}
    regions=[];used=set()
    for group in design['spatial_groups']:
        qid=group['reference_query']
        if qid not in complexes:raise ValueError('Unknown spatial reference complex')
        if complexes[qid].get('reference_state_check',{}).get('status')!='no_severe_overlap':
            raise ValueError('Spatial reference requires verified compatible receptor-state geometry; rebuild or review the contact ledger')
        atoms=complexes[qid]['ligand_atoms'];lookup={a['atom']:a for a in atoms}
        if len(lookup)!=len(atoms):raise ValueError('Ambiguous reference atom identity')
        if any(a not in lookup or lookup[a]['quality_issues'] for a in group['ligand_atoms']):
            raise ValueError('Missing or uncertain spatial reference atom')
        keys={(qid,a) for a in group['ligand_atoms']}
        if used & keys:raise ValueError('Reference atoms cannot define multiple spatial groups')
        used.update(keys)
        points=np.asarray([lookup[a]['point'] for a in group['ligand_atoms']],float)
        if points.ndim!=2 or points.shape[1]!=3 or not np.isfinite(points).all():raise ValueError('Invalid spatial reference coordinates')
        regions.append(dict(group,points=points.tolist(),coordinate_frame=ledger['coordinate_frame']))
    return regions


def grouped_terms(values,assignments,order,moved,design):
    scores={a:float(v) if i>=0 else 0. for a,v,i in zip(order,values,assignments)}
    family_scores={g['id']:max(scores[a] for a in g['anchor_ids']) for g in design.get('optional_groups',[])}
    regions=design.get('resolved_spatial_groups',[])
    if design.get('spatial_groups') and len(regions)!=len(design['spatial_groups']):
        raise ValueError('Spatial groups must be resolved during adoption')
    counts={g['id']:0 for g in regions};ambiguous=0
    if regions and len(moved):
        distances=np.column_stack([cdist(moved,np.asarray(g['points'])).min(axis=1) for g in regions])
        for row in distances:
            eligible=sorted((i for i,g in enumerate(regions) if row[i]<=g['radius']),key=lambda i:row[i])
            if not eligible:continue
            if len(eligible)>1 and row[eligible[1]]-row[eligible[0]]<=design.get('spatial_ambiguity',.5):
                ambiguous+=1;continue
            counts[regions[eligible[0]]['id']]+=1
    occupied=[g['id'] for g in regions if counts[g['id']]>=g['minimum_atoms']]
    reward=design['occupancy_rewards'][len(occupied)] if regions else 0.
    return dict(optional_group_scores=family_scores,spatial_atom_counts=counts,
                occupied_spatial_groups=occupied,occupied_group_count=len(occupied),
                ambiguous_spatial_atoms=ambiguous,occupancy_score=reward)
