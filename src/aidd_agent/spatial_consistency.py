"""Multiple explicit region definitions evaluated on one pose, without adoption."""
import argparse
from collections import Counter
import json
from pathlib import Path
import sqlite3

import numpy as np
from scipy.spatial.distance import cdist

from .expanded_wee1 import fingerprint
from .screening_selection import check_hashes


def occupancy(points, regions, margin):
    """Nearest-region assignment with ambiguity; no atom counts twice."""
    if regions is None:return None
    if not regions or len({g['id'] for g in regions})!=len(regions):
        raise ValueError('Nonempty unique region IDs required')
    counts={g['id']:0 for g in regions}
    for g in regions:
        p=np.asarray(g['points'],float)
        if p.ndim!=2 or p.shape[1]!=3 or not len(p) or not np.isfinite(p).all():
            raise ValueError('Invalid region coordinates')
        if not 0<g['radius']<=4 or type(g['minimum_atoms']) is not int or g['minimum_atoms']<1:
            raise ValueError('Invalid region thresholds')
    if not 0<=margin<=2:raise ValueError('Invalid ambiguity margin')
    points=np.asarray(points,float).reshape(-1,3)
    if not np.isfinite(points).all():raise ValueError('Nonfinite candidate points')
    if len(points):
        distances=np.column_stack([cdist(points,g['points']).min(axis=1) for g in regions])
        for row in distances:
            order=sorted((i for i,g in enumerate(regions) if row[i]<=g['radius']),key=lambda i:row[i])
            if not order:continue
            winner=order[0];g=regions[winner]
            if any(row[j]-row[winner]<=max(g.get('ambiguity_margin',margin),regions[j].get('ambiguity_margin',margin)) for j in order[1:]):continue
            counts[g['id']]+=1
    return {g['id']:counts[g['id']]>=g['minimum_atoms'] for g in regions}


def compare(points, definitions, margin=.5):
    if len(definitions)!=2:raise ValueError('Exactly two named definitions required')
    results={k:occupancy(points,v,margin) for k,v in definitions.items()}
    ids=set().union(*(set(v) for v in results.values() if v is not None))
    states={}
    for name in sorted(ids):
        values=[None if v is None else v.get(name) for v in results.values()]
        states[name]=('unknown' if None in values else 'agreement' if all(values)
                      else 'boundary' if any(values) else 'unoccupied')
    return dict(definitions=results,states=states,
        agreed_regions=sum(v=='agreement' for v in states.values()),
        boundary_regions=sum(v=='boundary' for v in states.values()))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--pilot',type=Path,required=True);p.add_argument('--batch',type=Path,required=True)
    p.add_argument('--definitions',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();report=json.loads((a.pilot/'report.json').read_text())
    check_hashes(report['output_hashes'])
    definition=json.loads(a.definitions.read_text())
    adopted=json.loads((a.pilot/'adopted-design/report.json').read_text())
    reference=adopted['reference']
    frame=reference.get('coordinate_frame',reference.get('query_id')) if isinstance(reference,dict) else reference
    if frame!=definition['coordinate_frame'] or adopted['target']['accession']!=definition['target']:
        raise ValueError('Region target or coordinate frame mismatch')
    from .gaussian_batch import ArtifactCatalogReader
    reader=ArtifactCatalogReader(a.batch/'artifacts/catalog.json')
    a.output.mkdir(parents=True,exist_ok=False);states=Counter();best={};n=0
    with sqlite3.connect((a.pilot/'candidates.sqlite').resolve().as_uri()+'?mode=ro',uri=True) as db, (a.output/'poses.jsonl').open('w') as stream:
        for mid,gid,payload in db.execute('SELECT mid,gid,payload FROM poses'):
            row=json.loads(payload);candidate=reader.get(gid)
            if candidate.molecule_id!=mid:raise ValueError('Pose molecule identity mismatch')
            matrix=np.asarray(row['transform']).reshape(4,4)
            moved=candidate.shape_points@matrix[:3,:3].T+matrix[:3,3]
            result=compare(moved,definition['definitions'],definition['ambiguity_margin'])
            states.update(f'{k}:{v}' for k,v in result['states'].items())
            best[mid]=max(best.get(mid,0),result['agreed_regions']);n+=1
            stream.write(json.dumps(dict(molecule_id=mid,global_id=gid,template_id=row['template_id'],**result))+'\n')
    result=dict(status='complete',production_changed=False,poses=n,region_states=dict(states),
        molecule_best_same_pose_agreement=dict(Counter(best.values())),
        limitations=['Only saved passing representatives are replayed; not an exhaustive region-aware search.',
            'Agreement is reference consistency, not independent confidence or activity.',
            'No new score coefficients or hard gate adopted.'],
        sources=fingerprint([a.definitions,a.pilot/'report.json',a.pilot/'candidates.sqlite']))
    (a.output/'report.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))


if __name__=='__main__':main()
