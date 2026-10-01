"""Read-only E094 size review and exact root-split replay; no chemistry recomputation."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import sqlite3

import numpy as np


def root_probe(x, minimum, cutoff):
    """Replay the v1 first proposal, including float32 arithmetic and tie policy."""
    x = np.asarray(x, dtype=np.float32)
    if x.ndim != 2 or x.shape[1] != 69 or not np.isfinite(x).all():
        raise ValueError('Expected finite 69-channel profiles')
    result = dict(population=len(x), active_channels=int((np.ptp(x, axis=0) > 1e-7).sum()))
    if len(x) < 2 * minimum:
        return dict(result, reason='population_floor')
    centered = x - x.mean(0)
    covariance = centered.T @ centered / len(x)
    eigenvalues, axes = np.linalg.eigh(covariance)
    axis = axes[:, -1]
    if axis[int(np.argmax(abs(axis)))] < 0:
        axis = -axis
    projection = centered @ axis
    ordered = np.argsort(projection, kind='stable')
    options = np.flatnonzero(np.diff(projection[ordered]) > 1e-7) + 1
    options = options[(options >= minimum) & (options <= len(x) - minimum)]
    result.update(first_pc_variance=float(eigenvalues[-1]), feasible_cuts=len(options))
    if not len(options):
        return dict(result, reason='no_distinct_projection_cut_meeting_floor')
    cut = int(options[np.argmin(abs(options - len(x) / 2))])
    delta = x[ordered[:cut]].mean(0) - x[ordered[cut:]].mean(0)
    contrast = float(np.sqrt(np.mean(delta ** 2)))
    chemical = np.concatenate([np.arange(i, i+11) for i in (0, 23, 46)])
    steric = np.concatenate([np.arange(i+11, i+23) for i in (0, 23, 46)])
    result.update(left_count=cut, right_count=len(x)-cut, centroid_contrast=contrast,
                  cutoff=cutoff, contrast_to_cutoff=contrast/cutoff,
                  chemistry_rms=float(np.sqrt(np.mean(delta[chemical] ** 2))),
                  steric_rms=float(np.sqrt(np.mean(delta[steric] ** 2))),
                  largest_channel_difference=float(abs(delta).max()),
                  reason='below_contrast_cutoff' if contrast < cutoff else 'root_split_would_pass')
    return result


def run(root, output, top):
    root, output = Path(root).resolve(), Path(output).resolve()
    if output == root or root in output.parents or output in root.parents:
        raise ValueError('Use a new diagnostic directory outside the completed build')
    paths = [root/'backbone/report.json', root/'properties/report.json',
             root/'profiles/report.json', root/'properties/blocks.csv']
    br, pr, fr = [json.loads(p.read_text()) for p in paths[:3]]
    if any(r.get('status') != 'complete' for r in (br, pr, fr)):
        raise ValueError('Completed builds required')
    rows = list(csv.DictReader(paths[3].open(newline='')))
    sizes = np.array([int(r['conformers']) for r in rows])
    if len(sizes) != pr['regular_work_blocks'] or int(sizes.sum()) != pr['regular_conformers']:
        raise ValueError('CSV/report population mismatch')
    result = dict(scope='Diagnostic replay only; no membership changes or recall claim',
        provenance={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        large_database_hashes_rechecked=False,
        regular_blocks=len(rows), regular_conformers=int(sizes.sum()),
        minimum=int(sizes.min()), maximum=int(sizes.max()), mean=float(sizes.mean()),
        size_quantiles=dict(zip(['p10','p50','p90','p95','p99'], map(float,np.quantile(sizes,[.1,.5,.9,.95,.99])))),
        largest_block_fraction=float(sizes.max()/sizes.sum()),
        largest_ten_fraction=float(np.sort(sizes)[-10:].sum()/sizes.sum()),
        below_minimum=int((sizes<br['minimum_regular_block_size']).sum()),
        split_rows=sum(r['split'].lower()=='true' for r in rows), probes=[],
        limitations=['Large databases rely on prior successful validation.',
          'Root replay uses all vectors of selected parents, not a sample; numerical near-ties can depend on NumPy/BLAS.',
          'Mean/std/max summaries cannot recover side-chain positions.',
          'Group RMS values describe the proposed split, not calibrated chemistry thresholds.'])
    output.mkdir(parents=True, exist_ok=False)
    (output/'report.json').write_text(json.dumps(dict(result,status='metadata_complete'),indent=2))
    with sqlite3.connect((root/'backbone/work_blocks.sqlite').as_uri()+'?mode=ro',uri=True) as db:
        parents=db.execute("SELECT block_id,n FROM block WHERE kind!='special_exhaustive' ORDER BY n DESC,block_id").fetchall()
    with sqlite3.connect((root/'profiles/profiles.sqlite').as_uri()+'?mode=ro',uri=True) as db:
        for parent,n in (parents if top==0 else parents[:top]):
            print(f'Replaying {parent}: {n:,} cached vectors',flush=True)
            matrix=np.empty((n,69),dtype=np.float32)
            count=0
            for count,(blob,) in enumerate(db.execute('SELECT vector FROM item WHERE parent=? ORDER BY cid',(parent,)),1):
                if count>n or blob is None or len(blob)!=69*4:
                    raise ValueError('Invalid cached profile coverage')
                matrix[count-1]=np.frombuffer(blob,dtype='<f4')
            if count!=n:raise ValueError('Cached profile count mismatch')
            probe=root_probe(matrix,br['minimum_regular_block_size'],pr['minimum_centroid_contrast'])
            probe['parent_block']=parent
            probe['maximum_children_policy']=pr['maximum_children_per_parent']
            result['probes'].append(probe)
            del matrix
            (output/'report.json').write_text(json.dumps(dict(result,status='running'),indent=2))
    result['status']='complete'
    (output/'report.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',required=True,type=Path)
    p.add_argument('--output',required=True,type=Path)
    p.add_argument('--top',type=int,default=10,help='Largest parent blocks to replay; 0 for all')
    a=p.parse_args()
    if a.top<0:p.error('--top must be nonnegative')
    run(a.root,a.output,a.top)
