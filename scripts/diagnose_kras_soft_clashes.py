"""Bounded collision audit of the exact E118 pilot; no gate changes or hit export."""
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
import json
from pathlib import Path
import time
import numpy as np
import kras_soft_rank as soft

THRESHOLDS = (1.0, 1.2, 1.5, 1.8, 2.0)


def distance_summary(points, matrices, state):
    moved = np.einsum('bij,nj->bni', matrices[:, :3, :3], points) + matrices[:, None, :3, 3]
    distances, partners = state['tree'].query(moved.reshape(-1, 3))
    distances = distances.reshape(len(matrices), -1)
    partners = partners.reshape(len(matrices), -1)
    return moved, distances, partners


def inspect_chunk(task):
    _, ids = task
    reader, _, states, chemistry_reader = soft.run.WORK
    output = {t: dict(minima=[], clash_counts=[], worst_partners={}, best=[]) for t in states}
    for gid in ids:
        candidate = reader.get(gid)
        targets = [t for t in states if soft.COUNTS.execute(
            'SELECT 1 FROM molecules WHERE target=? AND mid=?', (t, candidate.molecule_id)).fetchone()]
        if not targets: continue
        chemistry = chemistry_reader.get(gid)
        if (candidate.molecule_id,candidate.conformer_id) != (chemistry.molecule_id,chemistry.conformer_id):
            raise ValueError('Chemical identity mismatch')
        sets = soft.run.kw.chemical_sets(chemistry, candidate.feature_types)
        if not all(len(v) for v in sets.values()): continue
        for target in targets:
            state = states[target]; row = output[target]
            seeds, _ = soft.run.kw.prepare_seeds(candidate, state['query'], **soft.run.kw.POSE_PARAMETERS)
            matrices = np.array([s.transform_matrix for s in seeds]).reshape(-1, 4, 4)
            matrices = matrices[soft.run.kw.direct_contact_masks(candidate.shape_points, sets, state, matrices)['all']]
            if not len(matrices): continue
            moved, distances, partners = distance_summary(candidate.shape_points, matrices, state)
            minima = distances.min(axis=1)
            row['minima'].extend(minima.tolist())
            row['clash_counts'].extend((distances < 2.).sum(axis=1).tolist())
            labels = {tuple(np.round(a['xyz'], 5)): f"{a['resname']}{a['residue']}:{a['atom']}" for a in state['record']['receptor_atoms']}
            def label(index):
                return labels.get(tuple(np.round(state['geometry']['receptor'][index],5)), 'nonprotein_or_unmapped_environment')
            worst = partners[np.arange(len(matrices)), distances.argmin(axis=1)]
            counts = Counter(label(i) for i in worst)
            for key, value in counts.items(): row['worst_partners'][key] = row['worst_partners'].get(key,0)+value
            j = int(np.argmax(minima))
            atom_ids = np.argsort(distances[j])[:10]
            row['best'].append(dict(mid=candidate.molecule_id,cid=candidate.conformer_id,gid=gid,
                minimum_distance=float(minima[j]), atoms_below_2=int((distances[j]<2.).sum()),
                heavy_atoms=len(candidate.shape_points),transform=matrices[j].ravel().tolist(),
                closest_pairs=[dict(candidate_heavy_index=int(i),element_atomic_number=int(chemistry.atomic_numbers[i]),
                    distance=float(distances[j,i]),receptor_atom=label(partners[j,i]),
                    candidate_xyz=moved[j,i].tolist(), receptor_xyz=state['geometry']['receptor'][partners[j,i]].tolist()) for i in atom_ids]))
            row['best'] = sorted(row['best'],key=lambda r:(-r['minimum_distance'],r['gid']))[:30]
    return output


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--pilot',type=Path,required=True)
    p.add_argument('--direct',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--workers',type=int,default=20)
    args=p.parse_args(); kw=soft.run.kw
    report=kw.read(args.pilot/'report.json'); protocol=kw.read(args.pilot/'protocol.json')
    sample=protocol['sample']; common=protocol['common']; original=common['original']
    if report.get('status')!='complete' or report.get('scope')!='spread_pilot' or not 0<sample<=10000:
        p.error('A complete pilot of at most 10000 conformers is required')
    if args.workers<1: p.error('Positive worker count required')
    if args.output.exists(): p.error('Use a new diagnostic output directory')
    batch=Path(original['batch']);prepared=Path(original['prepared']);counts=args.direct/'counts.sqlite'
    for protected in (args.pilot,args.direct,batch,prepared):
        a,b=args.output.resolve(),protected.resolve()
        if a.is_relative_to(b) or b.is_relative_to(a): p.error('Output overlaps protected input')
    if soft.code_hash()!=common['code_hash'] or kw.sha(counts)!=common['counts_sha256']:
        p.error('Pilot code or direct pool changed')
    if soft.run.seal_inputs(batch,prepared,original['chunk_size'],original['sample_size'])!=original:
        p.error('Original input lineage changed')
    total=soft.run.catalog_count(kw.read(batch/'artifacts/catalog.json'))
    started=time.monotonic()
    aggregate={t:dict(minima=[],clash_counts=[],worst_partners=Counter(),best=[]) for t in kw.REFERENCES}
    with ProcessPoolExecutor(max_workers=args.workers,initializer=soft.initialize,initargs=(batch,prepared,counts)) as pool:
        for result in pool.map(inspect_chunk,soft.run.tasks(total,protocol['chunk_size'],sample)):
            for target,row in result.items():
                out=aggregate[target]
                for key in ('minima','clash_counts','best'): out[key].extend(row[key])
                out['worst_partners'].update(row['worst_partners'])
    targets={}
    for target,row in aggregate.items():
        values=np.array(row['minima']); clashes=np.array(row['clash_counts'])
        best={}
        for entry in sorted(row['best'],key=lambda r:(-r['minimum_distance'],r['gid'])):
            best.setdefault(entry['mid'],entry)
        state=kw.load_state(prepared,target)
        native=state['tree'].query(state['query']['shape_points'])[0]
        targets[target]=dict(direct_poses=len(values),
            minimum_distance_quantiles=dict(zip(('min','p25','median','p75','p95','max'),np.quantile(values,[0,.25,.5,.75,.95,1]).tolist())) if len(values) else {},
            diagnostic_pose_counts_by_distance={str(t):int((values>=t).sum()) for t in THRESHOLDS},
            atoms_below_2_quantiles=np.quantile(clashes,[0,.5,.95,1]).tolist() if len(clashes) else [],
            at_most_one_atom_below_2=int((clashes<=1).sum()),
            worst_partner_counts=dict(row['worst_partners'].most_common(20)),
            native_minimum_distance=float(native.min()),best_boundary_examples=list(best.values())[:20],
            matches_pilot_direct_count=len(values)==report['targets'][target]['direct_poses'],
            matches_pilot_geometry_count=int((values>=2.).sum())==report['targets'][target]['geometry_poses'])
    result=dict(status='complete',scope='diagnostic_only_not_candidates',targets=targets,
                seconds=time.monotonic()-started,pilot_report_sha256=kw.sha(args.pilot/'report.json'),
                limitations=['Distance cutoffs are sensitivity diagnostics, not adopted acceptance rules.',
                             'Worst-pair counts are pose counts, not distinct molecules.',
                             'Examples are capped closest-to-passing poses, not an unbiased sample.',
                             'No geometry refinement, affinity scoring or molecule export performed.'])
    if not all(t['matches_pilot_direct_count'] and t['matches_pilot_geometry_count'] for t in targets.values()):
        result['status']='failed_reconciliation'
    args.output.mkdir(parents=True);kw.save(args.output/'report.json',result)
    print(json.dumps({**result,'targets':{t:{k:v for k,v in r.items() if k!='best_boundary_examples'} for t,r in targets.items()}},indent=2))
    if result['status']!='complete':raise SystemExit(2)


if __name__=='__main__':main()
