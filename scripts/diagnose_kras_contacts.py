"""Compare strict KRAS gates with a diagnostic-only reference-match bypass."""
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
import json
from pathlib import Path
import time

import run_kras_weekend as run
import numpy as np

kw = run.kw


def diagnose(task):
    reader, feature_reader, states, chemical_reader = run.WORK
    output = {target: dict(counts=Counter(), matches=Counter(), molecules=[], examples=[])
              for target in states}
    for gid in task[1]:
        candidate = reader.get(gid)
        features = feature_reader.get(gid)
        chemistry = chemical_reader.get(gid)
        if (candidate.molecule_id, candidate.conformer_id) != (features.molecule_id, features.conformer_id):
            raise ValueError('Feature identity mismatch')
        if (candidate.molecule_id, candidate.conformer_id) != (chemistry.molecule_id, chemistry.conformer_id):
            raise ValueError('Chemical identity mismatch')
        sets = kw.chemical_sets(chemistry, candidate.feature_types)
        for target, state in states.items():
            row = output[target]
            strict = kw.evaluate(candidate, features, state, chemistry, sets)
            row['counts']['processed'] += 1
            row['counts']['strict_full'] += int(strict['flags']['full'])
            row['counts']['strict_contact'] += int(strict['flags']['contact'])
            n = strict['direct_region_pose_counts'].get('all', 0)
            row['counts']['direct_poses'] += n
            if not n:
                continue
            row['counts']['direct_conformers'] += 1
            original = kw.contact_mask

            def bypass(state, assignments, scores):
                indices = state['query']['anchor_feature_indices']
                lookup = {int(index): i for i, index in enumerate(indices)}
                required = [(95, 'HBA'), (69, 'HBD'), (62, 'HBD')]
                if target == 'KRAS_G12D':
                    required += [(12, 'HBD'), (68, 'HBA')]
                for assigned, score in zip(assignments, scores):
                    matched = []
                    for anchor in state['record']['anchors']:
                        index = lookup[anchor['feature_index']]
                        label = f"{anchor['target_residue']}_{anchor['feature_class']}_feature{anchor['feature_index']}"
                        if assigned[index] >= 0:
                            row['matches'][label + '_assigned'] += 1
                            if score[index] >= kw.POLICY['anchor']:
                                row['matches'][label + '_qualified'] += 1
                                matched.append(anchor)
                    for residue, kind in required:
                        ok = any(a['target_residue'] == residue and a['feature_class'] == kind for a in matched)
                        row['matches'][f'required_{residue}_{kind}_failed'] += int(not ok)
                    row['matches']['distinct_features_failed'] += int(len({a['feature_index'] for a in matched}) < kw.POLICY['distinct_features'])
                    row['matches']['distinct_residues_failed'] += int(len({a['target_residue'] for a in matched}) < kw.POLICY['residues'])
                return np.ones(len(assignments), dtype=bool)

            kw.contact_mask = bypass
            try:
                alternative = kw.evaluate(candidate, features, state, chemistry, sets)
            finally:
                kw.contact_mask = original
            for key in ('shape', 'color', 'geometry', 'full'):
                row['counts']['bypass_' + key] += int(alternative['flags'][key])
            if alternative['flags']['full']:
                row['molecules'].append(candidate.molecule_id)
            if alternative['best_contact_shape'] is not None:
                row['examples'].append(dict(gid=gid, molecule_id=candidate.molecule_id,
                    flags=alternative['flags'], best_contact_shape=alternative['best_contact_shape'],
                    passing_pose=alternative['passing_pose']))
                row['examples'] = sorted(row['examples'], key=lambda x: x['best_contact_shape']['shape'], reverse=True)[:10]
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pilot', type=Path, default=Path.home()/'kras-e116/pilot')
    parser.add_argument('--output', type=Path, default=Path.home()/'kras-e117-contact-diagnostic')
    parser.add_argument('--workers', type=int, default=20)
    args = parser.parse_args()
    protocol = kw.read(args.pilot/'protocol.json')
    if protocol['sample_size'] != 1000:
        parser.error('This protocol requires the original 1000-conformer pilot')
    if args.output.exists():
        parser.error('Use a new output directory; diagnostic output already exists')
    batch, prepared = Path(protocol['batch']), Path(protocol['prepared'])
    current = run.seal_inputs(batch, prepared, protocol['chunk_size'], protocol['sample_size'])
    if current != protocol:
        parser.error('Pilot input/code lineage mismatch')
    controls = {target: kw.self_control(kw.load_state(prepared, target)) for target in kw.REFERENCES}
    total = run.catalog_count(kw.read(batch/'artifacts/catalog.json'))
    args.output.mkdir(parents=True)
    kw.save(args.output/'protocol.json', dict(question='Reference matching gate ablation only',
        production_policy_unchanged=True, pilot_protocol=protocol,
        diagnostic_script_sha256=kw.sha(__file__)))
    merged = {target: dict(counts=Counter(), matches=Counter(), molecules=set(), examples=[]) for target in kw.REFERENCES}
    started = time.monotonic()
    with ProcessPoolExecutor(max_workers=args.workers, initializer=run.initialize, initargs=(batch, prepared)) as pool:
        for result in pool.map(diagnose, run.tasks(total, 32, 1000)):
            for target, row in result.items():
                dest = merged[target]
                dest['counts'].update(row['counts']); dest['matches'].update(row['matches'])
                dest['molecules'].update(row['molecules']); dest['examples'].extend(row['examples'])
    for row in merged.values():
        row['distinct_bypass_molecules'] = len(row.pop('molecules'))
        row['examples'] = sorted(row['examples'], key=lambda x: x['best_contact_shape']['shape'], reverse=True)[:20]
    report = dict(status='complete', scope='diagnostic_only', seconds=time.monotonic()-started,
        targets=merged, native_controls_passed=all(x['flags']['full'] for x in controls.values()),
        limitations=['Bypass hits are not production candidates.',
            'Failure counts overlap; feature scores combine positional and directional compatibility.',
            'This extra diagnostic work is not a production runtime estimate.'])
    kw.save(args.output/'report.json', report)
    print(json.dumps({**report, 'targets': {t: {k:v for k,v in r.items() if k!='examples'} for t,r in merged.items()}}, indent=2))


if __name__ == '__main__':
    main()
