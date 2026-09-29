"""Re-evaluate a sealed v2 pocket report into a fresh v3 recommendation."""
import argparse
import json
from pathlib import Path

import numpy as np

from aidd_agent.expanded_wee1 import fingerprint
from aidd_agent.pocket_states import render
from aidd_agent.receptor_advice import attach, summary
from aidd_agent.screening_selection import check_hashes


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--maximum-representatives', type=int, default=32)
    parser.add_argument('--coverage-fraction', type=float, default=.95)
    args = parser.parse_args()
    source = args.source.resolve(); output = args.output.resolve()
    if output.exists():
        raise FileExistsError('Use a fresh recommendation directory')
    report = json.loads(source.read_text())
    if report.get('descriptor_version') != 2:
        raise ValueError('Use a completed v2 source report with corrected fields')
    check_hashes(report['sources'])
    old = Path(report['artifact'])
    if str(old) not in report['sources']:
        raise ValueError('Pocket grids are not sealed by the source report')
    with np.load(old) as archive:
        arrays = {k:archive[k] for k in archive.files}
    if not np.allclose(arrays['distance'], report['pairwise']['distance']):
        raise ValueError('Report and saved matrix disagree')
    attach(report, dict(maximum_representatives=args.maximum_representatives,
                        coverage_fraction=args.coverage_fraction))
    arrays['legacy_distance'] = arrays['distance']
    arrays['distance'] = np.asarray(report['pairwise']['distance'])
    output.mkdir(parents=True)
    artifact = output/'pocket-grids.npz'
    np.savez_compressed(artifact, **arrays)
    report['artifact'] = str(artifact)
    report['outputs'] = dict(grids=str(artifact),html=str(output/'report.html'),
                             recommendation=str(output/'receptor-recommendation.txt'))
    report['sources'].update(fingerprint([source,artifact]))
    report['advice_code_hashes'] = fingerprint([Path(__file__),
        Path(__file__).resolve().parents[1]/'src/aidd_agent/receptor_advice.py'])
    report['review_scope'] = 'Saved corrected pocket descriptors; no coordinate recomputation or new literature retrieval'
    (output/'report.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    text = summary(report['receptor_advice'])
    (output/'receptor-recommendation.txt').write_text(text,encoding='utf-8')
    render(report,output/'report.html',arrays['grid'],arrays['masks'])
    print(text)
    print('Review: '+str(output/'report.html'))


if __name__ == '__main__':
    main()
