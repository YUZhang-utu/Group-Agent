"""Alternate bounded reference/batched kernels on deterministic numeric panels."""
import argparse
import json
import platform
from pathlib import Path
import time

import numpy as np

from aidd_agent.gaussian_overlay import pair_alignment_seeds
from aidd_agent.expanded_wee1 import fingerprint


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    if args.output.exists():raise ValueError('Use a fresh benchmark output')
    rng=np.random.default_rng(20260924);rows=[]
    for case in range(12):
        n=(8,20,40)[case%3]
        inputs=(rng.normal(size=(n,3)),rng.integers(1,4,n),rng.normal(size=(12,3)),rng.integers(1,4,12))
        measurements={'reference':[],'batched':[]}
        for repeat in range(5):
            results={}
            for backend in (('reference','batched') if repeat%2==0 else ('batched','reference')):
                started=time.perf_counter()
                results[backend]=pair_alignment_seeds(*inputs,max_seeds=512,backend=backend)
                measurements[backend].append(time.perf_counter()-started)
            if results['reference']!=results['batched']:raise ValueError('Ordered payload equivalence failed')
        medians={key:float(np.median(values)) for key,values in measurements.items()}
        rows.append(dict(case=case,features=n,seeds=len(results['reference']),seconds=measurements,
            medians=medians,speedup=medians['reference']/medians['batched'],exact_payloads=True))
    report=dict(status='passed',cases=rows,python=platform.python_version(),numpy=np.__version__,
        platform=platform.platform(),sources=fingerprint([Path(__file__),Path('src/aidd_agent/gaussian_overlay.py')]),
        scope='Exploratory synthetic bounded seed kernel only; excludes collision, scoring, I/O and workstation timings')
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status='passed',median_case_speedup=float(np.median([r['speedup'] for r in rows])))))


if __name__=='__main__':main()
