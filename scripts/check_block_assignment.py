"""Small deterministic equivalence and timing check; not a library benchmark."""
import argparse
import json
from pathlib import Path
import sys
import time

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'src'))
from aidd_agent.interaction_matching import _assignment_reference, _assignment_vectorized


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    rng = np.random.default_rng(106)
    cases = [(rng.integers(0, 4, (n,m))/3, rng.integers(0,4,n).astype(float))
             for n in (1,4,8,20) for m in (0,3,12,40) for _ in range(20)]
    for scores, weights in cases:
        np.testing.assert_array_equal(_assignment_reference(scores,weights),
                                      _assignment_vectorized(scores,weights))
    timings = []
    for n,m in ((8,12),(20,40),(40,60)):
        scores, weights = rng.random((n,m)), np.ones(n)
        measured = {}
        for name, function in [('python',_assignment_reference),('numpy',_assignment_vectorized)]:
            function(scores,weights)
            start = time.perf_counter()
            for _ in range(100):
                function(scores,weights)
            measured[name] = time.perf_counter()-start
        timings.append(dict(rows=n,columns=m,seconds=measured,
                            speed_ratio=measured['python']/measured['numpy']))
    result = dict(status='passed',exact_assignment_cases=len(cases),timings=timings,
                  scope='Synthetic solver microbenchmark; not workstation search acceleration or active recall')
    text = json.dumps(result,indent=2)
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(text,encoding='utf-8')
    print(text)


if __name__ == '__main__':
    main()
