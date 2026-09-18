"""Synthetic seed-generation microbenchmark; not end-to-end search latency."""
import argparse
import hashlib
import json
import platform
import time
from pathlib import Path

import numpy as np

from aidd_agent.gaussian_overlay import pair_alignment_seeds


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rng = np.random.default_rng(20260918)
    rows = []
    for case in range(12):
        count = (8, 16, 32)[case % 3]
        candidate = rng.normal(size=(count, 3))*2
        query = rng.normal(size=(5, 3))*2
        ct, qt = rng.integers(1, 4, count), rng.integers(1, 4, 5)
        timings = {"reference": [], "bounded": []}
        for repeat in range(3):
            for engine in (("reference", "bounded") if repeat % 2 == 0 else ("bounded", "reference")):
                started = time.perf_counter()
                result = pair_alignment_seeds(candidate, ct, query, qt,
                                               max_seeds=512 if engine == "bounded" else None)
                timings[engine].append(time.perf_counter()-started)
                if engine == "reference": full = result
                else: bounded = result
            if full[:512] != bounded:
                raise ValueError("Seed prefix changed")
        medians = {key: float(np.median(values)) for key, values in timings.items()}
        rows.append(dict(case=case, features=count, full_seeds=len(full), retained_seeds=len(bounded),
                         prefix_exact=True, seconds=timings, medians=medians,
                         speedup=medians["reference"]/medians["bounded"]))
    source = Path(__file__).resolve().parents[1] / "src/aidd_agent/gaussian_overlay.py"
    result = dict(status="passed", cases=rows, numpy=np.__version__, python=platform.python_version(),
                  platform=platform.platform(), source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                  scope="synthetic generator microbenchmark only; full query speedup unmeasured")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(dict(status="passed", cases=len(rows), min_speedup=min(r["speedup"] for r in rows),
                         max_speedup=max(r["speedup"] for r in rows))))


if __name__ == "__main__":
    main()
