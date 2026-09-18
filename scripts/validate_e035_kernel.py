"""Large synthetic equivalence stress test; never a real-library benchmark."""
import argparse
import json
from pathlib import Path
import time

import numpy as np

from aidd_agent.chemical_geometry import transform_directions
from aidd_agent.gaussian_overlay import apply_transform
from aidd_agent.interaction_matching import interaction_match
from aidd_agent.interaction_fast import match_batch


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cases", type=int, default=100000)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.cases < 4 or args.output.exists(): raise ValueError("Need >=4 cases and a new output")
    rng = np.random.default_rng(20260918); results = []; passed = True
    for group, (a, f) in enumerate(((3, 0), (3, 12), (5, 12), (5, 32))):
        count = args.cases // 4 + (group < args.cases % 4)
        qd = rng.normal(size=(a, 3)); qd /= np.linalg.norm(qd, axis=1, keepdims=True)
        q = (rng.normal(size=(a, 3)), rng.integers(1, 4, a), qd, rng.integers(0, 3, a, dtype=np.uint8), rng.uniform(.5, 2, a))
        times = dict(reference=0., batched=0.); max_error = 0.; mismatches = 0; ref_all = []; fast_all = []
        for start in range(0, count, 256):
            n = min(256, count-start)
            points = rng.normal(size=(n, f, 3)); types = rng.integers(1, 4, (n, f))
            directions = rng.normal(size=(n, f, 3)); directions /= np.linalg.norm(directions, axis=2, keepdims=True)
            kinds = rng.integers(0, 3, (n, f), dtype=np.uint8)
            matrices = np.tile(np.eye(4), (n, 1, 1)); theta = rng.uniform(-np.pi, np.pi, n)
            matrices[:, 0, 0] = matrices[:, 1, 1] = np.cos(theta)
            matrices[:, 1, 0] = np.sin(theta); matrices[:, 0, 1] = -np.sin(theta)
            matrices[:, :3, 3] = rng.normal(size=(n, 3))
            t = time.perf_counter()
            refs = [interaction_match(*q, apply_transform(points[i], matrices[i]), types[i],
                                      transform_directions(directions[i], matrices[i]), kinds[i]) for i in range(n)]
            times["reference"] += time.perf_counter() - t
            t = time.perf_counter()
            scores, assignments, anchors = match_batch(q, points, types, directions, kinds, matrices)
            times["batched"] += time.perf_counter() - t
            original = np.asarray([r["interaction_match_score"] for r in refs])
            original_anchors = np.asarray([r["anchor_scores"] for r in refs])
            max_error = max(max_error, float(np.max(np.abs(original-scores))), float(np.max(np.abs(original_anchors-anchors))))
            mismatches += int(np.sum(np.any(np.asarray([r["assignments"] for r in refs]) != assignments, axis=1)))
            ref_all.append(original.astype(np.float32)); fast_all.append(scores.astype(np.float32))
        ref, fast = np.concatenate(ref_all), np.concatenate(fast_all)
        ranks_equal = np.array_equal(np.argsort(-ref, kind="stable"), np.argsort(-fast, kind="stable"))
        ok = max_error <= 1e-6 and mismatches == 0 and ranks_equal; passed &= ok
        row = dict(anchors=a, features=f, cases=count, max_absolute_error=max_error,
                   assignment_mismatches=mismatches, float32_ranks_exact=bool(ranks_equal), passed=bool(ok), scoring_seconds=times)
        results.append(row); print(json.dumps(row), flush=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(dict(scope="synthetic kernel only; not real-library/pose/biology validation",
                                           cases=args.cases, passed=bool(passed), groups=results), indent=2), encoding="utf-8")
    if not passed: raise SystemExit(2)


if __name__ == "__main__": main()
