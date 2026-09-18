from types import SimpleNamespace
import sys
from pathlib import Path

import numpy as np
import pytest

from aidd_agent.gaussian_overlay import pair_alignment_seeds
from aidd_agent.gaussian_batch import _score_ids
from aidd_agent.workstation_suite import summarize


@pytest.mark.parametrize("seed", range(8))
def test_bounded_prefix_exact(seed):
    rng = np.random.default_rng(seed)
    candidate, query = rng.normal(size=(10, 3)), rng.normal(size=(4, 3))
    ct, qt = rng.integers(1, 3, 10), rng.integers(1, 3, 4)
    if seed == 0:
        candidate[:] = 0
    if seed == 1:
        candidate[1:] = candidate[0] + np.arange(9)[:, None]
    full = pair_alignment_seeds(candidate, ct, query, qt)
    for cap in (0, 1, 17, 512, len(full), len(full)+1):
        assert pair_alignment_seeds(candidate, ct, query, qt, max_seeds=cap) == full[:cap]
    with pytest.raises(ValueError, match="max_seeds"):
        pair_alignment_seeds(candidate, ct, query, qt, max_seeds=-1)


def test_gaussian_full_arrays_equal_with_bounded_seeds():
    rng = np.random.default_rng(41)
    query = dict(shape_points=rng.normal(size=(16, 3)), feature_points=rng.normal(size=(5, 3)),
                 feature_types=np.ones(5, dtype=int), anchored_weights=np.asarray([1, 2, 3, 1, 1.]),
                 anchor_feature_indices=np.arange(5))
    candidates = [SimpleNamespace(shape_points=rng.normal(size=(18, 3)), feature_points=rng.normal(size=(12, 3)),
                  feature_types=np.ones(12, dtype=int), molecule_id=str(i), conformer_id=str(i)) for i in range(2)]
    reader = SimpleNamespace(get=lambda i: candidates[i])
    kwargs = dict(sigma=1., cutoff=4.5, pair_tolerance=2., axial_samples=6, max_pair_seeds=17)
    reference = _score_ids(reader, query, np.arange(2), **kwargs)
    bounded = _score_ids(reader, query, np.arange(2), **kwargs, bounded_pair_seeds=True)
    assert reference["pair_seed_counts"].tolist() == [17, 17]
    for key in reference:
        np.testing.assert_array_equal(reference[key], bounded[key], err_msg=key)


def test_suite_excludes_profiles_and_partial_reuse():
    def record(name, seconds, eligible=True, profile=False):
        return dict(status="complete", scenario=name, profile=profile, report=dict(queries=[dict(
            query_id="q", execution_seconds=seconds, latency_eligible=eligible,
            gaussian_equivalence=dict(passed=True), annotation_equivalence=dict(passed=True))]))
    rows = summarize([record("old-reference", 4), record("old-reference", 6),
                      record("fine-bounded", 2), record("fine-bounded", 3),
                      record("fine-bounded", 100, profile=True),
                      record("old-bounded", 1, False), record("old-bounded", 1)])
    assert rows[0]["median_seconds"] == 5
    assert rows[1]["speedup_vs_old_reference"] == 2
    assert rows[2]["median_seconds"] is None


def test_suite_order_receipts_failure_and_resume(tmp_path, monkeypatch):
    from aidd_agent import workstation_suite as suite
    from aidd_agent import library_acceptance as ev
    batch, source = tmp_path / "batch", tmp_path / "source"
    for root, names in ((batch, ("artifacts/catalog.json", "chemical/catalog.json")),
                        (source, ("report.json", "protocol.json"))):
        for name in names:
            path = root / name; path.parent.mkdir(parents=True, exist_ok=True); ev.write(path, {})
    calls = []
    fail = [True]
    def fake_run(command, **kwargs):
        folder = Path(command[command.index("--output")+1]); folder.mkdir(exist_ok=True)
        calls.append((folder.name, "--resume" in command))
        if fail[0] and folder.name == "run-1-fine-reference":
            return SimpleNamespace(returncode=1)
        report = dict(status="complete", queries=[dict(query_id="q", execution_seconds=2,
                      latency_eligible="--profile-first-chunk" not in command,
                      gaussian_equivalence=dict(passed=True), annotation_equivalence=dict(passed=True))])
        ev.write(folder / "report.json", report)
        ev.write(folder / "RUN_STATUS.json", dict(status="complete", report_sha256=ev.sha(folder / "report.json")))
        return SimpleNamespace(returncode=0)
    monkeypatch.setattr(suite.subprocess, "run", fake_run)
    monkeypatch.setitem(sys.modules, "fcntl", SimpleNamespace(LOCK_EX=1, LOCK_NB=2, flock=lambda *_: None))
    args = SimpleNamespace(batch=batch, e034=source, output=tmp_path / "suite", workers=2, repeats=2,
                           profile=True, resume=False)
    assert suite.run(args) == 1
    assert len(calls) == 9
    assert calls[0][0] == "run-1-old-reference" and calls[4][0] == "run-2-old-bounded"
    assert calls[-1][0] == "profile-fine-bounded"
    assert ev.read(args.output / "RUN_STATUS.json")["status"] == "failed"
    fail[0] = False; args.resume = True; calls.clear()
    assert suite.run(args) == 0
    assert all(resume for _, resume in calls)
    assert len(ev.read(args.output / "report.json")["comparisons"]) == 4
    ev.write(source / "report.json", dict(changed=True))
    with pytest.raises(ValueError, match="Changed suite protocol"):
        suite.run(args)
