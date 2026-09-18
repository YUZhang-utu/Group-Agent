import json
from pathlib import Path
import platform
import os
import shutil
import sys
from types import SimpleNamespace

import numpy as np
import pytest

from aidd_agent.chemical_geometry import transform_directions
from aidd_agent.gaussian_overlay import apply_transform
from aidd_agent.gaussian_batch import _load_query, write_gaussian_query, ArtifactCatalogReader
from aidd_agent.chemical_companion import CHEM_META_DTYPE, ATOM_DTYPE, BOND_DTYPE, ChemicalCompanionReader
from aidd_agent import e035_validation as validation
from aidd_agent import library_acceptance as ev
from aidd_agent.interaction_matching import interaction_match, score_interaction_matches
from aidd_agent.interaction_fast import match_batch, InteractionFeatureReader
from aidd_agent.interaction_review import molecule_comparison
from aidd_agent.e035_validation import compare_sidecars, compare_score_arrays, stratified_ids, stress_chunk
from aidd_agent.expanded_wee1 import run_pose_stages
from test_expanded_wee1 import real_pose_fixture


def random_case(seed, batch=20, features=12, anchors=5):
    rng = np.random.default_rng(seed)
    directions = rng.normal(size=(batch, features, 3))
    directions /= np.linalg.norm(directions, axis=2, keepdims=True)
    qd = rng.normal(size=(anchors, 3)); qd /= np.linalg.norm(qd, axis=1, keepdims=True)
    q = (rng.normal(size=(anchors, 3)), rng.integers(1, 4, anchors), qd,
         rng.integers(0, 3, anchors, dtype=np.uint8), rng.uniform(.5, 2., anchors))
    cp = rng.normal(size=(batch, features, 3)); ct = rng.integers(1, 4, (batch, features))
    ck = rng.integers(0, 3, (batch, features), dtype=np.uint8)
    m = np.tile(np.eye(4), (batch, 1, 1))
    for matrix in m:
        u, _, v = np.linalg.svd(rng.normal(size=(3, 3))); r = u @ v
        if np.linalg.det(r) < 0: r[:, 0] *= -1
        matrix[:3, :3] = r; matrix[:3, 3] = rng.normal(size=3)
    return q, cp, ct, directions, ck, m


@pytest.mark.parametrize("features", [0, 1, 12, 32])
@pytest.mark.parametrize("anchors", [3, 5])
def test_batched_scalar_equivalence_signed_axial_and_empty(features, anchors):
    q, cp, ct, cd, ck, m = random_case(19, features=features, anchors=anchors)
    scores, assignments, anchors_out = match_batch(q, cp, ct, cd, ck, m)
    for i in range(len(cp)):
        r = interaction_match(*q, apply_transform(cp[i], m[i]), ct[i], transform_directions(cd[i], m[i]), ck[i])
        assert scores[i] == pytest.approx(r["interaction_match_score"], abs=1e-12)
        np.testing.assert_array_equal(assignments[i], r["assignments"])
        np.testing.assert_allclose(anchors_out[i], r["anchor_scores"], atol=1e-12, rtol=0)


def test_exact_tie_and_cutoff_assignments_stay_identical():
    q = (np.zeros((3, 3)), np.ones(3), np.zeros((3, 3)), np.zeros(3, dtype=np.uint8), np.ones(3))
    cp = np.zeros((2, 5, 3)); cp[1, :, 0] = 4.5
    cd = np.zeros_like(cp); types = np.ones((2, 5)); kinds = np.zeros((2, 5), dtype=np.uint8)
    matrices = np.tile(np.eye(4), (2, 1, 1))
    scores, assignments, _ = match_batch(q, cp, types, cd, kinds, matrices)
    for i in range(2):
        r = interaction_match(*q, cp[i], types[i], cd[i], kinds[i])
        np.testing.assert_array_equal(assignments[i], r["assignments"])
        assert scores[i] == r["interaction_match_score"]


@pytest.mark.parametrize("bad", ["nan", "homogeneous", "kind"])
def test_fast_validation_not_bypassed(bad):
    q, cp, ct, cd, ck, m = random_case(1)
    if bad == "nan": cp[0, 0, 0] = np.nan
    if bad == "homogeneous": m[0, 3, 0] = .1
    if bad == "kind": ck[0, 0] = 9
    with pytest.raises(ValueError): match_batch(q, cp, ct, cd, ck, m)


def test_real_readers_and_full_batched_sidecar_match_reference(tmp_path):
    root, batch, query, candidates = real_pose_fixture(tmp_path)
    run_pose_stages(root, batch, query, candidates, workers=1)
    optimized = root / "optimized.npz"
    score_interaction_matches(batch / "artifacts/catalog.json", batch / "chemical/catalog.json",
                              query / "gaussian-query.npz", root / "gaussian/refine/merged-scores.npz", optimized, engine="batched")
    assert compare_sidecars(root / "interaction-matches.npz", optimized)["passed"]
    reader = InteractionFeatureReader(batch / "artifacts/catalog.json", batch / "chemical/catalog.json")
    try:
        assert reader.get(0).molecule_id == "MOL-0"
        with pytest.raises(IndexError): reader.get(-1)
        with pytest.raises(IndexError): reader.get(3)
        with pytest.raises(ValueError): reader.get(.5)
    finally:
        reader.close()


def test_molecule_ranks_use_separate_representatives_and_no_group_backfill():
    ids = np.asarray([1, 2, 3, 4, 5, 6]); mols = np.asarray(["a", "a", "b", "c", "d", "e"])
    rows, summary, selected = molecule_comparison(ids, mols, np.asarray([10., 1, 9, 8, 2, 3]),
                                                 np.asarray([0., 9, 1, 2, 10, 3]), top=2, outside=3)
    a = next(r for r in rows if r["molecule_id"] == "a")
    assert a["gaussian_global_id"] == 1 and a["interaction_global_id"] == 2
    assert a["interaction_at_gaussian_pose"] == 0 and a["gaussian_at_interaction_pose"] == 1
    assert summary["different_representatives"] == 1
    assert [r["molecule_id"] for r in selected["gaussian_high_interaction_low"]] == ["b"]
    assert [r["molecule_id"] for r in selected["interaction_high_gaussian_low"]] == ["d"]
    assert [r["molecule_id"] for r in selected["both_high"]] == ["a"]


def test_rank_change_fails_even_with_small_numeric_error():
    result = compare_score_arrays(np.asarray([0, 1]), np.asarray(["a", "b"]),
                                  np.asarray([.1, .1]), np.asarray([.1, .10000001]))
    assert result["max_absolute_error"] < 1e-6
    assert not result["passed"]


def test_stratified_sample_distinct_reproducible_all_shards():
    catalog = dict(shards=[dict(name="a", global_id_start=0, conformers=2),
                           dict(name="b", global_id_start=2, conformers=98)])
    ids, counts = stratified_ids(catalog, 25, 12)
    assert len(np.unique(ids)) == 25 and len(ids) == 25
    assert min(counts.values()) >= 1
    assert np.array_equal(ids, stratified_ids(catalog, 25, 12)[0])
    assert np.array_equal(stratified_ids(catalog, 100, 12)[0], np.arange(100))
    with pytest.raises(ValueError): stratified_ids(catalog, 101, 12)


def test_scale_chunk_real_features_preserves_assignments_and_ids(tmp_path):
    root, batch, query_path, _ = real_pose_fixture(tmp_path)
    _, query = _load_query(query_path / "gaussian-query.npz")
    a = query["anchor_feature_indices"]
    q = tuple(query[k][a] for k in ("feature_points", "feature_types", "feature_directions", "feature_direction_kinds", "anchored_weights"))
    reader = InteractionFeatureReader(batch / "artifacts/catalog.json", batch / "chemical/catalog.json")
    try:
        result = stress_chunk(reader, np.arange(3), [q, q], root / "stress.npz", reference_readers=(
            ArtifactCatalogReader(batch / "artifacts/catalog.json"), ChemicalCompanionReader(batch / "chemical/catalog.json")))
    finally:
        reader.close()
    assert result["comparisons"] == 18 and result["assignment_mismatches"] == 0
    with np.load(root / "stress.npz") as archive:
        assert archive["global_ids"].tolist() == [0, 1, 2]


def test_scale_rejects_fast_reader_geometry_bug(tmp_path, monkeypatch):
    root, batch, _, _ = real_pose_fixture(tmp_path)
    reader = InteractionFeatureReader(batch / "artifacts/catalog.json", batch / "chemical/catalog.json")
    get = reader.get
    def corrupt(gid):
        result = get(gid); result.feature_points += 1
        return result
    monkeypatch.setattr(reader, "get", corrupt)
    try:
        with pytest.raises(ValueError, match="feature reader mismatch"):
            stress_chunk(reader, np.arange(3), [], root / "stress.npz", reference_readers=(
                ArtifactCatalogReader(batch / "artifacts/catalog.json"), ChemicalCompanionReader(batch / "chemical/catalog.json")))
    finally:
        reader.close()


def test_complete_e035_review_scale_and_resume_preserve_sources(tmp_path, monkeypatch):
    root, batch, query, candidates = real_pose_fixture(tmp_path)
    shard = batch / "chemical/shard-multi"
    meta = np.fromfile(shard / "chem-meta.bin", dtype=CHEM_META_DTYPE)
    meta["atoms"] = 4; meta["bonds"] = 3; meta.tofile(shard / "chem-meta.bin")
    atoms = np.zeros(4, dtype=ATOM_DTYPE); atoms["atomic_number"] = 6; atoms["formal_charge"][0] = 1
    atoms.tofile(shard / "atoms.bin")
    bonds = np.zeros(3, dtype=BOND_DTYPE); bonds["begin"] = [0, 1, 2]; bonds["end"] = [1, 2, 3]; bonds["order"] = 1
    bonds.tofile(shard / "bonds.bin")
    run_pose_stages(root, batch, query, candidates, workers=1)
    e034 = tmp_path / "e034"; e034.mkdir(); report = dict(status="complete", poses={})
    artifact, chemical = batch / "artifacts/catalog.json", batch / "chemical/catalog.json"
    for label, qid in (("8bju", "8BJU:QT9:A:601"), ("1x8b", "1X8B:824:A:901")):
        directory = e034 / label; qdir = directory / "query"; qdir.mkdir(parents=True)
        _, qa = _load_query(query / "gaussian-query.npz")
        write_gaussian_query(qdir / "gaussian-query.npz", **qa, source=dict(query_id=qid))
        (qdir / (qid.split(":")[0] + ".cif")).write_text("data_test\n")
        rigid = directory / "gaussian/refine/merged-scores.npz"; rigid.parent.mkdir(parents=True)
        shutil.copyfile(root / "gaussian/refine/merged-scores.npz", rigid)
        shutil.copyfile((root / "gaussian/refine/merged-scores.npz").with_suffix(".manifest.json"), rigid.with_suffix(".manifest.json"))
        side = directory / "interaction-matches.npz"
        manifest = score_interaction_matches(artifact, chemical, qdir / "gaussian-query.npz", rigid, side)
        report["poses"][qid] = dict(interaction=dict(manifest=manifest, timing=dict(status="measured_same_run", refine_wall_seconds=1.)),
                                    gaussian=dict(config=dict(query_sha256=ev.sha(qdir / "gaussian-query.npz"),
                                                             artifact_catalog_sha256=ev.sha(artifact)), final_result_sha256=ev.sha(rigid)))
    ev.write(e034 / "report.json", report)
    ev.write(e034 / "EXECUTION_COMPLETE.json", dict(status="complete", report_sha256=ev.sha(e034 / "report.json")))
    ev.write(e034 / "protocol.json", dict(host=platform.node(), cpu_count=os.cpu_count()))
    (batch / "run.lock").touch()
    before = validation.fingerprint(p for p in e034.rglob("*") if p.is_file())
    monkeypatch.setitem(sys.modules, "fcntl", SimpleNamespace(LOCK_EX=1, LOCK_NB=2, flock=lambda *_: None))
    monkeypatch.setattr(validation, "ensure_file_descriptor_limit", lambda _: None)
    args = SimpleNamespace(batch=batch, e034=e034, output=tmp_path / "e035", resume=False, repeats=3, scale=3)
    validation.run(args)
    result = ev.read(args.output / "report.json")
    assert result["scale"]["comparisons"] == 18
    assert all(q["equivalence_passed"] for q in result["queries"])
    review_dir = args.output / "8bju/review"
    exports = ev.read(review_dir / "poses.json")
    assert len(exports) == 3
    source_reader = validation.InteractionFeatureReader(artifact, chemical)
    from aidd_agent.gaussian_batch import ArtifactCatalogReader
    artifact_reader = ArtifactCatalogReader(artifact)
    for pose in exports:
        text = (review_dir / pose["file"]).read_text()
        assert "M  CHG" in text and "  4  3" in text.splitlines()[3]
        observed = np.asarray([[float(line[:10]), float(line[10:20]), float(line[20:30])] for line in text.splitlines()[4:8]])
        expected = apply_transform(artifact_reader.get(pose["global_id"]).shape_points, pose["transform"])
        np.testing.assert_allclose(observed, expected, atol=5.1e-5, rtol=0)
    source_reader.close()
    args.resume = True
    def never(*_, **__): pytest.fail("completed E035 stage recomputed")
    monkeypatch.setattr(validation, "benchmark_query", never)
    monkeypatch.setattr(validation, "stress_chunk", never)
    validation.run(args)
    assert before == validation.fingerprint(p for p in e034.rglob("*") if p.is_file())
    (review_dir / exports[0]["file"]).write_text("tampered")
    with pytest.raises(ValueError, match="changed completed output"):
        validation.run(args)
