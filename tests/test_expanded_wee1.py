import json
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np
import pytest

from aidd_agent import expanded_wee1 as e
from aidd_agent import library_acceptance as ev
from aidd_agent.chemical_companion import (
    ATOM_DTYPE, BOND_DTYPE, CHEM_META_DTYPE, FEATURE_DIRECTION_DTYPE, TORSION_DTYPE,
)
from aidd_agent.gaussian_batch import write_gaussian_query
from test_library_acceptance import library, ExactIndex
from test_gaussian_batch import _multi_artifact


def make_e033(batch, root):
    root.mkdir()
    _, acceptance = ev.accept_library(batch)
    params = dict(metadata_only=False)
    report = dict(acceptance=acceptance, parameters=params, calibration_status="gate_passed_on_panel")
    ev.write(root / "acceptance.json", acceptance)
    ev.write(root / "report.json", report)
    np.savez(root / "queries.npz", vectors=np.zeros((2, 60)))
    ev.write(root / "protocol.json", dict(parameters=params,
             artifact_catalog_sha256=ev.sha(batch / "artifacts/catalog.json"),
             faiss_manifest_sha256=ev.sha(batch / "faiss/manifest.json"),
             transform_sha256=ev.sha(batch / "faiss/transform.npz"), query_sha256=ev.sha(root / "queries.npz")))
    ev.write(root / "EVALUATION_COMPLETE.json", dict(report_sha256=ev.sha(root / "report.json"),
             acceptance_status="passed", calibration_status="gate_passed_on_panel"))
    return root


def test_acceptance_binds_current_library_and_rejects_corrupt_index(library, tmp_path):
    batch, _, _ = library
    evaluation = make_e033(batch, tmp_path / "e033")
    assert e.verify_acceptance(batch, evaluation)[1]["library_conformers"] == 12
    (batch / "faiss/index.faiss").write_bytes(b"altered index")
    with pytest.raises(ValueError, match="index checksum"):
        e.verify_acceptance(batch, evaluation)


@pytest.mark.parametrize("change", ["report", "transform", "acceptance", "protocol", "metadata_only"])
def test_e033_tamper_and_inadequate_acceptance_rejected(library, tmp_path, change):
    batch, _, _ = library
    evaluation = make_e033(batch, tmp_path / "e033")
    if change == "report":
        (evaluation / "report.json").write_text("{}")
    elif change == "transform":
        np.savez(batch / "faiss/transform.npz", mean=np.ones(60), std=np.ones(60))
    elif change == "acceptance":
        ev.write(evaluation / "acceptance.json", {})
    elif change == "protocol":
        protocol = ev.read(evaluation / "protocol.json")
        protocol["parameters"]["extra"] = 1
        ev.write(evaluation / "protocol.json", protocol)
    else:
        report = ev.read(evaluation / "report.json")
        report["acceptance"]["status"] = "metadata_passed_hashes_not_checked"
        ev.write(evaluation / "report.json", report)
        marker = ev.read(evaluation / "EVALUATION_COMPLETE.json")
        marker["report_sha256"] = ev.sha(evaluation / "report.json")
        ev.write(evaluation / "EVALUATION_COMPLETE.json", marker)
    with pytest.raises(ValueError):
        e.verify_acceptance(batch, evaluation)


def test_receipt_reuse_preserves_time_and_detects_changed_input_output(tmp_path):
    source = tmp_path / "source"; source.write_text("query")
    output = tmp_path / "output"
    def operation():
        output.write_text("result")
        return {"wall_seconds": 17.4}, [output]
    expected = e.checked_stage(tmp_path, "example", [source], operation)
    def never():
        pytest.fail("completed operation rerun")
    assert e.checked_stage(tmp_path, "example", [source], never) == expected
    source.write_text("other")
    with pytest.raises(ValueError, match="changed stage inputs"):
        e.checked_stage(tmp_path, "example", [source], never)
    source.write_text("query"); output.write_text("tampered")
    with pytest.raises(ValueError, match="changed completed output"):
        e.checked_stage(tmp_path, "example", [source], never)


def test_failed_operation_has_no_completed_receipt(tmp_path):
    def fail():
        raise RuntimeError("interrupted")
    with pytest.raises(RuntimeError):
        e.checked_stage(tmp_path, "stage", [], fail)
    assert not (tmp_path / "stage.stage.json").exists()


def test_gate_requires_both_queries_and_uses_fallback_only_if_passed():
    rows = [dict(nprobe=n, query=q, recall={"top1000_strict": r})
            for n, q, r in [(128, "a", .99), (128, "b", .94), (256, "a", .998), (256, "b", .99)]]
    assert e.choose_nprobe(rows, 2) == 256
    rows[1]["recall"]["top1000_strict"] = .95
    assert e.choose_nprobe(rows, 2) == 128
    rows[1]["recall"]["top1000_strict"] = .94
    rows[3]["recall"]["top1000_strict"] = .94
    assert e.choose_nprobe(rows, 2) is None
    assert e.choose_nprobe([rows[0], rows[0]], 2) is None


def test_timing_gate_uses_same_run_and_never_partial_resume():
    g = dict(stages=dict(refine=dict(wall_seconds_this_invocation=20., chunks_reused=0, chunks_computed=2)))
    assert e.timing_gate(g, 2.)["passed"] is True
    assert e.timing_gate(g, 3.)["passed"] is False
    g["stages"]["refine"]["chunks_reused"] = 1
    assert e.timing_gate(g, 1.)["passed"] is None
    assert e.timing_gate(g, 1.)["overhead_fraction"] is None
    with pytest.raises(ValueError):
        e.timing_gate(g, float("nan"))


def test_output_cannot_overwrite_inputs(tmp_path):
    batch, report, query = [tmp_path / x for x in ("batch", "e033", "query")]
    for output in (batch, batch / "new", report, report / "new", query, tmp_path):
        with pytest.raises(ValueError, match="separate"):
            e.check_output_location(output, batch, report, [query])


def test_retrieval_uses_query_specific_vectors_and_keeps_all_conformers(library, tmp_path, monkeypatch):
    batch, raw, _ = library
    catalog, acceptance = ev.accept_library(batch)
    index = ExactIndex(raw); index.nlist = 256
    monkeypatch.setitem(sys.modules, "faiss", SimpleNamespace(read_index=lambda _: index, omp_set_num_threads=lambda _: None))
    query_dirs = []
    for i in (0, 7):
        p = tmp_path / f"query-{i}"; p.mkdir(); np.save(p / "usrcat.npy", raw[i]); query_dirs.append(p)
    result, outputs = e.retrieve(batch, tmp_path / "retrieval", catalog, acceptance, query_dirs, 1)
    assert result["selected_nprobe"] == 128
    assert len(outputs) == 7
    for row in result["results"]:
        assert row["retained_conformers"] == 12
        assert row["recall"]["top1000_strict"] == 1
        with np.load(row["candidate_file"]) as c:
            assert len(c["global_ids"]) == 12
            # Different crystal descriptors have different nearest neighbors.
            assert c["global_ids"][0] == (0 if row["query"].startswith("8BJU") else 7)
            assert np.all(np.diff(c["squared_l2"]) >= 0)
            assert not any(key.startswith("cap_") for key in c.files)


def real_pose_fixture(tmp_path):
    batch = tmp_path / "batch"; artifacts = batch / "artifacts"; artifacts.mkdir(parents=True)
    old = _multi_artifact(artifacts, count=3)
    catalog = artifacts / "catalog.json"; catalog.write_bytes(old.read_bytes())
    chemical = batch / "chemical"; shard = chemical / "shard-multi"; shard.mkdir(parents=True)
    meta = np.zeros(3, dtype=CHEM_META_DTYPE)
    meta["global_id"] = np.arange(3); meta["feature_offset"] = [0, 2, 4]; meta["features"] = 2
    meta.tofile(shard / "chem-meta.bin")
    features = np.zeros(6, dtype=FEATURE_DIRECTION_DTYPE)
    features["direction"] = [1, 0, 0]; features["kind"] = 1
    features.tofile(shard / "feature-directions.bin")
    for name, dtype in [("atoms.bin", ATOM_DTYPE), ("bonds.bin", BOND_DTYPE), ("torsions.bin", TORSION_DTYPE),
                        ("feature-members.bin", "<u2"), ("torsion-members.bin", "<u2")]:
        np.empty(0, dtype=dtype).tofile(shard / name)
    for name in ("conformer_ids.bin", "molecule_ids.bin"):
        (shard / name).write_bytes((artifacts / "shard-multi" / name).read_bytes())
    ev.write(shard / "manifest.json", dict(conformers=3, global_id_start=0))
    ev.write(chemical / "catalog.json", dict(format="aidd-chemical-companion-catalog", version=1,
             artifact_v1_catalog_sha256=ev.sha(catalog), shards=[dict(name=shard.name, path=str(shard))]))
    run = tmp_path / "run"; query_dir = run / "query"; query_dir.mkdir(parents=True)
    write_gaussian_query(query_dir / "gaussian-query.npz", shape_points=[[0, 0, 0], [1, 0, 0], [0, 2, 0], [0, 0, 3]],
                         feature_points=[[0, 0, 0], [0, 0, 3]], feature_types=[1, 2], anchored_weights=[1.5, 1],
                         anchor_feature_indices=[0], feature_directions=[[1, 0, 0], [1, 0, 0]], feature_direction_kinds=[1, 1])
    candidates = run / "candidates.npz"
    np.savez(candidates, global_ids=np.asarray([0, 1, 2], dtype=np.int64), molecule_ids=np.asarray(["MOL-0", "MOL-1", "MOL-2"]))
    return run, batch, query_dir, candidates


def test_real_gaussian_and_e031_mmap_integration_with_completed_resume(tmp_path, monkeypatch):
    args = real_pose_fixture(tmp_path)
    result = e.run_pose_stages(*args, workers=1)
    assert result["interaction"]["counts"]["refined_conformers"] == 3
    assert result["interaction"]["counts"]["refined_molecules"] == 3
    assert result["interaction"]["counts"]["refine_molecule_retention_fraction"] == 1
    assert result["interaction"]["timing"]["status"] == "measured_same_run"
    assert set(result["interaction"]["manifest"]["analysis"]) == {
        "shape_only", "atomcentered_unweighted_joint", "atomcentered_anchored_joint"}
    def never(*args, **kwargs):
        pytest.fail("completed poses rerun")
    monkeypatch.setattr(e, "run_scaled_gaussian_reranking", never)
    monkeypatch.setattr(e, "score_interaction_matches", never)
    assert e.run_pose_stages(*args, workers=1) == result


def test_query_identity_checked_before_chemistry(tmp_path, monkeypatch):
    # A wrong manifest must fail even before parsing/constructing a molecule.
    source = tmp_path / "source"; source.mkdir()
    for name in ("8BJU.cif", "QT9.cif"):
        (source / name).write_text("not parsed")
    ev.write(source / "query_manifest.json", dict(query_id="1X8B:824:A:901", anchors=[{}]))
    with pytest.raises(ValueError, match="Expected locked query"):
        e.prepare_query(source, tmp_path / "output", e.QUERY_SPECS[0])


def test_failed_gate_report_distinguishes_execution_and_quality(tmp_path):
    report = dict(status="retrieval_gate_failed", retrieval=dict(results=[], selected_nprobe=None), poses={})
    e.render_report(report, tmp_path)
    text = (tmp_path / "report.md").read_text()
    assert "retrieval_gate_failed" in text
    assert "biological and docking quality remain unassessed" in text


def test_run_stops_before_pose_on_failed_gate_and_pins_resume(library, tmp_path, monkeypatch):
    batch, _, _ = library
    evaluation = make_e033(batch, tmp_path / "e033")
    (batch / "run.lock").touch()
    directories = []
    for spec in e.QUERY_SPECS:
        directory = tmp_path / spec[0]; directory.mkdir(); directories.append(directory)
        for name in (spec[2], spec[3], "query_manifest.json"):
            (directory / name).write_text("locked source fixture")
    monkeypatch.setitem(sys.modules, "fcntl", SimpleNamespace(LOCK_EX=1, LOCK_NB=2, flock=lambda *_: None))
    monkeypatch.setattr(e, "runtime_versions", lambda: dict(test="fixture"))
    monkeypatch.setattr(e, "ensure_file_descriptor_limit", lambda _: 4096)
    def prepare(directory, output, spec):
        output.mkdir(parents=True)
        np.save(output / "usrcat.npy", np.zeros(60, dtype=np.float32))
        return dict(query_id=spec[1]), [output / "usrcat.npy"]
    def retrieve(*args):
        return dict(results=[], selected_nprobe=None), []
    def no_poses(*args, **kwargs):
        pytest.fail("failed retrieval gate must not run Gaussian")
    monkeypatch.setattr(e, "prepare_query", prepare)
    monkeypatch.setattr(e, "retrieve", retrieve)
    monkeypatch.setattr(e, "run_pose_stages", no_poses)
    args = SimpleNamespace(batch=batch, e033=evaluation, output=tmp_path / "result", qt9=directories[0],
                           x8b=directories[1], workers=1, threads=1, resume=False)
    assert e.run(args) == 2
    assert ev.read(args.output / "RUN_STATUS.json")["status"] == "retrieval_gate_failed"
    args.resume = True
    monkeypatch.setattr(e, "prepare_query", no_poses)
    monkeypatch.setattr(e, "retrieve", no_poses)
    assert e.run(args) == 2
    args.workers = 2
    with pytest.raises(ValueError, match="Changed E034 inputs"):
        e.run(args)


@pytest.mark.parametrize("limits,raises,changes", [((1024, 65536), False, True),
                                                 ((65536, 65536), False, False),
                                                 ((1024, 4096), True, False)])
def test_many_shard_file_limit_checked_before_mapping(monkeypatch, limits, raises, changes):
    calls = []
    monkeypatch.setitem(sys.modules, "resource", SimpleNamespace(
        RLIMIT_NOFILE=7, RLIM_INFINITY=-1, getrlimit=lambda _: limits,
        setrlimit=lambda key, value: calls.append((key, value))))
    if raises:
        with pytest.raises(ValueError, match="hard limit"):
            e.ensure_file_descriptor_limit(365)
    else:
        assert e.ensure_file_descriptor_limit(365) == 9272
    assert bool(calls) == changes
    if changes:
        assert calls == [(7, (9272, 65536))]
