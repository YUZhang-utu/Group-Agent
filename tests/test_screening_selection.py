import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from aidd_agent import screening_selection as selection
from aidd_agent.expanded_wee1 import fingerprint
from aidd_agent.gaussian_batch import write_gaussian_query
from aidd_agent.library_acceptance import sha
from aidd_agent.prompt_workflow import create_plan, run_plan
from aidd_agent.chat_agent import ChatAgent, screening_summary, validate_route

O = selection.OBJECTIVE
Q = "8BJU:QT9:A:601"


def write(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data), encoding="utf-8")
    return path


def fixture_search(root):
    root.mkdir(parents=True, exist_ok=True)
    query_dir = root.parent / "crystal"
    mmcif = write(query_dir / "crystal.cif", {"fixture": True})
    ccd = write(query_dir / "ccd.cif", {"fixture": True})
    anchors = [dict(anchor_id=f"A{i}:HBA:{i}", feature_type="HBA", ligand_atom_indices=[i],
        atom_center=[i, 0, 0], weight=1.5, evidence=dict(interaction="direct_hydrogen_bond_geometry",
            protein_partner=dict(chain="A", residue="ASN", residue_number=str(100+i), atom_name="N"),
            distance_angstrom=3.0, angle_degrees=None, angle_status="not_evaluated_no_explicit_hydrogen")) for i in range(2)]
    manifest = write(query_dir / "query_manifest.json", dict(query_id=Q, anchors=anchors))
    query = query_dir / "gaussian-query.npz"
    source = dict(query_id=Q, anchor_mapping={a["anchor_id"]: i for i, a in enumerate(anchors)})
    for name, path in (("mmcif", mmcif), ("ccd", ccd), ("query_manifest", manifest)):
        source[name], source[name + "_sha256"] = str(path), sha(path)
    write_gaussian_query(query, shape_points=[[0,0,0]], feature_points=[[0,0,0],[1,0,0]], feature_types=[1,1],
                         anchored_weights=[1.5,1.5], anchor_feature_indices=[0,1], source=source)
    artifact = write(root.parent / "catalog.json", {"shards": []})
    chemical = write(root.parent / "chemical.json", {})
    rigid = root / "8bju/gaussian/refine/merged-scores.npz"
    rigid.parent.mkdir(parents=True)
    identities = dict(global_ids=np.array([1,2,3,4]), molecule_ids=np.array(["m1","m1","m2","m3"]),
                      conformer_ids=np.array(["c1","c2","c3","c4"]))
    transforms = np.tile(np.eye(4), (4,1,1))
    np.savez(rigid, **identities, **{O+"__objective":np.array([.9,.8,.7,.7]), O+"__transform":transforms})
    side = root / "8bju/interaction-matches.npz"
    np.savez(side, **identities, query_anchor_feature_indices=np.array([0,1]),
             **{O+"__anchor_scores":np.array([[.8,.1],[.1,.8],[.7,.7],[.7,.7]]),
                O+"__anchor_assignments":np.tile([0,1], (4,1))})
    write(side.with_suffix(".manifest.json"), dict(output_sha256=sha(side), inputs={
        key:dict(path=str(p),sha256=sha(p)) for key,p in
        (("query",query),("rigid_result",rigid),("artifact_catalog",artifact),("chemical_companion",chemical))}))
    schedule = root / "8bju/candidates.npz"
    np.savez(schedule, **identities)
    baseline = write(root.parent / "e034/report.json", dict(acceptance=dict(library_conformers=100, library_molecules=50)))
    write(root / "protocol.json", dict(sources=fingerprint([baseline, query.with_suffix(".manifest.json")])))
    config = dict(query=str(query), candidate_schedule=str(schedule), candidate_schedule_sha256=sha(schedule))
    path = write(root / "report.json", dict(status="complete", queries=[dict(query_id=Q,
        gaussian_equivalence=dict(passed=True), annotation_equivalence=dict(passed=True),
        gaussian=dict(config=config,final_result=str(rigid),final_result_sha256=sha(rigid)))]))
    write(root / "RUN_STATUS.json", dict(status="complete", report_sha256=sha(path)))
    return path


def test_evidence_and_same_pose_selection(tmp_path):
    search = fixture_search(tmp_path / "search")
    evidence = selection.review(search, tmp_path / "review")
    assert evidence["library"]["library_molecules"] == 50
    assert evidence["refined_molecule_union"] == 3
    anchors = [a["anchor_id"] for a in evidence["queries"][0]["anchors"]]
    assert evidence["queries"][0]["anchors"][0]["evidence"]["protein_partner"]["residue"] == "ASN"
    policy = dict(required_anchors=anchors, match_mode="all", minimum_score=.5)
    result = selection.preview(tmp_path / "review/report.json", tmp_path / "preview", policy)
    assert result["counts"]["matching_conformers"] == 2
    assert [r["molecule_id"] for r in result["representatives"]] == ["m2","m3"]
    assert not (tmp_path / "preview/selected-poses.sdf").exists()
    result = selection.preview(tmp_path / "review/report.json", tmp_path / "any", dict(policy,match_mode="any",max_molecules=1))
    assert result["counts"]["matching_molecules"] == 3
    assert result["counts"]["matching_conformers"] == 4
    assert result["representatives"][0]["global_id"] == 1
    assert result["counts"]["selected_molecules"] == 1


def test_zero_export_and_tampering(tmp_path):
    search = fixture_search(tmp_path / "search")
    evidence = selection.review(search, tmp_path / "review")
    policy = dict(required_anchors=[evidence["queries"][0]["anchors"][0]["anchor_id"]], match_mode="all", minimum_score=1.)
    selection.preview(tmp_path / "review/report.json", tmp_path / "preview", policy)
    result = selection.export(tmp_path / "preview/report.json", tmp_path / "export")
    assert result["counts"]["selected_molecules"] == 0
    assert Path(result["sdf"]).read_text() == ""
    assert json.loads(Path(result["ids"]).read_text()) == []
    Path(evidence["queries"][0]["sidecar"]).write_bytes(b"modified")
    with pytest.raises(ValueError, match="changed"):
        selection.preview(tmp_path / "review/report.json", tmp_path / "changed", policy)


@pytest.mark.parametrize("value", [0, -1, float("nan"), 1.1, True])
def test_invalid_threshold(value):
    with pytest.raises(ValueError):
        selection.validate_selection(dict(required_anchors=["anchor"],match_mode="all",minimum_score=value))


def test_export_preserves_ids_and_writes_only_selected_molecules(tmp_path, monkeypatch):
    search = fixture_search(tmp_path / "search")
    evidence = selection.review(search, tmp_path / "review")
    policy = dict(required_anchors=[a["anchor_id"] for a in evidence["queries"][0]["anchors"]], match_mode="all",minimum_score=.5)
    selection.preview(tmp_path / "review/report.json", tmp_path / "preview", policy)
    def get(gid):
        return SimpleNamespace(global_id=gid, molecule_id="m"+str(gid-1),conformer_id="c"+str(gid),
            shape_points=np.array([[0.,0.,0.]]), atomic_numbers=[6],bonds=[],formal_charges=[0])
    reader = SimpleNamespace(get=get)
    monkeypatch.setattr("aidd_agent.gaussian_batch.ArtifactCatalogReader", lambda _: reader)
    monkeypatch.setattr("aidd_agent.chemical_companion.ChemicalCompanionReader", lambda _: reader)
    monkeypatch.setattr("aidd_agent.expanded_wee1.ensure_file_descriptor_limit", lambda _: None)
    result = selection.export(tmp_path / "preview/report.json", tmp_path / "export")
    assert result["counts"]["selected_molecules"] == 2
    text = Path(result["sdf"]).read_text()
    assert text.count("$$$$") == 2 and "AIDD_MOLECULE_ID" in text
    assert result["docking_status"].startswith("not_run")


def attach_fixture(app, sid):
    ctx = app.context
    plan = create_plan(Path(ctx["db"]),ctx["user_id"],ctx["project_id"],"Synthetic search",local_plan=dict(
        version=1,summary="Synthetic search",clarifications=[],steps=[dict(id="search",action="search_3d",params=dict(query="wee1_qt9"))]))
    child = fixture_search(plan.parent / "execution/search/search")
    result = dict(status="search_completed",report=str(child))
    write(plan.parent / "execution/search.stage.json", dict(inputs={},outputs=fingerprint([p for p in child.parent.rglob("*") if p.is_file()]),result=result))
    report = write(plan.parent / "execution/report.json", dict(status="complete",steps=dict(search=dict(status="complete",action="search_3d",result=result))))
    write(report.parent / "RUN_STATUS.json",dict(status="complete",report_sha256=sha(report)))
    return app.attach(sid, str(plan))


def test_chat_review_preview_export_separate_tasks_and_owned_sources(tmp_path):
    app = ChatAgent(tmp_path,start=False)
    sid = app.new_session()
    try:
        source = attach_fixture(app, sid)
        app.ask(sid,"/evidence " + source,"deepseek")
        app.execute(app.task(sid))
        evidence_job = app.task(sid)
        assert evidence_job["status"] == "complete", app.snapshot(sid)
        summary = screening_summary(evidence_job)
        anchors = [a["anchor_id"] for a in summary["queries"][0]["anchors"]]
        app.router = lambda *args: dict(intent="select",message="Preview selection",request="",task_id=evidence_job["id"],
            selection=dict(required_anchors=anchors,match_mode="all",minimum_score=1.))
        app.ask(sid,"Keep both anchors at score 1; preview only.","deepseek")
        app.execute(app.task(sid))
        preview_job = app.task(sid)
        assert preview_job["status"] == "complete", app.snapshot(sid)
        assert screening_summary(preview_job)["counts"]["selected_molecules"] == 0
        assert not list(Path(preview_job["plan"]).parent.rglob("*.sdf"))
        other = app.new_session()
        with pytest.raises(ValueError,match="No matching"):
            app.ask(other,"/export " + preview_job["id"],"deepseek")
        app.ask(sid,"/export " + preview_job["id"],"deepseek")
        app.execute(app.task(sid))
        assert app.task(sid)["status"] == "complete", app.snapshot(sid)
        assert "not_run_requires_engine_specific_preparation" in app.ask(sid,"/results","deepseek")
        assert len(app.jobs(sid)) == 4
    finally:
        app.close()


def test_query_metadata_tampering_rejected(tmp_path):
    search = fixture_search(tmp_path / "search")
    manifest = tmp_path / "crystal/gaussian-query.manifest.json"
    data = json.loads(manifest.read_text()); data["source"]["anchor_mapping"] = {"forged":0}
    write(manifest,data)
    with pytest.raises(ValueError,match="locked search"):
        selection.review(search,tmp_path / "review")


def test_unknown_and_cross_query_anchor_policies_fail(tmp_path):
    search = fixture_search(tmp_path / "search")
    evidence = selection.review(search, tmp_path / "review")
    policy = dict(required_anchors=["missing"],match_mode="all",minimum_score=.5)
    with pytest.raises(ValueError,match="Unknown anchor"):
        selection.preview(tmp_path / "review/report.json",tmp_path / "preview",policy)
    other = json.loads(json.dumps(evidence["queries"][0]))
    other["query_id"] = "1X8B:824:A:901"
    for a in other["anchors"]: a["anchor_id"] = a["anchor_id"].replace(Q,other["query_id"])
    evidence["queries"].append(other)
    write(tmp_path / "review/report.json",evidence)
    policy["required_anchors"] = [evidence["queries"][0]["anchors"][0]["anchor_id"],other["anchors"][0]["anchor_id"]]
    with pytest.raises(ValueError,match="one crystal query"):
        selection.preview(tmp_path / "review/report.json",tmp_path / "preview",policy)


def test_source_tampering_is_reported_as_failed_task(tmp_path):
    app = ChatAgent(tmp_path,start=False)
    sid = app.new_session()
    try:
        source = attach_fixture(app,sid)
        source_job = app.task(sid,source)
        app.ask(sid,"/evidence " + source,"deepseek")
        data = json.loads(Path(source_job["report"]).read_text()); data["changed"] = True
        write(Path(source_job["report"]),data)
        app.execute(app.task(sid))
        assert app.task(sid)["status"] == "failed"
        assert "Source task must have a completed report receipt" in Path(app.task(sid)["report"]).read_text()
    finally: app.close()


def test_select_router_requires_policy_and_export_cannot_include_one():
    with pytest.raises(ValueError):
        validate_route(dict(intent="select",message="",task_id=None,request=""))
    with pytest.raises(ValueError):
        validate_route(dict(intent="export",message="",task_id=None,request="",selection={}))


@pytest.mark.parametrize("different", [False, True])
def test_legacy_mapping_requires_exact_query_reconstruction(tmp_path, monkeypatch, different):
    search = fixture_search(tmp_path / "search")
    query = tmp_path / "crystal/gaussian-query.npz"
    metadata_path = query.with_suffix(".manifest.json")
    metadata = json.loads(metadata_path.read_text())
    mapping = metadata["source"].pop("anchor_mapping")
    write(metadata_path,metadata)
    protocol_path = search.parent / "protocol.json"
    protocol = json.loads(protocol_path.read_text())
    protocol["sources"][str(metadata_path)] = sha(metadata_path)
    write(protocol_path,protocol)
    def rebuild(mmcif, ccd, manifest, output):
        arrays = selection.archive(query)
        if different: arrays["feature_points"][0,0] += .01
        np.savez(output, **arrays)
        return dict(source=dict(anchor_mapping=mapping))
    monkeypatch.setattr(selection,"prepare_gaussian_query",rebuild)
    if different:
        with pytest.raises(ValueError,match="reconstruction differs"):
            selection.review(search,tmp_path / "review")
    else:
        result = selection.review(search,tmp_path / "review")
        assert [a["score_column"] for a in result["queries"][0]["anchors"]] == [0,1]
