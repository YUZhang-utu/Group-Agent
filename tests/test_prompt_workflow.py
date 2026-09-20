import copy
import io
import json
from pathlib import Path
import sys
from types import SimpleNamespace
from urllib.error import HTTPError

import pytest

from aidd_agent import library_acceptance as ev
from aidd_agent.prompt_plan import validate_plan, strict_json, chat_plan
from aidd_agent.prompt_smoke import SMOKE_PLAN, OfflineServices, run_smoke
from aidd_agent.prompt_workflow import create_plan, run_plan, initialize_context
from aidd_agent.protein_data import fetch_protein, resolve_protein
from aidd_agent.registry import connect


def make_plan(tmp_path, steps=None):
    context = initialize_context(tmp_path / "home", "alice", "Prompt tests")
    payload = copy.deepcopy(SMOKE_PLAN)
    if steps is not None: payload["steps"] = steps
    fixture = tmp_path / "response.json"; ev.write(fixture, payload)
    path = create_plan(Path(context["db"]), context["user_id"], context["project_id"], "fixture prompt", response_file=fixture)
    return context, path


def run(context, path, **kwargs):
    return run_plan(Path(context["db"]), context["user_id"], context["project_id"], path, **kwargs)


@pytest.mark.parametrize("bad", ["command", "path", "future", "duplicate", "seed", "span", "query", "clarification"])
def test_invalid_model_actions_are_rejected(bad):
    plan = copy.deepcopy(SMOKE_PLAN)
    if bad == "command": plan["steps"][0]["action"] = "shell"
    if bad == "path": plan["steps"][0]["params"]["output"] = "../../victim"
    if bad == "future": plan["steps"][1]["params"]["protein_step"] = "prediction"
    if bad == "duplicate": plan["steps"][1]["id"] = "protein"
    if bad == "seed": plan["steps"][2]["params"]["seeds"] = [True]
    if bad == "span": plan["steps"][2]["params"]["end"] = 1
    if bad == "query": plan["steps"] = [dict(id="q", action="search_3d", params=dict(query="EGFR"))]
    if bad == "clarification": plan["clarifications"] = ["Which target?"]
    with pytest.raises(ValueError): validate_plan(plan)


def test_duplicate_json_and_nonfinite_rejected():
    for text in ('{"a":1,"a":2}', '{"a":NaN}'):
        with pytest.raises(ValueError): strict_json(text)


def test_compatible_chat_request_and_no_secret_in_payload(monkeypatch):
    monkeypatch.setenv("TEST_LLM_KEY", "secret-for-fixture")
    seen = []
    class Opener:
        def open(self, request, timeout):
            seen.append(request)
            return io.BytesIO(json.dumps(dict(choices=[dict(finish_reason="stop", message=dict(content=json.dumps(SMOKE_PLAN)))],
                                             id="fixture", usage=dict(total_tokens=123))).encode())
    plan, metadata = chat_plan("Prepare AF3 input for human WEE1", dict(base_url="https://example.org/v1", model="test",
                              api_key_env="TEST_LLM_KEY"), Opener())
    assert plan == SMOKE_PLAN
    request = seen[0]
    assert request.full_url == "https://example.org/v1/chat/completions"
    assert request.headers["Authorization"] == "Bearer secret-for-fixture"
    assert b"secret-for-fixture" not in request.data
    assert "secret-for-fixture" not in json.dumps(metadata)
    assert len(json.loads(request.data)["messages"]) == 2


def test_chat_errors_redact_body_and_reject_truncated_output(monkeypatch):
    class Error:
        def open(self, request, timeout):
            raise HTTPError(request.full_url, 401, "SECRET", {}, io.BytesIO(b"SECRET"))
    with pytest.raises(RuntimeError) as exc:
        chat_plan("hello", dict(base_url="https://example.org/v1", model="test"), Error())
    assert "401" in str(exc.value) and "SECRET" not in str(exc.value)
    class Truncated:
        def open(self, request, timeout):
            return io.BytesIO(json.dumps(dict(choices=[dict(finish_reason="length", message=dict(content="{}"))])).encode())
    with pytest.raises(ValueError, match="incomplete"):
        chat_plan("hello", dict(base_url="https://example.org/v1", model="test"), Truncated())
    with pytest.raises(ValueError, match="HTTPS"):
        chat_plan("hello", dict(base_url="http://public.example/v1", model="test"))


def test_protein_species_identity_and_ambiguity():
    record = fetch_protein("P30291", 9606, OfflineServices.fetch_json)
    assert record["length"] == 20
    with pytest.raises(ValueError, match="organism"):
        fetch_protein("P30291", 10090, OfflineServices.fetch_json)
    with pytest.raises(ValueError, match="unresolved"):
        resolve_protein("WEE1", 9606, lambda _: {"results": [{"primaryAccession": "P30291"}]*2})
    raw = OfflineServices.fetch_json("https://rest.uniprot.org/uniprotkb/P30291.json")
    raw["primaryAccession"] = "P00000"
    with pytest.raises(ValueError, match="accession mismatch"):
        fetch_protein("P30291", fetch=lambda _: raw)


def test_offline_smoke_and_context_idempotency(tmp_path):
    result = run_smoke(tmp_path / "smoke")
    assert result["status"] == "passed" and result["af3_inference"] == "not_run"
    a = initialize_context(tmp_path / "home")
    b = initialize_context(tmp_path / "home")
    assert a == b


def test_resume_tamper_and_project_authorization(tmp_path):
    context, path = make_plan(tmp_path)
    result, report = run(context, path, services=OfflineServices())
    assert result["status"] == "complete"
    with connect(Path(context["db"])) as connection:
        assert connection.execute("SELECT COUNT(*) FROM ai_recommendation").fetchone()[0] == 1
    other = initialize_context(tmp_path / "home", "bob", "Other")
    with pytest.raises(ValueError, match="escapes"):
        run(other, path, services=OfflineServices())
    (report.parent / "protein/protein.json").write_text("tampered")
    failed, _ = run(context, path, services=OfflineServices())
    assert failed["status"] == "failed" and failed["steps"]["structures"]["status"] == "not_run"
    ev.write(path, {})
    with pytest.raises(ValueError, match="ownership"):
        run(context, path, services=OfflineServices())


def test_af3_compute_gate_and_real_adapter_compilation(tmp_path):
    steps = copy.deepcopy(SMOKE_PLAN["steps"])
    steps.append(dict(id="fold", action="af3_run", params=dict(input_step="prediction")))
    context, path = make_plan(tmp_path, steps)
    runner = tmp_path / "installed_af3.py"; runner.write_text("# fixture runner")
    models = tmp_path / "models"; models.mkdir()
    database = tmp_path / "databases"; database.mkdir()
    profile = tmp_path / "af3.json"
    ev.write(profile, dict(backend="alphafold3", python=sys.executable, runner=str(runner),
                          model_parameters=str(models), databases=str(database), license_acknowledged=True))
    runtime = tmp_path / "runtime.json"; ev.write(runtime, dict(af3_profile=str(profile)))
    commands = []
    class FakeAF3(OfflineServices):
        @staticmethod
        def run_command(argv, log):
            commands.append(argv)
            output = Path(next(s.split("=", 1)[1] for s in argv if s.startswith("--output_dir=")))
            output.mkdir(parents=True)
            (output / "fixture_model.cif").write_text("data_fixture\n")
            ev.write(output / "fixture_summary_confidences.json", dict(ptm=.5, fixture=True))
            log.write_text("fixture execution, not AF3 inference")
    result, _ = run(context, path, runtime=runtime, services=FakeAF3())
    assert result["status"] == "blocked" and not commands
    result, report = run(context, path, runtime=runtime, services=FakeAF3(), allow_compute=True)
    assert result["status"] == "complete" and len(commands) == 1
    assert commands[0][:2] == [sys.executable, str(runner)]
    assert result["steps"]["fold"]["result"]["pose_quality"] == "not_independently_validated"
    run(context, path, runtime=runtime, services=FakeAF3(), allow_compute=True)
    assert len(commands) == 1
    runner.write_text("# changed runtime")
    with pytest.raises(ValueError, match="Changed workflow"):
        run(context, path, runtime=runtime, services=FakeAF3(), allow_compute=True)


def test_search_adapter_cannot_change_budget_or_e031_policy(tmp_path):
    context, path = make_plan(tmp_path, [dict(id="search", action="search_3d", params=dict(query="wee1_qt9"))])
    runtime = tmp_path / "runtime.json"
    ev.write(runtime, dict(search=dict(batch=str(tmp_path / "batch"), e034=str(tmp_path / "e034"),
             workers=2, coarse_chunk=500, refine_chunk=64, bounded_pair_seeds=False)))
    class FakeSearch(OfflineServices):
        @staticmethod
        def run_command(argv, log):
            assert argv[:3] == [sys.executable, "-m", "aidd_agent.fast_3d_search"]
            assert argv[argv.index("--query")+1] == "8bju"
            assert "--budget" not in argv and "--bounded-pair-seeds" not in argv
            output = Path(argv[argv.index("--output")+1]); output.mkdir()
            ev.write(output / "report.json", dict(status="complete"))
            ev.write(output / "RUN_STATUS.json", dict(status="complete", report_sha256=ev.sha(output / "report.json")))
    result, _ = run(context, path, runtime=runtime, services=FakeSearch(), allow_compute=True)
    assert result["status"] == "complete" and result["e031_changes_ranking"] is False
