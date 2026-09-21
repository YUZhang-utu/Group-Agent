import copy
import io
import json
from pathlib import Path, PurePosixPath

import pytest

from aidd_agent import prediction, prompt_workflow as workflow
from aidd_agent.llm_profiles import select_llm_profile
from aidd_agent.prompt_plan import chat_plan
from aidd_agent.prompt_smoke import SMOKE_PLAN, OfflineServices

ROOT = Path(__file__).resolve().parents[1]


def profile():
    data = json.loads((ROOT / "configs/models/alphafold3-apptainer.example.json").read_text())
    data["license_acknowledged"] = True
    return data


def test_container_command_preserves_mounts_and_xla():
    p = profile()
    argv = prediction.prediction_command(p, PurePosixPath("/jobs/input space/af3-input.json"), PurePosixPath("/jobs/output space"))
    assert argv[:3] == ["apptainer", "exec", "--nv"]
    assert argv[argv.index(p["image"])+1:][:2] == ["python", "/app/alphafold/run_alphafold.py"]
    binds = [argv[i+1] for i, x in enumerate(argv) if x == "--bind"]
    assert binds == ["/jobs/input space:/root/af_input:ro", "/jobs/output space:/root/af_output:rw",
                     p["model_parameters"]+":/root/models:ro", p["databases"]+":/root/public_databases:ro"]
    assert "--json_path=/root/af_input/af3-input.json" in argv
    assert "--flash_attention_implementation=xla" in argv
    with pytest.raises(ValueError, match="bind paths"):
        prediction.prediction_command(p, PurePosixPath("/jobs/bad,dir/input.json"), PurePosixPath("/jobs/out"))


def test_installation_checks_host_not_container_paths(tmp_path, monkeypatch):
    p = profile()
    for key in ("model_parameters", "databases"):
        d = tmp_path / key; d.mkdir(); p[key] = str(d)
    image = tmp_path / "af3.sif"; image.write_bytes(b"fixture image")
    p["image"] = str(image)
    monkeypatch.setattr(prediction.shutil, "which", lambda _: "/usr/bin/apptainer")
    prediction.validate_af3_installation(p)
    monkeypatch.setattr(prediction.shutil, "which", lambda _: None)
    with pytest.raises(ValueError, match="unavailable"):
        prediction.validate_af3_installation(p)
    image.unlink()
    with pytest.raises(ValueError, match="image"):
        prediction.validate_af3_installation(p)


def test_provider_precedence_and_missing_profile(tmp_path, monkeypatch):
    for key in ("AIDD_LLM_PROFILE", "AIDD_LLM_PROVIDER", "AIDD_LLM_CONFIG_DIR"):
        monkeypatch.delenv(key, raising=False)
    assert select_llm_profile().name == "gpt.json"
    monkeypatch.setenv("AIDD_LLM_PROFILE", str(ROOT / "configs/llm/gpt.json"))
    assert select_llm_profile("deepseek").name == "deepseek.json"
    assert select_llm_profile().name == "gpt.json"
    with pytest.raises(ValueError, match="OR"):
        select_llm_profile("gpt", ROOT / "configs/llm/gpt.json")
    with pytest.raises(ValueError, match="not found"):
        select_llm_profile("gpt", config_dir=tmp_path)
    monkeypatch.delenv("AIDD_LLM_PROFILE")
    monkeypatch.setenv("AIDD_LLM_PROVIDER", "deepseek")
    assert select_llm_profile().name == "deepseek.json"


@pytest.mark.parametrize("provider,key,url", [
    ("gpt", "OPENAI_API_KEY", "https://api.openai.com/v1/chat/completions"),
    ("deepseek", "DEEPSEEK_API_KEY", "https://api.deepseek.com/chat/completions")])
def test_provider_request_credential_isolation(provider, key, url, monkeypatch):
    monkeypatch.setenv("OPENAI_API_KEY", "gpt-fixture-secret")
    monkeypatch.setenv("DEEPSEEK_API_KEY", "deepseek-fixture-secret")
    cfg = json.loads(select_llm_profile(provider).read_text())
    class Opener:
        def open(self, req, timeout):
            assert req.full_url == url
            assert req.headers["Authorization"] == "Bearer " + ("gpt-fixture-secret" if key == "OPENAI_API_KEY" else "deepseek-fixture-secret")
            assert b"fixture-secret" not in req.data
            payload = json.loads(req.data)
            assert payload["model"] == cfg["model"]
            assert payload["response_format"] == {"type": "json_object"}
            return io.BytesIO(json.dumps({"choices": [{"finish_reason": "stop", "message": {"content": json.dumps(SMOKE_PLAN)}}]}).encode())
    plan, meta = chat_plan("Prepare WEE1 input", cfg, Opener())
    assert plan == SMOKE_PLAN and "fixture-secret" not in json.dumps(meta)


def test_container_failure_retry_reuse_and_image_change(tmp_path, monkeypatch):
    ctx = workflow.initialize_context(tmp_path / "home", "alice", "Container tests")
    payload = copy.deepcopy(SMOKE_PLAN)
    payload["steps"].append(dict(id="fold", action="af3_run", params=dict(input_step="prediction")))
    response = tmp_path / "response.json"; response.write_text(json.dumps(payload))
    plan = workflow.create_plan(Path(ctx["db"]), ctx["user_id"], ctx["project_id"], "fixture", response_file=response)
    p = profile()
    for key in ("model_parameters", "databases"):
        d = tmp_path / key; d.mkdir(); p[key] = str(d)
    image = tmp_path / "af3.sif"; image.write_bytes(b"fixture image"); p["image"] = str(image)
    pp = tmp_path / "profile.json"; pp.write_text(json.dumps(p))
    runtime = tmp_path / "runtime.json"; runtime.write_text(json.dumps(dict(af3_profile=str(pp))))
    monkeypatch.setattr(prediction.shutil, "which", lambda _: "/usr/bin/apptainer")
    # Generate Linux argv while the simulated host may be Windows.
    def command(cfg, input_path, output):
        assert input_path.is_file() and output.name.startswith("prediction-attempt-")
        cfg = dict(cfg, model_parameters="/models", databases="/db")
        return prediction.prediction_command(cfg, PurePosixPath("/task/af3-input.json"), PurePosixPath("/task/output"))
    monkeypatch.setattr(workflow, "prediction_command", command)
    calls = []
    class Fake(OfflineServices):
        @staticmethod
        def run_command(argv, log):
            calls.append(argv)
            assert argv[:3] == ["apptainer", "exec", "--nv"]
            output = sorted(log.parent.glob("prediction-attempt-*"))[-1]
            assert output.is_dir()
            if len(calls) == 1:
                (output / "partial.txt").write_text("failed fixture")
                raise RuntimeError("fixture failure")
            (output / "fixture_model.cif").write_text("data_fixture")
            (output / "fixture_summary_confidences.json").write_text('{"ptm": 0.5}')
            log.write_text("simulated AF3")
    def run(compute=True):
        return workflow.run_plan(Path(ctx["db"]), ctx["user_id"], ctx["project_id"], plan,
                                 runtime=runtime, allow_compute=compute, services=Fake())[0]
    assert run(False)["status"] == "blocked" and not calls
    assert run()["status"] == "failed"
    assert run()["status"] == "complete"
    assert (plan.parent / "execution/fold/prediction-attempt-001/partial.txt").is_file()
    assert run()["status"] == "complete" and len(calls) == 2
    image.write_bytes(b"changed image")
    with pytest.raises(ValueError, match="Changed workflow"):
        run()
