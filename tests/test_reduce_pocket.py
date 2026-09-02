import json
from pathlib import Path
from types import SimpleNamespace

from aidd_agent.reduce_pocket import load_reduce_profile, run_reduce


def test_reduce_profile_and_execution_manifest(tmp_path: Path):
    profile = tmp_path / "reduce.json"
    profile.write_text(json.dumps({"executable": "/opt/reduce", "het_dictionary": None,
                                   "build_arguments": ["-build"]}))
    pocket = tmp_path / "pocket.pdb"; pocket.write_text("ATOM\n")
    calls = []

    def runner(command, **kwargs):
        calls.append((command, kwargs))
        return SimpleNamespace(returncode=0, stdout="ATOM H\n", stderr="USER  MOD Single : A 10 HIS\n")

    manifest = run_reduce(profile, pocket, tmp_path / "out", runner=runner)
    assert calls[0][0] == ["/opt/reduce", "-build", str(pocket.resolve())]
    assert manifest["scope"].endswith("within_5A")
    assert Path(manifest["output"]["path"]).read_text() == "ATOM H\n"


def test_reduce_profile_rejects_non_flag_arguments(tmp_path: Path):
    profile = tmp_path / "bad.json"
    profile.write_text(json.dumps({"executable": "reduce", "build_arguments": ["input.pdb"]}))
    try:
        load_reduce_profile(profile)
    except ValueError as exc:
        assert "flags" in str(exc)
    else:
        raise AssertionError("unsafe argument was accepted")
