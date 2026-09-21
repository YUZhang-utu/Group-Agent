import json

import pytest

from aidd_agent.docking_discovery import discover


def test_explicit_missing_path_does_not_use_unrelated_path(tmp_path, monkeypatch):
    monkeypatch.setattr("aidd_agent.docking_discovery.shutil.which", lambda _: "/unrelated/tool")
    profile = tmp_path / "profile.json"
    profile.write_text(json.dumps({"schrodinger_root": str(tmp_path),
                                  "plants_executable": str(tmp_path / "missing")}))
    result = discover(profile)
    assert all(v is None for v in result["commands"].values())
    assert result["checks"]["PLANTS"]["status"] == "missing"


def test_home_expansion_and_execute_permission(tmp_path, monkeypatch):
    monkeypatch.setattr("aidd_agent.docking_discovery.Path.home", lambda: tmp_path)
    monkeypatch.setenv("HOME", str(tmp_path))
    monkeypatch.setenv("USERPROFILE", str(tmp_path))
    monkeypatch.setattr("aidd_agent.docking_discovery.shutil.which", lambda _: None)
    binary = tmp_path / "PLANTS1.2" / "PLANTS1.2_64bit"
    binary.parent.mkdir()
    binary.write_text("must never run")
    profile = tmp_path / "profile.json"
    profile.write_text(json.dumps({"plants_executable": "~/PLANTS1.2/PLANTS1.2_64bit"}))
    monkeypatch.setattr("aidd_agent.docking_discovery.os.access", lambda *_: True)
    assert discover(profile)["commands"]["PLANTS"] == str(binary)
    monkeypatch.setattr("aidd_agent.docking_discovery.os.access", lambda *_: False)
    assert discover(profile)["checks"]["PLANTS"]["status"] == "not_executable"


def test_profile_rejects_shell_fields(tmp_path):
    profile = tmp_path / "profile.json"
    profile.write_text('{"command": "arbitrary shell"}')
    with pytest.raises(ValueError, match="Unsupported"):
        discover(profile)
