from pathlib import Path
import re

import pytest

from aidd_agent.language_policy import contains_han
from aidd_agent.prompt_plan import ACTION_FIELDS, CAPABILITIES, validate_plan
from aidd_agent.workflow_skills import WORKFLOW_SKILLS, skills_for_plan


def test_skill_contract_covers_actions_and_local_references():
    root = Path(__file__).resolve().parents[1]
    actions = [action for skill in WORKFLOW_SKILLS for action in skill["actions"]]
    assert len(actions) == len(set(actions))
    assert set(actions) == set(ACTION_FIELDS)
    assert CAPABILITIES["workflow_skills"] == list(WORKFLOW_SKILLS)
    for skill in WORKFLOW_SKILLS:
        path = root / ".agents/skills" / skill["name"] / "SKILL.md"
        text = path.read_text(encoding="utf-8")
        assert text.startswith("---\nname: " + skill["name"] + "\ndescription: ")
        assert not contains_han(text)
        for link in re.findall(r"\]\(([^)]+)\)", text):
            assert (path.parent / link).resolve().is_file()


def test_english_prose_guard_and_plan_routing():
    plan = dict(version=1, summary="Retrieve protein and prepare AF3", clarifications=[], steps=[
        dict(id="protein", action="protein_fetch", params=dict(accession="P30291")),
        dict(id="fold", action="af3_prepare", params=dict(protein_step="protein", name="wee1"))])
    assert validate_plan(plan) == plan
    assert skills_for_plan(plan) == ["aidd-protein-preparation", "aidd-alphafold3", "aidd-prompt-workflows"]
    plan["summary"] = chr(0x4E2D)
    with pytest.raises(ValueError, match="English"): validate_plan(plan)
    assert contains_han(chr(0x20000))
    assert not contains_han("P30291; alpha=0.95; 1e-6")
    assert skills_for_plan(dict(steps=[])) == ["aidd-prompt-workflows"]
