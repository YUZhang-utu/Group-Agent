"""Compact planner routing, paired with repository-discoverable Codex skills."""

WORKFLOW_SKILLS = (
    {"name": "aidd-3d-search", "actions": ["consensus_budget", "budget_page", "search_3d", "review_screening", "classify_screening", "full_library_screen", "condition_funnel", "benchmark_funnel", "select_screening", "export_screening", "prepare_docking", "run_docking"],
     "description": "Calibrated WEE1 full-library retrieval and Gaussian refinement; E031 annotations only."},
    {"name": "aidd-protein-preparation", "actions": ["protein_fetch", "protein_resolve", "pdb_search", "pdb_fetch", "structure_survey", "structure_diversity", "structure_consensus", "consensus_recommend", "consensus_design", "consensus_funnel", "anchor_recommend", "anchor_design", "guided_funnel", "guided_select"],
     "description": "Verified UniProt target identity, sequences and RCSB structure evidence; no automatic receptor acceptance."},
    {"name": "aidd-alphafold3", "actions": ["af3_prepare", "af3_run"],
     "description": "Verified single-protein/CCD AF3 input and installed local execution; prediction is not biological validation."},
    {"name": "aidd-prompt-workflows", "actions": [],
     "description": "English structured planning, typed dependencies, Project scope and accurate blocked/failed states."},
)


def skills_for_plan(plan):
    actions = {step["action"] for step in plan["steps"]}
    return [skill["name"] for skill in WORKFLOW_SKILLS
            if not skill["actions"] or actions.intersection(skill["actions"])]
