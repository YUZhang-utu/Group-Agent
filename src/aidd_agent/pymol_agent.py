"""Skill-grounded planning for the existing desktop bridge."""
import hashlib
import json
from pathlib import Path
import time

from .prompt_plan import chat_plan
from .pymol_bridge import submit, write_json
from .pymol_program import compile_program, METHODS, SETTINGS, LABELS, REPRESENTATIONS

SYSTEM = '''Control the already-open task structures in desktop PyMOL.
Return JSON with exactly explanation and code (both strings). explanation must be English.
code contains only direct cmd.method calls with literal arguments, no imports, variables,
loops, functions, cmd.do, loading, saving, quitting, or arbitrary Python expressions.
An empty code is allowed when clarification is needed or the request is unsupported.
Use the provided scene: do not invent object IDs, chains, ligand identifiers or residues.
The user's request and the bounded API contract override upstream headless setup recipes.
Do not start/stop PyMOL, install packages, or use uv. Check CA-only structures before cartoons.
Atom selections are restricted to loaded task objects; default selection is all task atoms.
Use exact ligand_selection from the scene, especially for peptide ligands.
For multiple complexes, calculate neighborhoods and distances separately within each object.
Keep entire residues using byres for a ligand neighborhood. Never delete a ligand when
removing redundant protein chains. If several chains contact the ligand, ask which to keep.
For mixed coloring, assign different carbon colors then color N blue, O red, S yellow.
Generated selection/distance/color names must start ai_. delete supports exact ai_ names only.
Atom selections support letters, numbers, spaces, parentheses, . + * < > ! = - only.
Use quoted literal label expressions from the schema. Do not style measurement objects
with atom selections. distance mode 2 draws possible polar contacts, not established hydrogen
bonds or all interactions. Hydrophobic, pi, salt-bridge and water-mediated classification
requires separate scientific analysis. Do not invent these from distance lines alone.
PNG and PSE are saved automatically; session checkpoint and one-step undo are provided.
The model sees metadata and API receipts, not the rendered image. Do not claim visual review.
If a prior program failed and rollback succeeded, replace it with a complete corrected program.
'''


def skill_context():
    root=Path(__file__).parent/'skill_data'/'pymol'
    provenance=json.loads((root/'provenance.json').read_text(encoding='utf-8'))
    for name,digest in provenance['files'].items():
        if hashlib.sha256((root/name).read_bytes()).hexdigest()!=digest:
            raise ValueError('Bundled PyMOL skill integrity check failed: '+name)
    content='\n\n'.join((root/name).read_text(encoding='utf-8') for name in
        ('SKILL.md','references/PYMOL_REFERENCE.md','references/RECIPES.md'))
    return content,provenance


def validate_plan(value):
    if not isinstance(value,dict) or set(value)!={'explanation','code'}:
        raise ValueError('PyMOL plan requires explanation and code')
    if not isinstance(value['explanation'],str) or len(value['explanation'])>4000:
        raise ValueError('Invalid explanation')
    if not isinstance(value['code'],str):raise ValueError('Invalid code')
    if value['code'].strip():compile_program(value['code'])
    return value


def wait_result(queued,seconds=20):
    deadline=time.monotonic()+seconds
    path=Path(queued['result'])
    while time.monotonic()<deadline:
        if path.exists():return json.loads(path.read_text(encoding='utf-8'))
        time.sleep(.2)
    return dict(**queued,pending=True)


def run_agent(request,profile,root,catalog,executable=None,*,planner=None,dispatch=None,waiter=None):
    planner=planner or chat_plan;dispatch=dispatch or submit;waiter=waiter or wait_result
    root=Path(root)
    skill,provenance=skill_context()
    initial=dispatch(root,dict(operation='scene'),catalog,executable)
    observed=waiter(initial)
    if observed.get('status')!='complete':return observed
    scene=observed['scene']
    if not scene['objects']:raise ValueError('Open the task structures in PyMOL first')
    # Keep all object IDs and chain inventories; truncate only residue examples.
    scene=dict(scene,objects=[dict(row,protein_residues=row.get('protein_residues',[])[:20],
        residues_truncated=row.get('residue_count',0)>20) for row in scene['objects']])
    context=dict(request=request,scene=scene,previous_attempt=None)
    attempts=[]
    for attempt in range(2):
        try:
            plan,metadata=planner(json.dumps(context,ensure_ascii=False),profile,
                system_prompt=SYSTEM+'\nUpstream skill and references:\n'+skill,
                capabilities=dict(methods=METHODS,settings=sorted(SETTINGS),
                    labels=sorted(LABELS),representations=sorted(REPRESENTATIONS)),
                validator=validate_plan,max_prompt_chars=100000)
            validate_plan(plan)
        except (ValueError,SyntaxError) as exc:
            attempts.append(dict(status='invalid_plan',error=str(exc)))
            context['previous_attempt']=attempts[-1]
            continue
        if not plan['code'].strip():
            result=dict(status='clarification_required',message=plan['explanation'])
            attempts.append(result)
            break
        queued=dispatch(root,dict(operation='program',code=plan['code']),catalog,executable)
        result=waiter(queued)
        attempts.append(dict(plan=plan,model=metadata,result=result))
        if result.get('status')!='failed' or result.get('rollback')!='complete':break
        context['previous_attempt']=dict(code=plan['code'],error=result.get('error'),rollback='complete')
    else:
        result=dict(status='failed',message='No successful PyMOL program after two attempts')
    answer=dict(status=result['status'],result=result,attempts=attempts,skill=provenance,
        validation='API execution and atom counts only; rendered image was not inspected by the LLM.',
        undo='Use /view {"operation":"undo"} to restore the last successful program checkpoint.')
    write_json(root/(initial['operation_id']+'-agent.json'),answer)
    return answer
