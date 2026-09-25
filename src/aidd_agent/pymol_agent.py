"""Skill-grounded planning for the existing desktop bridge."""
import hashlib
import json
from pathlib import Path
import time

from .prompt_plan import chat_plan
from .pymol_bridge import submit, write_json
from .pymol_program import compile_program, METHODS, SETTINGS, LABELS, REPRESENTATIONS, MAX_CALLS, MAX_CODE_CHARS

SYSTEM = '''Control the already-open task structures in desktop PyMOL.
Return JSON with exactly explanation (English string) and calls (array).
Each call has exactly method (string) and arguments (object of named literal arguments).
Example: {"explanation":"Color protein.","calls":[{"method":"color",
"arguments":{"color":"cyan","selection":"v001 and polymer.protein"}}]}.
Do not generate Python source. All selections, object names and color names are JSON
strings, including identifiers starting with digits. No imports, loops or expressions.
Use an empty calls array when clarification is needed or the request is unsupported.
At most 512 calls and 64000 compiled characters are allowed. Combine common styling
across objects (protein cartoons and element colors); keep per-complex pocket/contact
operations separate. Do not drop requested complexes to meet the budget; ask to split
the request if necessary. The API schema lists parameter names and required counts.
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
For interaction display ALWAYS use the trusted typed_interactions method instead of
generating distance calls. It detects and deduplicates conservative typed hypotheses,
colors them, and replaces old contact lines. Arguments types and objects are optional
space-separated strings. Types: polar_contact salt_bridge pi_stacking cation_pi
halogen_bond hydrophobic. Default includes deduplicated hydrophobic contacts and processes enabled complexes
only. Explicit objects can select loaded complexes. Example call:
{"method":"typed_interactions","arguments":{"types":"polar_contact pi_stacking"}}.
It is a host adapter, not a native PyMOL API method. Never substitute all-pairs distances
when a type has zero results. Missing charges/aromaticity can prevent assignments.
Water bridges and metal coordination are not evaluated. Each type is deduplicated by residue pair.
Use typed_interactions as the LAST call after any representation/styling requests;
do not create distance objects after it. Its fixed per-type colors replace generic lines.
If a prior program failed and rollback succeeded, replace it with a complete corrected program.
For invalid_plan, nothing executed: inspect previous_attempt.plan and its error,
then return a complete corrected calls array, not just the repaired line.
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
    if not isinstance(value,dict) or set(value) not in ({'explanation','code'},{'explanation','calls'}):
        raise ValueError('PyMOL plan requires explanation and calls')
    if not isinstance(value['explanation'],str) or len(value['explanation'])>4000:
        raise ValueError('Invalid explanation')
    if 'calls' in value:
        calls=value['calls']
        if not isinstance(calls,list) or len(calls)>MAX_CALLS:raise ValueError(f'Use at most {MAX_CALLS} calls; combine shared styling across objects')
        lines=[]
        for call in calls:
            if not isinstance(call,dict) or set(call)!={'method','arguments'}:raise ValueError('Each call requires method and arguments')
            method=call['method'];args=call['arguments']
            if not isinstance(method,str) or method not in METHODS:raise ValueError('Unsupported method')
            if not isinstance(args,dict) or any(k not in METHODS[method][0] for k in args):raise ValueError('Unsupported argument names')
            lines.append('cmd.'+method+'('+', '.join(k+'='+repr(v) for k,v in args.items())+')')
        value=dict(explanation=value['explanation'],code='\n'.join(lines))
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
        captured={}
        def capture_and_validate(value):
            captured['plan']=value
            return validate_plan(value)
        try:
            plan,metadata=planner(json.dumps(context,ensure_ascii=False),profile,
                system_prompt=SYSTEM+'\nUpstream skill and references:\n'+skill,
                capabilities=dict(methods=METHODS,settings=sorted(SETTINGS),
                    labels=sorted(LABELS),representations=sorted(REPRESENTATIONS),
                    max_calls=MAX_CALLS,max_compiled_characters=MAX_CODE_CHARS),
                validator=capture_and_validate,max_prompt_chars=100000)
            plan=capture_and_validate(plan)
        except (ValueError,SyntaxError) as exc:
            failed=dict(status='invalid_plan',error=str(exc))
            if 'plan' in captured:
                raw=json.dumps(captured['plan'],ensure_ascii=False)
                failed['plan']=captured['plan'] if len(raw)<=70000 else dict(truncated_json=raw[:70000])
            if isinstance(exc,SyntaxError):failed.update(line=exc.lineno,column=exc.offset,source_line=exc.text)
            attempts.append(failed)
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
        result=dict(status='failed',message='No successful PyMOL program after two attempts',last_error=attempts[-1].get('error') or attempts[-1].get('result',{}).get('error'))
    answer=dict(status=result['status'],result=result,attempts=attempts,skill=provenance,
        validation='API execution and atom counts only; rendered image was not inspected by the LLM.',
        undo='Use /view {"operation":"undo"} to restore the last successful program checkpoint.')
    write_json(root/(initial['operation_id']+'-agent.json'),answer)
    return answer
