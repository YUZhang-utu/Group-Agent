"""Persistent conversations and a single background scientific job queue."""
from __future__ import annotations
import json
from contextlib import contextmanager
import os
from pathlib import Path
import signal
import sqlite3
import subprocess
import sys
import threading
import time
import uuid

from .llm_profiles import select_llm_profile
from .prompt_plan import CAPABILITIES, chat_plan
from .prompt_workflow import create_plan, initialize_context
from .language_policy import contains_han

WORKFLOWS = [
    {"name":"Structure-guided chat", "status":"survey, recommend, adopt/design and guided full-library adapters",
     "scope":"Mapped crystal recurrence, literature fallback, evidence-grounded LLM advice and editable coordinate-backed designs; mandatory geometry and receptor exclusion precede Gaussian. Exploratory, not calibrated new-target affinity."},
    {"name": "Protein and structure evidence", "status": "available", "scope": "Verified UniProt proteins and PDB retrieval; receptor selection still needs review."},
    {"name": "AF3 prediction", "status": "available", "scope": "Verified single protein and optional CCD ligands; installed native or Apptainer profile."},
    {"name": "3D screening", "status": "WEE1 templates; exhaustive coarse retrieval by default", "scope": "All conformer descriptors compared before budgeted Gaussian refinement. Approximate mode requires explicit selection. New exhaustive runs are not claimed equivalent to old ANN candidate sets."},
    {"name": "Uncapped library condition counts", "status": "full_count adapter; workstation validation pending", "scope": "Evaluate every conformer pose against classified features without descriptor Top-K or Gaussian Top-N. Long batch job, not all possible orientations or torsions. Reuse completed classification as query definition."},
    {"name": "Necessary-condition funnel", "status": "rule-based adapter; workstation speedup pending", "scope": "Use an explicit selection rule across all library conformers; reject only impossible feature type/assignment/pair geometry, then refine every survivor without Top-K. Per-feature counts are conditional on that rule."},
    {"name": "Screening evidence and human selection", "status": "available", "scope": "Review a completed search, preview explicit same-pose anchor conditions, then separately request SDF/ID export for docking preparation."},
    {"name": "New-target query preparation", "status": "structure_survey and adopted guided designs", "scope": "Requires mapped target-bound reference coordinates; exploratory thresholds are not independent retrieval calibration."},
    {"name": "Molecule aggregation and pocket QC", "status": "CLI components; chat integration pending", "scope": "Requires validated query poses, receptor/site definitions and aggregation artifacts."},
    {"name": "Flexible docking", "status": "Glide preparation/execution adapter; workstation validation pending", "scope": "Derive pocket from crystal ligand, prepare inputs, explicitly run licensed preparation and reference-gated docking. PLANTS execution and cross-docking validation remain pending."},
    {"name": "Biological enrichment and final selection", "status": "not validated", "scope": "Requires held-out active/decoy labels, diversity/property criteria and experimental evidence."},
]
ROUTER = """You are an AIDD conversational task coordinator. Return JSON only with exactly
intent, message, task_id, request, plus a selection object for select or a design object for design.
intent is run/status/results/resume/cancel/capabilities/clarify/evidence/classify/full_count/funnel/benchmark/coarse/select/export/prepare_docking/run_docking/recommend/adopt/design/guided.
For a pre-search protein/ligand PDB survey use run requesting protein verification and structure_survey.
For a resolution census and diverse ligand reference comparison use run requesting
structure_diversity, preserving the user's reference PDB. Do not claim polymer sequence
similarity is chemical similarity or that multi-reference library execution is available.
Use recommend on a completed structure_survey task for grounded advice. Use adopt when the user
explicitly accepts that proposal or delegates to it. Silence never means adoption.
Use design to edit an anchor_recommend task: include an additional design object with only
query_id, mandatory_anchors, alternative_groups, optional_anchors, evidence_ids or rationale.
Use exact IDs in evidence. Mandatory anchors are ALL, each alternative group is ANY in the same
pose, optional anchors are labels only. Do not invent coordinates or scientific evidence.
Use guided on an adopted anchor_design to run the full library with feature geometry and pocket
exclusion before Gaussian, preserving pose combinations and all scaffold-group members.
Use select on guided_funnel results to retain requested anchor combinations at the original
threshold. This keeps the adopted mandatory rules and exports pose/ID lists, not docking inputs.
These are separate chat tasks; progress questions must not launch them. If no usable crystal
coordinates exist, show literature and request the needed target-bound structure, not a made-up query.
Write message in English. task_id is an existing task ID from this session or null (latest).
request is a self-contained scientific request only for run, otherwise empty.
Use conversation context to resolve follow-up clarifications. Never invent biological sequences,
paths, results, unsupported tools or identifiers. Questions about progress/results NEVER run a new
job. Resume/cancel only when explicitly requested. For run retain all explicit user constraints.
Supported actions are supplied. Unsupported stages must explain missing capability using clarify,
not substitute WEE1 for another target. AF3 prediction cannot automatically create a calibrated
search query. Missing organism/construct/reference information requires clarification. Do not
claim a task is queued or complete: the local executor reports its actual status.
Use evidence to classify anchors and counts from a completed search task; no new search.
Use classify on a completed evidence task to derive broader crystal-based query features
(hydrophobic, aromatic, charge, water and metal hypotheses) and match stored 3D poses.
This requires no PLIP or docking. Halogen-specific features are not indexed.
Use coarse on a completed joint selection preview to audit 10000 spread IDs plus
saved boundary/positive rows without seeds or Gaussian. This is not a full-library count.
Use benchmark on a completed selection-preview task to test CPU and available GPU
performance and equivalence on a small sample before a full-library run.
Use funnel on a completed selection-preview task to apply its explicit same-pose
rule to the WHOLE library, with necessary-condition rejection before pose scoring.
The old preview candidate IDs and molecule cap do not constrain the full-library run.
Never invent a threshold or make every contact in a classification mandatory.
If no rule exists, explain that unconstrained individual-feature counts and a
combined feature rule answer different questions; do not launch full_count as a
substitute for a fast conditional funnel. full_count is an expensive baseline.
Use full_count on a completed classification task when the user wants every library
conformer evaluated for feature conditions, without a Top-K or Top-N membership limit.
This computes all feature scores first and reports existing diagnostic levels 0.25/0.5/0.75;
do not demand an extra threshold choice. It is a long pose-computation batch job.
Do not substitute search_3d exhaustive descriptor Top-K for this request.
Use select on its completed results for subsequent same-pose conditions.
Use prepare_docking after a completed export to derive the pocket from the crystal
ligand and generate reviewed Glide preparation/grid/redocking commands. No prepared
grid is needed in advance. Use run_docking only on an existing preparation task when
the user explicitly requests execution. This runs PrepWizard/LigPrep, grid generation,
reference redocking, then candidate docking only if the reference RMSD gate passes.
Do not run these during 3D feature classification. Local licenses/version compatibility
and preparation chemistry must be reviewed; PLANTS execution is not supported yet.
Use select for a selection preview from a completed evidence task. selection must contain
required_anchors (exact full IDs from that task), match_mode (all or any), minimum_score
(explicit user threshold in (0,1]), and optional max_molecules (explicit user cap).
Optional coarse_constraints requires all three explicit user fields: heavy_atom_ratio
([lower, upper] inclusive candidate/query heavy atom count ratio), maximum_extent_distance
(nonnegative normalized principal-extent distance), minimum_feature_coverage ((0,1]
query typed-feature count coverage). These add joint eligibility BEFORE pose generation.
Never invent these thresholds or claim a rejection percentage. Pocket-derived feature
IDs may be selected explicitly; do not claim receptor clashes are evaluated.
Never invent anchors, thresholds, caps or interaction evidence. Clarify missing conditions.
Select anchors from ONE query only per preview; all means the same pose satisfies all.
Use export only on a completed selection preview after a NEW explicit user request to export
or confirm that preview. This exports poses/IDs, does not run docking. Never combine selection
and export or run docking automatically. request must be empty for these three intents.
To inspect existing evidence or preview counts use results, avoiding duplicate review jobs.
JSON example: {"intent":"status","message":"Checking task status.","task_id":null,"request":""}.
"""


def validate_route(value):
    if not isinstance(value, dict) or set(value) not in ({"intent", "message", "task_id", "request"}, {"intent", "message", "task_id", "request", "selection"}, {"intent", "message", "task_id", "request", "design"}):
        raise ValueError("Invalid chat decision")
    if value["intent"] not in {"run", "status", "results", "resume", "cancel", "capabilities", "clarify", "evidence", "classify", "full_count", "funnel", "benchmark", "coarse", "select", "export", "prepare_docking", "run_docking", "recommend", "adopt", "design", "guided"}:
        raise ValueError("Unsupported chat intent")
    for field in ("message", "request"):
        if not isinstance(value[field], str) or len(value[field]) > 12000:
            raise ValueError("Invalid chat text")
    if contains_han(value["message"]): raise ValueError("Agent replies must be English")
    if value["task_id"] is not None and not isinstance(value["task_id"], str):
        raise ValueError("Invalid task reference")
    if (value["intent"] == "run") != bool(value["request"].strip()):
        raise ValueError("Execution requires a nonempty request; queries cannot contain one")
    if (value["intent"] == "select") != ("selection" in value):
        raise ValueError("Only select accepts a selection policy")
    if "selection" in value:
        from .screening_selection import validate_selection
        validate_selection(value["selection"])
    if (value['intent']=='design') != ('design' in value):raise ValueError('Only design accepts design edits')
    if 'design' in value and (not isinstance(value['design'],dict) or set(value['design'])-{'query_id','mandatory_anchors','alternative_groups','optional_anchors','evidence_ids','rationale'}):
        raise ValueError('Invalid design edits')
    return value


def read_json(path):
    try: return json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, ValueError): return None


def screening_summary(job):
    """Bounded actual evidence for conversation context and deterministic feedback."""
    from .project_context import ensure_within
    report = read_json(job["report"]) if job.get("report") else None
    if not report or not job.get("plan"): return {}
    for step in report.get("steps", {}).values():
        if step.get('action')=='structure_diversity' and step.get('status')=='complete':
            child=read_json(ensure_within(Path(step['result']['report']),Path(job['plan']).parent)) or {}
            summary={k:child[k] for k in ('kind','status','target','discovered_entries','census','quality_site_unique_ligands',
                'coverage_curve','outputs','failed_entries','limitations') if k in child}
            summary['proposed_references']=[{k:r[k] for k in ('query_id','resolution','heavy_atoms','r_free')}
                for r in child.get('proposed_references',[])]
            polymers=child.get('polymer_ligands',{})
            summary['polymer_ligands']=dict(unique_sequence_link_variants=polymers.get('unique_sequence_link_variants'),
                examples=[{k:r[k] for k in ('query_id','resolution','length','description','reference_readiness')}
                          for r in polymers.get('sequence_diversity_examples',[])],limitations=polymers.get('limitations',[]))
            return summary
        if step.get('action')=='guided_select' and step.get('status')=='complete':
            child=read_json(ensure_within(Path(step['result']['report']),Path(job['plan']).parent)) or {}
            return {k:child[k] for k in ('kind','status','matching_molecules','pose_records','selection','outputs','scope') if k in child}
        if step.get("status") != "complete" or step.get("action") not in {"review_screening", "classify_screening", "full_library_screen", "condition_funnel", "benchmark_funnel", "select_screening", "export_screening", "prepare_docking", "run_docking", "structure_survey", "anchor_recommend", "anchor_design", "guided_funnel"}: continue
        path = ensure_within(Path(step["result"]["report"]), Path(job["plan"]).parent)
        child = read_json(path) or {}
        if child.get('kind') in {'structure_survey','anchor_recommendation','anchor_design','guided_funnel'}:
            keys=('kind','status','readiness','target','discovered_structures','unanalyzed_structures',
                  'recurrence','recommendation','defaults','design','needs_input','literature','literature_advice',
                  'reference_control','stage_counts','counts','aggregation','outputs','failures')
            summary={k:child[k] for k in keys if k in child}
            summary['queries']=[dict(query_id=q['query_id'],anchors=[{k:a[k] for k in
                ('anchor_id','feature_class','target_residue') if k in a} for a in q['anchors']],
                evidence_id=q.get('evidence_id')) for q in child.get('queries',[])]
            return summary
        if child.get('kind') == 'coarse_audit':
            return {k:v for k,v in child.items() if k not in {'sources','code','sample_ids','sample_panel','scenarios'}}
        if child.get('kind') == 'hardware_benchmark':
            return {k:child[k] for k in ('kind','status','sample_count','sample_scope','scenarios','gpu','recommendation','sample_panel','reference_positive_count','positive_membership_validation','scalability_gate','scalability_scope') if k in child}
        if child.get('kind') in {'docking_preparation','docking_execution'}:
            return {k:child[k] for k in ('kind','status','engine','geometry','preparation_status','policy','reference_validation','pose_quality') if k in child}
        summary = {k:child[k] for k in ("kind", "library", "refined_molecule_union", "retrieved_molecule_union", "retrieved_conformer_union", "counts", "policy", "approval", "sdf", "ids", "preparation", "docking_status") if k in child}
        if child.get("kind") in ("full_library_conditions","condition_funnel"):
            summary.update(status=child['status'], interaction_table=str(path.parent/'interactions.html'),
                queries=[dict(query_id=q['query_id'], counts=q['counts'], full_coverage=q['full_coverage'],
                    anchors=[{k:a[k] for k in ('anchor_id','feature_class','diagnostic_counts')} for a in q['anchors']]) for q in child['queries']])
            if child.get('condition_policy'): summary['condition_policy'] = child['condition_policy']
        if child.get("kind") == "screening_evidence":
            if child.get('classification'):
                summary['classification']=child['classification'];summary['interaction_table']=str(path.parent/'interactions.html')
                summary['class_counts']={k:len(v) for k,v in child['classes'].items()}
                summary['unsupported_classes']=child['unsupported_classes']
            summary["queries"] = [dict(query_id=q["query_id"], counts=q["counts"], retrieval=q.get('retrieval',{}), anchors=[
                {k:a[k] for k in ("anchor_id", "feature_class", "ligand_atom_indices", "evidence_class", "evidence", "diagnostic_counts")} for a in q["anchors"]]) for q in child["queries"]]
            summary["interpretation"] = "Feature-match hypotheses, not validated candidate hydrogen bonds; diagnostic score thresholds are not probabilities or accepted cutoffs."
        return summary
    return {}


def screening_feedback(summary):
    if summary.get('kind')=='structure_diversity':
        return json.dumps(summary,indent=2)+'\nReview the reference proposal. Polymer chemical preparation and multi-reference execution require further implementation.'
    if summary.get('kind')=='guided_selection':return json.dumps(summary,indent=2)
    if summary.get('kind') in {'structure_survey','anchor_recommendation','anchor_design','guided_funnel'}:
        kind=summary['kind']
        next_step={'structure_survey':'Use /recommend for evidence-grounded LLM advice.',
            'anchor_recommendation':'Edit the proposed anchor IDs in chat, or use /adopt to accept the proposal.',
            'anchor_design':'Use /guided to start the full-library scan with this sealed design.',
            'guided_funnel':'Review the pose combinations and scaffold groups in the reported outputs.'}[kind]
        if summary.get('readiness')=='needs_structure_input':
            next_step='Provide a target-bound PDB ID in chat for a new survey. Literature alone does not provide 3D coordinates. Use /recommend to summarize available literature.'
        return json.dumps(summary,indent=2)+'\n'+next_step
    lines = [summary["kind"].replace("_", " ").capitalize()]
    if summary['kind'] in {'coarse_audit','hardware_benchmark','docking_preparation','docking_execution'}:return json.dumps(summary,indent=2)
    if summary.get("library"):
        lines.append(f"Library: {summary['library']['library_conformers']:,} conformers / {summary['library']['library_molecules']:,} source-grouped molecules.")
    if summary['kind'] in ('full_library_conditions','condition_funnel'):
        lines.append('No descriptor Top-K or Gaussian Top-N membership truncation. Necessary-condition rejects may have no calculated pose; survivors use heuristic rigid poses.')
        if summary.get('condition_policy'): lines.append('Selection rule: '+json.dumps(summary['condition_policy']))
        for q in summary['queries']:
            c=q['counts']
            lines.append(f"{q['query_id']}: {c['evaluated_conformers']:,}/{c['total_conformers']:,} conformers evaluated; {c['evaluated_molecules']:,} molecules; full coverage: {q['full_coverage']}")
            if summary['kind']=='condition_funnel':
                lines.append('Funnel stages: '+json.dumps(c))
                lines.append('Feature counts below are conditional on passing the selected rule.')
            for a in q['anchors']:
                lines.append(a['anchor_id']+' | '+ '; '.join(f"score >= {d['minimum_score']}: {d['molecules']} molecules" for d in a['diagnostic_counts']))
        lines.append('Classification table: '+summary['interaction_table'])
        lines.append('Preview all/any conditions on this task, then explicitly export. Diagnostic thresholds are exploratory.')
        return '\n'.join(lines)
    for q in summary.get("queries", []):
        c = q["counts"]
        lines.append(f"{q['query_id']}: retrieved {c['retrieved_conformers']:,} conformers / {c['retrieved_molecules']:,} molecules; refined {c['refined_conformers']:,} conformers / {c['refined_molecules']:,} molecules.")
        retrieval=q.get('retrieval',{})
        if retrieval.get('mode')=='exhaustive':lines.append(f"Exhaustive coarse coverage: {retrieval['compared_conformers']:,}/{retrieval['total_conformers']:,} conformers compared.")
        else:lines.append('This source search does not establish exhaustive descriptor coverage.')
        for a in q["anchors"]:
            evd = a["evidence"]
            partner = evd.get("protein_partner", {})
            lines.append(f"Anchor {a['anchor_id']} | {a['feature_class']} | ligand atom indices {a['ligand_atom_indices']}")
            lines.append("Crystal partner: " + json.dumps(partner) + f"; distance {evd.get('distance_angstrom')} A; angle {evd.get('angle_degrees')}; {evd.get('angle_status', 'angle status unreported')}.")
            lines.append("Diagnostic molecule counts: " + "; ".join(f"score >= {d['minimum_score']}: {d['molecules']}" for d in a["diagnostic_counts"]))
    if "refined_molecule_union" in summary:
        lines.append(f"Across queries: {summary['retrieved_molecule_union']:,} retrieved / {summary['refined_molecule_union']:,} refined distinct source-grouped molecules. Do not sum per-query molecule counts.")
        lines.append("HBA = ligand acceptor; HBD = ligand donor. These are crystal-derived feature hypotheses, not verified candidate-protein contacts.")
        if summary.get('classification'):
            lines.append('Class counts: '+json.dumps(summary['class_counts']))
            lines.append('Classification table: '+summary['interaction_table'])
            lines.append('Unassessed classes: '+', '.join(summary['unsupported_classes']))
        else: lines.append('Other interaction classes are unassessed. Use /classify REVIEW_TASK_ID for broader crystal-derived feature hypotheses.')
        lines.append("Next: specify exact anchor IDs from one query, all/any mode and a minimum score in (0,1]. Diagnostic thresholds are illustrative, not validated cutoffs. Request a preview first.")
    if "counts" in summary:
        lines.extend(f"{k}: {v:,}" for k,v in summary["counts"].items())
        lines.append("Selection policy: " + json.dumps(summary["policy"]))
    for key in ("approval", "sdf", "ids", "preparation", "docking_status"):
        if key in summary: lines.append(key + ": " + str(summary[key]))
    if summary["kind"] == "selection_preview":
        lines.append("To export this exact preview, send a separate confirmation or /export TASK_ID. Original search ranks and E031 annotations remain unchanged.")
    return "\n".join(lines)


class ChatAgent:
    def __init__(self, storage, runtime=None, config_dir=None, allow_compute=False, *, planner=create_plan, router=None, start=True):
        self.root = Path(storage).resolve() / "chat"
        self.root.mkdir(parents=True, exist_ok=True)
        self.db = self.root / "conversations.sqlite3"
        self.context = initialize_context(Path(storage))
        self.runtime = Path(runtime).resolve() if runtime else None
        self.config_dir = config_dir
        self.allow_compute = allow_compute
        self.planner, self.router = planner, router
        self.stop = threading.Event()
        self.lock = threading.RLock()
        self.process = None
        self.worker = None
        with self.connect() as db:
            db.executescript("""CREATE TABLE IF NOT EXISTS sessions(id TEXT PRIMARY KEY,title TEXT,created REAL);
            CREATE TABLE IF NOT EXISTS messages(id INTEGER PRIMARY KEY,session TEXT,role TEXT,text TEXT,created REAL);
            CREATE TABLE IF NOT EXISTS jobs(id TEXT PRIMARY KEY,session TEXT,request TEXT,provider TEXT,
            profile TEXT,status TEXT,plan TEXT,report TEXT,log TEXT,error TEXT,cancel INTEGER DEFAULT 0,created REAL);""")
            db.execute("UPDATE jobs SET status='interrupted',error='Server stopped during execution; explicitly resume this task.' WHERE status IN ('planning','running')")
        if start:
            self.worker = threading.Thread(target=self.work, daemon=True)
            self.worker.start()

    @contextmanager
    def connect(self):
        db = sqlite3.connect(self.db, timeout=30)
        db.row_factory = sqlite3.Row
        try:
            with db: yield db
        finally:
            db.close()

    def session(self, sid):
        with self.connect() as db:
            if db.execute("SELECT 1 FROM sessions WHERE id=?", (sid,)).fetchone() is None:
                raise ValueError("Unknown session")

    def new_session(self):
        sid = uuid.uuid4().hex
        with self.connect() as db:
            db.execute("INSERT INTO sessions VALUES(?,?,?)", (sid, "New conversation", time.time()))
        return sid

    def message(self, sid, role, text):
        with self.connect() as db:
            db.execute("INSERT INTO messages(session,role,text,created) VALUES(?,?,?,?)", (sid, role, text, time.time()))

    def jobs(self, sid):
        with self.connect() as db:
            return [dict(r) for r in db.execute("SELECT * FROM jobs WHERE session=? ORDER BY created", (sid,))]

    def task(self, sid, jid=None):
        jobs = self.jobs(sid)
        if jid is None and jobs: return jobs[-1]
        for job in jobs:
            if job["id"] == jid: return job
        raise ValueError("No matching task in this conversation")

    def attach(self, sid, plan_path):
        """Attach an existing owned Project run for inspection without rerunning it."""
        from .prompt_workflow import project_root
        from .project_context import ensure_within
        from .library_acceptance import sha
        self.session(sid)
        ctx = self.context
        path = ensure_within(Path(plan_path).resolve(), project_root(Path(ctx["db"]),ctx["user_id"],ctx["project_id"]))
        envelope = read_json(path)
        seal = read_json(path.parent / "plan-seal.json")
        if not envelope or envelope.get("user") != ctx["user_id"] or envelope.get("project") != ctx["project_id"]:
            raise ValueError("Existing plan must belong to the active Project")
        if seal != dict(plan_sha256=sha(path),user=ctx["user_id"],project=ctx["project_id"]):
            raise ValueError("Existing plan seal mismatch")
        report_path = path.parent / "execution/report.json"
        report = read_json(report_path)
        marker = read_json(report_path.parent / "RUN_STATUS.json")
        status = "interrupted"
        if report and marker:
            if marker.get("report_sha256") != sha(report_path): raise ValueError("Existing report receipt mismatch")
            status = report.get("status")
            if status not in {"complete","failed","blocked"}: status = "interrupted"
        jid = uuid.uuid4().hex[:16]
        with self.connect() as db:
            db.execute("INSERT INTO jobs(id,session,request,provider,profile,status,plan,report,created) VALUES(?,?,?,?,?,?,?,?,?)",
                       (jid,sid,envelope["plan"]["summary"],"imported","",status,str(path),str(report_path),time.time()))
        self.message(sid,"assistant",f"Attached existing task {jid}: {status}. No computation was launched. Report: {report_path}")
        return jid

    def snapshot(self, sid=None):
        with self.connect() as db:
            sessions = [dict(r) for r in db.execute("SELECT * FROM sessions ORDER BY created DESC")]
            messages = []
            if sid:
                self.session(sid)
                messages = [dict(r) for r in db.execute("SELECT role,text FROM messages WHERE session=? ORDER BY id", (sid,))]
        tasks = []
        for job in self.jobs(sid) if sid else []:
            report = read_json(job["report"]) if job["report"] else None
            tail = ""
            if job["log"]:
                try:
                    with Path(job["log"]).open("rb") as f:
                        f.seek(0, 2); f.seek(max(0, f.tell()-6000)); tail = f.read().decode("utf-8", "replace")
                except OSError: pass
            if job["plan"]:
                execution = Path(job["plan"]).parent / "execution"
                for child_log in sorted(execution.glob("*/execution.log"))[-3:]:
                    from .project_context import ensure_within
                    try:
                        ensure_within(child_log, execution)
                        with child_log.open("rb") as f:
                            f.seek(0, 2); f.seek(max(0, f.tell()-4000))
                            tail += "\n"+str(child_log)+"\n"+f.read().decode("utf-8", "replace")
                    except (OSError, ValueError): pass
            tasks.append(dict(id=job["id"], status=job["status"], provider=job["provider"], request=job["request"],
                              plan=job["plan"], report=job["report"], log=job["log"], error=job["error"],
                              details=report, progress=tail))
        return dict(sessions=sessions, messages=messages, tasks=tasks, workflows=WORKFLOWS,
                    compute_enabled=self.allow_compute)

    def ask(self, sid, text, provider):
        self.session(sid)
        if provider not in {"gpt", "deepseek"}: raise ValueError("Unknown provider")
        if not isinstance(text, str) or not 1 <= len(text.strip()) <= 8000: raise ValueError("Message must contain 1-8000 characters")
        # Serialize decisions, not scientific execution. Session messages stay ordered.
        with self.lock:
            self.message(sid, "user", text)
            with self.connect() as db:
                db.execute("UPDATE sessions SET title=? WHERE id=? AND title='New conversation'", (text[:70], sid))
            parts = text.strip().split()
            command = parts[0].lower()
            if command in {"/status", "/results", "/capabilities", "/resume", "/cancel", "/evidence", "/classify", "/full_count", "/funnel", "/benchmark", "/coarse", "/export", "/prepare_docking", "/run_docking", "/recommend", "/adopt", "/guided"}:
                if len(parts) > 2: raise ValueError("Use /command followed by an optional task ID")
                decision = dict(intent=command[1:], message="", task_id=parts[1] if len(parts)==2 else None, request="")
            else:
                profile = select_llm_profile(provider, config_dir=self.config_dir)
                snapshot = self.snapshot(sid)
                context = dict(conversation=snapshot["messages"][-12:], tasks=[{k: j[k] for k in ("id","status","request")} for j in snapshot["tasks"][-10:]])
                # Bounded dialogue context; no files, sequences or execution logs sent to provider.
                for item in context["conversation"]: item["text"] = item["text"][:1000]
                for item in context["tasks"]: item["request"] = item["request"][:300]
                # Scientific summaries are deliberately sent for anchor-aware follow-ups;
                # no library structures, raw arrays, logs or credentials are uploaded.
                for item in context["tasks"][-3:]:
                    item["screening_evidence"] = screening_summary(self.task(sid, item["id"]))
                    for q in item['screening_evidence'].get('queries',[]):
                        q['total_anchors']=len(q['anchors'])
                        q['anchors']=[{k:a[k] for k in ('anchor_id','feature_class')} for a in q['anchors'][:64]]
                        q['context_note']='At most 64 anchor IDs per query shown; ask for exact IDs from the full table if missing.'
                context["latest_message"] = text
                prompt = json.dumps(context, ensure_ascii=False)
                if len(prompt) > 20000: context["conversation"] = context["conversation"][-4:]; prompt = json.dumps(context, ensure_ascii=False)
                while len(prompt) > 19000 and context["tasks"]:
                    context["tasks"].pop(0)
                    prompt = json.dumps(context, ensure_ascii=False)
                decision = (self.router(prompt, profile) if self.router else chat_plan(prompt, read_json(profile),
                    system_prompt=ROUTER, capabilities=dict(actions=CAPABILITIES, workflows=WORKFLOWS), validator=validate_route)[0])
                validate_route(decision)
            intent = decision["intent"]
            if intent in {'recommend','adopt','design','guided'}:
                expected={'recommend':'structure_survey','adopt':'anchor_recommend','design':'anchor_recommend','guided':'anchor_design'}[intent]
                if decision['task_id'] is None:
                    candidates=[j for j in self.jobs(sid) if j['status']=='complete' and j.get('report') and
                        any(s.get('action')==expected and s.get('status')=='complete' for s in
                            (read_json(j['report']) or {}).get('steps',{}).values())]
                    if not candidates:raise ValueError('Complete '+expected+' in this conversation first')
                    job=candidates[-1]
                else:job=self.task(sid,decision['task_id'])
                if job['status']!='complete' or not job.get('plan'):raise ValueError('Choose a completed source task')
                details=read_json(job['report']) or {}
                if sum(s.get('action')==expected and s.get('status')=='complete' for s in details.get('steps',{}).values())!=1:
                    raise ValueError('Wrong source task for '+intent)
                action={'recommend':'anchor_recommend','adopt':'anchor_design','design':'anchor_design','guided':'guided_funnel'}[intent]
                params=dict(source_run=Path(job['plan']).parent.name)
                if intent=='recommend':params['provider']=provider
                if intent=='design':params['design']=decision['design']
                plan=dict(version=1,summary='Structure-guided workflow: '+intent,clarifications=[],
                          steps=[dict(id='guided',action=action,params=params)])
                ctx=self.context
                path=create_plan(Path(ctx['db']),ctx['user_id'],ctx['project_id'],text,local_plan=plan)
                jid=uuid.uuid4().hex[:16]
                with self.connect() as db:
                    db.execute('INSERT INTO jobs(id,session,request,provider,profile,status,plan,report,log,created) VALUES(?,?,?,?,?,?,?,?,?,?)',
                        (jid,sid,plan['summary'],provider,'','queued',str(path),str(path.parent/'execution/report.json'),str(self.root/(jid+'.log')),time.time()))
                answer=f"Task {jid} queued: {intent}. Source task: {job['id']}."
                self.message(sid,'assistant',answer)
                return answer
            if intent in {"evidence", "classify", "full_count", "funnel", "benchmark", "coarse", "select", "export", "prepare_docking", "run_docking"}:
                if intent in {'benchmark','funnel','coarse'} and decision['task_id'] is None:
                    candidates = [j for j in self.jobs(sid) if j['status']=='complete' and j.get('report') and
                        any(s.get('action')=='select_screening' and s.get('status')=='complete'
                            for s in (read_json(j.get('report')) or {}).get('steps',{}).values())]
                    if not candidates: raise ValueError('Create a selection preview in this conversation first')
                    job = candidates[-1]
                else:
                    job = self.task(sid, decision["task_id"])
                if job["status"] != "complete" or not job["plan"]:
                    raise ValueError("Choose a completed source task from this conversation")
                expected = {"evidence":{"search_3d"}, "classify":{"review_screening"}, "full_count":{"classify_screening"}, "funnel":{"select_screening"}, "coarse":{"select_screening"}, "benchmark":{"select_screening"}, "select":{"review_screening","classify_screening","full_library_screen","condition_funnel","guided_funnel"}, "export":{"select_screening"}, "prepare_docking":{"export_screening"}, "run_docking":{"prepare_docking"}}[intent]
                details = read_json(job["report"]) or {}
                if sum(s.get("action") in expected and s.get("status") == "complete" for s in details.get("steps", {}).values()) != 1:
                    raise ValueError("Choose a completed " + '/'.join(sorted(expected)) + " task")
                action = {"evidence":"review_screening", "classify":"classify_screening", "full_count":"full_library_screen", "funnel":"condition_funnel", "coarse":"benchmark_funnel", "benchmark":"benchmark_funnel", "select":"select_screening", "export":"export_screening", "prepare_docking":"prepare_docking", "run_docking":"run_docking"}[intent]
                params = dict(source_run=Path(job["plan"]).parent.name)
                if intent == "select": params.update(decision["selection"])
                if intent == 'select' and any(s.get('action')=='guided_funnel' for s in details.get('steps',{}).values()):
                    action='guided_select'
                if intent == "coarse": params["coarse_only"] = True
                local_plan = dict(version=1, summary={"evidence":"Review crystal anchors and search counts", "classify":"Classify crystal-derived 3D feature hypotheses", "full_count":"Evaluate all library conformers against classified features without Top-K or Top-N", "funnel":"Check all library conformers with necessary conditions before uncapped precise matching", "coarse":"Audit joint coarse selectivity without generating poses", "benchmark":"Compare equivalent CPU and available CUDA pose scoring on a bounded spread-out pilot",
                    "select":"Preview explicit same-pose anchor selection", "export":"Export the confirmed selection for docking preparation", "prepare_docking":"Prepare crystal-derived pocket and Glide workflow", "run_docking":"Execute prepared Glide validation and gated docking"}[intent],
                    clarifications=[], steps=[dict(id="screening", action=action, params=params)])
                ctx = self.context
                plan = create_plan(Path(ctx["db"]), ctx["user_id"], ctx["project_id"], text, local_plan=local_plan)
                jid = uuid.uuid4().hex[:16]
                with self.connect() as db:
                    db.execute("INSERT INTO jobs(id,session,request,provider,profile,status,plan,report,log,created) VALUES(?,?,?,?,?,?,?,?,?,?)",
                        (jid,sid,local_plan["summary"],provider,"","queued",str(plan),str(plan.parent / "execution/report.json"),str(self.root / (jid+".log")),time.time()))
                answer = f"Task {jid} queued: {local_plan['summary']}. Source task: {job['id']}. No new search or docking job is launched."
                if intent=='benchmark':answer=f"Task {jid} queued for a bounded hardware/equivalence pilot, not full-library screening. Source task: {job['id']}."
                if intent=='funnel':answer=f"Task {jid} queued for a full-library necessary-condition funnel using the selected rule, with no Top-K or Top-N cap. Source task: {job['id']}."
                if intent=='full_count':answer=f"Task {jid} queued for full-library pose and feature evaluation with no Top-K/Top-N membership limit. This can be a long batch job. Source task: {job['id']}."
                if intent=='run_docking':answer=f"Task {jid} queued for preparation, reference redocking and gated candidate docking. Source task: {job['id']}."
            elif intent == "run":
                jid = uuid.uuid4().hex[:16]
                log = str(self.root / (jid + ".log"))
                with self.connect() as db:
                    db.execute("INSERT INTO jobs(id,session,request,provider,profile,status,log,created) VALUES(?,?,?,?,?,?,?,?)",
                               (jid,sid,decision["request"],provider,str(profile),"queued",log,time.time()))
                answer = f"Task {jid} queued. I will show its plan, progress and result paths here."
            elif intent in {"status", "results", "resume", "cancel"}:
                job = self.task(sid, decision["task_id"])
                if intent == "cancel":
                    if job["status"] not in {"queued", "planning", "running"}:
                        raise ValueError("Task is not active; there is no computation to cancel")
                    with self.connect() as db:
                        db.execute("UPDATE jobs SET cancel=1 WHERE id=?", (job["id"],))
                    answer = f"Cancellation requested for {job['id']}. Existing artifacts will be retained."
                elif intent == "resume":
                    if job["status"] not in {"failed","blocked","interrupted","cancelled","complete"}:
                        raise ValueError("Task is already queued or active")
                    with self.connect() as db:
                        db.execute("UPDATE jobs SET status='queued',cancel=0,error=NULL WHERE id=?", (job["id"],))
                    answer = f"Task {job['id']} queued for resume using its existing plan and receipts."
                else:
                    view = next(t for t in self.snapshot(sid)["tasks"] if t["id"] == job["id"])
                    lines = [f"Task {job['id']}: {job['status']}."]
                    if job["error"]: lines.append(job["error"])
                    details = view["details"] or {}
                    for name, step in details.get("steps", {}).items():
                        lines.append(f"{name}: {step.get('status', 'unknown')}")
                        if step.get("error"): lines.append(step["error"])
                        if intent == "results":
                            for key, value in step.get("result", {}).items():
                                if key in {"model", "report", "source_url", "status", "accession", "pose_quality", "biological_quality"}:
                                    lines.append(f"{key}: {value}")
                            result = step.get("result", {})
                            confidence = result.get("confidence", {})
                            for key in ("ptm", "ranking_score", "fraction_disordered"):
                                if key in confidence: lines.append(f"{key}: {confidence[key]}")
                            nested = result.get("report")
                            if nested and job["plan"]:
                                from .project_context import ensure_within
                                try:
                                    nested_path = ensure_within(Path(nested).resolve(), Path(job["plan"]).parent)
                                    child = read_json(nested_path) or {}
                                    for query in child.get("queries", []):
                                        if "query_id" in query:
                                            lines.append(f"{query['query_id']}: execution {query.get('execution_seconds', 'unreported')} s; Gaussian {query.get('gaussian_seconds', 'unreported')} s; E031 {query.get('annotation_seconds', 'unreported')} s")
                                except ValueError: lines.append("Nested report is outside this task; not opened.")
                    if job["report"]: lines.append("Report: "+job["report"])
                    if job["plan"]: lines.append("Plan: "+job["plan"])
                    if job["log"]: lines.append("Log: "+job["log"])
                    if not details: lines.append("No completed scientific report is available yet.")
                    if intent == "results":
                        summary = screening_summary(job)
                        if summary: lines.append(screening_feedback(summary))
                    answer = "\n".join(lines)
            elif intent == "capabilities": answer = "\n\n".join(w["name"]+": "+w["status"]+". "+w["scope"] for w in WORKFLOWS)
            else: answer = decision["message"]
            self.message(sid, "assistant", answer)
            return answer

    def update(self, jid, **fields):
        with self.connect() as db:
            db.execute("UPDATE jobs SET "+",".join(k+"=?" for k in fields)+" WHERE id=?", (*fields.values(),jid))

    def work(self):
        while not self.stop.wait(.25):
            with self.connect() as db:
                job = db.execute("SELECT * FROM jobs WHERE status='queued' ORDER BY created LIMIT 1").fetchone()
            if job:
                try: self.execute(dict(job))
                except Exception as exc:
                    # Exception bodies from remote services must not expose credentials.
                    self.update(job["id"], status="failed", error=str(exc) if isinstance(exc, ValueError) else type(exc).__name__)
                    self.message(job["session"], "assistant", f"Task {job['id']} failed. Check its task panel and local execution log.")

    @staticmethod
    def terminate(process):
        if process.poll() is not None: return
        if os.name == "posix":
            try: os.killpg(process.pid, signal.SIGTERM)
            except ProcessLookupError: return
        else: subprocess.run(["taskkill", "/PID", str(process.pid), "/T", "/F"], capture_output=True)
        try: process.wait(timeout=5)
        except subprocess.TimeoutExpired:
            if os.name == "posix": os.killpg(process.pid, signal.SIGKILL)
            else: process.kill()
            process.wait(timeout=5)

    def execute(self, job):
        jid, sid = job["id"], job["session"]
        if job["cancel"]: self.update(jid,status="cancelled"); return
        self.update(jid,status="planning")
        ctx = self.context
        if not job["plan"]:
            plan = self.planner(Path(ctx["db"]),ctx["user_id"],ctx["project_id"],job["request"],llm_profile=Path(job["profile"]))
            job["plan"] = str(plan)
            job["report"] = str(plan.parent / "execution/report.json")
            self.update(jid,plan=job["plan"],report=job["report"])
            envelope = read_json(plan)
            self.message(sid,"assistant",f"Task {jid}: "+envelope["plan"]["summary"])
        if self.task(sid,jid)["cancel"] or self.stop.is_set():
            self.update(jid,status="interrupted" if self.stop.is_set() else "cancelled"); return
        cmd = [sys.executable,"-u","-m","aidd_agent.prompt_workflow","run","--db",ctx["db"],"--user",ctx["user_id"],
               "--project",ctx["project_id"],"--plan",job["plan"]]
        if self.runtime: cmd += ["--runtime",str(self.runtime)]
        if self.allow_compute: cmd += ["--allow-compute"]
        options = {"start_new_session":True} if os.name == "posix" else {"creationflags":subprocess.CREATE_NO_WINDOW}
        if self.config_dir:
            options['env']=dict(os.environ,AIDD_LLM_CONFIG_DIR=str(self.config_dir))
        self.update(jid,status="running")
        with Path(job["log"]).open("a",encoding="utf-8") as log:
            with subprocess.Popen(cmd,stdout=log,stderr=subprocess.STDOUT,**options) as process:
                self.process = process
                while process.poll() is None:
                    if self.stop.wait(.25) or self.task(sid,jid)["cancel"]:
                        self.terminate(process)
                        self.update(jid,status="interrupted" if self.stop.is_set() else "cancelled")
                        self.process = None
                        return
                self.process = None
                report = read_json(job["report"])
                status = report.get("status") if report else "failed"
                if status not in {"complete","blocked","failed"} or (process.returncode != 0 and status == "complete"): status="failed"
                self.update(jid,status=status,error=None if status=="complete" else "Inspect the task report and log; resolve configuration before retrying.")
                self.message(sid,"assistant",f"Task {jid}: {status}. Report: {job['report']}")
                if status == "complete":
                    summary = screening_summary(job)
                    if summary: self.message(sid,"assistant",screening_feedback(summary))
                if status == "blocked" and report and report.get("clarifications"):
                    self.message(sid,"assistant"," ".join(report["clarifications"]))

    def close(self):
        self.stop.set()
        if self.worker: self.worker.join(timeout=8)
        if self.process and self.process.poll() is None: self.terminate(self.process)
