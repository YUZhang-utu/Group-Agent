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
    {"name": "Protein and structure evidence", "status": "available", "scope": "Verified UniProt proteins and PDB retrieval; receptor selection still needs review."},
    {"name": "AF3 prediction", "status": "available", "scope": "Verified single protein and optional CCD ligands; installed native or Apptainer profile."},
    {"name": "3D screening", "status": "calibrated WEE1 only", "scope": "USRCAT retrieval, Gaussian rigid refinement and annotation-only E031 for QT9/824."},
    {"name": "New-target query preparation", "status": "not exposed in chat", "scope": "Requires reference ligand/site, query preparation and independent retrieval calibration."},
    {"name": "Molecule aggregation and pocket QC", "status": "CLI components; chat integration pending", "scope": "Requires validated query poses, receptor/site definitions and aggregation artifacts."},
    {"name": "Flexible docking and rescoring", "status": "not implemented as a validated executor", "scope": "Schrodinger/Maestro and PLANTS reported installed; CLI/license/grid discovery, receptor/ligand preparation and redocking/cross-docking remain pending."},
    {"name": "Biological enrichment and final selection", "status": "not validated", "scope": "Requires held-out active/decoy labels, diversity/property criteria and experimental evidence."},
]
ROUTER = """You are an AIDD conversational task coordinator. Return JSON only with exactly
intent, message, task_id, request. intent is run/status/results/resume/cancel/capabilities/clarify.
Write message in English. task_id is an existing task ID from this session or null (latest).
request is a self-contained scientific request only for run, otherwise empty.
Use conversation context to resolve follow-up clarifications. Never invent biological sequences,
paths, results, unsupported tools or identifiers. Questions about progress/results NEVER run a new
job. Resume/cancel only when explicitly requested. For run retain all explicit user constraints.
Supported actions are supplied. Unsupported stages must explain missing capability using clarify,
not substitute WEE1 for another target. AF3 prediction cannot automatically create a calibrated
search query. Missing organism/construct/reference information requires clarification. Do not
claim a task is queued or complete: the local executor reports its actual status.
JSON example: {"intent":"status","message":"Checking task status.","task_id":null,"request":""}.
"""


def validate_route(value):
    if not isinstance(value, dict) or set(value) != {"intent", "message", "task_id", "request"}:
        raise ValueError("Invalid chat decision")
    if value["intent"] not in {"run", "status", "results", "resume", "cancel", "capabilities", "clarify"}:
        raise ValueError("Unsupported chat intent")
    for field in ("message", "request"):
        if not isinstance(value[field], str) or len(value[field]) > 12000:
            raise ValueError("Invalid chat text")
    if contains_han(value["message"]): raise ValueError("Agent replies must be English")
    if value["task_id"] is not None and not isinstance(value["task_id"], str):
        raise ValueError("Invalid task reference")
    if (value["intent"] == "run") != bool(value["request"].strip()):
        raise ValueError("Execution requires a nonempty request; queries cannot contain one")
    return value


def read_json(path):
    try: return json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, ValueError): return None


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
            if command in {"/status", "/results", "/capabilities", "/resume", "/cancel"}:
                if len(parts) > 2: raise ValueError("Use /command followed by an optional task ID")
                decision = dict(intent=command[1:], message="", task_id=parts[1] if len(parts)==2 else None, request="")
            else:
                profile = select_llm_profile(provider, config_dir=self.config_dir)
                snapshot = self.snapshot(sid)
                context = dict(conversation=snapshot["messages"][-12:], tasks=[{k: j[k] for k in ("id","status","request")} for j in snapshot["tasks"][-10:]])
                # Bounded dialogue context; no files, sequences or execution logs sent to provider.
                for item in context["conversation"]: item["text"] = item["text"][:1000]
                for item in context["tasks"]: item["request"] = item["request"][:300]
                context["latest_message"] = text
                prompt = json.dumps(context, ensure_ascii=False)
                if len(prompt) > 20000: context["conversation"] = context["conversation"][-4:]; prompt = json.dumps(context, ensure_ascii=False)
                decision = (self.router(prompt, profile) if self.router else chat_plan(prompt, read_json(profile),
                    system_prompt=ROUTER, capabilities=dict(actions=CAPABILITIES, workflows=WORKFLOWS), validator=validate_route)[0])
                validate_route(decision)
            intent = decision["intent"]
            if intent == "run":
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
                if status == "blocked" and report and report.get("clarifications"):
                    self.message(sid,"assistant"," ".join(report["clarifications"]))

    def close(self):
        self.stop.set()
        if self.worker: self.worker.join(timeout=8)
        if self.process and self.process.poll() is None: self.terminate(self.process)
