"""E037: fixed factorial comparisons, reverse-order repeats and one report."""
from __future__ import annotations

import argparse
import os
from pathlib import Path
import subprocess
import sys
import time

import numpy as np

from . import library_acceptance as ev
from .expanded_wee1 import fingerprint
from .gaussian_batch import _atomic_json

SCENARIOS = (
    ("old-reference", 2000, 250, False),
    ("fine-reference", 500, 64, False),
    ("fine-bounded", 500, 64, True),
    ("old-bounded", 2000, 250, True),
)


def summarize(runs):
    groups = {}
    for run in runs:
        if run["status"] != "complete" or run.get("profile"):
            continue
        for query in run["report"]["queries"]:
            key = (run["scenario"], query["query_id"])
            row = groups.setdefault(key, dict(scenario=key[0], query_id=key[1], seconds=[], eligible=True))
            row["seconds"].append(query["execution_seconds"])
            row["eligible"] &= bool(query["latency_eligible"] and query["gaussian_equivalence"]["passed"]
                                    and query["annotation_equivalence"]["passed"])
    rows = []
    for row in groups.values():
        row["eligible"] &= len(row["seconds"]) >= 2
        row["median_seconds"] = float(np.median(row["seconds"])) if row["eligible"] else None
        rows.append(row)
    baselines = {r["query_id"]: r["median_seconds"] for r in rows if r["scenario"] == "old-reference"}
    for row in rows:
        baseline = baselines.get(row["query_id"])
        row["speedup_vs_old_reference"] = baseline / row["median_seconds"] if baseline and row["median_seconds"] else None
    return rows


def run(args):
    import fcntl
    output = args.output.resolve()
    batch, source = args.batch.resolve(), args.e034.resolve()
    ev.require(args.repeats >= 2 and args.workers > 0, "At least 2 repeats and positive workers required")
    for protected in (batch, source):
        ev.require(not output.is_relative_to(protected) and not protected.is_relative_to(output), "Output overlaps inputs")
    ev.require(not output.exists() or args.resume, "Suite exists; use --resume or a new output")
    output.mkdir(parents=True, exist_ok=True)
    with (output / "suite.lock").open("a") as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        protocol = dict(format="aidd-e037-suite", workers=args.workers, repeats=args.repeats,
                        profile=args.profile, batch=str(batch), e034=str(source),
                        sources=fingerprint([source / "report.json", source / "protocol.json",
                                             batch / "artifacts/catalog.json", batch / "chemical/catalog.json"]),
                        code=fingerprint(sorted(Path(__file__).parent.glob("*.py"))))
        path = output / "protocol.json"
        if path.exists():
            ev.require(ev.read(path) == protocol, "Changed suite protocol; use fresh output")
        elif any(output.glob("run-*")):
            raise ValueError("Existing suite runs without protocol; preserve directory and use a new output")
        _atomic_json(path, protocol)
        tasks = []
        for repeat in range(args.repeats):
            scenarios = SCENARIOS if repeat % 2 == 0 else tuple(reversed(SCENARIOS))
            tasks.extend((repeat, scenario, False) for scenario in scenarios)
        if args.profile:
            tasks.append((0, SCENARIOS[2], True))
        runs = []
        _atomic_json(output / "RUN_STATUS.json", dict(status="running", total=len(tasks)))
        for repeat, (name, coarse, refine, bounded), profile in tasks:
            folder = output / (f"profile-{name}" if profile else f"run-{repeat+1}-{name}")
            command = [sys.executable, "-u", "-m", "aidd_agent.fast_3d_search",
                       "--batch", str(batch), "--e034", str(source), "--output", str(folder),
                       "--workers", str(args.workers), "--coarse-chunk", str(coarse),
                       "--refine-chunk", str(refine), "--fresh-retrieval"]
            if bounded: command.append("--bounded-pair-seeds")
            if profile: command.append("--profile-first-chunk")
            if folder.exists(): command.append("--resume")
            ev.log(f"Suite {len(runs)+1}/{len(tasks)}: {folder.name}")
            started = time.perf_counter()
            log = output / (folder.name + ".log")
            env = dict(os.environ)
            for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
                env[key] = "1"
            with log.open("a", encoding="utf-8") as stream:
                completed = subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT, env=env)
            record = dict(scenario=name, repeat=repeat+1, profile=profile, output=str(folder),
                          log=str(log), process_wall_seconds=time.perf_counter()-started,
                          status="failed", returncode=completed.returncode)
            if completed.returncode == 0:
                try:
                    marker = ev.read(folder / "RUN_STATUS.json")
                    ev.require(marker["status"] == "complete" and marker["report_sha256"] == ev.sha(folder / "report.json"),
                               "Child report receipt mismatch")
                    report = ev.read(folder / "report.json")
                    ev.require(report["status"] == "complete", "Child incomplete")
                    record.update(status="complete", report=report, report_sha256=marker["report_sha256"])
                except (OSError, ValueError, KeyError) as exc:
                    record["error"] = str(exc)
            runs.append(record)
            _atomic_json(output / "progress.json", dict(completed=len(runs), total=len(tasks), runs=runs))
            ev.log(f"{folder.name}: {record['status']}; log={log}")
        result = dict(status="complete" if all(r["status"] == "complete" for r in runs) else "failed",
                      runs=runs, comparisons=summarize(runs),
                      interpretation="exploratory paired scheduling/kernel timings; reversed order, no SLA; profiles excluded",
                      e031_changes_ranking=False)
        _atomic_json(output / "report.json", result)
        lines = ["# E037 workstation suite", "", f"Status: {result['status']}", "",
                 "| Scenario | Query | Execution median (s) | Speedup vs old/reference | Eligible |",
                 "|---|---|---:|---:|---|"]
        for row in result["comparisons"]:
            seconds = "unavailable" if row["median_seconds"] is None else f"{row['median_seconds']:.3f}"
            speed = "unavailable" if row["speedup_vs_old_reference"] is None else f"{row['speedup_vs_old_reference']:.2f}"
            lines.append(f"| {row['scenario']} | {row['query_id']} | {seconds} | {speed} | {row['eligible']} |")
        lines += ["", "All scenarios preserve candidate and seed budgets; E031 remains annotation-only.",
                  "Profile runs and partial chunk reuse are excluded from latency comparison.",
                  "Startup/integrity/index load is in child reports; execution excludes final equivalence checks.",
                  "Inspect failed child logs before interpreting results; no whole-library exact scans or million stress rerun."]
        (output / "report.md").write_text("\n".join(lines)+"\n", encoding="utf-8")
        _atomic_json(output / "RUN_STATUS.json", dict(status=result["status"], report_sha256=ev.sha(output / "report.json")))
        return 0 if result["status"] == "complete" else 1


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("batch", "e034", "output"):
        p.add_argument("--"+name, type=Path, required=True)
    p.add_argument("--workers", type=int, default=16)
    p.add_argument("--repeats", type=int, default=2)
    p.add_argument("--profile", action="store_true")
    p.add_argument("--resume", action="store_true")
    return run(p.parse_args())


if __name__ == "__main__":
    raise SystemExit(main())
