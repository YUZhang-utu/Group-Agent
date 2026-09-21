"""Convenient workstation launcher using a persisted active Project context."""
import argparse
import json
from pathlib import Path

from aidd_agent.llm_profiles import select_llm_profile
from aidd_agent.prompt_workflow import initialize_context, create_plan, run_plan


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--storage-root", type=Path, required=True)
    parser.add_argument("--username", default="workstation")
    parser.add_argument("--project-name", default="Prompt AIDD")
    provider = parser.add_mutually_exclusive_group()
    provider.add_argument("--llm-profile", type=Path)
    provider.add_argument("--provider", choices=("gpt", "deepseek"))
    parser.add_argument("--llm-config-dir", type=Path)
    parser.add_argument("--runtime", type=Path)
    text = parser.add_mutually_exclusive_group(required=True)
    text.add_argument("--prompt"); text.add_argument("--prompt-file", type=Path)
    parser.add_argument("--plan-only", action="store_true")
    parser.add_argument("--allow-compute", action="store_true")
    args = parser.parse_args()
    args.llm_profile = select_llm_profile(args.provider, args.llm_profile, args.llm_config_dir)
    context = initialize_context(args.storage_root, args.username, args.project_name)
    prompt = args.prompt if args.prompt is not None else args.prompt_file.read_text(encoding="utf-8")
    plan = create_plan(Path(context["db"]), context["user_id"], context["project_id"], prompt, llm_profile=args.llm_profile)
    print(json.dumps(dict(plan=str(plan), context=context), ensure_ascii=False))
    if args.plan_only: return 0
    report, path = run_plan(Path(context["db"]), context["user_id"], context["project_id"], plan,
                            runtime=args.runtime, allow_compute=args.allow_compute)
    print(json.dumps(dict(status=report["status"], report=str(path)), ensure_ascii=False))
    return 0 if report["status"] == "complete" else 2


if __name__ == "__main__": raise SystemExit(main())
