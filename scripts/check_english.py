"""Check maintained Git-visible UTF-8 content for Han text; not a language classifier.

Includes tracked and nonignored untracked files. Excludes Git history, ignored raw
artifacts and external environments. This guard complements English authoring.
"""
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from aidd_agent.language_policy import contains_han


def check(root=ROOT):
    names = subprocess.check_output(
        ["git", "ls-files", "--cached", "--others", "--exclude-standard", "-z"], cwd=root
    ).decode("utf-8").split("\0")
    failures = []
    checked = 0
    for name in sorted(set(names) - {""}):
        if contains_han(name): failures.append(f"Non-English filename: {name!a}")
        path = root / name
        if not path.is_file(): continue
        raw = path.read_bytes()
        if b"\0" in raw: continue
        try: text = raw.decode("utf-8-sig")
        except UnicodeError: continue
        checked += 1
        for number, line in enumerate(text.splitlines(), 1):
            if contains_han(line): failures.append(f"{name}:{number}: Han text detected")
    return checked, failures


if __name__ == "__main__":
    count, failures = check()
    print("\n".join(failures) if failures else f"English-content guard passed: {count} maintained text files; no Han text.")
    raise SystemExit(bool(failures))
