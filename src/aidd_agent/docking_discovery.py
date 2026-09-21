"""Read-only docking installation discovery; never execute discovered programs."""

import argparse
import json
import os
import shutil
from pathlib import Path


def discover(profile=None):
    config = {}
    if profile:
        config = json.loads(Path(profile).expanduser().read_text(encoding="utf-8"))
        if not isinstance(config, dict) or set(config) - {"schrodinger_root", "plants_executable"}:
            raise ValueError("Unsupported docking profile fields")
        if any(not isinstance(v, str) or not v.strip() for v in config.values()):
            raise ValueError("Docking paths must be non-empty strings")
    root = config.get("schrodinger_root") or os.environ.get("SCHRODINGER")
    root = str(Path(root).expanduser()) if root else None
    explicit_plants = config.get("plants_executable") or os.environ.get("AIDD_PLANTS_EXECUTABLE")
    commands, checks = {}, {}
    for name in ("maestro", "glide", "ligprep", "prepwizard", "structconvert", "PLANTS", "plants"):
        if name in ("PLANTS", "plants"):
            candidates = ([Path(explicit_plants).expanduser()] if explicit_plants else
                          [Path.home() / "PLANTS1.2" / "PLANTS1.2_64bit"])
            found = None if explicit_plants else (shutil.which(name) or shutil.which("PLANTS1.2_64bit"))
        else:
            candidates = ([Path(root) / name, Path(root) / "utilities" / name] if root else [])
            found = None if root else shutil.which(name)
        selected = Path(found) if found else next((p for p in candidates if p.is_file()), None)
        executable = bool(selected and os.access(selected, os.X_OK))
        commands[name] = str(selected) if executable else None
        checks[name] = {"status": "executable" if executable else
                        "not_executable" if selected else "missing",
                        "candidates": [str(p) for p in candidates],
                        "selected": str(selected) if selected else None}
    return {"schrodinger_root": root, "commands": commands, "checks": checks,
            "scope": "Read-only installation discovery; licenses, grids, docking execution and quality not tested"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile")
    args = parser.parse_args()
    print(json.dumps(discover(args.profile), indent=2))


if __name__ == "__main__":
    main()
