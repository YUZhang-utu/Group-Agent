from __future__ import annotations

import importlib.util
import json
import platform
import shutil
import sys


def environment_report() -> dict:
    modules = {name: importlib.util.find_spec(name) is not None
               for name in ("numpy", "scipy", "rdkit", "Bio", "gemmi", "pymol")}
    return {
        "python": sys.version.split()[0], "platform": platform.platform(),
        "modules": modules,
        "executables": {name: shutil.which(name) for name in ("git", "ssh", "plink", "pymol")},
        "ready": modules["numpy"] and modules["scipy"] and modules["rdkit"] and modules["Bio"],
    }


def main() -> int:
    report = environment_report()
    print(json.dumps(report, indent=2))
    return 0 if report["ready"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
