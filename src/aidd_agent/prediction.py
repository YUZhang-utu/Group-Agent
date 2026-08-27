from __future__ import annotations

import json
from pathlib import Path
import re
from typing import Any, Sequence

from .project_context import ensure_within

SAFE_NAME = re.compile(r"^[A-Za-z0-9_.-]+$")
AA_SEQUENCE = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]+$")


def load_model_profile(path: Path) -> dict[str, Any]:
    profile = json.loads(path.read_text(encoding="utf-8"))
    required = {"backend", "python", "runner", "model_parameters", "databases"}
    missing = sorted(required - profile.keys())
    if missing:
        raise ValueError(f"Model profile is missing: {', '.join(missing)}")
    if profile["backend"] not in {"alphafold2", "alphafold3"}:
        raise ValueError("Unsupported prediction backend")
    if profile["backend"] == "alphafold3" and profile.get("license_acknowledged") is not True:
        raise ValueError("AlphaFold 3 model-parameter terms must be acknowledged in the local profile")
    return profile


def write_alphafold3_input(name: str, sequences: Sequence[str], output: Path,
                           project_root: Path, *, seeds: Sequence[int] = (1,)) -> Path:
    if not SAFE_NAME.fullmatch(name):
        raise ValueError("Unsafe job name")
    cleaned = [sequence.replace(" ", "").replace("\n", "").upper() for sequence in sequences]
    if not cleaned or any(not AA_SEQUENCE.fullmatch(sequence) for sequence in cleaned):
        raise ValueError("Protein sequences must use the 20 standard amino-acid letters")
    if not seeds or any(int(seed) < 0 for seed in seeds):
        raise ValueError("At least one non-negative model seed is required")
    target = ensure_within(output, project_root)
    entities = [{"protein": {"id": chr(65 + index), "sequence": sequence}}
                for index, sequence in enumerate(cleaned)]
    payload = {"name": name, "modelSeeds": [int(seed) for seed in seeds],
               "sequences": entities, "dialect": "alphafold3", "version": 1}
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return target


def prediction_command(profile: dict[str, Any], input_path: Path,
                       output_dir: Path) -> list[str]:
    base = [str(profile["python"]), str(profile["runner"])]
    if profile["backend"] == "alphafold3":
        return base + [f"--json_path={input_path}", f"--output_dir={output_dir}",
                       f"--model_dir={profile['model_parameters']}",
                       f"--db_dir={profile['databases']}"]
    return base + [f"--fasta_paths={input_path}", f"--output_dir={output_dir}",
                   f"--data_dir={profile['databases']}",
                   f"--model_preset={profile.get('model_preset', 'monomer')}"]

