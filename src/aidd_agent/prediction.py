from __future__ import annotations

import json
import csv
import hashlib
from pathlib import Path
import re
from typing import Any, Sequence

from .project_context import ensure_within

SAFE_NAME = re.compile(r"^[A-Za-z0-9_.-]+$")
AA_SEQUENCE = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]+$")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def inspect_alphafold3_output(
    output_dir: Path, *, construct_name: str, model_version: str,
    chain_ids: Sequence[str], input_path: Path | None = None,
    runtime: dict[str, Any] | None = None,
) -> dict[str, Any]:
    root = output_dir.resolve()
    if not root.is_dir():
        raise FileNotFoundError(f"AlphaFold 3 output directory does not exist: {root}")
    final_models = sorted(
        (path for path in root.rglob("*_model.cif")
         if "seed-" not in path.name and not any(
             part.startswith("seed-") for part in path.relative_to(root).parts)),
        key=lambda path: (len(path.relative_to(root).parts), str(path)),
    )
    if len(final_models) != 1:
        raise ValueError(
            f"Expected one final AlphaFold 3 model in {root}, found {len(final_models)}")
    structure = final_models[0]
    ranking_files = sorted(root.rglob("*_ranking_scores.csv"))
    ranking_score = None
    ranking_path = ranking_files[0] if len(ranking_files) == 1 else None
    if ranking_path:
        with ranking_path.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        values = [float(row["ranking_score"]) for row in rows if row.get("ranking_score")]
        ranking_score = max(values) if values else None
    summary_files = sorted(
        path for path in root.rglob("*_summary_confidences.json")
        if "seed-" not in path.name and not any(
            part.startswith("seed-") for part in path.relative_to(root).parts))
    confidence: dict[str, Any] = {}
    summary_path = summary_files[0] if len(summary_files) == 1 else None
    if summary_path:
        confidence = json.loads(summary_path.read_text(encoding="utf-8"))
    provenance: dict[str, Any] = {
        "output_dir": str(root),
        "ranking_path": str(ranking_path) if ranking_path else None,
        "ranking_sha256": _sha256(ranking_path) if ranking_path else None,
        "confidence_path": str(summary_path) if summary_path else None,
        "confidence_sha256": _sha256(summary_path) if summary_path else None,
        "runtime": runtime or {},
    }
    if input_path:
        resolved_input = input_path.resolve()
        if not resolved_input.is_file():
            raise FileNotFoundError(f"AlphaFold 3 input does not exist: {resolved_input}")
        provenance.update({"input_path": str(resolved_input),
                           "input_sha256": _sha256(resolved_input)})
    return {
        "schema_version": 1,
        "backend": "alphafold3",
        "model_name": "AlphaFold 3",
        "model_version": model_version,
        "construct_name": construct_name,
        "chain_ids": list(chain_ids),
        "structure_path": str(structure),
        "structure_sha256": _sha256(structure),
        "ranking_score": ranking_score,
        "confidence": confidence,
        "provenance": provenance,
    }


def load_model_profile(path: Path) -> dict[str, Any]:
    profile = json.loads(path.read_text(encoding="utf-8"))
    backend = profile.get("backend")
    common = {"backend"}
    required_by_backend = {
        "alphafold2": {"python", "runner", "model_parameters", "databases"},
        "alphafold3": {"python", "runner", "model_parameters", "databases"},
        "boltz2": {"executable", "cache"},
        "chai1": {"executable"},
    }
    if backend not in required_by_backend:
        raise ValueError("Unsupported prediction backend")
    required = common | required_by_backend[backend]
    missing = sorted(required - profile.keys())
    if missing:
        raise ValueError(f"Model profile is missing: {', '.join(missing)}")
    if profile["backend"] == "alphafold3" and profile.get("license_acknowledged") is not True:
        raise ValueError("AlphaFold 3 model-parameter terms must be acknowledged in the local profile")
    return profile


def _clean_sequences(sequences: Sequence[str]) -> list[str]:
    cleaned = [sequence.replace(" ", "").replace("\n", "").upper()
               for sequence in sequences]
    if not cleaned or any(not AA_SEQUENCE.fullmatch(sequence) for sequence in cleaned):
        raise ValueError("Protein sequences must use the 20 standard amino-acid letters")
    if len(cleaned) > 26:
        raise ValueError("At most 26 protein chains are supported")
    return cleaned


def write_alphafold3_input(name: str, sequences: Sequence[str], output: Path,
                           project_root: Path, *, seeds: Sequence[int] = (1,)) -> Path:
    if not SAFE_NAME.fullmatch(name):
        raise ValueError("Unsafe job name")
    cleaned = _clean_sequences(sequences)
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


def write_prediction_input(backend: str, name: str, sequences: Sequence[str],
                           output: Path, project_root: Path,
                           *, seeds: Sequence[int] = (1,)) -> Path:
    if backend == "alphafold3":
        return write_alphafold3_input(name, sequences, output, project_root, seeds=seeds)
    if backend not in {"boltz2", "chai1", "alphafold2"}:
        raise ValueError("Unsupported prediction backend")
    if not SAFE_NAME.fullmatch(name):
        raise ValueError("Unsafe job name")
    cleaned = _clean_sequences(sequences)
    target = ensure_within(output, project_root)
    target.parent.mkdir(parents=True, exist_ok=True)
    if backend == "boltz2":
        lines = ["version: 1", "sequences:"]
        for index, sequence in enumerate(cleaned):
            lines.extend(["  - protein:", f"      id: {chr(65 + index)}",
                          f"      sequence: {sequence}"])
        target.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    else:
        lines = []
        for index, sequence in enumerate(cleaned):
            chain = chr(65 + index)
            header = f">protein|{name}_{chain}" if backend == "chai1" else f">{name}_{chain}"
            lines.extend([header, sequence])
        target.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    return target


def prediction_command(profile: dict[str, Any], input_path: Path,
                       output_dir: Path) -> list[str]:
    backend = profile["backend"]
    if backend == "boltz2":
        command = [str(profile["executable"]), "predict", str(input_path),
                   "--out_dir", str(output_dir), "--cache", str(profile["cache"])]
        if profile.get("use_msa_server", False):
            command.append("--use_msa_server")
        if profile.get("use_potentials", True):
            command.append("--use_potentials")
        return command
    if backend == "chai1":
        command = [str(profile["executable"]), "fold", str(input_path), str(output_dir)]
        if profile.get("use_msa_server", False):
            command.append("--use-msa-server")
        if profile.get("use_templates_server", False):
            command.append("--use-templates-server")
        return command
    base = [str(profile["python"]), str(profile["runner"])]
    if backend == "alphafold3":
        return base + [f"--json_path={input_path}", f"--output_dir={output_dir}",
                       f"--model_dir={profile['model_parameters']}",
                       f"--db_dir={profile['databases']}"]
    return base + [f"--fasta_paths={input_path}", f"--output_dir={output_dir}",
                   f"--data_dir={profile['databases']}",
                   f"--model_preset={profile.get('model_preset', 'monomer')}"]
