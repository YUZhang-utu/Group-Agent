from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
from pathlib import Path
import re
from typing import Iterator

MOLECULE_MARKER = "@<TRIPOS>MOLECULE"
SECTION_PREFIX = "@<TRIPOS>"
CONFORMER_RE = re.compile(r"^(?P<base>.+?)[_-]conf[_-]?(?P<index>\d+)$", re.IGNORECASE)


class Mol2Error(ValueError):
    """Raised when a MOL2 record cannot be safely registered."""


@dataclass(frozen=True)
class Mol2Record:
    source_path: Path
    record_index: int
    name: str
    molecule_name: str
    conformer_index: int
    atom_count: int
    bond_count: int
    raw_text: str
    content_sha256: str
    topology_sha256: str
    warnings: tuple[str, ...]


def split_conformer_name(name: str) -> tuple[str, int]:
    """Split only a terminal `_confN`/`_conf_N` suffix.

    Restricting the match to the end avoids silently changing legitimate names
    that merely contain the word `conf`.
    """
    cleaned = name.strip()
    match = CONFORMER_RE.fullmatch(cleaned)
    if not match:
        return cleaned, 0
    return match.group("base"), int(match.group("index"))


def iter_mol2_blocks(path: Path) -> Iterator[tuple[int, str]]:
    """Yield individual MOL2 molecule blocks without loading the whole file."""
    current: list[str] = []
    index = 0
    with path.open("r", encoding="utf-8", errors="replace", newline=None) as handle:
        for line in handle:
            if line.strip().upper() == MOLECULE_MARKER:
                if current:
                    yield index, "".join(current)
                    index += 1
                current = [line]
            elif current:
                current.append(line)
    if current:
        yield index, "".join(current)
    elif path.stat().st_size:
        raise Mol2Error(f"No {MOLECULE_MARKER} section in {path}")


def _sections(lines: list[str]) -> dict[str, list[str]]:
    sections: dict[str, list[str]] = {}
    current: str | None = None
    for line in lines:
        stripped = line.strip()
        if stripped.upper().startswith(SECTION_PREFIX):
            current = stripped.upper()
            sections.setdefault(current, [])
        elif current is not None:
            sections[current].append(line)
    return sections


def parse_mol2_block(path: Path, record_index: int, raw_text: str) -> Mol2Record:
    lines = raw_text.splitlines()
    sections = _sections(lines)
    molecule_lines = sections.get(MOLECULE_MARKER, [])
    if len(molecule_lines) < 2:
        raise Mol2Error(f"Missing name/count line in {path} record {record_index}")

    name = molecule_lines[0].strip()
    if not name:
        raise Mol2Error(f"Empty molecule name in {path} record {record_index}")

    counts = molecule_lines[1].split()
    try:
        declared_atoms = int(counts[0])
        declared_bonds = int(counts[1])
    except (IndexError, ValueError) as exc:
        raise Mol2Error(f"Invalid count line in {path} record {record_index}") from exc

    atom_lines = [line for line in sections.get("@<TRIPOS>ATOM", []) if line.strip()]
    bond_lines = [line for line in sections.get("@<TRIPOS>BOND", []) if line.strip()]
    warnings: list[str] = []
    if not atom_lines:
        raise Mol2Error(f"Missing ATOM section in {path} record {record_index}")
    if declared_atoms != len(atom_lines):
        raise Mol2Error(
            f"Atom count mismatch in {path} record {record_index}: "
            f"declared {declared_atoms}, found {len(atom_lines)}"
        )
    if declared_bonds != len(bond_lines):
        raise Mol2Error(
            f"Bond count mismatch in {path} record {record_index}: "
            f"declared {declared_bonds}, found {len(bond_lines)}"
        )
    if declared_bonds == 0:
        warnings.append("record_has_no_bonds")

    # A coordinate-independent fingerprint useful for detecting repeated
    # conformers. This is not a chemical identity key and is deliberately kept
    # separate from molecule-name grouping.
    atom_types = []
    for line in atom_lines:
        fields = line.split()
        if len(fields) < 6:
            raise Mol2Error(f"Malformed atom line in {path} record {record_index}")
        atom_types.append(fields[5])
    bond_types = []
    for line in bond_lines:
        fields = line.split()
        if len(fields) < 4:
            raise Mol2Error(f"Malformed bond line in {path} record {record_index}")
        bond_types.append(f"{fields[1]}-{fields[2]}:{fields[3]}")
    topology = "|".join(atom_types) + "||" + "|".join(bond_types)

    molecule_name, conformer_index = split_conformer_name(name)
    return Mol2Record(
        source_path=path.resolve(),
        record_index=record_index,
        name=name,
        molecule_name=molecule_name,
        conformer_index=conformer_index,
        atom_count=declared_atoms,
        bond_count=declared_bonds,
        raw_text=raw_text,
        content_sha256=sha256(raw_text.encode("utf-8")).hexdigest(),
        topology_sha256=sha256(topology.encode("utf-8")).hexdigest(),
        warnings=tuple(warnings),
    )


def iter_mol2_records(path: Path) -> Iterator[Mol2Record]:
    for record_index, raw_text in iter_mol2_blocks(path):
        yield parse_mol2_block(path, record_index, raw_text)


def discover_mol2_files(input_path: Path, recursive: bool = True) -> list[Path]:
    if input_path.is_file():
        if input_path.suffix.lower() != ".mol2":
            raise Mol2Error(f"Input file is not MOL2: {input_path}")
        return [input_path]
    if not input_path.is_dir():
        raise Mol2Error(f"Input path does not exist: {input_path}")
    iterator = input_path.rglob("*") if recursive else input_path.glob("*")
    return sorted(path for path in iterator if path.is_file() and path.suffix.lower() == ".mol2")
