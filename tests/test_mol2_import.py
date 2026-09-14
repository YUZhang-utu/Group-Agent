from pathlib import Path
import sqlite3

from aidd_agent.importer import import_library
from aidd_agent.mol2 import iter_mol2_records, split_conformer_name
from aidd_agent import registry
from aidd_agent.registry import connect, initialize, library_summary, register_record


def block(name: str, offset: float = 0.0) -> str:
    return f"""@<TRIPOS>MOLECULE
{name}
2 1 0 0 0
SMALL
USER_CHARGES

@<TRIPOS>ATOM
1 C1 {offset:.1f} 0.0 0.0 C.3 1 MOL 0.0
2 N1 {offset + 1:.1f} 0.0 0.0 N.3 1 MOL 0.0
@<TRIPOS>BOND
1 1 2 1
"""


def test_terminal_conformer_name_parsing() -> None:
    assert split_conformer_name("macro_A_conf0") == ("macro_A", 0)
    assert split_conformer_name("macro_A_conf_12") == ("macro_A", 12)
    assert split_conformer_name("macro_confidence") == ("macro_confidence", 0)


def test_reads_multiple_records_from_one_file(tmp_path: Path) -> None:
    source = tmp_path / "batch.mol2"
    source.write_text(block("macro_A_conf0") + block("macro_A_conf1", 2.0), encoding="utf-8")
    records = list(iter_mol2_records(source))
    assert len(records) == 2
    assert [record.molecule_name for record in records] == ["macro_A", "macro_A"]
    assert [record.conformer_index for record in records] == [0, 1]


def test_import_is_idempotent_and_groups_conformers(tmp_path: Path) -> None:
    source = tmp_path / "batch.mol2"
    source.write_text(
        block("macro_A_conf0") + block("macro_A_conf1", 2.0) + block("macro_B", 4.0),
        encoding="utf-8",
    )
    database = tmp_path / "aidd.sqlite3"
    first = import_library(database, "macrocycles-v1", source)
    second = import_library(database, "macrocycles-v1", source)
    assert (first.inserted, first.invalid) == (3, 0)
    assert (second.duplicates, second.invalid) == (3, 0)
    with connect(database) as connection:
        summary = library_summary(connection, "macrocycles-v1")
    assert summary["molecule_count"] == 2
    assert summary["conformer_count"] == 3
    assert summary["conformers_per_molecule"] == {"1": 1, "2": 1}


def test_conflicting_same_index_is_reported(tmp_path: Path) -> None:
    source = tmp_path / "batch.mol2"
    source.write_text(block("macro_A_conf0") + block("macro_A_conf0", 3.0), encoding="utf-8")
    report = import_library(tmp_path / "aidd.sqlite3", "macrocycles-v1", source)
    assert report.inserted == 1
    assert report.invalid == 1
    assert "Conformer index conflict" in report.issues[0].error


def test_random_primary_id_collisions_are_retried(tmp_path: Path, monkeypatch) -> None:
    database = tmp_path / "aidd.sqlite3"
    source = tmp_path / "batch.mol2"
    source.write_text(block("macro_A_conf0"), encoding="utf-8")
    record = next(iter_mol2_records(source))
    initialize(database)

    generated = iter([
        "MOL-AAAAAAAAAAAA", "MOL-BBBBBBBBBBBB",
        "CNF-AAAAAAAAAAAA", "CNF-BBBBBBBBBBBB",
    ])
    monkeypatch.setattr(registry, "stable_id", lambda _prefix: next(generated))

    with connect(database) as connection:
        connection.execute(
            "INSERT INTO library(id, name, created_at) VALUES (?, ?, ?)",
            ("LIB-TEST", "test", "2026-09-14"),
        )
        connection.execute(
            "INSERT INTO molecule(id, library_id, source_name, created_at) VALUES (?, ?, ?, ?)",
            ("MOL-AAAAAAAAAAAA", "LIB-TEST", "existing", "2026-09-14"),
        )
        connection.execute(
            "INSERT INTO conformer VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
            ("CNF-AAAAAAAAAAAA", "MOL-AAAAAAAAAAAA", 0, "existing_conf0", 1, 0,
             "old-content", "old-topology", "old.mol2", 0, "[]", "2026-09-14"),
        )

        assert register_record(connection, "LIB-TEST", record) == "inserted"
        molecule = connection.execute(
            "SELECT id FROM molecule WHERE source_name = ?", ("macro_A",)
        ).fetchone()
        conformer = connection.execute(
            "SELECT id FROM conformer WHERE molecule_id = ?", (molecule["id"],)
        ).fetchone()

    assert molecule["id"] == "MOL-BBBBBBBBBBBB"
    assert conformer["id"] == "CNF-BBBBBBBBBBBB"
