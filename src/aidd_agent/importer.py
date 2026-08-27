from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Callable

from .mol2 import Mol2Error, discover_mol2_files, iter_mol2_records
from .registry import connect, ensure_library, initialize, register_record


@dataclass
class ImportIssue:
    path: str
    record_index: int | None
    error: str


@dataclass
class ImportReport:
    library: str
    files_seen: int = 0
    records_seen: int = 0
    inserted: int = 0
    duplicates: int = 0
    invalid: int = 0
    dry_run: bool = False
    issues: list[ImportIssue] | None = None

    def to_dict(self) -> dict[str, object]:
        result = asdict(self)
        result["issues"] = [asdict(issue) for issue in (self.issues or [])]
        return result


def import_library(
    db_path: Path,
    library_name: str,
    input_path: Path,
    *,
    recursive: bool = True,
    dry_run: bool = False,
    progress_every: int = 0,
    progress: Callable[[ImportReport, Path], None] | None = None,
) -> ImportReport:
    report = ImportReport(library=library_name, dry_run=dry_run, issues=[])
    files = discover_mol2_files(input_path, recursive=recursive)
    report.files_seen = len(files)
    if not files:
        raise Mol2Error(f"No MOL2 files found under {input_path}")

    if dry_run:
        for path in files:
            try:
                for record in iter_mol2_records(path):
                    report.records_seen += 1
                    if progress and progress_every and report.records_seen % progress_every == 0:
                        progress(report, path)
                    if record.warnings:
                        report.issues.append(ImportIssue(str(path), record.record_index, ",".join(record.warnings)))
            except (Mol2Error, OSError) as exc:
                report.invalid += 1
                report.issues.append(ImportIssue(str(path), None, str(exc)))
        return report

    initialize(db_path)
    with connect(db_path) as connection:
        library_id = ensure_library(connection, library_name)
        for path in files:
            try:
                for record in iter_mol2_records(path):
                    report.records_seen += 1
                    try:
                        outcome = register_record(connection, library_id, record)
                    except ValueError as exc:
                        report.invalid += 1
                        report.issues.append(ImportIssue(str(path), record.record_index, str(exc)))
                        continue
                    if outcome == "inserted":
                        report.inserted += 1
                    else:
                        report.duplicates += 1
                    if record.warnings:
                        report.issues.append(ImportIssue(str(path), record.record_index, ",".join(record.warnings)))
                    if progress and progress_every and report.records_seen % progress_every == 0:
                        progress(report, path)
            except (Mol2Error, OSError) as exc:
                report.invalid += 1
                report.issues.append(ImportIssue(str(path), None, str(exc)))
    return report
