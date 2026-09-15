"""Read-only extraction of an ambiguous source conformer and its registry peer.

Streams source files; never changes the registry, source records or batch state.
Outputs JSON and bounded text differences to stdout for workstation diagnosis.
"""
from __future__ import annotations

import argparse
import difflib
import json
from pathlib import Path
import sqlite3

from aidd_agent.mol2 import iter_mol2_records


def describe(record):
    return dict(source_path=str(record.source_path), record_index=record.record_index,
                name=record.name, content_sha256=record.content_sha256,
                topology_sha256=record.topology_sha256,
                atom_count=record.atom_count, bond_count=record.bond_count)


def diagnose(dbpath, source, molecule, index):
    # mode=ro rejects a missing DB instead of creating a misleading empty one.
    with sqlite3.connect(dbpath.resolve().as_uri() + '?mode=ro', uri=True) as db:
        db.row_factory = sqlite3.Row
        rows = db.execute(
            'SELECT c.*, m.library_id FROM conformer c JOIN molecule m ON m.id=c.molecule_id '
            'WHERE m.source_name=? AND c.conformer_index=?', (molecule, index)).fetchall()
    existing = []
    for row in rows:
        item = dict(row)
        record = None
        try:
            record = next((r for r in iter_mol2_records(Path(row['source_path']))
                           if r.record_index == row['source_record_index']), None)
            item['source_record_found'] = record is not None
            item['source_matches_registered_hash'] = (
                record.content_sha256 == row['content_sha256'] if record else False)
        except (OSError, ValueError) as exc:
            item['source_read_error'] = str(exc)
        existing.append((item, record))
    incoming = [r for r in iter_mol2_records(source)
                if r.molecule_name == molecule and r.conformer_index == index]
    comparisons = []
    # Include same-file peers because the failed transaction may have rolled
    # back the only earlier occurrence, leaving no committed registry row.
    peers = [(f'registry:{item["id"]}', record) for item, record in existing if record]
    peers += [(f'incoming:{r.record_index}', r) for r in incoming]
    for record in incoming:
        for label, peer in peers:
            if peer.source_path == record.source_path and peer.record_index == record.record_index:
                continue
            diff = list(difflib.unified_diff(peer.raw_text.splitlines(), record.raw_text.splitlines(),
                                            fromfile=label, tofile=f'incoming:{record.record_index}', n=2))
            comparisons.append(dict(incoming_record_index=record.record_index, peer=label,
                                    same_raw_text=peer.raw_text == record.raw_text,
                                    same_whitespace_tokens=peer.raw_text.split() == record.raw_text.split(),
                                    same_ordered_topology_hash=peer.topology_sha256 == record.topology_sha256,
                                    diff=diff[:100], diff_truncated=len(diff) > 100))
    return dict(molecule=molecule, conformer_index=index, record_indices='zero-based',
                registered=[item for item, _ in existing], incoming=[describe(r) for r in incoming],
                comparisons=comparisons,
                note='Topology hash is not a chemical identity or stereochemistry test. No data modified.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--registry', type=Path, required=True)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--molecule', required=True)
    parser.add_argument('--conformer-index', type=int, required=True)
    args = parser.parse_args()
    print(json.dumps(diagnose(args.registry, args.source, args.molecule, args.conformer_index),
                     indent=2, ensure_ascii=False))


if __name__ == '__main__':
    main()
