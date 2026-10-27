"""Select scored KRAS molecules and export hash-checked original MOL2 batches."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import sqlite3
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from aidd_agent.mol2 import iter_mol2_blocks, parse_mol2_block, load_rdkit_mol2
from aidd_agent.kras_weekend import read, save, sha

TARGETS = ('KRAS_G12D', 'KRAS_G12V')


def readonly(path):
    return sqlite3.connect(Path(path).resolve().as_uri() + '?mode=ro', uri=True)


def selected_rows(db, target, count, fraction=.8):
    """Score quota followed by coarse-graph bucket round-robin, never random IDs."""
    first = int(count * fraction)
    seen = set()
    for mid, score, payload in db.execute(
            'SELECT mid,score,payload FROM molecules WHERE target=? ORDER BY score DESC,mid LIMIT ?', (target, first)):
        seen.add(mid)
        yield mid, score, payload, 'score'
    # Ranking within each bucket happens after excluding the leading score quota.
    # Stable score/mid ordering makes selection independent of worker completion.
    query = '''WITH leading AS (
        SELECT mid FROM molecules WHERE target=? ORDER BY score DESC,mid LIMIT ?),
        remaining AS (
        SELECT mid,score,payload,bucket FROM molecules WHERE target=?
        AND mid NOT IN (SELECT mid FROM leading)),
        ordered AS (
        SELECT *,ROW_NUMBER() OVER (PARTITION BY bucket ORDER BY score DESC,mid) AS depth FROM remaining)
        SELECT mid,score,payload FROM ordered ORDER BY depth,score DESC,mid LIMIT ?'''
    for mid, score, payload in db.execute(query, (target, first, target, count - len(seen))):
        if mid in seen:
            raise ValueError('Duplicate selection')
        seen.add(mid)
        yield mid, score, payload, 'diversity_proxy'


def rename_record(raw, alias):
    lines = raw.splitlines(keepends=True)
    if len(lines) < 3 or lines[0].strip() != '@<TRIPOS>MOLECULE':
        raise ValueError('Invalid MOL2 header')
    lines[1] = alias + '\n'
    return ''.join(lines).rstrip('\n') + '\n'


def source_pose_checker(batch,prepared):
    """Recheck the exact source coordinates rather than only quantized artifacts."""
    import numpy as np
    from aidd_agent.gaussian_batch import ArtifactCatalogReader
    from aidd_agent.chemical_companion import ChemicalCompanionReader
    from aidd_agent import kras_weekend as kw
    reader=ArtifactCatalogReader(batch/'artifacts/catalog.json')
    chemistry=ChemicalCompanionReader(batch/'chemical/catalog.json')
    states={target:kw.load_state(prepared,target) for target in TARGETS}
    def check(molecule,target,pose):
        candidate=reader.get(pose['gid']);chem=chemistry.get(pose['gid'])
        if any((r.molecule_id,r.conformer_id)!=(pose['mid'],pose['cid']) for r in (candidate,chem)):
            raise ValueError('Export artifact/chemical identity mismatch')
        heavy=[a.GetIdx() for a in molecule.GetAtoms() if a.GetAtomicNum()>1]
        numbers=[molecule.GetAtomWithIdx(i).GetAtomicNum() for i in heavy]
        points=molecule.GetConformer().GetPositions()[heavy]
        if numbers!=list(chem.atomic_numbers) or points.shape!=candidate.shape_points.shape:
            raise ValueError('Export source heavy-atom order mismatch')
        if np.max(np.abs(points-candidate.shape_points))>.02:
            raise ValueError('Export source/artifact coordinate mismatch')
        matrix=np.array(pose['transform']).reshape(4,4);state=states[target]
        minimum=float(state['tree'].query(points@matrix[:3,:3].T+matrix[:3,3])[0].min())
        sets=kw.chemical_sets(chem,candidate.feature_types)
        direct=bool(kw.direct_contact_masks(points,sets,state,matrix[None])['all'][0])
        return minimum>=2. and direct,minimum
    return check


class BatchWriter:
    def __init__(self, root, size):
        self.root = root; self.size = size; self.count = 0; self.stream = None; self.batches = []
        root.mkdir()

    def add(self, raw):
        if self.count % self.size == 0:
            self.close()
            path = self.root / f'batch_{len(self.batches) + 1:03d}.mol2'
            self.stream = path.open('w', encoding='utf-8', newline='')
            self.batches.append(dict(path=str(path), molecules=0))
        self.stream.write(raw); self.count += 1; self.batches[-1]['molecules'] += 1
        return len(self.batches), self.batches[-1]['molecules']

    def close(self):
        if self.stream is not None:
            self.stream.close(); self.stream = None
            self.batches[-1]['sha256'] = sha(Path(self.batches[-1]['path']))


def export(ranking, output, count=1000000, batch_size=100000, allow_pilot=False, allow_budgeted=False):
    ranking = Path(ranking).resolve(); output = Path(output).resolve()
    report = read(ranking / 'report.json'); protocol = read(ranking / 'protocol.json')
    budgeted = (allow_budgeted and report.get('scope') == 'budgeted_direct_pool_refinement'
                and report.get('budget_finished') is True and report.get('status') in ('complete','partial'))
    if report['status'] != 'complete' and not budgeted:
        raise ValueError('Complete scoring report required')
    if not report.get('full_catalog') and not allow_pilot and not budgeted:
        raise ValueError('Pilot cannot be exported as a full-library selection')
    if report.get('ranking_sha256') != sha(ranking / 'ranking.sqlite'):
        raise ValueError('Ranking database checksum mismatch')
    if count < 1 or batch_size < 1:
        raise ValueError('Positive molecule and batch counts required')
    original = protocol['common']['original']; batch = Path(original['batch']).resolve()
    prepared = Path(original['prepared']).resolve()
    for protected in (ranking, batch, prepared):
        if output.is_relative_to(protected) or protected.is_relative_to(output):
            raise ValueError('Output overlaps protected input')
    if output.exists():
        raise ValueError('Use a new export directory; existing files are never overwritten')
    output.mkdir(parents=True)
    save(output / 'source_ranking_report.json', report)
    save(output / 'source_ranking_protocol.json', protocol)
    save(output / 'report.json', dict(status='running', requested_per_target=count))
    save(output / 'protocol.json', dict(ranking_protocol_sha256=sha(ranking/'protocol.json'),
        ranking_report_sha256=sha(ranking/'report.json'), requested_per_target=count, batch_size=batch_size,
        selection='80% score; 20% round-robin coarse graph buckets', coordinates='original_source',
        representative='highest soft-score conformer; not an energy minimum', pilot=report.get('scope')=='spread_pilot',
        source_scope=report.get('scope'),source_processed=report.get('processed'),budgeted=budgeted,
        export_code_sha256=sha(Path(__file__))))
    started = time.monotonic()
    db = sqlite3.connect(output / 'selection.sqlite')
    db.executescript('''CREATE TABLE selected(target TEXT,mid TEXT,cid TEXT,score REAL,reason TEXT,payload TEXT,
        source TEXT,record_index INTEGER,source_sha TEXT, PRIMARY KEY(target,mid));
        CREATE INDEX source_records ON selected(source,record_index);''')
    source_db = readonly(ranking/'ranking.sqlite'); registry = readonly(batch/'registry.sqlite3')
    try:
        for target in TARGETS:
            eligible = source_db.execute('SELECT COUNT(*) FROM molecules WHERE target=?', (target,)).fetchone()[0]
            if eligible != report['targets'][target]['eligible_molecules']:
                raise ValueError('Ranking/report molecule count mismatch')
            with db:
                for mid, score, payload, reason in selected_rows(source_db, target, count):
                    pose = json.loads(payload)[0]
                    if pose['mid'] != mid or pose['score'] != score:
                        raise ValueError('Stored pose identity/score mismatch')
                    record = registry.execute('SELECT source_path,source_record_index,content_sha256 FROM conformer WHERE id=? AND molecule_id=?',
                                              (pose['cid'], mid)).fetchone()
                    if record is None:
                        raise ValueError('Selected conformer absent from original registry: ' + mid)
                    db.execute('INSERT INTO selected VALUES (?,?,?,?,?,?,?,?,?)',
                               (target, mid, pose['cid'], score, reason, payload, *record))
            print(json.dumps(dict(target=target, selected=db.execute('SELECT COUNT(*) FROM selected WHERE target=?', (target,)).fetchone()[0])), flush=True)
    finally:
        source_db.close(); registry.close()
    writers = {t: BatchWriter(output/t, batch_size) for t in TARGETS}
    manifests = {t: (output/t/'molecules.csv').open('w', encoding='utf-8', newline='') for t in TARGETS}
    fields = ['alias','molecule_id','conformer_id','score','selection_reason','batch','record_in_batch',
              'source_path','source_record_index','source_sha256','shape','color','feature_match','envelope',
              'minimum_distance','transform','retained_conformer_ids']
    csvs = {t: csv.writer(f) for t, f in manifests.items()}
    for writer in csvs.values():
        writer.writerow(fields)
    rejected = {t: 0 for t in TARGETS}; seen = 0
    recheck=source_pose_checker(batch,prepared) if report.get('recheck_source_pose_required') else None
    qc_file = (output/'rejected.csv').open('w', encoding='utf-8', newline='')
    qc = csv.writer(qc_file); qc.writerow(['target','mid','cid','reason'])
    try:
        sources = [r[0] for r in db.execute('SELECT DISTINCT source FROM selected ORDER BY source')]
        for source in sources:
            pending = {}
            for row in db.execute('SELECT * FROM selected WHERE source=? ORDER BY record_index,target', (source,)):
                pending.setdefault(row[7], []).append(row)
            for index, raw in iter_mol2_blocks(Path(source)):
                if index not in pending:
                    continue
                records = pending.pop(index)
                digest = hashlib.sha256(raw.encode('utf-8')).hexdigest()
                if any(r[8] != digest for r in records):
                    raise ValueError('Original MOL2 record hash mismatch: ' + source + ':' + str(index))
                parse_mol2_block(Path(source), index, raw)
                reason = None
                try:
                    molecule, mode = load_rdkit_mol2(raw, records[0][2])
                    if mode != 'strict':
                        reason = 'non_strict_sanitization'
                    elif any(not math.isfinite(v) for point in molecule.GetConformer().GetPositions() for v in point):
                        reason = 'nonfinite_coordinates'
                except Exception as exc:
                    reason = 'parse_failure:' + type(exc).__name__
                for target, mid, cid, score, selected_by, payload, _, _, _ in records:
                    seen += 1
                    poses = json.loads(payload); pose = poses[0]
                    target_reason=reason
                    if target_reason is None and recheck:
                        valid,_=recheck(molecule,target,pose)
                        if not valid:target_reason='exact_source_pose_failed_contacts_or_geometry'
                    if target_reason:
                        rejected[target] += 1; qc.writerow([target, mid, cid, target_reason]); continue
                    alias = mid + '__' + cid
                    batch_id, record_id = writers[target].add(rename_record(raw, alias))
                    csvs[target].writerow([alias,mid,cid,score,selected_by,batch_id,record_id,source,index,digest,
                        pose['shape'],pose['color'],pose['feature_match'],pose['envelope'],pose['minimum_distance'],
                        json.dumps(pose['transform']),json.dumps([p['cid'] for p in poses])])
                if not pending:
                    break
            if pending:
                raise ValueError('Source missing selected records: ' + source)
            for f in manifests.values(): f.flush()
            qc_file.flush()
            print(json.dumps(dict(source_complete=source, checked_target_records=seen,
                                  exported={t:w.count for t,w in writers.items()})), flush=True)
    finally:
        for w in writers.values(): w.close()
        for f in manifests.values(): f.close()
        qc_file.close(); db.close()
    targets = {t: dict(exported_molecules=w.count, requested=count, shortfall=count-w.count,
                      rejected=rejected[t], batches=w.batches,
                      manifest_sha256=sha(output/t/'molecules.csv')) for t,w in writers.items()}
    result = dict(status='complete', scope='exploratory_docking_input_not_activity_hits', targets=targets,
                  exact_source_pose_rechecked=bool(recheck),
                  source_scope=report.get('scope'), source_processed=report.get('processed'),
                  source_scheduled=report.get('scheduled'), source_full_catalog=report.get('full_catalog'),
                  elapsed_seconds=time.monotonic()-started, quota_met=all(w.count==count for w in writers.values()),
                  limitations=['Molecule IDs are library identities, not canonical chemical deduplication.',
                               'Original coordinates, charges and atom/bond records preserved; only title changed.',
                               'No protonation, tautomer enumeration, minimization or docking performed.',
                               'Strict parse failures are reported, not replaced with unevaluated molecules.',
                               'One scored conformer per molecule; original coordinates are not receptor-aligned.'])
    save(output/'report.json',result)
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--ranking', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--count', type=int, default=1000000)
    p.add_argument('--batch-size', type=int, default=100000)
    p.add_argument('--allow-pilot', action='store_true', help='Only for a separately labelled diagnostic export')
    p.add_argument('--allow-budgeted', action='store_true', help='Export an explicitly labelled completed compute-budget snapshot')
    args = p.parse_args()
    print(json.dumps(export(args.ranking,args.output,args.count,args.batch_size,args.allow_pilot,args.allow_budgeted),indent=2),flush=True)


if __name__ == '__main__':
    main()
