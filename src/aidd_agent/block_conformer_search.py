"""Contact-first refinement restricted to frozen block members, retaining conformers."""
from collections import defaultdict
import csv
import hashlib
import json
import math
from pathlib import Path
import sqlite3

import numpy as np

from .block_sampling import check_outputs, readonly, snapshot
from .expanded_wee1 import fingerprint
from .final_work_blocks import read, sha
from .joint_spatial_profiles import save
from .screening_selection import check_hashes


def profile_paths(profile):
    profile = Path(profile)
    return {k:(Path(v) if Path(v).is_absolute() else profile.parent/v).resolve()
            for k,v in read(profile).items()}


def members(profile, selection, sample_report, output):
    """Bind membership to the exact databases used for the original score panel."""
    paths = profile_paths(profile)
    check_outputs(Path(sample_report).parent)
    signature = read(Path(sample_report).parent/'signature.json')
    check_hashes(signature['sources'])
    for old in signature['database_snapshots']:
        if snapshot(old['path']) != old:
            raise ValueError('Original block population database changed')
    scheme = selection['scheme']
    needed = [paths['profiles']/'profiles.sqlite', paths['descriptors']/'descriptors.sqlite']
    if scheme != 'E094':
        needed.append(paths[scheme]/'property_blocks.sqlite')
    recorded = {r['path'] for r in signature['database_snapshots']}
    if not all(str(p) in recorded for p in needed):
        raise ValueError('Configured block profile differs from original scoring panel')
    with sqlite3.connect(output) as db:
        db.executescript('''CREATE TABLE member(cid TEXT PRIMARY KEY, mid TEXT, file_id INTEGER,
            idx INTEGER, block_id TEXT, artifact_cid TEXT UNIQUE, gid INTEGER UNIQUE, raw_sha TEXT);
            CREATE TABLE files(id INTEGER PRIMARY KEY,path TEXT);
            CREATE INDEX member_file ON member(file_id,idx);''')
        with readonly(paths['profiles']/'profiles.sqlite') as source:
            if scheme != 'E094':
                source.execute('ATTACH DATABASE ? AS partition', ((paths[scheme]/'property_blocks.sqlite').as_uri()+'?mode=ro',))
            db.executemany('INSERT INTO files VALUES(?,?)', source.execute('SELECT id,path FROM files'))
            for block in selection['block_ids']:
                if scheme == 'E094':
                    rows = source.execute('SELECT cid,mid,file_id,idx,parent FROM item WHERE parent=?', (block,))
                    db.executemany('INSERT INTO member(cid,mid,file_id,idx,block_id) VALUES(?,?,?,?,?)', rows)
                else:
                    rows = source.execute('SELECT i.cid,i.mid,i.file_id,i.idx,m.block_id FROM partition.membership m JOIN item i USING(cid) WHERE m.block_id=?',(block,))
                    db.executemany('INSERT INTO member(cid,mid,file_id,idx,block_id) VALUES(?,?,?,?,?)',rows)
                n = db.execute('SELECT count(*) FROM member WHERE block_id=?',(block,)).fetchone()[0]
                with readonly(Path(sample_report).parent/'samples.sqlite') as old:
                    expected = old.execute('SELECT n FROM population WHERE scheme=? AND block_id=?',(scheme,block)).fetchone()
                if expected != (n,) or n == 0:
                    raise ValueError('Selected block population differs from original sampled panel')
        db.commit()
    return paths


def map_artifacts(batch, database, descriptors):
    """Map source record hashes to artifact IDs; never match by molecule name alone."""
    from . import library_acceptance as ev
    artifact = read(batch/'artifacts/catalog.json')
    chemical = read(batch/'chemical/catalog.json')
    if chemical['artifact_v1_catalog_sha256'] != sha(batch/'artifacts/catalog.json') or chemical['library_id'] != artifact['library_id']:
        raise ValueError('Chemical/artifact lineage mismatch')
    for kind, catalog in [('artifacts', artifact), ('chemical', chemical)]:
        for row in catalog['shards']:
            directory = ev.shard_path(row,batch/kind/'catalog.json')
            if sha(directory/'manifest.json') != row['manifest_sha256']:
                raise ValueError('Changed library shard manifest')
            ev.verify_manifest(directory, full=False)
    with sqlite3.connect(database) as db, readonly(batch/'registry.sqlite3') as registry, readonly(descriptors/'descriptors.sqlite') as desc:
        # One indexed source-file query at a time avoids a full in-memory library map.
        for fid, path in db.execute('SELECT id,path FROM files').fetchall():
            wanted = {idx:cid for cid,idx in db.execute('SELECT cid,idx FROM member WHERE file_id=?',(fid,))}
            if not wanted:
                continue
            print(f'Mapping {len(wanted):,} selected source conformers: {path}',flush=True)
            for aid,idx,digest in registry.execute('SELECT id,source_record_index,content_sha256 FROM conformer WHERE source_path=?',(path,)):
                cid = wanted.get(idx)
                if cid is None:
                    continue
                row = desc.execute('SELECT payload FROM descriptor WHERE cid=?',(cid,)).fetchone()
                if row is None or json.loads(row[0])['provenance']['content_sha256'] != digest:
                    raise ValueError('Source descriptor / library record hash mismatch')
                if db.execute('SELECT artifact_cid FROM member WHERE cid=?',(cid,)).fetchone()[0] is not None:
                    raise ValueError('Ambiguous source-to-artifact mapping')
                db.execute('UPDATE member SET artifact_cid=?,raw_sha=? WHERE cid=?',(aid,digest,cid))
        wanted = {r[0] for r in db.execute('SELECT artifact_cid FROM member WHERE artifact_cid IS NOT NULL')}
        for row in artifact['shards']:
            directory = ev.shard_path(row,batch/'artifacts/catalog.json')
            manifest = read(directory/'manifest.json')
            count, start = manifest['conformers'], manifest['global_id_start']
            ids = np.memmap(directory/'conformer_ids.bin', dtype='S16', mode='r', shape=(count,))
            db.executemany('UPDATE member SET gid=? WHERE artifact_cid=?',
                           ((start+i,bytes(v).decode()) for i,v in enumerate(ids) if bytes(v).decode() in wanted))
        missing = db.execute('SELECT count(*) FROM member WHERE gid IS NULL').fetchone()[0]
        total = db.execute('SELECT count(*) FROM member').fetchone()[0]
        db.commit()
        if missing:
            raise ValueError(f'Search library lacks {missing}/{total} selected source conformers. Build matching search artifacts before continuing; no outside-block fallback.')
        return np.array([r[0] for r in db.execute('SELECT gid FROM member ORDER BY gid')],dtype=np.int64)


def rank_conformers(output, design, database, count):
    """Fuse per-template ranks by conformer, with no molecule collapse."""
    with sqlite3.connect(database) as members_db, sqlite3.connect(output/'ranking.sqlite') as db:
        db.executescript('''DROP TABLE IF EXISTS poses; DROP TABLE IF EXISTS template_ranks;
            DROP TABLE IF EXISTS ranking;
            CREATE TABLE poses(gid INTEGER,template TEXT,contact REAL,gaussian REAL,payload TEXT,
                               PRIMARY KEY(gid,template));''')
        for ti, template in enumerate(design['templates']):
            for path in sorted((output/'chunks'/f'{ti:02d}').glob('*.poses.jsonl')):
                check_hashes(read(path.with_name(path.name.replace('.poses.jsonl','.receipt.json')))['files'])
                with path.open(encoding='utf-8') as stream:
                    for line in stream:
                        row = json.loads(line)
                        if type(row['global_id']) is not int or members_db.execute('SELECT 1 FROM member WHERE gid=?',(row['global_id'],)).fetchone() is None:
                            raise ValueError('Search returned a conformer outside the selected blocks')
                        if not all(math.isfinite(row[k]) for k in ('optional_score','gaussian_same_pose','composite_score')):
                            raise ValueError('Nonfinite search score')
                        transform = np.asarray(row['transform']).reshape(4,4)
                        if not np.isfinite(transform).all():
                            raise ValueError('Invalid search pose transform')
                        db.execute('''INSERT INTO poses VALUES(?,?,?,?,?) ON CONFLICT(gid,template) DO UPDATE SET
                            contact=excluded.contact,gaussian=excluded.gaussian,payload=excluded.payload
                            WHERE (excluded.contact,excluded.gaussian)>(poses.contact,poses.gaussian)''',
                                   (row['global_id'],template['query_id'],row['optional_score'],row['gaussian_same_pose'],line))
        db.execute('CREATE TABLE template_ranks AS SELECT gid,template,ROW_NUMBER() OVER(PARTITION BY template ORDER BY contact DESC,gaussian DESC,gid) rank FROM poses')
        db.execute('''CREATE TABLE ranking AS SELECT ROW_NUMBER() OVER(ORDER BY positive DESC,rrf DESC,gid) rank,gid,rrf,support FROM
            (SELECT t.gid,MAX(p.contact>0) positive,SUM(1.0/(60+t.rank)) rrf,COUNT(*) support
             FROM template_ranks t JOIN poses p USING(gid,template) GROUP BY t.gid)''')
        db.execute('CREATE UNIQUE INDEX ranked_gid ON ranking(gid)')
        db.commit()
        rows = []
        for rank,gid,rrf,support in db.execute('SELECT * FROM ranking WHERE rank<=? ORDER BY rank',(count,)):
            payload = db.execute('SELECT p.payload FROM poses p JOIN template_ranks t USING(gid,template) WHERE gid=? ORDER BY t.rank,template LIMIT 1',(gid,)).fetchone()[0]
            rows.append(dict(rank=rank,gid=gid,rrf=rrf,template_support=support,pose=json.loads(payload)))
        return rows, db.execute('SELECT count(*) FROM ranking').fetchone()[0]


def export_selected(output, database, rows, scheme):
    """Preserve original MOL2 chemistry and unique source conformers for PLANTS."""
    from .mol2 import iter_mol2_blocks
    output.mkdir(parents=True,exist_ok=True)
    identity = dict(membership_sha256=sha(database),scheme=scheme,
                    selection_sha256=hashlib.sha256(json.dumps(rows,sort_keys=True).encode()).hexdigest())
    if (output/'report.json').exists():
        previous = check_outputs(output)
        if previous.get('identity') != identity:
            raise ValueError('Export selection changed')
        return previous
    wanted = defaultdict(dict)
    with readonly(database) as db:
        for rank,row in enumerate(rows,1):
            cid,mid,fid,idx,block,aid,gid,digest = db.execute('SELECT * FROM member WHERE gid=?',(row['gid'],)).fetchone()
            path = db.execute('SELECT path FROM files WHERE id=?',(fid,)).fetchone()[0]
            row.update(cid=cid,molecule_id=mid,block_id=block,artifact_cid=aid,alias=f'c{rank:09d}',source_path=path,source_index=idx,source_sha256=digest)
            wanted[path][idx] = row
    chunks = []
    pending, metadata = [], []
    (output/'ligands').mkdir(exist_ok=True)
    def flush():
        if not pending:
            return
        name = f'ligands/{len(chunks):06d}.mol2'
        (output/name).write_text(''.join(pending),encoding='utf-8',newline='')
        chunks.append(dict(path=name,sha256=sha(output/name),count=len(pending),records=list(metadata)))
        pending.clear()
        metadata.clear()
    for path, records in wanted.items():
        for index, raw in iter_mol2_blocks(Path(path)):
            row = records.pop(index,None)
            if row is None:
                continue
            if hashlib.sha256(raw.encode()).hexdigest() != row['source_sha256']:
                raise ValueError('Original selected MOL2 record changed')
            lines = raw.splitlines(keepends=True)
            lines[1] = row['alias']+'\n'
            pending.append(''.join(lines).rstrip('\n')+'\n')
            metadata.append(dict(cid=row['cid'],mid=row['molecule_id'],alias=row['alias'],source_index=index,source_sha256=row['source_sha256']))
            if len(pending) == 100:
                flush()
            if not records:
                break
        if records:
            raise ValueError('Selected source records are missing')
    flush()
    dbpath = output/'samples.sqlite'
    if dbpath.exists():
        dbpath.unlink()  # Named derived export only; sealed input membership is untouched.
    with sqlite3.connect(dbpath) as db:
        db.executescript('''CREATE TABLE sample(scheme TEXT,block_id TEXT,cid TEXT,PRIMARY KEY(scheme,block_id,cid));
            CREATE TABLE selected(cid TEXT PRIMARY KEY,alias TEXT UNIQUE,mid TEXT);
            CREATE TABLE population(scheme TEXT,block_id TEXT,n INTEGER);''')
        db.executemany('INSERT INTO sample VALUES(?,?,?)',((scheme,r['block_id'],r['cid']) for r in rows))
        db.executemany('INSERT INTO selected VALUES(?,?,?)',((r['cid'],r['alias'],r['molecule_id']) for r in rows))
        # These are enriched follow-up panels; n denotes selected panel size, not uniform sampling.
        db.execute('INSERT INTO population SELECT scheme,block_id,count(*) FROM sample GROUP BY scheme,block_id')
    save(output/'exports.json',dict(chunks=chunks))
    with (output/'candidates.jsonl').open('w',encoding='utf-8') as stream:
        for row in rows:
            stream.write(json.dumps(row)+'\n')
    names = ['samples.sqlite','exports.json','candidates.jsonl',*[c['path'] for c in chunks]]
    report = dict(kind='block_search_export',status='complete',identity=identity,unique_conformers=len(rows),
                  unique_molecules=len({r['molecule_id'] for r in rows}),candidate_unit='conformer',
                  output_hashes={n:sha(output/n) for n in names},
                  limitations=['Search-enriched panel, not uniform block sampling or population estimation.',
                               'Original MOL2 chemistry retained; protonation review remains separate.'])
    save(output/'report.json',report)
    return report


def docking_candidates(output, selected, docking):
    """Save per-conformer docking priority, without claiming uniform block estimates."""
    with readonly(selected/'samples.sqlite') as db:
        identity = {cid:(mid,bid) for cid,mid,bid in db.execute('SELECT cid,mid,block_id FROM selected JOIN sample USING(cid)')}
    rows = []
    with (docking/'scores.csv').open(encoding='utf-8',newline='') as stream:
        for row in csv.DictReader(stream):
            if row['cid'] not in identity:
                raise ValueError('Docking score references an unselected conformer')
            mid,bid = identity[row['cid']]
            row.update(molecule_id=mid,block_id=bid)
            if row['status'] == 'ok':
                row['score'] = float(row['score'])
                if not math.isfinite(row['score']):
                    raise ValueError('Nonfinite docking score')
            rows.append(row)
    rows.sort(key=lambda r:(r['status'] != 'ok',r['score'] if r['status']=='ok' else 0,r['cid'],r['receptor']))
    for i,row in enumerate(rows,1):
        row['rank'] = i if row['status']=='ok' else ''
    if rows:
        with (output/'docking_candidates.csv').open('w',encoding='utf-8',newline='') as stream:
            writer = csv.DictWriter(stream,fieldnames=list(rows[0]))
            writer.writeheader();writer.writerows(rows)
    return [r for r in rows if r['status']=='ok'][:10]


def inherited_plants_profile(full, configured, receptor):
    """Reuse the scored panel's reviewed receptor/site and engine version."""
    from .block_plants import profile_read
    plants = Path(full['source_report'])
    original = read(plants)
    if sha(plants) != full['identity']['source_sha256']:
        raise ValueError('Original PLANTS source changed')
    signature = read(plants.parent/'signature.json')
    prepared = Path(original['prepared_report'])
    if sha(prepared) != signature['prepared_sha256']:
        raise ValueError('Original PLANTS preparation changed')
    check_outputs(prepared.parent)
    saved = read(prepared.parent/'profile.json')
    saved['receptors'] = [r for r in saved['receptors'] if r['id']==receptor]
    if len(saved['receptors']) != 1:
        raise ValueError('Selected receptor is absent from the original scored panel')
    configured = profile_read(configured)
    if sha(Path(configured['executable'])) != signature['executable_sha256']:
        raise ValueError('Configured PLANTS engine differs from the original panel; review a new protocol')
    saved['executable'] = configured['executable']
    saved['workers'] = configured['workers']
    saved['timeout_seconds'] = configured['timeout_seconds']
    saved['receptors'][0]['mol2'] = str(prepared.parent/'receptors'/(receptor+'.mol2'))
    return saved, [plants,plants.parent/'signature.json',prepared,prepared.parent/'profile.json',
                   Path(saved['receptors'][0]['mol2']),Path(saved['executable'])]


def run(request, output, cfg):
    from .prompt_workflow import Blocked
    from .budget_screen import run_refinement
    from .consensus_design import adopt
    from .contact_policy import require_protein_contacts, migrate
    from .block_plants import prepare, run as dock
    search, block = cfg.get('search',{}), cfg.get('block_evaluation',{})
    if not search.get('batch') or not block.get('sampling_profile') or not block.get('plants_profile'):
        raise Blocked('Configure search.batch, block_evaluation.sampling_profile and plants_profile')
    batch = Path(search['batch'])
    analysis = read(request['analysis'])
    full = read(analysis['source_report'])
    profile, inherited_inputs = inherited_plants_profile(full,block['plants_profile'],request['selection']['receptor'])
    if sha(Path(analysis['source_report'])) != analysis['ranking']['source_report_sha256']:
        raise ValueError('Original full score report changed')
    if full['identity']['profile_inputs'].get(full['sample_report']) != sha(Path(full['sample_report'])):
        raise ValueError('Original sample report changed after EquiScore')
    sample_signature = read(Path(full['sample_report']).parent/'signature.json')
    check_hashes(sample_signature['sources'])
    for record in sample_signature['database_snapshots']:
        if snapshot(record['path']) != record:
            raise ValueError('Original block population database changed')
    paths = [Path(request['analysis']),Path(request['query']),Path(full['sample_report']),
             Path(block['sampling_profile']),Path(block['plants_profile']),
             batch/'artifacts/catalog.json',batch/'chemical/catalog.json',*inherited_inputs]
    protocol = dict(request=request,inputs=fingerprint(paths),code=fingerprint(sorted(Path(__file__).parent.glob('*.py'))),
                    library_registry=snapshot(batch/'registry.sqlite3'),
                    ranking='contact-first per-template ranks fused by conformer RRF; positive-contact tier first')
    if (output/'protocol.json').exists() and read(output/'protocol.json') != protocol:
        raise ValueError('Search inputs changed; use a new campaign')
    save(output/'protocol.json',protocol)
    if (output/'report.json').exists() and read(output/'report.json').get('status') == 'complete':
        return check_outputs(output)
    design_dir = output/'query'
    if not (design_dir/'report.json').exists():
        source = read(request['query'])
        check_hashes(source['sources'])
        if source['kind'] == 'consensus_recommendation':
            if not (output/'protein-only-proposal/report.json').exists():
                migrate(Path(request['query']),output/'protein-only-proposal')
            adopt(output/'protein-only-proposal/report.json',design_dir)
        else:
            require_protein_contacts(source['design'],source['anchors'])
            design_dir.mkdir(exist_ok=True)
            save(design_dir/'report.json',source)
    design = read(design_dir/'report.json')
    check_hashes(design['sources'])
    if design.get('failed_template_self_controls') or design.get('readiness') == 'needs_template_state_review':
        raise Blocked('The selected query requires template-state review before block search')
    # MDM2 is this campaign's explicit target; never use the legacy WEE1 shortcuts.
    if design['target'].get('accession') != 'Q00987':
        raise ValueError('This continuation requires the confirmed human MDM2 query (Q00987)')
    database = output/'members.sqlite'
    marker = output/'members-seal.json'
    if marker.exists():
        check_hashes(read(marker))
        with readonly(database) as db:
            ids = np.array([r[0] for r in db.execute('SELECT gid FROM member ORDER BY gid')],dtype=np.int64)
    else:
        if database.exists():
            database.unlink()  # Rebuild only an unfinished derived mapping.
        save(output/'progress.json',dict(stage='mapping_block_members',selection=request['selection']))
        print('Binding selected block membership to the original panel and search library.',flush=True)
        locations = members(block['sampling_profile'],request['selection'],full['sample_report'],database)
        ids = map_artifacts(batch,database,locations['descriptors'])
        save(marker,fingerprint([database]))
    save(output/'progress.json',dict(stage='3d_refinement',selected_block_conformers=len(ids),requested_conformers=request['conformers']))
    print(f'Refining {len(ids):,} selected-block conformers; retaining up to {request["conformers"]:,} conformers, not molecules.',flush=True)
    run_refinement(batch,output,design,ids,search['workers'],search['refine_chunk'])
    rows, ranked = rank_conformers(output,design,database,request['conformers'])
    if not rows:
        raise ValueError('No qualified 3D search poses; docking was not submitted')
    exported = export_selected(output/'selected',database,rows,request['selection']['scheme'])
    check_hashes(protocol['inputs'])
    for record in sample_signature['database_snapshots']:
        if snapshot(record['path']) != record:
            raise ValueError('Original block population changed during search')
    outputs = ['protocol.json','members.sqlite','members-seal.json','ranking.sqlite','selected/report.json','query/report.json']
    results = dict(kind='block_search_dock',status='complete',selection=request['selection'],query=str(design_dir/'report.json'),
                   requested_conformers=request['conformers'],searched_conformers=len(ids),ranked_conformers=ranked,
                   exported_conformers=exported['unique_conformers'],unique_molecules=exported['unique_molecules'],
                   shortfall=request['conformers']-len(rows),receptors=1,planned_pairs=len(rows),
                   candidate_unit='conformer',docking_requested=request['dock'],biological_validation='not_run')
    if request['dock']:
        save(output/'plants-profile.json',profile)
        prepare(output/'selected/report.json',output/'plants-profile.json',output/'prepare')
        save(output/'progress.json',dict(stage='docking',planned_pairs=len(rows)))
        result = dock(output/'prepare/report.json',output/'docking')
        results.update(status=result['status'],scored=result['scored'],failed_jobs=result['failed_jobs'],
                       docking_report=str(output/'docking/report.json'),scores=str(output/'docking/scores.csv'))
        results['top_candidates'] = docking_candidates(output,output/'selected',output/'docking')
        outputs.extend(['plants-profile.json','prepare/report.json','docking/report.json','docking/scores.csv','docking_candidates.csv'])
    results['outputs'] = {n:str(output/n) for n in outputs}
    results['output_hashes'] = {n:sha(output/n) for n in outputs}
    results['limitations'] = ['Block priorities and contact-first retrieval are exploratory, not affinity or experimental enrichment.',
                             'Follow-up conformers are search-selected; do not report population-weighted block estimates.',
                             'Shortfalls are retained explicitly; no outside-block conformers are added.']
    save(output/'report.json',results)
    if results['status'] != 'complete':
        raise RuntimeError('Follow-up docking incomplete; inspect the saved docking report before resuming')
    return results
