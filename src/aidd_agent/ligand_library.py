"""Trusted ligand-led full-catalog Gaussian search and original MOL2 handoff."""
import hashlib
import json
from pathlib import Path
import sqlite3
import time
import uuid

import numpy as np

from .expanded_wee1 import fingerprint
from .gaussian_batch import _atomic_json
from .library_acceptance import read, sha
from .project_context import ensure_within
from .screening_selection import check_hashes

OBJECTIVE = 'atomcentered_unweighted_joint'


def validate(args):
    if not isinstance(args, dict) or args.get('operation') not in {'inspect', 'start'}:
        raise ValueError('ligand_library operation must be inspect or start')
    allowed = {'operation', 'task_id'} if args['operation'] == 'inspect' else {
        'operation', 'task_id', 'query_id', 'attachment_id', 'export_molecules'}
    if set(args) - allowed:
        raise ValueError('Unknown ligand_library fields')
    if args['operation'] == 'inspect':
        if not isinstance(args.get('task_id'), str): raise ValueError('inspect requires task_id')
        return
    if ('task_id' in args) == ('attachment_id' in args):
        raise ValueError('Choose one owned task query OR one uploaded 3D ligand')
    if 'task_id' in args and (not isinstance(args['task_id'], str) or not isinstance(args.get('query_id'), str)):
        raise ValueError('Choose an exact query_id returned by inspect for the source task')
    if 'attachment_id' in args and ('query_id' in args or not isinstance(args['attachment_id'], str)):
        raise ValueError('Uploaded ligand uses attachment_id only')
    count = args.get('export_molecules')
    if type(count) is not int or not 1 <= count <= 1000000:
        raise ValueError('Specify export_molecules from 1 to 1000000; this caps delivery only')


def queries(app, sid, task_id):
    from .prompt_workflow import project_root
    ctx = app.context_for(sid)
    project = project_root(Path(ctx['db']), ctx['user_id'], ctx['project_id'])
    job = app.task(sid, task_id)
    if job['status'] != 'complete' or not job.get('report'):
        raise ValueError('Choose a completed structure/query task')
    report_path = ensure_within(Path(job['report']), project)
    result = []
    for step in read(report_path).get('steps', {}).values():
        if step.get('status') != 'complete' or step.get('action') not in {
                'structure_consensus', 'pocket_consensus', 'target_query_prepare', 'consensus_design', 'research_task'}:
            continue
        path = ensure_within(Path(step['result']['report']), project)
        child = read(path)
        if child.get('sources'):check_hashes(child['sources'])
        elif child.get('kind')!='molecule_prepare':raise ValueError('Prepared query has no source provenance')
        rows = child.get('prepared_complexes', child.get('templates', []))
        if child.get('kind') == 'target_query': rows = [child['query']]
        if child.get('kind') == 'molecule_prepare':
            check_hashes(child['outputs'])
            rows = child['prepared_queries']
        for row in rows:
            if not row.get('query_npz') or not row.get('query_id'): continue
            query = ensure_within(Path(row['query_npz']), project)
            manifest = query.with_suffix('.manifest.json')
            result.append(dict(query_id=row['query_id'], query_npz=str(query),
                               sources=fingerprint([report_path, path, query, manifest])))
    if len({r['query_id'] for r in result}) != len(result):
        raise ValueError('Ambiguous query IDs; select a task with one unambiguous prepared query set')
    return result


def handle(app, sid, args):
    from .chat_uploads import attachment
    from .prompt_workflow import project_root, create_plan
    validate(args); app.session(sid)
    if args['operation'] == 'inspect': return dict(queries=queries(app, sid, args['task_id']))
    ctx = app.context_for(sid); project = project_root(Path(ctx['db']), ctx['user_id'], ctx['project_id'])
    request = dict(export_molecules=args['export_molecules'])
    if 'attachment_id' in args:
        record = attachment(app, sid, args['attachment_id'])
        if Path(record['path']).suffix.lower() not in {'.sdf', '.mol2'}:
            raise ValueError('Use a single 3D ligand in SDF or MOL2; 2D/SMILES need a separate conformer preparation step')
        request.update(ligand=record['path'], query_id='upload-'+record['id'],
                       input_hashes={record['path']: record['sha256']})
    else:
        candidates = [r for r in queries(app, sid, args['task_id']) if r['query_id'] == args['query_id']]
        if len(candidates) != 1: raise ValueError('query_id is absent from the selected task')
        row = candidates[0]
        request.update(query_npz=row['query_npz'], query_id=row['query_id'], input_hashes=row['sources'])
    with app.lock:
        identifier = uuid.uuid4().hex[:24]
        folder = project/'ligand-library'/identifier; folder.mkdir(parents=True)
        path = folder/'request.json'; _atomic_json(path, request)
        plan = dict(version=1, summary='Full-library ligand 3D search and MOL2 export', clarifications=[],
                    steps=[dict(id='search_export', action='ligand_library_search',
                                params=dict(request_id=identifier, request_sha256=sha(path)))])
        plan_path = create_plan(Path(ctx['db']), ctx['user_id'], ctx['project_id'], json.dumps(args), local_plan=plan)
        jid = uuid.uuid4().hex[:16]
        with app.connect() as db:
            db.execute('INSERT INTO jobs(id,session,request,provider,profile,status,plan,report,log,created) VALUES(?,?,?,?,?,?,?,?,?,?)',
                       (jid, sid, plan['summary'], 'local-coordinator', '', 'queued', str(plan_path),
                        str(plan_path.parent/'execution/report.json'), str(app.root/(jid+'.log')), time.time()))
        return dict(status='queued', task_id=jid, full_library_requested=True,
                    query_id=request['query_id'], export_molecules=request['export_molecules'])


def request_file(execution, params):
    project = execution.parent.parent.parent
    path = ensure_within(project/'ligand-library'/params['request_id']/'request.json', project)
    if sha(path) != params['request_sha256']: raise ValueError('Ligand search request changed')
    return path


def prepare_ligand(path, output, query_id):
    from rdkit import RDConfig
    from rdkit.Chem import ChemicalFeatures
    from .target_query import load_ligand
    from .conformer_artifacts import FEATURE_TYPES
    from .chemical_companion import _feature_direction
    from .gaussian_batch import write_gaussian_query
    mol = load_ligand(path)
    xyz = mol.GetConformer().GetPositions()
    shape = xyz[[a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]]
    if not mol.GetConformer().Is3D() or np.linalg.norm(np.std(shape, axis=0)) < 1e-5:
        raise ValueError('Upload a nondegenerate 3D conformer; no coordinates are generated implicitly')
    factory = ChemicalFeatures.BuildFeatureFactory(str(Path(RDConfig.RDDataDir)/'BaseFeatures.fdef'))
    features = [f for f in factory.GetFeaturesForMol(mol) if f.GetFamily() in FEATURE_TYPES]
    if not features: raise ValueError('Ligand has no supported 3D chemical features')
    directions = [_feature_direction(mol, f.GetFamily(), tuple(f.GetAtomIds())) for f in features]
    write_gaussian_query(output, shape_points=shape,
        feature_points=[xyz[list(f.GetAtomIds())].mean(axis=0) for f in features],
        feature_types=[FEATURE_TYPES[f.GetFamily()] for f in features],
        feature_directions=[d[0] for d in directions], feature_direction_kinds=[d[1] for d in directions],
        anchored_weights=np.ones(len(features)), anchor_feature_indices=[],
        source=dict(query_id=query_id, ligand=str(path), provenance='uploaded_3d_ligand_no_receptor'))


def rank_chunks(stage, total, chunk_size, ranking, query_id):
    """Keep best actual pose per molecule; ties use global conformer ID."""
    with sqlite3.connect(ranking) as db:
        db.executescript('DROP TABLE IF EXISTS best; DROP TABLE IF EXISTS ranking; DROP TABLE IF EXISTS poses; DROP TABLE IF EXISTS template_ranks; '
            'CREATE TABLE best(mid TEXT PRIMARY KEY,gid INTEGER,score REAL,payload TEXT); '
            'CREATE TABLE ranking(rank INTEGER PRIMARY KEY,mid TEXT,rrf REAL,template_support INTEGER); '
            'CREATE TABLE poses(mid TEXT,template TEXT,payload TEXT); '
            'CREATE TABLE template_ranks(mid TEXT,template TEXT,rank INTEGER);')
        for index, start in enumerate(range(0, total, chunk_size)):
            with np.load(stage/'chunks'/f'chunk-{index:06d}.npz', allow_pickle=False) as arrays:
                if not np.array_equal(arrays['global_ids'], np.arange(start, min(start+chunk_size, total))):
                    raise ValueError('Full-library chunk coverage mismatch')
                for i, gid in enumerate(arrays['global_ids']):
                    mid = str(arrays['molecule_ids'][i]); score = float(arrays[OBJECTIVE+'__objective'][i])
                    transform = arrays[OBJECTIVE+'__transform'][i]
                    if not np.isfinite(score) or not np.isfinite(transform).all():
                        raise ValueError('Nonfinite candidate score or transform')
                    payload = json.dumps(dict(global_id=int(gid), molecule_id=mid, template=query_id,
                        gaussian_score=score, objective=OBJECTIVE, transform=transform.tolist()))
                    db.execute('INSERT INTO best VALUES(?,?,?,?) ON CONFLICT(mid) DO UPDATE SET gid=excluded.gid,score=excluded.score,payload=excluded.payload '
                        'WHERE excluded.score>best.score OR (excluded.score=best.score AND excluded.gid<best.gid)',
                        (mid, int(gid), score, payload))
            db.commit()
        for rank, (mid, payload) in enumerate(db.execute('SELECT mid,payload FROM best ORDER BY score DESC,mid'), 1):
            db.execute('INSERT INTO ranking VALUES(?,?,?,1)', (rank, mid, 1/(60+rank)))
            db.execute('INSERT INTO poses VALUES(?,?,?)', (mid, query_id, payload))
            db.execute('INSERT INTO template_ranks VALUES(?,?,?)', (mid, query_id, rank))
        db.execute('CREATE INDEX pose_lookup ON poses(mid,template)')
        db.execute('CREATE INDEX template_lookup ON template_ranks(mid,template)')
        return db.execute('SELECT COUNT(*) FROM ranking').fetchone()[0]


def execute(params, directory, execution, cfg, allow_compute):
    from .workflow_errors import Blocked
    from .gaussian_batch import _execute_stage
    from .full_library_screen import catalog_count
    from .budget_export import export
    from .prompt_workflow import file_lock
    if not allow_compute: raise Blocked('Enable --allow-compute for full-library scoring and MOL2 export')
    search = cfg.get('search', {})
    if not search.get('batch'): raise Blocked('Configure search.batch for the local molecular library')
    path = request_file(execution, params); request = read(path); check_hashes(request['input_hashes'])
    batch = Path(search['batch']).resolve(); catalog = batch/'artifacts/catalog.json'
    catalog_data = read(catalog); chemistry = read(batch/'chemical/catalog.json')
    ranges = lambda data: sorted((row['global_id_start'], row['conformers']) for row in data['shards'])
    if (chemistry.get('library_id') != catalog_data.get('library_id') or
            ranges(chemistry) != ranges(catalog_data) or chemistry.get('artifact_v1_catalog_sha256') != sha(catalog)):
        raise ValueError('Artifact and chemistry catalogs identify different libraries or counts')
    query = Path(request['query_npz']) if request.get('query_npz') else directory/'query.npz'
    if request.get('ligand') and (not query.exists() or not query.with_suffix('.manifest.json').exists()):
        prepare_ligand(request['ligand'], query, request['query_id'])
    workers, chunk_size = search.get('workers', 1), search.get('refine_chunk', 64)
    if type(workers) is not int or workers < 1 or type(chunk_size) is not int or chunk_size < 1:
        raise ValueError('Positive integer workers/refine_chunk required')
    sources = {**request['input_hashes'], **fingerprint([path, query, query.with_suffix('.manifest.json'),
        catalog, batch/'chemical/catalog.json', batch/'registry.sqlite3'])}
    from rdkit import rdBase
    protocol = dict(sources=sources, code=fingerprint(sorted(Path(__file__).parent.glob('*.py'))),
        numpy=np.__version__, rdkit=rdBase.rdkitVersion,
        workers=workers, chunk_size=chunk_size, objective=OBJECTIVE,
        parameters=dict(sigma=1., cutoff=4.5, pair_tolerance=2., axial_samples=6,
                        max_pair_seeds=512, bounded_pair_seeds=True, backend='numpy'))
    run = directory/'search'; run.mkdir(exist_ok=True)
    with file_lock(run/'run.lock'):
        saved = run/'protocol.json'
        if saved.exists() and read(saved) != protocol: raise ValueError('Search input/configuration changed; start a fresh task')
        _atomic_json(saved, protocol)
        total = catalog_count(catalog_data)
        if total < 1: raise ValueError('The configured library has no conformers')
        ids = np.arange(total, dtype=np.int64)
        config_hash = hashlib.sha256(json.dumps(protocol, sort_keys=True).encode()).hexdigest()
        _, stage = _execute_stage(run/'scores', 'full_library', ids, artifact_catalog=catalog,
            query_path=query, parameters=protocol['parameters'], config_hash=config_hash,
            workers=workers, chunk_size=chunk_size, resume=True, progress_every=1, merge=False)
        check_hashes(sources)
        if (run/'report.json').exists():
            ranked = read(run/'report.json'); check_hashes(ranked['outputs'])
            count = ranked['ranked_molecules']
        else:
            count = rank_chunks(run/'scores', total, chunk_size, run/'ranking.sqlite', request['query_id'])
            from .budget_screen import readonly
            with readonly(batch/'registry.sqlite3') as registry:
                expected = registry.execute('SELECT COUNT(*) FROM molecule WHERE library_id=?',
                                            (catalog_data['library_id'],)).fetchone()[0]
            if count != expected: raise ValueError('Whole-library molecule coverage differs from the registry')
            _atomic_json(run/'report.json', dict(status='complete', full_coverage=True,
                processed_conformers=total, ranked_molecules=count, sources=sources,
                ranking='Best Gaussian shape/color score per molecule; single-template RRF is monotonic in rank',
                outputs=fingerprint([run/'ranking.sqlite'])))
        page = directory/'mol2'; delivery = export(batch, run, page, 1, request['export_molecules'])
        result = dict(kind='ligand_library_search', status='complete', query_id=request['query_id'],
            full_coverage=True, processed_conformers=total, ranked_molecules=count,
            exported_molecules=delivery['exported_molecules'], shortfall=delivery['shortfall'],
            report=str(directory/'report.json'), ranking=str(run/'ranking.sqlite'),
            original_mol2_directory=str(page/'original'), posed_mol2_directory=str(page/'posed'),
            names=str(page/'molecules.csv'), evidence=str(page/'molecules.jsonl'), stage=stage,
            sources=sources, export_manifest=str(page/'report.json'),
            outputs=fingerprint([page/'report.json', page/'molecules.csv', page/'molecules.jsonl', run/'ranking.sqlite']),
            limitations=['Ligand 3D shape/color similarity, not receptor contact screening or predicted affinity',
                         'All stored conformers scored with finite rigid-pose seeds; no new internal conformers generated',
                         'Output cap selects candidates for review; it is not an activity threshold'])
        _atomic_json(directory/'report.json', result)
        return result
