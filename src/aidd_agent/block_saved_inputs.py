"""Discover saved CLI evidence and adopt the original multi-cocrystal query."""
import hashlib
from pathlib import Path

from .final_work_blocks import read, sha
from .joint_spatial_profiles import save
from .project_context import ensure_within


def query_header(path):
    path = Path(path).resolve()
    report = read(path)
    if (report.get('kind') != 'consensus_design' or report.get('status') != 'complete'
            or report.get('target', {}).get('accession') != 'Q00987'
            or not report.get('templates')):
        raise ValueError('Select the saved MDM2 adopted-design/report.json from the original multi-cocrystal search')
    return path, report


def discover(app, sid, project):
    """Read bounded known result locations; never infer absence from an empty session."""
    paths = set()
    for pattern in ('*/execution/equiscore-analysis/report.json',
                    '*/execution/*/equiscore-analysis/report.json',
                    '*/execution/*/blocks/report.json'):
        paths.update((project/'runs').glob(pattern))
    with app.connect() as db:
        jobs = list(db.execute("SELECT plan,report FROM jobs WHERE status='complete'"))
    for job in jobs:
        try:
            path = ensure_within(Path(job['report']), project)
            for step in read(path).get('steps', {}).values():
                if step.get('status') == 'complete' and step.get('result', {}).get('report'):
                    paths.add(ensure_within(Path(step['result']['report']), project))
        except (OSError, ValueError, TypeError, KeyError):
            continue
    roots = [project]
    if app.runtime:
        cfg = read(app.runtime)
        batch = cfg.get('search', {}).get('batch')
        if batch:
            batch = Path(batch)
            roots.append((batch if batch.is_absolute() else Path(app.runtime).parent/batch).resolve().parent)
    for root in roots:
        for pattern in ('*/adopted-design/report.json', '*/search/adopted-design/report.json'):
            paths.update(root.glob(pattern))
    folder = project/'block-input-discovery'/sid
    folder.mkdir(parents=True, exist_ok=True)
    candidates = []
    for path in sorted(paths)[:300]:
        try:
            if path.stat().st_size > 32*1024*1024:
                continue
            r = read(path)
            kind = r.get('kind')
            if kind == 'consensus_design':
                _, r = query_header(path)
                category = 'query'
                details = dict(target=r['target'],query_ids=[q['query_id'] for q in r['templates']],
                               templates=len(r['templates']),readiness=r.get('readiness'),
                               failed_template_self_controls=r.get('failed_template_self_controls',[]))
            elif kind == 'block_analyze':
                from .block_campaign import analysis_header
                _, r = analysis_header(path, project)
                if r.get('adoption'):
                    continue  # Original saved analysis, not copies created by this importer.
                category, details = 'scores', dict(ranking=r['ranking'])
            else:
                continue
            digest = sha(path)
            identifier = hashlib.sha256((str(path.resolve())+digest).encode()).hexdigest()[:24]
            row = dict(candidate_id=identifier,category=category,report=str(path.resolve()),
                       report_sha256=digest,session=sid,details=details)
            save(folder/(identifier+'.json'),row)
            candidates.append(row)
        except (OSError, ValueError, KeyError, TypeError):
            continue
    return dict(candidates=candidates,computation_launched=False,
                scope='Active Project results and adopted query packages next to the configured search library; bounded known layouts, not an exhaustive disk scan.',
                next_step='Import matching saved scores and the original query package by candidate_id. If absent, request its saved adopted-design/report.json path, not a nonexistent Chat task ID. Do not create a replacement query.')


def resolve(args, sid, project, category):
    if args.get('report'):
        return Path(args['report']).resolve()
    entry = ensure_within(project/'block-input-discovery'/sid/(args['candidate_id']+'.json'),project)
    row = read(entry)
    if row.get('session') != sid or row.get('category') != category:
        raise ValueError('Select the matching saved input discovered in this conversation')
    path = Path(row['report'])
    if sha(path) != row['report_sha256']:
        raise ValueError('Discovered input changed; discover its final version again')
    return path


def adopt_query(request, output):
    from .expanded_wee1 import fingerprint
    from .screening_selection import check_hashes
    from .gaussian_batch import _load_query
    path, report = query_header(request['saved_query'])
    if sha(path) != request['saved_query_sha256']:
        raise ValueError('Original query changed after import submission')
    check_hashes(report['sources'])
    files = [path,Path(report['consensus_npz'])]
    files.extend(Path(q['query_npz']) for q in report['templates'])
    for query in files[1:]:
        if not query.is_absolute():
            raise ValueError('Saved query assets must retain their original absolute paths')
        _load_query(query)
    files.extend(p.with_suffix('.manifest.json') for p in list(files[1:]))
    references = []
    for template in report['templates']:
        aligned = Path(template['query_npz'])
        aligned_manifest, _ = _load_query(aligned)
        native = aligned.parent/'native.npz'
        native_manifest, _ = _load_query(native)
        source = native_manifest['source']
        if source.get('query_id') != template['query_id'] or aligned_manifest['source'].get('query_id') != template['query_id']:
            raise ValueError('Co-crystal ligand identity disagrees with the saved template')
        original = {source[k]:source[k+'_sha256'] for k in ('mmcif','ccd','query_manifest')}
        check_hashes(original)
        fields = template['query_id'].split(':')
        if len(fields) != 4:
            raise ValueError('Co-crystal reference requires PDB:CCD:chain:residue identity')
        references.append(dict(query_id=template['query_id'],pdb_id=fields[0],ligand_ccd=fields[1],
            ligand_chain=fields[2],ligand_residue=fields[3],mmcif=source['mmcif'],ccd=source['ccd'],
            native_query=str(native),aligned_query=str(aligned),
            protein_alignment=aligned_manifest['source'].get('protein_alignment'),
            coordinate_scope='Observed co-crystal ligand coordinates; original saved alignment retained'))
        files.extend([native,native.with_suffix('.manifest.json'),*map(Path,original)])
    sources = {**report['sources'],**fingerprint(files)}
    save(output/'reference_ligands.json',dict(kind='saved_cocrystal_ligand_references',references=references))
    result = dict(report,sources=sources,query_adoption=dict(source=str(path),source_sha256=sha(path),
                  template_count=len(report['templates']),query_ids=[q['query_id'] for q in report['templates']],
                  policy='Original multi-cocrystal query reused verbatim; no new template selection, recommendation or search.' ))
    result['outputs'] = dict(report.get('outputs',{}),reference_ligands=str(output/'reference_ligands.json'))
    check_hashes(sources)
    save(output/'report.json',result)
    return result
