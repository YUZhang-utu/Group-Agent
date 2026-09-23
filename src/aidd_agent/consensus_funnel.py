"""Full-library multi-template execution with explicit union and pose provenance."""
import copy
import json
from pathlib import Path
import sqlite3

from . import library_acceptance as ev
from .expanded_wee1 import fingerprint
from .gaussian_batch import _atomic_json
from .screening_selection import check_hashes


def run(source, output, search):
    from .preselection_full import run as scan
    source=Path(source);out=Path(output);out.mkdir(parents=True,exist_ok=True)
    adopted=ev.read(source);check_hashes(adopted['sources'])
    if adopted['readiness']!='ready_for_consensus_funnel':raise ValueError('Adopt a consensus design first')
    batch=Path(search['batch']);_,acceptance=ev.accept_library(batch,full=False)
    design=adopted['design'];sources={**adopted['sources'],**fingerprint([source,batch/'artifacts/catalog.json',batch/'chemical/catalog.json'])}
    reports=[]
    result=dict(kind='consensus_funnel',status='running',full_library=True,full_coverage=False,
        anchor_order=adopted['anchor_order'],policy=dict(minimum_score=design['minimum_score']),design=design,
        library={k:acceptance[k] for k in ('library_conformers','library_molecules')},template_reports=reports,
        scope='Each selected template scans every catalog conformer; same-pose scores and anchor masks; no candidate Top-K',
        scheduling='Sequential template scans, parallel conformer chunks; repeated catalog I/O is included in runtime',
        independent_active_validation=adopted['independent_active_validation'])
    from .active_controls import run as controls
    ev.log('Checking independent binding controls separately from full-library coverage')
    controls_report=controls(adopted,out/'active-controls')
    result['independent_active_validation']=controls_report
    for i,template in enumerate(adopted['templates']):
        folder=out/f'template-{i:03d}';inputs=folder/'inputs';inputs.mkdir(parents=True,exist_ok=True)
        q=dict(template,anchors=adopted['anchors'],consensus_npz=adopted['consensus_npz'],
            artifact_catalog=str((batch/'artifacts/catalog.json').resolve()),chemical_companion=str((batch/'chemical/catalog.json').resolve()))
        classified=inputs/'classification.json'
        _atomic_json(classified,dict(kind='screening_evidence',status='complete',queries=[q],sources=sources,
            library=result['library'],classification='decoupled-pocket-consensus-v1'))
        selection=inputs/'selection.json'
        _atomic_json(selection,dict(kind='selection_preview',query=q,evidence_report=str(classified.resolve()),guided_design=design,
            policy=dict(required_anchors=adopted['anchor_order'],match_mode='any',minimum_score=design['minimum_score'],coarse_constraints=design['coarse_constraints']),
            sources={**sources,**fingerprint([classified])}))
        _atomic_json(out/'report.json',result)
        child=scan(selection,folder/'full-library',workers=search['workers'],chunk_size=search['refine_chunk'])
        _atomic_json(folder/'full-library/report.json',child)
        reports.append(dict(template_id=template['query_id'],report=str((folder/'full-library/report.json').resolve()),
                            full_coverage=child['full_coverage'],counts=child['counts']))
    dbpath=out/'candidates.sqlite'
    with sqlite3.connect(dbpath) as db:
        db.executescript('DROP TABLE IF EXISTS poses; DROP TABLE IF EXISTS members; CREATE TABLE poses(mid TEXT,template TEXT,mask INTEGER,quality REAL,payload TEXT,spatial TEXT NOT NULL,PRIMARY KEY(mid,template,mask,spatial)); CREATE TABLE members(mid TEXT PRIMARY KEY,cluster TEXT);')
        for meta in reports:
            child=ev.read(meta['report']);check_hashes(child['output_hashes'])
            with sqlite3.connect(Path(child['outputs']['candidates.sqlite']).resolve().as_uri()+'?mode=ro',uri=True) as incoming:
                for mid,mask,quality,payload in incoming.execute('SELECT mid,mask,quality,payload FROM poses'):
                    row=json.loads(payload);row['template_id']=meta['template_id'];row['coordinate_frame']=adopted['reference']['coordinate_frame']
                    signature=json.dumps(sorted(row.get('occupied_spatial_groups',[])),separators=(',',':'))
                    db.execute('INSERT INTO poses VALUES(?,?,?,?,?,?)',(mid,meta['template_id'],mask,quality,json.dumps(row),signature))
                for mid,cluster in incoming.execute('SELECT mid,cluster FROM members'):db.execute('INSERT OR IGNORE INTO members VALUES(?,?)',(mid,cluster))
        count=db.execute('SELECT COUNT(*) FROM members').fetchone()[0];poses=db.execute('SELECT COUNT(*) FROM poses').fetchone()[0]
    result.update(status='complete',full_coverage=all(r['full_coverage'] for r in reports),matching_molecules=count,pose_records=poses,
        outputs={'candidates.sqlite':str(dbpath.resolve())},output_hashes=fingerprint([dbpath]),sources=sources,
        acceptance='Full catalog coverage is not independent active recall, exhaustive orientations, docking or affinity validation')
    _atomic_json(out/'report.json',result);return result
