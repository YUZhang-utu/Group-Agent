"""Read owned task structures and contact evidence for Chat and desktop PyMOL."""
from collections import Counter
import hashlib
import json
from pathlib import Path

from .project_context import ensure_within


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def task_report(job, project):
    if job.get('status')!='complete' or not job.get('report') or not job.get('plan'):
        raise ValueError('Choose a completed structure task')
    return read(ensure_within(Path(job['report']),project))


def consensus_report(child, project):
    # Follow only known report links, bounded to the owned project.
    for _ in range(4):
        if child.get('kind')=='structure_consensus': return child
        if not child.get('survey'):
            # Older adopted designs kept source hashes but did not include a survey link.
            matches=[]
            for name,digest in child.get('sources',{}).items():
                if Path(name).name!='report.json':continue
                path=ensure_within(Path(name),project)
                candidate=read(path)
                if candidate.get('kind')=='structure_consensus':
                    if hashlib.sha256(path.read_bytes()).hexdigest()!=digest:raise ValueError('Consensus source changed')
                    matches.append(candidate)
            if len(matches)==1:return matches[0]
            if len(matches)>1:raise ValueError('Ambiguous consensus source; select the consensus task directly')
            return child
        child=read(ensure_within(Path(child['survey']),project))
    raise ValueError('Too many structure report links')


def structures(job, project, collection='diverse'):
    report=task_report(job,project); result=[]
    def add(path,label,ligand=None,transform=None,expected=None,kind='experimental'):
        path=ensure_within(Path(path),project)
        digest=hashlib.sha256(path.read_bytes()).hexdigest()
        if expected and digest!=expected: raise ValueError('Structure provenance mismatch')
        result.append(dict(id=f'v{len(result)+1:03d}',path=str(path),label=label,
                           ligand=ligand or {},transform=transform,sha256=digest,kind=kind))
    for step in report.get('steps',{}).values():
        if step.get('status')!='complete': continue
        payload=step.get('result',{});action=step.get('action')
        if action=='af3_run':
            add(payload['model'],'AF3 predicted complex',expected=payload.get('structure_sha256'),kind='prediction')
        elif action=='pdb_fetch' and payload.get('path'):
            add(payload['path'],payload.get('pdb_id','PDB structure'),expected=payload.get('sha256'))
        elif payload.get('report'):
            path=ensure_within(Path(payload['report']),project); child=consensus_report(read(path),project)
            if child.get('kind')=='structure_consensus':
                admitted={r['query_id']:r for r in child.get('cohort',{}).get('admitted',[])}
                chosen=set(child.get('proposed_template_ids',[]))
                for row in child.get('prepared_complexes',[]):
                    qid=row['query_id']
                    if collection=='diverse' and chosen and qid not in chosen: continue
                    native=Path(row['query_npz']).parent/'native.manifest.json'
                    source=read(ensure_within(native,project))['source']
                    pdb,ccd,chain,residue=qid.split(':')
                    add(source['mmcif'],qid,dict(chain=chain,resi=residue,resn=ccd),
                        admitted[qid]['admission']['transform'],source.get('mmcif_sha256'))
            elif child.get('kind')=='structure_diversity':
                for row in child.get('proposed_references',[]):
                    qid=row['query_id'];pdb,ccd,chain,residue=qid.split(':')
                    add(path.parent/'structures'/f'{pdb}.cif',qid,dict(chain=chain,resi=residue,resn=ccd))
    if not result: raise ValueError('No supported structures in this task; choose AF3, PDB download, diversity or consensus results')
    if len(result)>100: raise ValueError('Viewer catalog exceeds 100 structures; use a smaller reference proposal')
    return result


def confidence(job,project):
    from .af3_analysis import explain_confidence
    report=task_report(job,project)
    steps=[s for s in report.get('steps',{}).values() if s.get('action')=='af3_run' and s.get('status')=='complete']
    if len(steps)!=1: raise ValueError('Choose a task with one completed AF3 prediction')
    return explain_confidence(steps[0]['result'].get('confidence',{}))


def contacts(job,project,query_id=None,residue=None,contact_class=None,limit=50):
    report=task_report(job,project)
    for step in report.get('steps',{}).values():
        payload=step.get('result',{})
        if step.get('status')!='complete' or not payload.get('report'): continue
        child=consensus_report(read(ensure_within(Path(payload['report']),project)),project)
        ledger_path=child.get('contact_evidence',{}).get('path')
        if not ledger_path: continue
        path=ensure_within(Path(ledger_path),project)
        expected=child.get('sources',{}).get(str(path))
        if expected and hashlib.sha256(path.read_bytes()).hexdigest()!=expected:
            raise ValueError('Contact ledger hash mismatch')
        ledger=read(path);rows=[];queries=[]
        for complex_row in ledger['complexes']:
            qid=complex_row['query_id'];queries.append(qid)
            if query_id and qid!=query_id: continue
            for pair in complex_row['pairs']:
                if residue and residue not in {str(pair.get('target_residue')),str(pair.get('protein_author_residue'))}:continue
                if contact_class and contact_class not in pair['classes']:continue
                rows.append(dict(query_id=qid,**pair))
        classes=Counter(c for row in rows for c in row['classes'])
        return dict(total_matching_pairs=len(rows),class_counts=dict(classes),pairs=rows[:limit],
                    truncated=len(rows)>limit,available_queries=queries,full_ledger=str(path),
                    coordinate_frame=ledger.get('coordinate_frame'),policy=ledger.get('policy'),
                    meanings={'donor_acceptor_proximity':'Nearby complementary assigned donor/acceptor atoms; angles and protonation are not validated.',
                              'heavy_atom_proximity':'Observed heavy atoms within the recorded cutoff; not an interaction energy.',
                              'carbon_sulfur_proximity':'Carbon/sulfur proximity; not a confirmed hydrophobic interaction.',
                              'halogen_proximity_not_halogen_bond':'Halogen proximity only; not a geometry-validated halogen bond.'})
    raise ValueError('No atom-contact ledger in this task. Use a completed consensus task, or open AF3/PDB in PyMOL and request polar contacts')
