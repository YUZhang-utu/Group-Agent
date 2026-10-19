"""E112: native-receptor paired audit or bounded fresh per-state pose search."""
import argparse
import csv
import hashlib
from contextlib import closing
import json
from pathlib import Path
import time
import numpy as np
from scipy.spatial import cKDTree
import audit_strict_cached_poses as base


def rigid(value):
    m = np.asarray(value, float).reshape(4, 4)
    base.require(np.isfinite(m).all() and np.allclose(m[3], [0,0,0,1]) and
                 np.allclose(m[:3,:3] @ m[:3,:3].T, np.eye(3), atol=1e-5) and
                 np.linalg.det(m[:3,:3]) > 0, 'Invalid rigid transform')
    return m


def move(xyz, matrix):
    return np.asarray(xyz) @ matrix[:3,:3].T + matrix[:3,3]


def paired_receptor(native_points, native_ligand, aligned_ligand, matrix):
    matrix = rigid(matrix)
    base.require(np.allclose(move(native_ligand, matrix), aligned_ligand, atol=1e-5),
                 'Native/aligned ligand coordinate mismatch')
    native_tree = cKDTree(base.points(native_points))
    receptor = move(native_points, matrix)
    before = native_tree.query(native_ligand)[0]
    after = cKDTree(receptor).query(aligned_ligand)[0]
    base.require(np.allclose(before, after, atol=1e-5), 'Rigid receptor/ligand distance invariance failed')
    return receptor, float(before.min())


def filtered_contacts(rows, original_anchors, receptor_atoms, chain):
    allowed = {(int(a['target_residue']), a['protein_atom'], a['feature_class']) for a in original_anchors}
    canonical = {(a['auth_asym_id'], str(a['auth_seq_id'])):int(a['canonical_residue']) for a in receptor_atoms}
    available = {(int(a['canonical_residue']),a['auth_atom_id']) for a in receptor_atoms}
    selected = []
    seen = set()
    for row in rows:
        partner = row['evidence'].get('protein_partner', {})
        if partner.get('chain') != chain:
            continue
        residue = canonical.get((chain,str(partner.get('residue_number'))))
        atom = partner.get('atom_name','')
        kind = row.get('interaction_class',row['feature_class'])
        signature = (residue,atom,kind)
        identity = (signature,int(row['feature_index']))
        if signature not in allowed or identity in seen:
            continue
        if not all((residue,name) in available for name in atom.split('/')):
            continue
        seen.add(identity)
        selected.append(dict(anchor_id=row['anchor_id'], feature_index=int(row['feature_index']),
                             target_residue=residue, protein_atom=atom, feature_class=kind,
                             evidence=row['evidence']))
    return selected


def native_contacts(state, features, transform):
    from aidd_agent.interaction_fast import match_batch
    anchors = state['anchors']
    if not anchors:
        return False, [], []
    indices = sorted({a['feature_index'] for a in anchors})
    data = state['aligned']
    arrays = tuple(data[k][indices] for k in ('feature_points','feature_types','feature_directions',
                                             'feature_direction_kinds','anchored_weights'))
    arrays = (*arrays[:-1], np.ones(len(indices)))
    _, assignments, scores = match_batch(arrays, features.feature_points[None], features.feature_types[None],
                                         features.feature_directions[None], features.feature_kinds[None],
                                         rigid(transform)[None])
    cols = {index:i for i,index in enumerate(indices)}
    hits = [a for a in anchors if assignments[0,cols[a['feature_index']]] >= 0 and
            scores[0,cols[a['feature_index']]] >= base.POLICY['anchor_threshold']]
    residues = sorted({a['target_residue'] for a in hits})
    ile = any(a['target_residue']==61 and a['protein_atom']=='CD1' and
              a['feature_class']=='hydrophobic' for a in hits)
    return bool(ile and len(residues)>=3), residues, [a['anchor_id'] for a in hits]


def build_state(template, query):
    from aidd_agent.gaussian_batch import _load_query
    from aidd_agent.screening_selection import anchor_evidence
    from aidd_agent.classified_features import crystal_environment, extra_hypotheses
    from aidd_agent.consensus_admission import read_structure
    path = Path(template['query_npz'])
    manifest, aligned = _load_query(path)
    source = manifest['source']
    base.require(source['query_id']==template['query_id'], 'Aligned template identity mismatch')
    admission = source['protein_alignment']
    matrix = rigid(admission['transform'])
    chain = admission['target_chain']
    native_path = path.parent/'native.npz'
    native_manifest, native = _load_query(native_path)
    ns = native_manifest['source']
    base.require(ns['query_id']==template['query_id'], 'Native template identity mismatch')
    base.require(np.array_equal(native['feature_types'],aligned['feature_types']) and
                 np.allclose(move(native['feature_points'],matrix),aligned['feature_points'],atol=1e-5) and
                 np.allclose(native['feature_directions']@matrix[:3,:3].T,aligned['feature_directions'],atol=1e-5),
                 'Native/aligned feature mismatch')
    rows, paths, _ = anchor_evidence(native_path)
    cif = Path(ns['mmcif'])
    structure = read_structure(cif,query['target'])
    base.require(chain in structure['chains'], 'Native target chain not verified')
    atoms = [a for a in structure['chains'][chain] if not a['label_alt_id'] and
             float(a['occupancy'])>=0.9 and not a.get('pdbx_PDB_ins_code')]
    base.require(atoms, 'No unambiguous fully occupied receptor atoms')
    protein_atoms, _, _ = crystal_environment(cif,template['query_id'])
    atom_keys = {(a['auth_asym_id'],str(a['auth_seq_id']),a['auth_atom_id']) for a in atoms}
    protein_atoms = [a for a in protein_atoms if (a['chain'],str(a['residue_number']),a['atom_name']) in atom_keys]
    rows += extra_hypotheses(template['query_id'],native,protein_atoms)
    anchors = filtered_contacts(rows,query['anchors'],atoms,chain)
    receptor, native_minimum = paired_receptor(np.array([a['xyz'] for a in atoms]),native['shape_points'],
                                               aligned['shape_points'],matrix)
    tree = cKDTree(receptor)
    contact, regions = base.make_regions(aligned['shape_points'],tree)
    support = any(a['target_residue']==61 and a['protein_atom']=='CD1' and a['feature_class']=='hydrophobic'
                  for a in anchors) and len({a['target_residue'] for a in anchors})>=3
    resources = set(map(Path,paths)) | {path,path.with_suffix('.manifest.json'),native_path,native_path.with_suffix('.manifest.json')}
    for resource in resources:
        digest = query['sources'].get(str(resource)) or query['sources'].get(str(resource.resolve()))
        base.require(digest is not None,'Native resource missing from original source seals: '+str(resource))
        base.verify(resource,digest)
    record = dict(template=template['query_id'],chain=chain,native_mmcif=str(cif),transform=matrix.tolist(),
                  receptor_atoms=len(atoms),native_minimum_distance=native_minimum,anchors=anchors,
                  original_protein_alignment=admission, receptor_scope='Verified target chain only; no added symmetry mates or waters',
                  supported_for_requested_contact_rule=support,
                  excluded_atom_policy='Exclude alternate locations, insertion codes and occupancy <0.9',
                  sources={str(p):base.sha(p) for p in resources},
                  region_points=[r.tolist() for r in regions] if regions else [])
    return dict(aligned=aligned,tree=tree,contact=contact,regions=regions,anchors=anchors,record=record)


def run(args):
    from types import SimpleNamespace
    from aidd_agent.gaussian_batch import ArtifactCatalogReader
    from aidd_agent.gaussian_overlay import score_overlay
    from aidd_agent.interaction_fast import InteractionFeatureReader
    from aidd_agent.block_sampling import snapshot
    start = time.monotonic()
    previous, prepared, reference, out = [p.resolve() for p in (args.previous,args.prepared,args.reference,args.output)]
    base.require(args.seconds>0,'Positive time budget required')
    base.require(all(not out.is_relative_to(p) and not p.is_relative_to(out) for p in [previous,prepared,reference]),
                 'Output overlaps inputs')
    out.mkdir(parents=True,exist_ok=False)
    base.save(out/'report.json',dict(status='incomplete'))
    old = base.read(previous/'report.json')
    old_protocol = base.read(previous/'protocol.json')
    base.require(old['status']=='complete','Need complete E111 audit')
    for name,digest in old['output_hashes'].items(): base.verify(previous/name,digest)
    prep, ref = base.read(prepared/'report.json'),base.read(reference/'report.json')
    base.require(prep['status']==ref['status']=='complete','Incomplete original input')
    base.require(base.sha(prepared/'report.json')==old_protocol['prepared_report_sha256']==ref['prepared_report_sha256'],
                 'Prepared lineage mismatch')
    base.require(base.sha(reference/'report.json')==old_protocol['reference_report_sha256'],'Reference lineage mismatch')
    base.verify(prepared/'query.json',old_protocol['query_sha256'])
    for name in ['scores.sqlite','schedule.sqlite']: base.verify(reference/name,ref['output_hashes'][name])
    query = base.read(prepared/'query.json')
    batch = Path(prep['batch'])
    base.require(snapshot(batch/'registry.sqlite3')==prep['registry_snapshot'],'Registry changed')
    code_checks = {}
    scientific_helpers = ['gaussian_batch.py', 'gaussian_overlay.py', 'interaction_fast.py',
                          'interaction_matching.py', 'consensus_admission.py', 'classified_features.py',
                          'screening_selection.py', 'guided_filters.py', 'contact_groups.py',
                          'full_library_screen.py', 'chemical_geometry.py', 'conformer_artifacts.py',
                          'chemical_companion.py']
    for name in scientific_helpers:
        digest = prep['code'][name]
        source_path = Path(__file__).resolve().parents[1]/'src'/'aidd_agent'/Path(name).name
        base.require(source_path.is_file(), 'Missing pinned scientific helper: '+name)
        raw = source_path.read_bytes()
        normalized = raw.replace(b'\r\n', b'\n')
        base.require(digest in [hashlib.sha256(raw).hexdigest(), hashlib.sha256(normalized).hexdigest()],
                     'Pinned scientific helper changed: '+name)
        code_checks[name] = digest
    states = {t['query_id']:build_state(t,query) for t in query['templates']}
    common_path=Path(query['design']['receptor_npz'])
    base.verify(common_path,query['sources'][str(common_path)])
    with np.load(common_path,allow_pickle=False) as z: common=cKDTree(base.points(z['points']))
    base.save(out/'native_states.json',[s['record'] for s in states.values()])
    controls=[]
    for key,s in states.items():
        data=s['aligned']
        features=SimpleNamespace(feature_points=data['feature_points'],feature_types=data['feature_types'],
                                 feature_directions=data['feature_directions'],feature_kinds=data['feature_direction_kinds'])
        ok,residues,hits=native_contacts(s,features,np.eye(4))
        metrics=base.geometry_metrics(data['shape_points'],s['tree'],s['contact'],s['regions'])
        flags=base.gates(ok,1,1,metrics)
        controls.append(dict(template=key,native_contact_supported=s['record']['supported_for_requested_contact_rule'],
                             residues=residues,common_minimum_distance=float(common.query(data['shape_points'])[0].min()),
                             native_minimum_distance=metrics['minimum_distance'],**flags,accepted=base.accepted(flags)))
    base.save(out/'reference_controls.json',controls)
    cross=[]
    for target,s in states.items():
        for ligand,c in states.items():
            if target==ligand: continue
            a,b=s['aligned'],c['aligned']
            score=score_overlay(a['shape_points'],b['shape_points'],a['feature_points'],a['feature_types'],
                                b['feature_points'],b['feature_types'],np.eye(4).ravel(),cutoff=None)
            cross.append(dict(reference=target,ligand=ligand,shape=score['shape']['tanimoto'],
                              color=score['color']['tanimoto'],
                              shape_color_pass=score['shape']['tanimoto']>=0.5 and score['color']['tanimoto']>=0.3))
    base.write_csv(out/'cross_template_scores.csv',cross)
    settings=dict(policy=base.POLICY,previous=str(previous),prepared=str(prepared),reference=str(reference),
                  previous_report_sha256=base.sha(previous/'report.json'),query_sha256=old_protocol['query_sha256'],
                  script_sha256=base.sha(__file__),helper_sha256=base.sha(Path(base.__file__)),seconds=args.seconds,
                  mode='fresh_bounded_pose_pilot' if args.fresh else 'paired_cached_poses_no_new_alignment',
                  sample=args.sample,workers=args.workers,pinned_code_checks=code_checks)
    base.save(out/'protocol.json',settings)
    base.require(any(s['record']['supported_for_requested_contact_rule'] for s in states.values()),
                 'No native state supports requested contacts; inspect native_states.json and reference_controls.json')
    with (previous/'sample.csv').open(newline='',encoding='utf-8') as f: selected=list(csv.DictReader(f))
    base.require(len(selected)==old['processed_conformers'] and len({r['gid'] for r in selected})==len(selected),
                 'E111 sample count/identity mismatch')
    with base.ro(reference/'schedule.sqlite') as db:
        for row in selected:
            found=db.execute('SELECT cid,mid FROM selected WHERE gid=?',(int(row['gid']),)).fetchone()
            base.require(found==(row['cid'],row['mid']),'Sample differs from original schedule')
    if args.fresh:
        return fresh_pilot(args,selected,states,query,batch,out,controls,start)
    with (previous/'poses.csv').open(newline='',encoding='utf-8') as f:
        old_scores={(int(r['gid']),r['template']):(float(r['shape_tanimoto']),float(r['color_tanimoto'])) for r in csv.DictReader(f)}
    reader=ArtifactCatalogReader(batch/'artifacts/catalog.json')
    features=InteractionFeatureReader(batch/'artifacts/catalog.json',batch/'chemical/catalog.json')
    setup=time.monotonic()-start; tick=time.monotonic()
    poses=[]; results=[]
    try:
        with base.ro(reference/'scores.sqlite') as db:
            for item in selected:
                if time.monotonic()-tick>args.seconds: break
                gid=int(item['gid']); saved=list(db.execute('SELECT template,payload FROM poses WHERE gid=? ORDER BY template',(gid,)))
                cumulative={k:False for k in base.GATES}
                solo={k:False for k in base.GATES}
                if saved:
                    candidate=reader.get(gid); cf=features.get(gid)
                    base.require(candidate.conformer_id==cf.conformer_id==item['cid'] and candidate.molecule_id==cf.molecule_id,
                                 'Candidate identity mismatch')
                    for key,payload in saved:
                        p=json.loads(payload); m=rigid(p['transform']);s=states[key]
                        base.require(p['global_id']==gid and p['conformer_id']==candidate.conformer_id and
                                     p['molecule_id']==candidate.molecule_id,'Saved pose identity mismatch')
                        ok,residues,hits=native_contacts(s,cf,m)
                        shape,color=old_scores[(gid,key)]
                        base.require(abs((shape+color)/2-p['gaussian_same_pose'])<1e-4,'Changed cached shape/color')
                        moved=move(candidate.shape_points,m)
                        metric=base.geometry_metrics(moved,s['tree'],s['contact'],s['regions'])
                        flags=base.gates(ok,shape,color,metric); prefix=True
                        for k in base.GATES:
                            solo[k]|=flags[k];prefix=prefix and flags[k];cumulative[k]|=prefix
                        poses.append(dict(gid=gid,mid=item['mid'],template=key,residues=json.dumps(residues),
                                          native_contact_supported=s['record']['supported_for_requested_contact_rule'],
                                          shape_tanimoto=shape,color_tanimoto=color,
                                          common_minimum_distance=float(common.query(moved)[0].min()),
                                          native_minimum_distance=metric['minimum_distance'],
                                          envelope_coverage=metric['envelope_coverage'],
                                          region_coverages=json.dumps(metric['region_coverages']),
                                          **flags,accepted=base.accepted(flags)))
                results.append(dict(gid=gid,mid=item['mid'],old_pass=bool(saved),**cumulative,
                                    **{'standalone_'+k:v for k,v in solo.items()}))
                if len(results)%100==0: print('Audited '+str(len(results))+'/'+str(len(selected)),flush=True)
    finally:
        features.close()
        if hasattr(reader,'close'):reader.close()
    base.write_csv(out/'poses.csv',poses);base.write_csv(out/'conformers.csv',results)
    complete=len(results)==len(selected)
    report=dict(status='complete' if complete else 'time_budget_exhausted',kind='native_state_paired_cached_audit',
                processed_conformers=len(results),selected_conformers=len(selected),
                old_pass_conformers=sum(r['old_pass'] for r in results),
                cumulative_conformer_counts={k:sum(r[k] for r in results) for k in base.GATES},
                standalone_conformer_counts={k:sum(r['standalone_'+k] for r in results) for k in base.GATES},
                strict_pass_fraction=sum(r['regions'] for r in results)/len(results) if complete and results else None,
                reference_controls_passed=sum(r['accepted'] for r in controls),reference_controls=len(controls),
                supported_native_contact_states=sum(s['record']['supported_for_requested_contact_rule'] for s in states.values()),
                setup_seconds=setup,audit_seconds=time.monotonic()-tick,total_seconds=time.monotonic()-start,
                limitations=['No new pose search; same cached shape/color cannot be rescued by changing receptor.',
                             'Native contact hypotheses are recomputed; not validated energetic interactions.',
                             'Only adopted residue/atom/type combinations retained; unsupported native states are explicit.',
                             'Reference envelope coverage is not canonical pocket-volume evaluation.',
                             'Fixed-pose results do not measure fresh state-specific search recall.',
                             'Cross-template calibration is exploratory; no thresholds changed.'])
    report['output_hashes']={p.name:base.sha(p) for p in out.iterdir() if p.is_file() and p.name!='report.json'}
    base.save(out/'report.json',report);print(json.dumps(report,indent=2))
    return report


_FRESH = None


def fresh_init(batch, states, coarse, deadline):
    global _FRESH
    from aidd_agent.gaussian_batch import ArtifactCatalogReader
    from aidd_agent.interaction_fast import InteractionFeatureReader
    _FRESH = (ArtifactCatalogReader(Path(batch)/'artifacts/catalog.json'),
              InteractionFeatureReader(Path(batch)/'artifacts/catalog.json',Path(batch)/'chemical/catalog.json'),
              states, coarse, deadline)


def fresh_one(item):
    from aidd_agent.gaussian_batch import prepare_seeds
    from aidd_agent.gaussian_overlay import score_overlay
    from aidd_agent.full_library_screen import POSE_PARAMETERS
    from aidd_agent.guided_filters import adaptive_stage
    reader, feature_reader, states, coarse, deadline = _FRESH
    gid=int(item['gid'])
    cumulative={k:False for k in base.GATES}
    solo={k:False for k in base.GATES}
    result=dict(gid=gid,mid=item['mid'],complete=False,tested_seeds=0,coarse_pass=False,
                best_shape_after_contact=None,best_color_after_contact=None,passing_poses=[])
    def finish(done):
        return dict(result,complete=done,**cumulative,**{'standalone_'+k:v for k,v in solo.items()})
    if time.monotonic()>deadline: return finish(False)
    candidate=reader.get(gid); features=feature_reader.get(gid)
    base.require(candidate.conformer_id==features.conformer_id==item['cid'] and
                 candidate.molecule_id==features.molecule_id,'Fresh candidate identity mismatch')
    if adaptive_stage(candidate,coarse)<3: return finish(True)
    result['coarse_pass']=True
    for key,state in states.items():
        if not state['record']['supported_for_requested_contact_rule']: continue
        if time.monotonic()>deadline: return finish(False)
        seeds,_=prepare_seeds(candidate,state['aligned'],**POSE_PARAMETERS)
        stored=False
        for seed in seeds:
            if time.monotonic()>deadline: return finish(False)
            m=rigid(seed.transform_matrix)
            result['tested_seeds']+=1
            ok,residues,hits=native_contacts(state,features,m)
            if not ok: continue
            cumulative['contacts']=solo['contacts']=True
            data=state['aligned']
            scores=score_overlay(data['shape_points'],candidate.shape_points,data['feature_points'],data['feature_types'],
                                 candidate.feature_points,candidate.feature_types,m.ravel(),cutoff=None)
            shape,color=scores['shape']['tanimoto'],scores['color']['tanimoto']
            for name,value in [('best_shape_after_contact',shape),('best_color_after_contact',color)]:
                result[name]=max(result[name],value) if result[name] is not None else value
            metrics=base.geometry_metrics(move(candidate.shape_points,m),state['tree'],state['contact'],state['regions'])
            flags=base.gates(ok,shape,color,metrics);prefix=True
            for k in base.GATES:
                solo[k]|=flags[k];prefix=prefix and flags[k];cumulative[k]|=prefix
            if base.accepted(flags) and not stored:
                result['passing_poses'].append(dict(gid=gid,mid=item['mid'],template=key,residues=residues,
                                                    shape=shape,color=color,transform=m.ravel().tolist(),**metrics))
                stored=True
    return finish(True)


def fresh_pilot(args, selected, states, query, batch, out, controls, started):
    from concurrent.futures import ProcessPoolExecutor
    import multiprocessing
    base.require(args.sample>0 and 1<=args.workers<=32,'Invalid sample/workers')
    selected=selected[:args.sample]
    base.require(selected,'Empty sample')
    base.write_csv(out/'sample.csv',selected)
    setup=time.monotonic()-started;tick=time.monotonic();deadline=tick+args.seconds
    outcomes=[]
    with ProcessPoolExecutor(max_workers=args.workers,mp_context=multiprocessing.get_context('spawn'),
                             initializer=fresh_init,initargs=(str(batch),states,query['design']['adaptive_coarse'],deadline)) as pool:
        for outcome in pool.map(fresh_one,selected,chunksize=1):
            outcomes.append(outcome)
            if len(outcomes)%10==0: print('Returned '+str(len(outcomes))+'/'+str(len(selected)),flush=True)
    processed=[r for r in outcomes if r['complete']]
    complete=len(processed)==len(selected)
    passes=[p for r in outcomes for p in r['passing_poses']]
    base.save(out/'passing_poses.json',passes)
    base.write_csv(out/'conformers.csv',[{k:v for k,v in r.items() if k!='passing_poses'} for r in outcomes])
    report=dict(status='complete' if complete else 'time_budget_exhausted',kind='native_state_fresh_bounded_pose_pilot',
                selected_conformers=len(selected),complete_conformers=len(processed),
                incomplete_conformers=len(outcomes)-len(processed),
                coarse_pass_complete_conformers=sum(r['coarse_pass'] for r in processed),
                cumulative_conformer_counts={k:sum(r[k] for r in processed) for k in base.GATES},
                diagnostics_conditional_on_contact_counts={k:sum(r['standalone_'+k] for r in processed) for k in base.GATES},
                strict_pass_fraction=sum(r['regions'] for r in processed)/len(processed) if complete else None,
                tested_seeds=sum(r['tested_seeds'] for r in outcomes),
                supported_native_contact_states=sum(s['record']['supported_for_requested_contact_rule'] for s in states.values()),
                template_controls_passed=sum(r['accepted'] for r in controls),template_controls=len(controls),
                workers=args.workers,setup_seconds=setup,search_seconds=time.monotonic()-tick,total_seconds=time.monotonic()-started,
                limitations=['Fresh bounded pair/principal-axis seed search, not exhaustive pose search or docking.',
                             'Same-pose gates; no score/interaction mixing across states.',
                             'Only source-supported selected native contact modes; unsupported states are reported, not silently fabricated.',
                             'Shape/color and geometry thresholds unchanged and exploratory.',
                             'Contact matching is a geometric hypothesis, not binding-energy validation.',
                             'Coverage is a reference contact-envelope proxy, not canonical pocket volumetry.',
                             'Partial conformers excluded from complete counts; no overall fraction on timeout.',
                             'Conformer sample, not a molecule-complete population.'])
    report['output_hashes']={p.name:base.sha(p) for p in out.iterdir() if p.is_file() and p.name!='report.json'}
    base.save(out/'report.json',report);print(json.dumps(report,indent=2));return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['previous','prepared','reference','output']: p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--seconds',type=float,default=600)
    p.add_argument('--fresh',action='store_true')
    p.add_argument('--sample',type=int,default=200)
    p.add_argument('--workers',type=int,default=8)
    run(p.parse_args())


if __name__=='__main__':main()
