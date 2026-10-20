"""E114 local shape diagnostics from sealed E113 recovered coordinates."""
import argparse
from pathlib import Path
import time
import numpy as np
from scipy.optimize import minimize
from scipy.spatial import cKDTree
from scipy.spatial.transform import Rotation
import audit_native_state_poses as native
base=native.base
from aidd_agent.gaussian_overlay import gaussian_overlap, tanimoto


def primitives(query,candidate,sigma=1.):
    a=gaussian_overlap(query,query,sigma=sigma,cutoff=None)
    b=gaussian_overlap(candidate,candidate,sigma=sigma,cutoff=None)
    c=gaussian_overlap(query,candidate,sigma=sigma,cutoff=None)
    bound=np.sqrt(a*b)
    return dict(query_self=a,candidate_self=b,cross=c,shape=tanimoto(c,a,b),
                norm_upper_bound=tanimoto(bound,a,b))


def delta_transform(parameters,center):
    rotation=Rotation.from_rotvec(parameters[:3]).as_matrix()
    matrix=np.eye(4);matrix[:3,:3]=rotation
    matrix[:3,3]=center+parameters[3:]-rotation@center
    return matrix


def refine(query,candidate,deadline,max_calls=180):
    values=primitives(query,candidate)
    center=candidate.mean(axis=0);best=[values['shape'],np.eye(4)];calls=0
    class Stop(Exception):pass
    def objective(x):
        nonlocal calls
        if time.monotonic()>deadline:raise Stop('deadline')
        if calls>=max_calls:raise Stop('evaluation_limit')
        calls+=1
        matrix=delta_transform(x,center)
        cross=gaussian_overlap(query,native.move(candidate,matrix),cutoff=None)
        score=tanimoto(cross,values['query_self'],values['candidate_self'])
        if score>best[0]:best[:]=[score,matrix]
        return -score
    try:
        fit=minimize(objective,np.zeros(6),method='L-BFGS-B',bounds=[(-.5,.5)]*3+[(-2.,2.)]*3,
                     options=dict(maxfun=max_calls,maxiter=30,ftol=1e-10))
        status=str(fit.message);success=bool(fit.success)
    except Stop as error:
        status=str(error);success=False
    return dict(shape=best[0],delta=best[1].tolist(),calls=calls,success=success,status=status)


def point_coverage(query,candidate):
    return dict(reference_near_candidate=float(np.mean(cKDTree(candidate).query(query)[0]<=1.5)),
                candidate_near_reference=float(np.mean(cKDTree(query).query(candidate)[0]<=1.5)))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--receipt',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--seconds',type=float,default=120)
    args=parser.parse_args();started=time.monotonic()
    base.require(args.seconds>0,'Positive budget required')
    source=args.receipt.resolve();out=args.output.resolve()
    base.require(not out.is_relative_to(source) and not source.is_relative_to(out),'Overlapping input/output')
    report=base.read(source/'report.json')
    base.require(report['status']=='complete','Complete E113 receipt required')
    for name,digest in report['output_hashes'].items():base.verify(source/name,digest)
    rows=base.read(source/'pose_diagnostics.json')
    states={s['template']:s for s in base.read(source/'native_states.json')}
    references={}
    for row in rows:
        if row['kind']=='self_control' and row['complete']:
            pose=row['best']['shape']
            base.require(abs(pose['shape']-1)<1e-8,'Self reference not recovered exactly')
            references[row['id']]=base.points(pose['aligned_heavy_atom_coordinates'])
    panel=[]
    for row in rows:
        best=row.get('best',{})
        objective=('coverage_with_geometry' if row['kind']=='library' else 'shape')
        if row['kind']!='self_control' and objective in best:
            panel.append((row,objective,'geometry_near_miss' if row['kind']=='library' else 'leave_self_out'))
    library=[r for r in rows if r['kind']=='library' and 'shape' in r.get('best',{})]
    for row in sorted(library,key=lambda r:(-r['best']['shape']['shape'],r['id']))[:20]:
        panel.append((row,'shape','top_shape'))
    out.mkdir(parents=True,exist_ok=False)
    base.save(out/'report.json',dict(status='incomplete'))
    output=[];proposals=[];tick=time.monotonic();deadline=tick+args.seconds
    for row,objective,category in panel:
        if time.monotonic()>deadline:break
        pose=row['best'][objective];query=references[pose['template']]
        candidate=base.points(pose['aligned_heavy_atom_coordinates'])
        before=primitives(query,candidate)
        base.require(abs(before['shape']-pose['shape'])<1e-6,'Reconstructed reference score mismatch')
        optimized=refine(query,candidate,deadline)
        delta=np.array(optimized['delta']);moved=native.move(candidate,delta)
        original_transform=native.rigid(pose['transform'])
        result=dict(id=row['id'],mid=row['mid'],kind=row['kind'],category=category,template=pose['template'],
                    objective=objective,reference_atoms=len(query),candidate_atoms=len(candidate),
                    atom_ratio=len(candidate)/len(query),original_shape=pose['shape'],color=pose['color'],
                    contact_coverage=pose['coverage'],original_geometry_pass=pose['geometry_pass'],
                    original_primitives=before,original_point_coverage=point_coverage(query,candidate),
                    optimized=optimized,optimized_point_coverage=point_coverage(query,moved),
                    sigma_sensitivity={str(s):primitives(query,candidate,s)['shape'] for s in [.75,1.,1.5]},
                    original_metrics={k:pose[k] for k in ['minimum_distance','envelope_coverage','region_coverages']},
                    optimized_geometry_validated=False)
        regions=states[pose['template']]['region_points']
        if regions:
            contact=np.concatenate(regions)
            result['reference_contact_atom_count']=len(contact)
            result['optimized_contact_envelope']=float(np.mean(cKDTree(moved).query(contact)[0]<=1.5))
        output.append(result)
        if row['kind'] in ('library','leave_self_out'):
            proposals.append(dict(kind=row['kind'],id=row['id'],mid=row['mid'],objective=objective,
                                  template=pose['template'],category=category,
                                  original_transform=pose['transform'],
                                  transform=(delta@original_transform).ravel().tolist(),
                                  original_shape=pose['shape'],optimized_shape=optimized['shape']))
    base.save(out/'diagnostics.json',output)
    base.save(out/'proposals.json',proposals)
    base.save(out/'protocol.json',dict(receipt_report_sha256=base.sha(source/'report.json'),
        script_sha256=base.sha(__file__),seconds=args.seconds,panel_size=len(panel),
        reconstruction_tolerance=1e-6,rotation_component_bound=.5,translation_component_bound=2.,max_calls=180))
    result=dict(status='complete' if len(output)==len(panel) and all(r['optimized']['status']!='deadline' for r in output)
                else 'time_budget_exhausted',selected_poses=len(panel),processed_poses=len(output),
                library_proposals=sum(p['kind']=='library' for p in proposals),
                control_proposals=sum(p['kind']!='library' for p in proposals),reference_clouds=len(references),
                shape_crossings=sum(r['original_shape']<.5<=r['optimized']['shape'] for r in output),
                setup_seconds=tick-started,optimization_seconds=time.monotonic()-tick,
                total_seconds=time.monotonic()-started,
                limitations=['Outcome-selected saved poses, not population pass rate.',
                             'Rigid local optimization only; no new conformers or exhaustive search.',
                             'Optimized receptor contacts, clashes and color not checked locally.',
                             'Self-norm bound is necessary only; a high bound does not imply attainable fit.',
                             'Point coverage and kernel-width comparisons are diagnostics, not molecular volume.'])
    result['output_hashes']={p.name:base.sha(p) for p in out.iterdir() if p.name!='report.json'}
    base.save(out/'report.json',result);print(base.json.dumps(result,indent=2))


if __name__=='__main__':main()
