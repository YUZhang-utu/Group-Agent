"""E122 exploratory receptor-derived grid scoring; no molecular rejection proof."""
import numpy as np
from scipy.ndimage import map_coordinates
from scipy.spatial import cKDTree
import kras_soft_rank as soft

PROFILES = dict(lenient=(.10,.60,.15), balanced=(.05,.40,.25), strict=(0.,.25,.35))
SETTINGS = dict(step=.75,padding=8.,conformers=3,exact_poses=8,
                ray_distances=[2.,4.,6.,8.],ray_atom_radius=1.7,buried_rays=3,
                overlap_distance=2.2,severe_distance=1.2)
DIRECTIONS=np.vstack((np.eye(3),-np.eye(3)))


def burial(points,tree):
    points=np.asarray(points).reshape(-1,3)
    hits=np.zeros((len(points),6),dtype=bool)
    for distance in SETTINGS['ray_distances']:
        probes=points[:,None,:]+DIRECTIONS[None,:,:]*distance
        hits |= tree.query(probes.reshape(-1,3))[0].reshape(-1,6)<SETTINGS['ray_atom_radius']
    return hits.sum(axis=1)>=SETTINGS['buried_rays']


def make_grid(state,step=None):
    step=step or SETTINGS['step'];points=state['query']['shape_points']
    origin=np.floor((points.min(axis=0)-SETTINGS['padding'])/step)*step
    shape=np.ceil((points.max(axis=0)+SETTINGS['padding']-origin)/step).astype(int)+1
    xyz=np.indices(shape).reshape(3,-1).T*step+origin
    distances=np.empty(len(xyz),dtype=np.float32); buried=np.empty(len(xyz),dtype=np.float32)
    for start in range(0,len(xyz),20000):
        p=xyz[start:start+20000];distances[start:start+len(p)]=state['tree'].query(p)[0]
        buried[start:start+len(p)]=burial(p,state['tree'])
    return dict(origin=origin,step=step,distance=distances.reshape(shape),buried=buried.reshape(shape))


def query_grid(points,grid):
    shape=points.shape[:-1];coords=((points.reshape(-1,3)-grid['origin'])/grid['step']).T
    inside=np.all((coords>=0)&(coords<=np.array(grid['distance'].shape)[:,None]-1),axis=0)
    distance=map_coordinates(grid['distance'],coords,order=1,mode='constant',cval=0.,prefilter=False)
    buried=map_coordinates(grid['buried'],coords,order=1,mode='constant',cval=0.,prefilter=False)
    # Outside-box poses cannot be certified by an absent grid value.
    return distance.reshape(shape),buried.reshape(shape),inside.reshape(shape)


def metrics(distances,buried):
    d=np.asarray(distances)
    overlap=np.mean(np.maximum(0.,SETTINGS['overlap_distance']-d)**2,axis=-1)
    severe=np.mean(d<SETTINGS['severe_distance'],axis=-1)
    core=np.mean(buried,axis=-1)
    return dict(score=core-overlap-2*severe,mean_overlap=overlap,
                severe_fraction=severe,buried_fraction=core,minimum_distance=d.min(axis=-1))


def profile_flags(row):
    return {name:bool(row['severe_fraction']<=a+1e-12 and row['mean_overlap']<=b and
                      row['buried_fraction']>=c) for name,(a,b,c) in PROFILES.items()}


def evaluate(candidate,chem,state,grid):
    kw=soft.run.kw;sets=kw.chemical_sets(chem,candidate.feature_types)
    seeds,_=kw.prepare_seeds(candidate,state['query'],**kw.POSE_PARAMETERS)
    matrices=np.array([s.transform_matrix for s in seeds]).reshape(-1,4,4)
    matrices=matrices[kw.direct_contact_masks(candidate.shape_points,sets,state,matrices)['all']]
    result=dict(direct_poses=len(matrices),exact_poses=0,best=None,
                profiles={p:False for p in PROFILES},profile_poses={},outside_grid_poses=0)
    if not len(matrices): return result
    moved=np.einsum('bij,nj->bni',matrices[:,:3,:3],candidate.shape_points)+matrices[:,None,:3,3]
    distance,core,inside=query_grid(moved,grid);approx=metrics(distance,core)
    valid=inside.all(axis=1);result['outside_grid_poses']=int((~valid).sum())
    choices=np.flatnonzero(valid)
    choices=choices[np.argsort(-approx['score'][choices],kind='stable')[:SETTINGS['exact_poses']]]
    for i in choices:
        p=moved[i];d=state['tree'].query(p)[0];row={k:float(v) for k,v in metrics(d,burial(p,state['tree'])).items()}
        row.update(gid=int(candidate.global_id),cid=candidate.conformer_id,
                   transform=matrices[i].ravel().tolist(),grid_distance_max_error=float(np.max(np.abs(d-distance[i]))))
        flags=profile_flags(row)
        for name,passed in flags.items():
            if passed and (name not in result['profile_poses'] or row['score']>result['profile_poses'][name]['score']):
                result['profile_poses'][name]=row
        result['profiles']={name:result['profiles'][name] or flags[name] for name in PROFILES}
        result['exact_poses']+=1
        if result['best'] is None or row['score']>result['best']['score']:result['best']=row
    return result


def controls(state,grid):
    from types import SimpleNamespace
    kw=soft.run.kw;q=state['query'];p=q['shape_points'];chem=SimpleNamespace(**state['record']['chemistry'])
    sets=kw.chemical_sets(chem,q['feature_types'])
    contacts=bool(kw.direct_contact_masks(p,sets,state,np.eye(4)[None])['all'][0])
    d=state['tree'].query(p)[0];row={k:float(v) for k,v in metrics(d,burial(p,state['tree'])).items()}
    a,_,_=query_grid(p,grid)
    shifted=p+(state['geometry']['receptor'][0]-p[0])
    negative_min=float(state['tree'].query(shifted)[0].min())
    return dict(native_metrics=row,native_profiles=profile_flags(row),native_contacts=contacts,
                grid_distance_max_error=float(np.max(np.abs(d-a))),
                forced_overlap_detected=negative_min<SETTINGS['severe_distance'])
