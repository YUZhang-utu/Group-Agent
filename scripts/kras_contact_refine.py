"""E120 contact-constrained rigid refinement; internal molecular geometry is fixed."""
from pathlib import Path
import hashlib
import numpy as np
from scipy.optimize import least_squares
import kras_soft_rank as soft

SETTINGS = dict(version='E120-rigid-contact-v1', minimum_distance=2.0,
                optimization_distance=2.4, contact_margin=.05, angle_bound=.65,
                translation_bound=4., max_nfev=80, starts=3)


def digest():
    return hashlib.sha256(Path(__file__).read_bytes().replace(b'\r\n',b'\n')).hexdigest()


def rotation_derivatives(angles):
    x,y,z=angles; cx,cy,cz=np.cos(angles); sx,sy,sz=np.sin(angles)
    rx=np.array([[1,0,0],[0,cx,-sx],[0,sx,cx]])
    ry=np.array([[cy,0,sy],[0,1,0],[-sy,0,cy]])
    rz=np.array([[cz,-sz,0],[sz,cz,0],[0,0,1]])
    dx=np.array([[0,0,0],[0,-sx,-cx],[0,cx,-sx]])
    dy=np.array([[-sy,0,cy],[0,0,0],[-cy,0,-sy]])
    dz=np.array([[-sz,-cz,0],[cz,-sz,0],[0,0,0]])
    return rz@ry@rx, (rz@ry@dx,rz@dy@rx,dz@ry@rx)


def contact_groups(sets,state):
    atoms=state['record']['receptor_atoms']
    def option(kind,residue,names,limit):
        partners=np.array([a['xyz'] for a in atoms if a['residue']==residue and a['atom'] in names])
        if not len(sets[kind]) or not len(partners):
            raise ValueError('Missing mandatory contact atoms')
        return sets[kind],partners,limit
    p=soft.run.kw.POLICY['polar_distance'];h=soft.run.kw.POLICY['hydrophobic_distance']
    groups=[('His95',[option('acceptor',95,{'ND1','NE2'},p)]),
            ('Asp69',[option('donor',69,{'OD1','OD2'},p)]),
            ('Glu62',[option('donor',62,{'OE1','OE2'},p)]),
            ('Met72',[option('carbon',72,{'CB','CG','CE'},h)])]
    polar=option('acceptor',68,{'NE','NH1','NH2'},p)
    if state['record']['target']=='KRAS_G12D':
        groups.extend([('Arg68_polar',[polar]),('Asp12',[option('donor',12,{'OD1','OD2'},p)])])
    else:
        groups.append(('Arg68',[polar,option('carbon',68,{'CB','CG','CD'},h)]))
    # G12V Val12 exclusion is covered by the whole-receptor floor and rechecked
    # using the unchanged direct_contact_masks before acceptance.
    return groups


def contact_margins(points,groups):
    values=[]; gradients=[]
    for name,options in groups:
        best=None
        for indices,partners,limit in options:
            delta=points[indices,None,:]-partners[None,:,:]
            distance=np.linalg.norm(delta,axis=2)
            i,j=np.unravel_index(np.argmin(distance),distance.shape)
            margin=float(distance[i,j]-limit)
            if best is None or margin<best[0]:
                best=(margin,int(indices[i]),delta[i,j]/max(float(distance[i,j]),1e-10))
        values.append(best[0]);gradients.append((best[1],best[2]))
    return np.array(values),gradients


class Objective:
    def __init__(self,points,initial,state,sets):
        self.points=np.asarray(points);self.initial=np.array(initial).reshape(4,4)
        self.base=self.points@self.initial[:3,:3].T+self.initial[:3,3]
        self.center=self.base.mean(axis=0);self.relative=self.base-self.center
        self.state=state;self.sets=sets;self.groups=contact_groups(sets,state)
        self.best=None;self.calls=0

    def calculate(self,x,remember=True):
        rotation,derivatives=rotation_derivatives(x[:3])
        points=self.relative@rotation.T+self.center+x[3:]
        jac=np.zeros((len(points),3,6))
        for k,d in enumerate(derivatives):jac[:,:,k]=self.relative@d.T
        jac[:,:,3:]=np.eye(3)
        distances,indices=self.state['tree'].query(points)
        delta=points-self.state['geometry']['receptor'][indices]
        unit=delta/np.maximum(distances[:,None],1e-10)
        clash=np.maximum(SETTINGS['optimization_distance']-distances,0)
        clash_jac=-np.einsum('ni,nij->nj',unit,jac)
        clash_jac[clash==0]=0
        margins,contact_derivatives=contact_margins(points,self.groups)
        residual=np.maximum(margins+SETTINGS['contact_margin'],0)
        contact_jac=np.array([direction@jac[index] for index,direction in contact_derivatives])
        contact_jac[residual==0]=0
        values=np.concatenate([clash,residual]);gradient=np.vstack([clash_jac,contact_jac])
        increment=np.eye(4);increment[:3,:3]=rotation
        increment[:3,3]=self.center+x[3:]-rotation@self.center
        matrix=increment@self.initial
        minimum=float(distances.min());cost=float(values@values)
        if remember:
            self.calls+=1
            if minimum>=SETTINGS['minimum_distance'] and np.all(margins<=0):
                valid=soft.run.kw.direct_contact_masks(self.points,self.sets,self.state,matrix[None])['all'][0]
                if valid and (self.best is None or cost<self.best['cost']):
                    self.best=dict(matrix=matrix.copy(),minimum_distance=minimum,cost=cost,
                                   contact_margins={g[0]:float(v) for g,v in zip(self.groups,margins)})
        return values,gradient,matrix,minimum,margins


def optimize_pose(points,initial,state,sets,max_nfev=None):
    objective=Objective(points,initial,state,sets)
    start=objective.calculate(np.zeros(6))
    bounds=np.array([SETTINGS['angle_bound']]*3+[SETTINGS['translation_bound']]*3)
    fit=least_squares(lambda x:objective.calculate(x)[0],np.zeros(6),
        jac=lambda x:objective.calculate(x,remember=False)[1],bounds=(-bounds,bounds),
        max_nfev=max_nfev or SETTINGS['max_nfev'],ftol=1e-5,xtol=1e-5,gtol=1e-5)
    final=objective.calculate(fit.x)
    evidence=dict(initial_minimum_distance=start[3],endpoint_minimum_distance=final[3],
        endpoint_contact_margins={g[0]:float(v) for g,v in zip(objective.groups,final[4])},
        optimizer_success=bool(fit.success),nfev=int(fit.nfev),objective_calls=objective.calls,
        start_transform=np.asarray(initial).ravel().tolist(),endpoint_transform=final[2].ravel().tolist())
    if objective.best is None:return None,evidence
    return objective.best['matrix'],evidence


def score_matrix(candidate,features,chemistry,state,matrix):
    """Use unchanged E118 scoring on one already recovered pose, without new seeds."""
    from types import SimpleNamespace
    original=soft.run.kw.prepare_seeds
    try:
        soft.run.kw.prepare_seeds=lambda *args,**kwargs:([SimpleNamespace(transform_matrix=matrix.ravel())],0)
        return soft.score_candidate(candidate,features,chemistry,state)[0]
    finally:
        soft.run.kw.prepare_seeds=original


def refine_candidate(candidate,features,chemistry,state,starts=3,max_nfev=80,initials=None):
    sets=soft.run.kw.chemical_sets(chemistry,candidate.feature_types)
    stats=dict(seeds=0,direct_poses=0,optimized_starts=0,geometry_conformers=0)
    if not all(len(v) for v in sets.values()):return None,stats
    if initials is None:
        seeds,_=soft.run.kw.prepare_seeds(candidate,state['query'],**soft.run.kw.POSE_PARAMETERS)
        matrices=np.array([s.transform_matrix for s in seeds]).reshape(-1,4,4)
    else:matrices=np.asarray(initials).reshape(-1,4,4)
    stats['seeds']=len(matrices)
    matrices=matrices[soft.run.kw.direct_contact_masks(candidate.shape_points,sets,state,matrices)['all']]
    stats['direct_poses']=len(matrices)
    if not len(matrices):return None,stats
    moved=np.einsum('bij,nj->bni',matrices[:,:3,:3],candidate.shape_points)+matrices[:,None,:3,3]
    distances=state['tree'].query(moved.reshape(-1,3))[0].reshape(len(matrices),-1)
    cost=(np.maximum(SETTINGS['optimization_distance']-distances,0)**2).sum(axis=1)
    chosen=np.argsort(cost,kind='stable')[:starts]
    best=None
    for j in chosen:
        matrix,evidence=optimize_pose(candidate.shape_points,matrices[j],state,sets,max_nfev)
        stats['optimized_starts']+=1
        if matrix is None:continue
        pose=score_matrix(candidate,features,chemistry,state,matrix)
        if pose is None:raise ValueError('Optimizer feasible incumbent failed independent contact/geometry check')
        pose['refinement']=evidence
        if best is None or pose['score']>best['score']:best=pose
    # A feasible unrefined seed must not be lost due to the optimization cap.
    for j in np.flatnonzero(distances.min(axis=1)>=SETTINGS['minimum_distance']):
        pose=score_matrix(candidate,features,chemistry,state,matrices[j])
        if pose is not None and (best is None or pose['score']>best['score']):
            best=pose;best['refinement']=dict(status='already_feasible',nfev=0)
    stats['geometry_conformers']=int(best is not None)
    return best,stats
