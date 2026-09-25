"""Conservative live-scene contact hypotheses; not interaction energies or PLIP."""
from collections import Counter, defaultdict
import csv
from pathlib import Path

import numpy as np

COLORS=dict(polar_contact='yellow',salt_bridge='magenta',pi_stacking='cyan',
            cation_pi='orange',hydrophobic='gray70',halogen_bond='green')
DEFAULT_TYPES=['polar_contact','salt_bridge','pi_stacking','cation_pi','halogen_bond']
RULES=dict(polar_max_A=3.5,salt_max_A=4.0,pi_max_A=5.5,pi_angle_deviation_deg=30,
           pi_offset_max_A=2.0,cation_pi_max_A=5.0,cation_pi_face_angle_max_deg=30,
           hydrophobic_max_A=4.0,halogen_max_A=3.5,halogen_donor_min_deg=150,
           halogen_acceptor_min_deg=90,halogen_acceptor_max_deg=150)


def residue(a):return (a.get('segi',''),a['chain'],a['resi'],a['resn'])
def distance(a,b):return float(np.linalg.norm(np.asarray(a)-np.asarray(b)))
def angle(a,b):
    norm=np.linalg.norm(a)*np.linalg.norm(b)
    return None if norm<1e-8 else float(np.degrees(np.arccos(np.clip(np.dot(a,b)/norm,-1,1))))


def aromatic_rings(atoms,edges,aromatic):
    graph=defaultdict(set)
    for a,b in edges:
        if a in aromatic and b in aromatic:graph[a].add(b);graph[b].add(a)
    cycles=set()
    def walk(path):
        for nxt in graph[path[-1]]:
            if nxt==path[0] and len(path) in (5,6):
                cycles.add(tuple(sorted(path)))
            elif nxt>path[0] and nxt not in path and len(path)<6:walk(path+[nxt])
    for start in sorted(graph):walk([start])
    result=[]
    for ids in sorted(cycles):
        points=np.array([atoms[i]['coord'] for i in ids]);center=points.mean(axis=0)
        _,s,vh=np.linalg.svd(points-center,full_matrices=False)
        if s[1]<.2 or np.max(np.abs((points-center)@vh[-1]))>.15:continue
        result.append(dict(ids=ids,center=center,normal=vh[-1]))
    return result


def detect(atoms,edges,ligand,protein,aromatic,polar_pairs):
    """All atom IDs are live object-local indices. No cross-complex pairs."""
    neighbors=defaultdict(set)
    for a,b in edges:neighbors[a].add(b);neighbors[b].add(a)
    rows=[]
    def add(kind,lids,pids,lp=None,pp=None,**geometry):
        lids=tuple(sorted(lids));pids=tuple(sorted(pids))
        lp=np.mean([atoms[i]['coord'] for i in lids],axis=0) if lp is None else lp
        pp=np.mean([atoms[i]['coord'] for i in pids],axis=0) if pp is None else pp
        rows.append(dict(kind=kind,ligand_atoms=list(lids),protein_atoms=list(pids),
            ligand_residue=list(residue(atoms[lids[0]])),protein_residue=list(residue(atoms[pids[0]])),
            ligand_point=np.asarray(lp).tolist(),protein_point=np.asarray(pp).tolist(),
            distance_A=distance(lp,pp),geometry=geometry))
    # Use only explicit formal charges. Missing charges do not imply neutral chemistry.
    for l in sorted(ligand):
        for p in sorted(protein):
            d=distance(atoms[l]['coord'],atoms[p]['coord'])
            if 2<=d<=4 and atoms[l]['charge']*atoms[p]['charge']<0:
                add('salt_bridge',[l],[p],charge_basis='explicit loaded formal charges')
    salt_pairs={(l,p) for r in rows for l in r['ligand_atoms'] for p in r['protein_atoms']}
    for l,p in sorted(set(polar_pairs)):
        if l not in ligand or p not in protein or (l,p) in salt_pairs:continue
        if 2.2<=distance(atoms[l]['coord'],atoms[p]['coord'])<=3.5:
            add('polar_contact',[l],[p],basis='PyMOL donor-to-acceptor directional match; protonation unvalidated')
    lr=aromatic_rings(atoms,edges,aromatic & ligand)
    pr=aromatic_rings(atoms,edges,aromatic & protein)
    for l in lr:
        for p in pr:
            delta=l['center']-p['center'];d=np.linalg.norm(delta)
            theta=angle(l['normal'],p['normal']);theta=min(theta,180-theta)
            offsets=[np.linalg.norm(delta-np.dot(delta,n)*n) for n in (l['normal'],p['normal'])]
            if 3<=d<=5.5 and (theta<=30 or theta>=60) and min(offsets)<=2:
                add('pi_stacking',l['ids'],p['ids'],l['center'],p['center'],angle_deg=theta,offset_A=float(min(offsets)))
    for rings,charged,ring_is_ligand in ((lr,protein,True),(pr,ligand,False)):
        for ring in rings:
            for i in sorted(charged):
                if atoms[i]['charge']<=0:continue
                delta=np.array(atoms[i]['coord'])-ring['center'];d=np.linalg.norm(delta)
                theta=angle(delta,ring['normal'])
                if theta is None:continue
                theta=min(theta,180-theta)
                if 2.5<=d<=5 and theta<=30:
                    add('cation_pi',ring['ids'] if ring_is_ligand else [i],
                        [i] if ring_is_ligand else ring['ids'],face_angle_deg=theta,charge_basis='explicit loaded formal charge')
    # Halogen donor must be bonded to carbon; acceptor typing supplied by PyMOL.
    for l in sorted(ligand):
        if atoms[l]['elem'] not in {'Cl','Br','I'}:continue
        carbons=[i for i in neighbors[l] if atoms[i]['elem']=='C']
        if len(carbons)!=1:continue
        x=np.array(atoms[l]['coord']);c=np.array(atoms[carbons[0]]['coord'])
        for p in sorted(protein):
            if not atoms[p].get('acceptor') or atoms[p]['elem'] not in {'N','O','S'}:continue
            a=np.array(atoms[p]['coord']);d=distance(x,a)
            donor_angle=angle(c-x,a-x)
            heavy=[i for i in neighbors[p] if atoms[i]['elem']!='H']
            if len(heavy)!=1 or donor_angle is None:continue
            acceptor_angle=angle(x-a,np.array(atoms[heavy[0]]['coord'])-a)
            if 2.5<=d<=3.5 and donor_angle>=150 and acceptor_angle is not None and 90<=acceptor_angle<=150:
                add('halogen_bond',[l],[p],donor_angle_deg=donor_angle,acceptor_angle_deg=acceptor_angle)
    pi_pairs={(l,p) for r in rows if r['kind']=='pi_stacking' for l in r['ligand_atoms'] for p in r['protein_atoms']}
    def nonpolar(i):
        return atoms[i]['elem']=='C' and atoms[i]['charge']==0 and bool(neighbors[i]) and all(atoms[n]['elem'] in {'C','H'} for n in neighbors[i])
    for l in sorted(i for i in ligand if nonpolar(i)):
        for p in sorted(i for i in protein if nonpolar(i)):
            if (l,p) not in pi_pairs and 3<=distance(atoms[l]['coord'],atoms[p]['coord'])<=4:
                add('hydrophobic',[l],[p],basis='nonpolar carbon neighborhood and distance')
    # Display representatives, not every atom pair. Full accepted rows remain exported.
    representatives={}
    for row in sorted(rows,key=lambda r:(r['distance_A'],r['ligand_atoms'],r['protein_atoms'])):
        key=(row['kind'],tuple(row['ligand_residue']),tuple(row['protein_residue']))
        representatives.setdefault(key,row)
    return rows,list(representatives.values())


def display(cmd,message,catalog,root,checkpoint=True):
    from .pymol_bridge import selection,write_json
    from .pymol_program import _UNDO
    root=Path(root);view=message['view'];types=view.get('types',DEFAULT_TYPES)
    enabled=set(cmd.get_names('objects',enabled_only=1))
    ids=view.get('objects') or [r['id'] for r in catalog if r['id'] in enabled]
    if not ids:raise ValueError('Enable a complex or specify its object ID first')
    if set(ids)-{r['id'] for r in catalog} or set(ids)-set(cmd.get_names('objects')):raise ValueError('Unknown or unloaded complex object')
    selected={r['id']:r for r in catalog if r['id'] in ids}
    results={}
    # Complete all detection before touching display state.
    for obj,row in selected.items():
        def indices(sel):return {a.index for a in cmd.get_model(sel,state=1).atom}
        lig=selection(obj,row,dict(target='ligand'))
        pro=f'({obj} and polymer.protein) and not ({lig})'
        ligand=indices(lig);protein=indices(f'({pro}) within 7 of ({lig})')
        if not ligand or not protein:raise ValueError('No ligand or nearby protein atoms for '+obj)
        model=cmd.get_model(obj,state=1);acceptors=indices(obj+' and acceptors')
        atoms={a.index:dict(coord=list(a.coord),elem=a.symbol,charge=getattr(a,'formal_charge',0),
            name=a.name,segi=getattr(a,'segi',''),chain=a.chain,resi=a.resi,resn=a.resn,acceptor=a.index in acceptors) for a in model.atom}
        # Exclude alternate atoms unless blank or A; record this explicit display policy.
        admitted={a.index for a in model.atom if getattr(a,'alt','') in ('','A') and getattr(a,'q',1)>0}
        ligand &= admitted;protein &= admitted
        edges=[(model.atom[b.index[0]].index,model.atom[b.index[1]].index) for b in model.bond]
        polar=[]
        for left,right,reverse in ((lig+' and donors',pro+' and acceptors',False),
                                   (pro+' and donors',lig+' and acceptors',True)):
            for a,b in cmd.find_pairs(left,right,state1=1,state2=1,cutoff=3.5,mode=1,angle=45):
                if a[0]!=obj or b[0]!=obj:continue
                polar.append((b[1],a[1]) if reverse else (a[1],b[1]))
        raw,reps=detect(atoms,edges,ligand,protein,indices(obj+' and aromatic'),polar)
        results[obj]=dict(accepted=raw,representatives=reps,
            chemistry=dict(ligand_formal_charge_atoms=sum(atoms[i]['charge']!=0 for i in ligand),
                           protein_formal_charge_atoms=sum(atoms[i]['charge']!=0 for i in protein)),
            atom_metadata={str(i):{k:v for k,v in a.items() if k!='coord'} for i,a in atoms.items() if i in ligand|protein})
    before=cmd.get_session();artifacts={};prefix=root/message['id']
    if checkpoint:
        cmd.save(str(prefix)+'-before.pse');artifacts['checkpoint']=str(prefix)+'-before.pse'
    try:
        # Replace previous proximity displays, including generated AI distance objects.
        for name in cmd.get_names('objects'):
            if name.startswith('ai_') and cmd.get_type(name)=='object:measurement':cmd.disable(name)
            if any(name.startswith(r['id']+'_ix_') or name in (r['id']+'_contacts',r['id']+'_polar_contacts') for r in catalog):cmd.disable(name)
        for obj in ids:
            for kind in COLORS:cmd.delete(obj+'_ix_'+kind)
            points=obj+'_ix_points';cmd.delete(points)
            for n,r in enumerate(results[obj]['representatives']):
                if r['kind'] not in types:continue
                endpoints=[]
                for side in ('ligand','protein'):
                    name=f'p{n}_{side}'
                    cmd.pseudoatom(points,name=name,pos=r[side+'_point'],state=1)
                    endpoints.append(f'{points} and name {name}')
                name=obj+'_ix_'+r['kind']
                cmd.distance(name,*endpoints,cutoff=8,mode=0,state=1)
                cmd.set('dash_color',COLORS[r['kind']],name);cmd.hide('labels',name)
            cmd.hide('everything',points)
        report=dict(status='complete',operation='typed_interactions',objects=ids,legend={k:COLORS[k] for k in types},
            displayed_counts={o:dict(Counter(r['kind'] for r in v['representatives'] if r['kind'] in types)) for o,v in results.items()},
            results=results,rules=RULES,scope='Conservative geometric hypotheses using live state 1 and PyMOL chemical typing. Not PLIP or confirmed physical interactions.',
            display_policy='One nearest representative per type and ligand/protein residue pair; hydrophobic contacts opt-in. Alternate locations blank/A only.',
            limitations=['Missing formal charges suppress salt/cation-pi assignments; no protonation inference.',
                'Missing aromatic bond typing suppresses ring assignments.',
                'Water bridges and metal coordination not evaluated; separate chemistry is required.'])
        path=Path(str(prefix)+'-interactions.json');write_json(path,report);artifacts['interactions_json']=str(path)
        table=Path(str(prefix)+'-interactions.csv')
        with table.open('w',newline='',encoding='utf-8') as handle:
            writer=csv.writer(handle);writer.writerow(['object','type','ligand_residue','protein_residue','ligand_atoms','protein_atoms','distance_A','displayed'])
            for obj,v in results.items():
                for r in v['accepted']:writer.writerow([obj,r['kind'],r['ligand_residue'],r['protein_residue'],r['ligand_atoms'],r['protein_atoms'],r['distance_A'],r in v['representatives'] and r['kind'] in types])
        artifacts['interactions_csv']=str(table)
        if checkpoint:
            cmd.png(str(prefix)+'.png',width=1400,height=1000,ray=0,quiet=1);cmd.save(str(prefix)+'.pse')
            artifacts.update(snapshot=str(prefix)+'.png',save_session=str(prefix)+'.pse')
            _UNDO[str(root.resolve())]=before
        return {k:v for k,v in report.items() if k!='results'}|dict(artifacts=artifacts)
    except Exception:
        cmd.set_session(before)
        raise
