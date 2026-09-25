"""Local, file-queued PyMOL commands. No model-provided Python or selections."""
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re
import subprocess
import threading
import time
import uuid

OPERATIONS = {'open','cartoon','sticks','surface','hide_surface','polar_contacts',
              'contacts','hide_contacts','zoom','color','label_residues','hide_labels',
              'show','hide','snapshot','save_session','status','align','rotate','background','transparency',
              'chains','remove_chain','pocket_view','interaction_overview'}
TARGETS = {'all','protein','ligand','pocket','water'}
COLORS = {'cyan','green','yellow','orange','magenta','white','gray','red','blue','marine','salmon'}


def write_json(path, value):
    path = Path(path)
    temp = path.with_suffix('.tmp')
    temp.write_text(json.dumps(value, indent=2, allow_nan=False), encoding='utf-8')
    temp.replace(path)


def validate_view(value):
    if not isinstance(value, dict) or set(value)-{'operation','objects','target','color','cutoff','chain','residue','reference','collection','axis','angle','opacity','scheme','radius'}:
        raise ValueError('Unsupported viewer fields')
    if value.get('operation') not in OPERATIONS:
        raise ValueError('Unsupported PyMOL operation')
    if value.get('target','all') not in TARGETS:
        raise ValueError('Unsupported viewer target')
    if value.get('scheme','solid') not in {'solid','element','chain','rainbow'}:
        raise ValueError('Unsupported color scheme')
    radius=value.get('radius',5)
    if type(radius) not in (int,float) or not math.isfinite(radius) or not 1<=radius<=12:
        raise ValueError('Pocket radius must be 1..12 angstroms')
    if value['operation']=='remove_chain' and (not value.get('objects') or not value.get('chain')):
        raise ValueError('Specify objects and the exact protein chain to remove; use chains first')
    if value.get('collection','diverse') not in {'diverse','all_admitted'}:raise ValueError('Unknown structure collection')
    if value.get('axis','y') not in {'x','y','z'}:raise ValueError('Invalid rotation axis')
    for field,low,high,default in [('angle',-360,360,30),('opacity',0,1,.5)]:
        v=value.get(field,default)
        if type(v) not in (int,float) or not math.isfinite(v) or not low<=v<=high:raise ValueError('Invalid '+field)
    if 'color' in value and value['color'] not in COLORS:
        raise ValueError('Unsupported color')
    ids = value.get('objects', [])
    if not isinstance(ids,list) or len(ids)>100 or any(not isinstance(i,str) or not re.fullmatch(r'v\d{3}',i) for i in ids):
        raise ValueError('Use object IDs from the viewer catalog')
    if 'reference' in value and not re.fullmatch(r'v\d{3}',str(value['reference'])):
        raise ValueError('Invalid reference object')
    for field, pattern in [('chain',r'[A-Za-z0-9]{1,8}'),('residue',r'-?\d+[A-Za-z]?')]:
        if field in value and (not isinstance(value[field],str) or not re.fullmatch(pattern,value[field])):
            raise ValueError('Invalid '+field)
    cutoff = value.get('cutoff',3.6)
    if type(cutoff) not in (int,float) or not math.isfinite(cutoff) or not 1<=cutoff<=8:
        raise ValueError('Contact cutoff must be 1..8 angstroms')
    return value


def selection(obj, row, view):
    target = view.get('target','all')
    ligand = row.get('ligand', {})
    lig = f'({obj} and organic)'
    if ligand:
        # These identifiers are derived from the sealed artifact, never arbitrary selections.
        for key, pattern in [('chain',r'[A-Za-z0-9]{1,8}'),('resi',r'-?\d+[A-Za-z]?'),('resn',r'[A-Za-z0-9_]{1,12}')]:
            if key in ligand and not re.fullmatch(pattern,str(ligand[key])):
                raise ValueError('Unsupported structure identifier: '+key)
        terms = [obj] + [f'{key} {ligand[key]}' for key in ('chain','resi','resn') if key in ligand]
        lig = '('+' and '.join(terms)+')'
    protein = f'({obj} and polymer.protein)'
    sel = {'all':obj,'protein':protein,'ligand':lig,'water':f'({obj} and solvent)',
           'pocket':f'({protein} and byres ({protein} within {view.get("radius",5)} of {lig}))'}[target]
    if view.get('chain'): sel = f'({sel} and chain {view["chain"]})'
    if view.get('residue'): sel = f'({sel} and resi {view["residue"]})'
    return sel


def element_colors(cmd, sel, carbon='green'):
    cmd.color(carbon,sel)
    for element,color in [('N','blue'),('O','red'),('S','yellow'),('P','orange'),('H','white'),('F','cyan'),('Cl','green'),('Br','salmon'),('I','violet')]:
        cmd.color(color,f'({sel}) and elem {element}')


def residue_rows(cmd,sel):
    return [dict(chain=c,residue=i,resname=n) for c,i,n in sorted({
        (a.chain,a.resi,a.resn) for a in cmd.get_model(sel,state=1).atom})]


def execute_command(cmd, message, catalog, root):
    """Run against PyMOL's cmd API; called only by the desktop bridge."""
    view = validate_view(message['view']); op = view['operation']
    rows = {r['id']:r for r in catalog}
    ids = view.get('objects') or list(rows)
    if set(ids)-rows.keys(): raise ValueError('Unknown viewer object')
    output = dict(status='complete',operation=op,objects=ids,artifacts={})
    if op == 'open':
        for obj in ids:
            row = rows[obj]; path = Path(row['path']).resolve()
            if hashlib.sha256(path.read_bytes()).hexdigest()!=row['sha256']:
                raise ValueError('Structure changed since viewer request')
            cmd.delete(obj)
            cmd.load(str(path),obj)
            if row.get('transform'):
                import numpy as np
                matrix = np.asarray(row['transform'], dtype=float)
                if matrix.shape!=(4,4) or not np.isfinite(matrix).all(): raise ValueError('Invalid alignment')
                xyz = cmd.get_coords(obj, state=1)
                cmd.load_coords(xyz @ matrix[:3,:3].T + matrix[:3,3],obj,state=1)
            cmd.hide('everything',obj)
            cmd.show('cartoon',selection(obj,row,dict(target='protein')))
            cmd.show('sticks',selection(obj,row,dict(target='ligand')))
            cmd.color(['cyan','green','orange','magenta','yellow'][list(rows).index(obj)%5],obj)
            element_colors(cmd,selection(obj,row,dict(target='ligand')))
        cmd.zoom(' or '.join(ids))
        output['alignment'] = 'Recorded transforms applied where available; other structures remain in their original frame'
        return output
    if op == 'status':
        output['loaded_objects']=cmd.get_names('objects'); return output
    loaded = set(cmd.get_names('objects'))
    if set(ids)-loaded: raise ValueError('Open the selected structures before controlling them')
    if op=='interaction_overview':
        results={}
        for name,operation in [('pocket','pocket_view'),('nearby','contacts'),('polar','polar_contacts')]:
            child=dict(view,operation=operation)
            child.pop('cutoff',None)
            result=execute_command(cmd,dict(id=message['id']+'-'+name,view=child),catalog,root)
            results[name]=result
            output['artifacts'].update({name+'_'+k:v for k,v in result['artifacts'].items()})
        output['results']=results
        output['scope']='Pocket residues (default 5 A), heavy-atom proximity (4.5 A), polar candidates (3.6 A). Categories overlap; do not sum counts.'
        output['not_evaluated']=['validated hydrogen bonds','salt bridges','pi stacking','cation-pi','halogen bonds','water bridges','metal coordination','interaction energies']
        return output
    if op=='chains':
        output['chains']={}
        for obj in ids:
            protein=selection(obj,rows[obj],dict(target='protein'))
            nearby={a['chain'] for a in residue_rows(cmd,selection(obj,rows[obj],dict(target='pocket',radius=view.get('radius',5))))}
            output['chains'][obj]=[dict(chain=c,near_ligand=c in nearby) for c in sorted({a['chain'] for a in residue_rows(cmd,protein)})]
        output['limitation']='Near-ligand geometry is not proof of biological importance or redundancy. Specify a chain explicitly before removal.'
        return output
    if op=='remove_chain':
        selections={obj:selection(obj,rows[obj],dict(target='protein',chain=view['chain'])) for obj in ids}
        if any(not cmd.count_atoms(s) for s in selections.values()):raise ValueError('Requested protein chain is absent in one or more objects')
        output['removed_atoms']={obj:cmd.count_atoms(s) for obj,s in selections.items()}
        for obj,sel in selections.items():
            cmd.remove(sel)
            cmd.delete(obj+'_contacts');cmd.delete(obj+'_polar_contacts')
        output['scope']='Protein atoms removed from the displayed copy only; source files and ligand are preserved. Reopen the task to restore.'
        return output
    if op == 'align':
        ref = view.get('reference')
        if ref not in rows or ref not in loaded: raise ValueError('Choose a loaded reference object ID')
        output['alignments']={obj:cmd.align(f'{obj} and polymer.protein and name CA',
                                          f'{ref} and polymer.protein and name CA') for obj in ids if obj!=ref}
        output['scope']='Display-only CA alignment; does not establish same-pocket admission or change source coordinates'
        return output
    all_contacts=[]
    for obj in ids:
        row=rows[obj]; sel=selection(obj,row,view)
        if op in {'cartoon','sticks','surface'}: cmd.show(op,sel)
        elif op=='hide_surface': cmd.hide('surface',sel)
        elif op=='transparency':cmd.set('transparency',1-view.get('opacity',.5),sel)
        elif op=='color':
            scheme=view.get('scheme','solid')
            if scheme=='element':element_colors(cmd,sel,view.get('color','green'))
            elif scheme=='rainbow':cmd.spectrum('count','rainbow',sel)
            elif scheme=='chain':
                for i,c in enumerate(sorted({a['chain'] for a in residue_rows(cmd,sel)})):
                    if not re.fullmatch(r'[A-Za-z0-9]{1,8}',c):raise ValueError('Unsupported chain identifier')
                    cmd.color(['cyan','green','orange','magenta','yellow','salmon'][i%6],f'({sel}) and chain {c}')
            else:cmd.color(view.get('color','cyan'),sel)
        elif op=='pocket_view':
            lig=selection(obj,row,dict(target='ligand'))
            pocket=selection(obj,row,dict(target='pocket',radius=view.get('radius',5),**{k:view[k] for k in ('chain','residue') if k in view}))
            if not cmd.count_atoms(lig):raise ValueError('No ligand atoms for '+obj)
            cmd.hide('sticks',selection(obj,row,dict(target='protein')))
            cmd.show('cartoon',selection(obj,row,dict(target='protein')))
            cmd.show('sticks',pocket);cmd.show('sticks',lig)
            element_colors(cmd,pocket,'cyan');element_colors(cmd,lig,'green')
            cmd.label(f'({pocket}) and name CA','"%s%s/%s" % (resn,resi,chain)')
            output.setdefault('residues',{})[obj]=residue_rows(cmd,pocket)
        elif op=='show': cmd.enable(obj)
        elif op=='hide': cmd.disable(obj)
        elif op=='label_residues': cmd.label(f'({sel} and name CA)','"%s%s/%s" % (resn,resi,chain)')
        elif op=='hide_labels': cmd.hide('labels',sel)
        elif op=='hide_contacts':
            cmd.delete(obj+'_contacts');cmd.delete(obj+'_polar_contacts')
        elif op in {'polar_contacts','contacts'}:
            lig=selection(obj,row,dict(target='ligand'))
            protein=selection(obj,row,dict(target='protein',**{k:view[k] for k in ('chain','residue') if k in view}))
            if not cmd.count_atoms(lig) or not cmd.count_atoms(protein):
                raise ValueError('No ligand/protein atoms for '+obj)
            cutoff=view.get('cutoff',3.6 if op=='polar_contacts' else 4.5)
            mode=1 if op=='polar_contacts' else 0
            left=lig+' and not hydro'; right=protein+' and not hydro'
            if mode:
                left+=' and (donors or acceptors)'; right+=' and (donors or acceptors)'
            pairs=cmd.find_pairs(left,right,state1=1,state2=1,cutoff=cutoff,mode=mode,angle=45)
            contact_object=obj+('_polar_contacts' if mode else '_contacts')
            cmd.delete(contact_object)
            for a,b in pairs:
                sa=f'{a[0]} and index {a[1]}'; sb=f'{b[0]} and index {b[1]}'
                distance=cmd.get_distance(sa,sb,state=1)
                cmd.distance(contact_object,sa,sb)
                aa=cmd.get_model(sa,state=1).atom[0]; bb=cmd.get_model(sb,state=1).atom[0]
                all_contacts.append(dict(object=obj,query_id=row['label'],ligand_atom=aa.name,
                    ligand_chain=aa.chain,ligand_residue=aa.resi,protein_atom=bb.name,
                    protein_chain=bb.chain,protein_residue=bb.resi,protein_resname=bb.resn,
                    distance_angstrom=distance,kind='pymol_polar_candidate' if mode else 'heavy_atom_proximity'))
            cmd.show('sticks',selection(obj,row,dict(target='pocket')))
            if pairs:
                cmd.set('dash_color','yellow' if mode else 'gray',contact_object)
                cmd.hide('labels',contact_object)
    if op=='pocket_view':
        cmd.zoom(' or '.join(selection(obj,rows[obj],dict(target='ligand'))+' or '+selection(obj,rows[obj],dict(target='pocket',radius=view.get('radius',5))) for obj in ids))
        path=root/(message['id']+'-residues.json');write_json(path,output.get('residues',{}))
        output['artifacts']['residues_json']=str(path)
    if op=='zoom': cmd.zoom(' or '.join(selection(obj,rows[obj],view) for obj in ids))
    if op=='rotate':cmd.turn(view.get('axis','y'),view.get('angle',30))
    if op=='background':cmd.bg_color(view.get('color','white'))
    if op in {'snapshot','save_session'}:
        suffix='.png' if op=='snapshot' else '.pse'
        path=root/(message['id']+suffix)
        if op=='snapshot': cmd.png(str(path),width=1600,height=1000,dpi=150,ray=0,quiet=1)
        else: cmd.save(str(path))
        output['artifacts'][op]=str(path)
    if op in {'polar_contacts','contacts'}:
        path=root/(message['id']+'-contacts.json')
        write_json(path,all_contacts)
        csv_path=path.with_suffix('.csv')
        if all_contacts:
            with csv_path.open('w',newline='',encoding='utf-8') as handle:
                writer=csv.DictWriter(handle,fieldnames=list(all_contacts[0]));writer.writeheader();writer.writerows(all_contacts)
            output['artifacts']['contacts_csv']=str(csv_path)
        output.update(contact_count=len(all_contacts),examples=all_contacts[:30],
                      limitation='PyMOL geometry/atom typing only; protonation and chemistry are unvalidated. Not all physical interactions or confirmed hydrogen bonds.')
        output['artifacts']['contacts_json']=str(path)
    return output


def start_bridge(queue):
    """Called by a locally generated PyMOL startup script, inside the GUI process."""
    from pymol import cmd
    root=Path(queue).resolve()
    def heartbeat():
        while True:
            write_json(root/'heartbeat.json',dict(pid=os.getpid(),time=time.time()))
            time.sleep(1)
    def loop():
        managed=set()
        while True:
            for path in sorted(root.glob('*.request.json')):
                result_path=path.with_name(path.name.replace('.request.json','.result.json'))
                if result_path.exists(): continue
                message=json.loads(path.read_text(encoding='utf-8'))
                try:
                    if not re.fullmatch(r'\d+-[a-f0-9]{12}',message['id']): raise ValueError('Invalid operation ID')
                    if message['view']['operation']=='open':
                        for obj in managed:
                            cmd.delete(obj)
                            cmd.delete(obj+'_contacts')
                            cmd.delete(obj+'_polar_contacts')
                        managed={r['id'] for r in message['catalog']}
                    result=execute_command(cmd,message,message['catalog'],root)
                except Exception as exc:
                    result=dict(status='failed',error=str(exc))
                write_json(result_path,dict(id=message['id'],**result))
            time.sleep(.3)
    threading.Thread(target=heartbeat,daemon=True).start()
    threading.Thread(target=loop,daemon=True).start()


def bridge_status(root):
    root=Path(root)
    heartbeat=root/'heartbeat.json'
    try: alive=time.time()-json.loads(heartbeat.read_text())['time']<10
    except (OSError,ValueError,KeyError): alive=False
    receipts=sorted(root.glob('*.result.json'))
    return dict(connected=alive,latest=json.loads(receipts[-1].read_text(encoding='utf-8')) if receipts else None)


def submit(root, view, catalog, executable=None):
    validate_view(view);root=Path(root).resolve();root.mkdir(parents=True,exist_ok=True)
    if not bridge_status(root)['connected']:
        if not executable: raise ValueError('Configure runtime.pymol.executable with the installed desktop PyMOL executable')
        launch=root/'startup.py'
        launch.write_text('import sys\nsys.path.insert(0, '+repr(str(Path(__file__).resolve().parents[1]))+')\n'
                          'from aidd_agent.pymol_bridge import start_bridge\nstart_bridge('+repr(str(root))+')\n',encoding='utf-8')
        # The user explicitly requests a visible PyMOL window. No shell is involved.
        with (root/'pymol.log').open('ab') as log:
            process=subprocess.Popen([str(executable),'-r',str(launch)],stdout=log,stderr=log,shell=False)
        deadline=time.monotonic()+8
        while not bridge_status(root)['connected'] and time.monotonic()<deadline:
            if process.poll() is not None: break
            time.sleep(.2)
        if not bridge_status(root)['connected']:
            raise RuntimeError('PyMOL did not connect. Check desktop DISPLAY, installation and '+str(root/'pymol.log'))
    identifier=str(time.time_ns())+'-'+uuid.uuid4().hex[:12]
    path=root/(identifier+'.request.json')
    write_json(path,dict(id=identifier,view=view,catalog=catalog))
    return dict(status='queued',operation_id=identifier,result=str(root/(identifier+'.result.json')),
                message='Queued in desktop PyMOL; check /pymol_status for execution result')
