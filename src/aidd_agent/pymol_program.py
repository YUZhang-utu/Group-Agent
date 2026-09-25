"""Interpret bounded literal PyMOL API programs without Python eval/exec."""
import ast
import json
import math
from pathlib import Path
import re

# Values are positional API parameter names, with required argument counts.
METHODS = {
    'show': (['representation','selection'],1),
    'hide': (['representation','selection'],1),
    'color': (['color','selection'],1),
    'spectrum': (['expression','palette','selection','minimum','maximum'],2),
    'set': (['name','value','selection'],2),
    'set_color': (['name','rgb'],2),
    'select': (['name','selection'],2),
    'remove': (['selection'],1),
    'delete': (['name'],1),
    'label': (['selection','expression'],2),
    'zoom': (['selection','buffer'],0),
    'orient': (['selection'],0),
    'center': (['selection'],0),
    'turn': (['axis','angle'],2),
    'move': (['axis','distance'],2),
    'bg_color': (['color'],1),
    'count_atoms': (['selection'],0),
    'get_distance': (['atom1','atom2'],2),
    'distance': (['name','selection1','selection2','cutoff','mode'],3),
    'align': (['mobile','target','cutoff','cycles'],2),
    'super': (['mobile','target','cutoff','cycles'],2),
    'dss': (['selection'],0),
    'rebuild': (['selection','representation'],0),
}
SELECTION_FIELDS={'selection','selection1','selection2','atom1','atom2','mobile','target'}
SETTINGS={'transparency','cartoon_transparency','stick_transparency','sphere_transparency',
          'stick_radius','sphere_scale','cartoon_fancy_helices','cartoon_flat_sheets',
          'cartoon_smooth_loops','cartoon_trace_atoms','cartoon_color','surface_color',
          'dash_color','dash_width','dash_gap','label_color','label_size',
          'orthoscopic','depth_cue','ray_opaque_background','two_sided_lighting',
          'ambient','direct','specular','shininess','surface_quality','cartoon_sampling'}
LABELS={'"%s%s/%s" % (resn,resi,chain)','resn','resi','name','chain','""'}
REPRESENTATIONS={'everything','cartoon','sticks','lines','surface','spheres','ribbon','labels','dots','nonbonded','mesh'}
_UNDO={}


def selection_text(value):
    if not isinstance(value,str) or not 1<=len(value)<=1200 or not re.fullmatch(r'[A-Za-z0-9_ .+*()<>!=\-]+',value):
        raise ValueError('Unsupported selection syntax')
    level=0
    for char in value:
        level += (char=='(')-(char==')')
        if level<0:raise ValueError('Unbalanced selection')
    if level:raise ValueError('Unbalanced selection')
    return value


def compile_program(code):
    if not isinstance(code,str) or not 1<=len(code)<=16000:raise ValueError('Program must contain 1..16000 characters')
    tree=ast.parse(code)
    if not 1<=len(tree.body)<=80:raise ValueError('Use 1..80 direct cmd calls')
    calls=[]
    for statement in tree.body:
        call=statement.value if isinstance(statement,ast.Expr) else None
        if not isinstance(call,ast.Call) or not isinstance(call.func,ast.Attribute) or not isinstance(call.func.value,ast.Name) or call.func.value.id!='cmd':
            raise ValueError('Only direct cmd.method calls with literal arguments are supported')
        method=call.func.attr
        if method not in METHODS:raise ValueError('Unsupported PyMOL API method: '+method)
        fields,required=METHODS[method]
        if len(call.args)>len(fields):raise ValueError('Too many arguments')
        try:
            values={fields[i]:ast.literal_eval(v) for i,v in enumerate(call.args)}
            for keyword in call.keywords:
                if keyword.arg not in fields or keyword.arg in values:raise ValueError('Invalid or repeated keyword')
                values[keyword.arg]=ast.literal_eval(keyword.value)
        except (ValueError,TypeError):raise ValueError('Arguments must be declared literal values') from None
        if any(field not in values for field in fields[:required]):raise ValueError('Missing argument for '+method)
        for key,value in values.items():
            if key in SELECTION_FIELDS:selection_text(value)
            elif key=='rgb':
                if not isinstance(value,(list,tuple)) or len(value)!=3 or any(type(v) not in (int,float) or not math.isfinite(v) or not 0<=v<=1 for v in value):raise ValueError('RGB requires three numbers in [0,1]')
            elif key=='expression':
                if value not in (LABELS if method=='label' else {'b','q','count'}):raise ValueError('Unsupported expression')
            elif key=='representation':
                if value not in REPRESENTATIONS:raise ValueError('Unsupported representation')
            elif key=='name':
                if method=='set':
                    if value not in SETTINGS:raise ValueError('Unsupported display setting')
                elif not isinstance(value,str) or not re.fullmatch(r'ai_[A-Za-z0-9_]{1,50}',value):raise ValueError('Generated names must start with ai_')
            elif key in {'color','palette'}:
                if not isinstance(value,str) or not re.fullmatch(r'[A-Za-z][A-Za-z0-9_]{0,60}|0x[0-9a-fA-F]{6}',value):raise ValueError('Invalid color/palette')
            elif key=='axis':
                if value not in {'x','y','z'}:raise ValueError('Invalid axis')
            elif key=='value' and isinstance(value,str):
                if not re.fullmatch(r'[A-Za-z0-9_ .\-]{1,60}',value):raise ValueError('Invalid setting value')
            elif type(value) not in (int,float) or not math.isfinite(value) or abs(value)>10000:
                raise ValueError('Invalid numeric argument')
        if method=='distance' and (values.get('mode',0) not in (0,2) or not 0<values.get('cutoff',3.6)<=8):raise ValueError('Distance mode must be 0 or 2, cutoff <=8 A')
        if method in {'align','super'} and (type(values.get('cycles',5)) is not int or not 0<=values.get('cycles',5)<=5):raise ValueError('Alignment cycles must be 0..5')
        if method=='set' and not values['name'].endswith('_color'):
            v=values['value'];name=values['name']
            low,high=(-4,2) if name=='surface_quality' else (1,30) if name=='cartoon_sampling' else (0,1) if 'transparency' in name else (0,100)
            if type(v) not in (int,float) or not low<=v<=high:raise ValueError('Display setting outside supported range')
        calls.append((method,values))
    return calls


def scene(cmd,catalog):
    from .pymol_bridge import selection
    loaded=set(cmd.get_names('objects'));objects=[]
    for row in catalog:
        obj=row['id']
        if obj not in loaded:continue
        atoms=cmd.get_model(obj,state=1).atom
        protein=cmd.get_model(f'{obj} and polymer.protein',state=1).atom
        residues=sorted({(a.chain,a.resi,a.resn) for a in protein})
        lig=selection(obj,row,dict(target='ligand'))
        nearby=cmd.get_model(f'({obj} and polymer.protein) within 5 of {lig}',state=1).atom
        objects.append(dict(id=obj,label=row['label'],kind=row.get('kind','unknown'),atoms=len(atoms),
            ligand_selection=lig,ligand_atoms=cmd.count_atoms(lig),
            protein_chains=sorted({a.chain for a in protein}),
            chains_within_5A_of_ligand=sorted({a.chain for a in nearby}),
            protein_residues=[dict(chain=c,resi=i,resn=n) for c,i,n in residues[:200]],
            residue_count=len(residues),residues_truncated=len(residues)>200,
            ca_only=bool(protein) and all(a.name=='CA' for a in protein)))
    return dict(objects=objects,selections=[n for n in cmd.get_names('selections') if n.startswith('ai_')],
                scope='Current displayed state 1; no coordinates sent to the model. Residue lists may be truncated.')


def run_program(cmd,message,catalog,root):
    from .pymol_bridge import write_json
    view=message['view'];root=Path(root);key=str(root.resolve())
    if view['operation']=='undo':
        if key not in _UNDO:raise ValueError('No program checkpoint in this PyMOL process')
        cmd.set_session(_UNDO.pop(key))
        return dict(status='complete',operation='undo',scene=scene(cmd,catalog),artifacts={})
    if view['operation']=='scene':return dict(status='complete',operation='scene',scene=scene(cmd,catalog),artifacts={})
    calls=compile_program(view['code'])
    loaded=set(cmd.get_names('objects'));ids=[r['id'] for r in catalog if r['id'] in loaded]
    if not ids:raise ValueError('Open the task structures first')
    scope='('+' or '.join(ids)+')'
    before=cmd.get_session()
    prefix=root/message['id'];checkpoint=Path(str(prefix)+'-before.pse')
    cmd.save(str(checkpoint))
    codepath=Path(str(prefix)+'.py');codepath.write_text(view['code'],encoding='utf-8')
    output=dict(status='complete',operation='program',calls=[],artifacts=dict(code=str(codepath),checkpoint=str(checkpoint)))
    try:
        for method,raw in calls:
            values=dict(raw)
            if 'selection' in METHODS[method][0]:values.setdefault('selection','all')
            for field in SELECTION_FIELDS.intersection(values):
                values[field]=scope+' and ('+values[field]+')'
            if method=='set' and 'selection' not in raw:
                # Known global display settings only; all other settings stay object-scoped.
                if values['name'] in {'orthoscopic','depth_cue','ray_opaque_background','two_sided_lighting','ambient','direct','specular','shininess'}:values.pop('selection',None)
            if method in {'align','super'}:values.update(mobile_state=1,target_state=1)
            if method in {'distance','get_distance'}:values['state']=1
            if method=='distance':values.setdefault('cutoff',3.6)
            selected={f:cmd.count_atoms(values[f]) for f in SELECTION_FIELDS.intersection(values)}
            if method!='count_atoms' and selected and not all(selected.values()):raise ValueError('Empty selection in '+method+': '+json.dumps(raw))
            result=getattr(cmd,method)(**values)
            output['calls'].append(dict(method=method,arguments=raw,selected_atoms=selected,result=result))
        png=Path(str(prefix)+'.png');pse=Path(str(prefix)+'.pse')
        cmd.png(str(png),width=1400,height=1000,dpi=150,ray=0,quiet=1)
        cmd.save(str(pse))
        output['artifacts'].update(snapshot=str(png),save_session=str(pse))
        output['scene']=scene(cmd,catalog)
        _UNDO[key]=before
    except Exception as exc:
        try:
            cmd.set_session(before);rollback='complete'
        except Exception as rollback_exc:rollback='failed: '+str(rollback_exc)
        output.update(status='failed',error=str(exc),rollback=rollback)
    write_json(Path(str(prefix)+'-execution.json'),output)
    return output
