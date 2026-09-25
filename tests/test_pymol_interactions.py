import json
from types import SimpleNamespace

import numpy as np
import pytest

from aidd_agent.pymol_interactions import detect, display, DEFAULT_TYPES, aromatic_indices
from aidd_agent.pymol_program import compile_program
from aidd_agent.pymol_bridge import validate_view, execute_command


def atom(xyz,elem='C',charge=0,resi='1',acceptor=False):
    return dict(coord=xyz,elem=elem,charge=charge,chain='A',resi=resi,resn='LIG',acceptor=acceptor)


def ring(z=0,x=0,tilt=False):
    return [([x+np.cos(t),0,z+np.sin(t)] if tilt else [x+np.cos(t),np.sin(t),z])
            for t in np.arange(6)*np.pi/3]


def test_ring_geometry_rejects_side_by_side_and_keeps_stack():
    def run(x):
        atoms={i:atom(p,resi=str(i//6)) for i,p in enumerate(ring()+ring(z=3.5,x=x))}
        edges=[(i,(i//6)*6+(i+1)%6) for i in range(12)]
        return detect(atoms,edges,set(range(6)),set(range(6,12)),set(range(12)),[])[0]
    assert any(r['kind']=='pi_stacking' for r in run(0))
    assert not any(r['kind']=='pi_stacking' for r in run(3))


def test_charge_polar_overlap_and_missing_charge():
    atoms={1:atom([0,0,0],'N',1),2:atom([3,0,0],'O',-1,'2')}
    raw,reps=detect(atoms,[],{1},{2},set(),[(1,2),(1,2)])
    assert [r['kind'] for r in raw]==['salt_bridge']
    atoms[1]['charge']=0
    raw,reps=detect(atoms,[],{1},{2},set(),[(1,2),(1,2)])
    assert [r['kind'] for r in raw]==['polar_contact']


def test_hydrophobic_dedup_and_polar_carbon_exclusion():
    atoms={1:atom([0,0,0]),2:atom([0,1,0]),3:atom([3.5,0,0],resi='2'),4:atom([3.5,1,0],resi='2')}
    raw,reps=detect(atoms,[(1,2),(3,4)],{1,2},{3,4},set(),[])
    assert len(raw)==4 and len(reps)==1
    atoms[2]['elem']='O'
    assert not detect(atoms,[(1,2),(3,4)],{1,2},{3,4},set(),[])[0]
    assert 'hydrophobic' in DEFAULT_TYPES


def test_cation_pi_requires_face_geometry():
    atoms={i:atom(p) for i,p in enumerate(ring())}
    atoms[6]=atom([0,0,4],'N',1,'2')
    edges=[(i,(i+1)%6) for i in range(6)]
    assert any(r['kind']=='cation_pi' for r in detect(atoms,edges,set(range(6)),{6},set(range(6)),[])[0])
    atoms[6]['coord']=[4,0,0]
    assert not detect(atoms,edges,set(range(6)),{6},set(range(6)),[])[0]


class Scene:
    def __init__(self):self.calls=[]
    def get_names(self,*args,**kwargs):return ['v001','v002','ai_old'] if not kwargs else ['v001']
    def get_type(self,name):return 'object:measurement'
    def get_model(self,sel,state=1):
        def a(i,element,xyz):return SimpleNamespace(index=i,symbol=element,coord=xyz,formal_charge=0,
            name=element+str(i),segi='',chain='A',resi=str(i),resn='LIG' if i==1 else 'ASN',alt='',q=1)
        atoms=[a(1,'N',[0,0,0]),a(2,'O',[3,0,0])]
        if 'aromatic' in sel:raise RuntimeError(' Error: invalid selection name aromatic')
        elif 'polymer' in sel:atoms=atoms[1:]
        elif 'organic' in sel:atoms=atoms[:1]
        return SimpleNamespace(atom=atoms,bond=[])
    def find_pairs(self,left,right,**kwargs):
        return [(('v001',1),('v001',2))] if 'organic' in left else []
    def get_session(self):return {}
    def save(self,path):
        from pathlib import Path
        Path(path).write_text('session')
    def png(self,path,**kwargs):self.save(path)
    def __getattr__(self,name):return lambda *args,**kwargs:self.calls.append((name,args,kwargs))


def test_display_only_enabled_complex_replaces_lines_and_exports(tmp_path):
    cmd=Scene();catalog=[dict(id='v001',label='one'),dict(id='v002',label='two')]
    result=execute_command(cmd,dict(id='1-abcdefabcdef',view=dict(operation='interaction_overview')),catalog,tmp_path)
    assert result['objects']==['v001']
    assert result['displayed_counts']=={'v001':{'polar_contact':1}}
    assert ('disable',('ai_old',),{}) in cmd.calls
    assert any(n=='set' and args==('dash_color','yellow','v001_ix_polar_contact') for n,args,kwargs in cmd.calls)
    assert len([1 for n,*_ in cmd.calls if n=='distance'])==1
    with open(result['artifacts']['interactions_json']) as f:assert json.load(f)['results']['v001']['accepted']


def test_type_contract():
    validate_view(dict(operation='typed_interactions',types=['pi_stacking']))
    compile_program('cmd.typed_interactions(types="pi_stacking salt_bridge", objects="v001")')
    with pytest.raises(ValueError):compile_program('cmd.typed_interactions(types="made_up")')
    with pytest.raises(ValueError):validate_view(dict(operation='typed_interactions',types=['all']))


def test_halogen_requires_both_angles():
    atoms={1:atom([0,0,0],'Cl'),2:atom([-1,0,0]),
           3:atom([3,0,0],'O',resi='2',acceptor=True),4:atom([3.5,.866,0],resi='2')}
    raw,_=detect(atoms,[(1,2),(3,4)],{1,2},{3,4},set(),[])
    assert any(r['kind']=='halogen_bond' for r in raw)
    atoms[2]['coord']=[0,1,0]
    assert not any(r['kind']=='halogen_bond' for r in detect(atoms,[(1,2),(3,4)],{1,2},{3,4},set(),[])[0])


def test_zero_typed_hits_does_not_fall_back_to_distance(tmp_path):
    cmd=Scene()
    result=display(cmd,dict(id='1-abcdefabcdef',view=dict(operation='typed_interactions',types=['salt_bridge'])),
                   [dict(id='v001',label='one')],tmp_path)
    assert result['displayed_counts']=={'v001':{}}
    assert not any(n=='distance' for n,*_ in cmd.calls)
    assert not any(n=='hide' and args[-1]=='v001_ix_points' for n,args,kwargs in cmd.calls)


def test_aromatic_bond_indices_are_object_indices_not_positions():
    model=SimpleNamespace(atom=[SimpleNamespace(index=i) for i in (5,12,99)],
                          bond=[SimpleNamespace(index=[0,1],order=4),SimpleNamespace(index=[1,2],order=1)])
    assert aromatic_indices(model)=={5,12}


def test_empty_pymol_error_keeps_stage_and_selection(tmp_path):
    class Broken(Scene):
        def get_model(self,sel,state=1):raise RuntimeError(' Error: ')
    with pytest.raises(RuntimeError,match='detect interactions for v001: cmd.get_model') as error:
        display(Broken(),dict(id='1-abcdefabcdef',view=dict(operation='typed_interactions')),
                [dict(id='v001',label='one')],tmp_path)
    assert 'organic' in str(error.value)


def test_live_pymol_aromatic_model_when_available():
    pymol2=pytest.importorskip('pymol2',reason='Real PyMOL runtime is only available on the workstation')
    with pymol2.PyMOL() as pm:
        pm.cmd.fragment('benzene','v001')
        model=pm.cmd.get_model('v001',state=1)
        assert len(aromatic_indices(model))==6


def test_deposited_aromatic_flags_recover_kekule_rings_without_changing_coordinates():
    from aidd_agent.pymol_chemistry import apply_components
    atoms={i:dict(atom(p,resi=str(i//6)),name='C'+str(i%6)) for i,p in enumerate(ring()+ring(z=3.5))}
    edges=[(i,(i//6)*6+(i+1)%6) for i in range(12)]
    component=dict(atoms={f'C{i}':dict(element='C',aromatic=True) for i in range(6)},
                   bonds=[[f'C{i}',f'C{(i+1)%6}'] for i in range(6)])
    metadata=dict(status='available',components=dict(LIG=component))
    assert not any(r['kind']=='pi_stacking' for r in detect(atoms,edges,set(range(6)),set(range(6,12)),set(),[])[0])
    aromatic,info=apply_components(atoms,edges,metadata,set())
    assert len(aromatic)==12 and info['component_typed_atoms']==12
    assert any(r['kind']=='pi_stacking' for r in detect(atoms,edges,set(range(6)),set(range(6,12)),aromatic,[])[0])
    # Removed bonds and element mismatches must not reconstruct the original ring.
    aromatic,_=apply_components(atoms,edges[1:],metadata,set())
    assert 0 not in aromatic and 1 not in aromatic
    atoms[0]['elem']='N'
    aromatic,_=apply_components(atoms,edges,metadata,set())
    assert 0 not in aromatic


def test_component_metadata_reads_deposited_flags_and_checks_hash(tmp_path):
    pytest.importorskip('gemmi')
    from aidd_agent.pymol_chemistry import component_metadata
    path=tmp_path/'complex.cif'
    path.write_text('''data_example
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.pdbx_aromatic_flag
LIG C1 C Y
LIG C2 C Y
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
LIG C1 C2
''')
    result=component_metadata(path)
    assert result['components']['LIG']['atoms']['C1']['aromatic']
    assert result['components']['LIG']['bonds']==[['C1','C2']]
    with pytest.raises(ValueError,match='provenance'):component_metadata(path,'incorrect')


def test_display_uses_source_typing_to_draw_cyan_pi_object(tmp_path):
    class RingScene(Scene):
        def get_model(self,sel,state=1):
            atoms=[SimpleNamespace(index=i+1,name=f'C{i%6}',symbol='C',coord=p,
                formal_charge=0,segi='',chain='A',resi=str(i//6),resn='LIG',alt='',q=1)
                for i,p in enumerate(ring()+ring(z=3.5))]
            bonds=[SimpleNamespace(index=[i,(i//6)*6+(i+1)%6],order=1+i%2) for i in range(12)]
            if 'acceptors' in sel:return SimpleNamespace(atom=[],bond=[])
            if 'polymer' in sel:return SimpleNamespace(atom=atoms[6:],bond=[])
            if 'organic' in sel:return SimpleNamespace(atom=atoms[:6],bond=[])
            return SimpleNamespace(atom=atoms,bond=bonds)
        def find_pairs(self,*args,**kwargs):return []
    metadata=dict(status='available',components=dict(LIG=dict(
        atoms={f'C{i}':dict(element='C',aromatic=True) for i in range(6)},
        bonds=[[f'C{i}',f'C{(i+1)%6}'] for i in range(6)])))
    cmd=RingScene()
    result=display(cmd,dict(id='replay',view=dict(operation='typed_interactions')),
                   [dict(id='v001',chemical_components=metadata)],tmp_path)
    assert result['displayed_counts']['v001']==dict(pi_stacking=1)
    assert result['chemistry']['v001']['ligand_rings']==1
    assert ('set',('dash_color','cyan','v001_ix_pi_stacking'),{}) in cmd.calls
    assert result['type_counts']['v001']['hydrophobic']['requested']
