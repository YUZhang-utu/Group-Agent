"""Target identity, adopted receptor cleaning and native-frame SPORES preparation."""
from io import StringIO
from pathlib import Path
import subprocess

import numpy as np
import pytest

from aidd_agent.plants_receptors import extract, prepare, site_geometry, verify_heavy_coordinates
from aidd_agent.protein_data import resolve_pdb_target
from aidd_agent.joint_spatial_profiles import save
from aidd_agent.final_work_blocks import sha, read


def test_pdb_target_mapping_is_evidence_bound_and_ambiguous_entities_pause():
    def fetch(url):
        if '/entry/' in url:return dict(rcsb_entry_container_identifiers=dict(polymer_entity_ids=['1','2']))
        if '/polymer_entity/' in url:
            return dict(rcsb_polymer_entity_container_identifiers=dict(reference_sequence_identifiers=[dict(database_name='UniProt',database_accession='Q00987' if url.endswith('/1') else 'P04637')]))
        return dict(primaryAccession='Q00987',sequence=dict(value='AAA',length=3),organism=dict(taxonId=9606))
    result=resolve_pdb_target('6Y4Q',fetch=fetch)
    assert result['status']=='needs_target_identity' and len(result['candidates'])==2
    result=resolve_pdb_target('6Y4Q',entity_id='1',fetch=fetch)
    assert result['protein']['accession']=='Q00987'


def receptor_fixture(root):
    from Bio.PDB import PDBParser, MMCIFIO
    lines=[]
    atoms=[('ATOM','CA','ALA','A',1,(0,0,0),'C'),('ATOM','CA','ALA','A',2,(1,0,0),'C'),
           ('HETATM','C1','LIG','B',201,(4,0,0),'C'),('HETATM','C2','LIG','B',201,(4,2,0),'C'),
           ('HETATM','O','HOH','A',300,(3,1,0),'O'),('ATOM','CA','ALA','C',1,(30,0,0),'C')]
    for i,(kind,name,res,chain,num,xyz,element) in enumerate(atoms,1):
        x,y,z=xyz;lines.append(f'{kind:<6}{i:5d} {name:^4s} {res:3s} {chain}{num:4d}    {x:8.3f}{y:8.3f}{z:8.3f}{1.:6.2f}{20.:6.2f}          {element:>2s}  ')
    structure=PDBParser(QUIET=True).get_structure('test',StringIO('\n'.join(lines)+'\nEND\n'))
    source=root/'1ABC.cif'; writer=MMCIFIO();writer.set_structure(structure);writer.save(str(source))
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict
    data=MMCIF2Dict(str(source));data['_atom_site.auth_comp_id']=data['_atom_site.label_comp_id'];data['_atom_site.auth_atom_id']=data['_atom_site.label_atom_id']
    writer.set_dict(data);writer.save(str(source))
    row=dict(id='1ABC:A',pdb_id='1ABC',target_chain='A',structure_path=str(source),queries=[dict(query_id='1ABC:LIG:B:201')])
    return source,row


def fake_spores(argv,*,cwd,stdout,stderr,shell,timeout):
    assert argv[-4:]==['--mode','complete','protein-input.pdb','protein.mol2'] and shell is False
    text=(Path(cwd)/'protein-input.pdb').read_text()
    assert 'LIG' not in text and 'HOH' not in text and '30.000' not in text
    (Path(cwd)/'protein.mol2').write_text('@<TRIPOS>MOLECULE\nprotein\n2 0 0 0 0\nPROTEIN\nUSER_CHARGES\n@<TRIPOS>ATOM\n1 CA 0 0 0 C.3 1 ALA 0\n2 CA 1 0 0 C.3 2 ALA 0\n')
    return subprocess.CompletedProcess(argv,0)


def test_cleaning_spores_centroid_and_profile_are_automatic_after_adoption(tmp_path,monkeypatch):
    source,row=receptor_fixture(tmp_path)
    pockets=tmp_path/'pockets.json';save(pockets,dict(status='complete',clusters=[dict(id='p1',representative=row['id'])],structures=[row],sources={str(source):sha(source)}))
    adoption=tmp_path/'adoption.json';save(adoption,dict(status='complete',kind='pocket_adoption',pocket_report=str(pockets),selected_cluster_ids=['p1'],sources={str(pockets):sha(pockets),str(source):sha(source)}))
    exe=tmp_path/'SPORES_64bit';exe.write_text('fixture')
    profile=tmp_path/'tools.json';save(profile,dict(spores_executable=str(exe),spores_mode='complete',plants_executable=str(exe)))
    monkeypatch.setattr('aidd_agent.plants_receptors.subprocess.run',fake_spores)
    result=prepare(adoption,profile,tmp_path/'prepared')
    assert result['receptors']==1
    cfg=read(tmp_path/'prepared/plants-profile.json');r=cfg['receptors'][0]
    assert r['center']==[4.,1.,0.] and r['radius']==6
    assert prepare(adoption,profile,tmp_path/'prepared')==result
    with pytest.raises(ValueError,match='Explicitly adopted'):prepare(pockets,profile,tmp_path/'bad')
    with pytest.raises(ValueError,match='correspondence|moved'):verify_heavy_coordinates([('C',[10,0,0]),('C',[11,0,0])],r['mol2'])


def test_site_radius_contains_native_ligand():
    xyz=np.array([[0,0,0],[2,0,0],[7,3,0.]])
    g=site_geometry(xyz)
    assert np.allclose(g['center'],xyz.mean(0))
    assert np.linalg.norm(xyz-g['center'],axis=1).max()+5<=g['radius']
