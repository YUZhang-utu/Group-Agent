from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from aidd_agent.fast_3d_search import exhaustive_candidates
from aidd_agent.classified_features import extra_hypotheses, write_tables


def test_exhaustive_visits_every_shard_and_matches_dense_topk(tmp_path):
    rng=np.random.default_rng(45)
    vectors=rng.normal(size=(19,60)).astype('float32');vectors[-1]=0;vectors[-2]=0
    ids=np.array([('m'+str(i)).encode() for i in range(19)],dtype='S16')
    rows=[]
    for start,stop in ((0,7),(7,12),(12,19)):
        path=tmp_path/str(start);path.mkdir()
        vectors[start:stop].tofile(path/'usrcat.f32.bin');ids[start:stop].tofile(path/'molecule_ids.bin')
        rows.append(dict(name=str(start),path=str(path),global_id_start=start,conformers=stop-start))
    corpus=SimpleNamespace(rows=rows,total=19,fetch=lambda selected:(ids[selected],None))
    result=exhaustive_candidates(corpus,np.zeros(60,dtype='float32'),np.zeros(60,dtype='float32'),
        np.ones(60,dtype='float32'),tmp_path/'candidates.npz',budget=6,chunk_size=3)
    d=np.einsum('ij,ij->i',vectors,vectors);expected=np.lexsort((np.arange(19),d))[:6]
    with np.load(tmp_path/'candidates.npz') as data:
        assert np.array_equal(data['global_ids'],expected)
        assert np.array_equal(data['squared_l2'],d[expected])
    assert result['compared_conformers']==result['total_conformers']==19
    assert result['coverage_fraction']==1.
    corpus.rows[1]['global_id_start']=8
    with pytest.raises(ValueError,match='gaps or overlap'):
        exhaustive_candidates(corpus,np.zeros(60),np.zeros(60),np.ones(60),tmp_path/'bad.npz')


def test_crystal_hypotheses_have_distinct_categories_and_no_candidate_contact_claim():
    query=dict(feature_points=np.array([[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]]),feature_types=np.array([5,3,2]),
               feature_directions=np.zeros((3,3)),feature_direction_kinds=np.zeros(3))
    def atom(res,name,xyz,num):
        return dict(residue=res,atom_name=name,xyz=np.array(xyz),chain='A',residue_number=num,donor=True,acceptor=True)
    atoms=[atom('ALA','CB',[3,0,0],'1'),atom('ASP','OD1',[0,3,0],'2'),atom('ASP','OD2',[0,3,1],'2')]
    water=atom('HOH','O',[0,2,0],'3');metal=atom('ZN','ZN',[0,0,2],'4')
    rows=extra_hypotheses('TEST:LIG:A:1',query,atoms,[water],[metal])
    assert {'hydrophobic','salt_bridge','water_bridge','metal_coordination'} <= {r['feature_class'] for r in rows}
    assert all('not a candidate-protein contact' in r['interpretation'] for r in rows)
    assert all(r['evidence_class']=='crystal_derived_geometric_hypothesis' for r in rows)


def test_classification_html_escapes_evidence_and_shows_unassessed_categories(tmp_path):
    anchor=dict(anchor_id='<script>x</script>',feature_index=0,interaction_class='hydrophobic',ligand_atom_indices=[2],
        diagnostic_counts=[dict(minimum_score=.5,molecules=7)],evidence_class='hypothesis',interpretation='not validated',
        evidence=dict(protein_partner=dict(chain='A',residue='ALA',residue_number='1',atom_name='CB'),distance_angstrom=3.))
    write_tables(dict(queries=[dict(query_id='test',anchors=[anchor])]),tmp_path)
    text=(tmp_path/'interactions.html').read_text()
    assert '<script>' not in text and '&lt;script&gt;' in text
    assert 'Not indexed' in text and 'No hypothesis detected' in text
    assert 'molecules_at_score_0_5' in text


def test_classification_preview_export_lineage(tmp_path,monkeypatch):
    import json
    from test_screening_selection import fixture_search,write
    from aidd_agent import classified_features as cf, screening_selection as ss
    from aidd_agent.library_acceptance import sha
    search=fixture_search(tmp_path/'search')
    qpath=tmp_path/'crystal/gaussian-query.manifest.json'
    q=json.loads(qpath.read_text());q['source']['feature_atom_indices']=[[0],[1]];write(qpath,q)
    protocol=search.parent/'protocol.json';p=json.loads(protocol.read_text());p['sources'][str(qpath)]=sha(qpath);write(protocol,p)
    ss.review(search,tmp_path/'review')
    monkeypatch.setattr(cf,'crystal_environment',lambda *a:([],[],[]))
    def score(artifact,chemical,query,ids,mols,confs,objectives,transforms,**kw):
        n=len(ids);k=len(query['anchor_feature_indices']);o=objectives[0]
        return {o:np.zeros(n)},{o:np.full((n,k),-1)},{o:np.zeros((n,k))}
    monkeypatch.setattr(cf,'score_batched',score)
    classified=cf.classify(tmp_path/'review/report.json',tmp_path/'classified')
    assert (tmp_path/'classified/interactions.html').is_file()
    policy=dict(required_anchors=[classified['queries'][0]['anchors'][0]['anchor_id']],match_mode='all',minimum_score=.5)
    ss.preview(tmp_path/'classified/report.json',tmp_path/'preview',policy)
    result=ss.export(tmp_path/'preview/report.json',tmp_path/'export')
    assert result['counts']['selected_molecules']==0


def test_pocket_geometry_is_derived_and_translation_covariant():
    from aidd_agent.docking_preparation import pocket_geometry
    points=np.array([[-2.,0,0],[2,4,6]])
    base=pocket_geometry(points);moved=pocket_geometry(points+10)
    assert base['center_angstrom']==[0.,2.,3.]
    assert np.allclose(np.array(base['center_angstrom'])+10,moved['center_angstrom'])
    assert base['radius_angstrom']==moved['radius_angstrom']
    with pytest.raises(ValueError):pocket_geometry([[float('nan'),0,0]])


@pytest.mark.parametrize('passed',[True,False])
def test_docking_reference_gate_blocks_candidate_submission(tmp_path,monkeypatch,passed):
    from test_screening_selection import write
    from aidd_agent import docking_preparation as dp
    from aidd_agent.expanded_wee1 import fingerprint
    executable=tmp_path/'engine';executable.write_text('fixture');executable.chmod(0o700)
    inputs=tmp_path/'inputs';inputs.mkdir()
    native=inputs/'reference-native.sdf';native.write_text('fixture')
    source=tmp_path/'source';source.write_text('source')
    commands=[dict(name='reference_conversion',argv=[str(executable),'reference'],expected='reference-poses.sdf'),
              dict(name='candidate_docking',argv=[str(executable),'candidate'],expected='candidates_pv.maegz')]
    plan=write(inputs/'report.json',dict(kind='docking_preparation',sources=fingerprint([source]),outputs=fingerprint([native]),commands=commands,redocking_rmsd_threshold_angstrom=2.))
    calls=[]
    def run(argv,cwd,**kwargs):
        calls.append(argv[-1]);(cwd/('reference-poses.sdf' if argv[-1]=='reference' else 'candidates_pv.maegz')).write_text('fixture')
        return SimpleNamespace(returncode=0)
    monkeypatch.setattr(dp.subprocess,'run',run)
    monkeypatch.setattr(dp,'reference_rmsd',lambda *a:dict(best_in_place_heavy_atom_rmsd=1. if passed else 4.))
    if passed:
        result=dp.run(plan,tmp_path/'execution');assert result['status']=='complete';assert calls==['reference','candidate']
    else:
        with pytest.raises(ValueError,match='gate failed'):dp.run(plan,tmp_path/'execution')
        assert calls==['reference']
