"""Native controls and tiny-library integration for the shared KRAS runner."""
from pathlib import Path
import sys
from types import SimpleNamespace
import numpy as np
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import run_kras_weekend as run
from aidd_agent import kras_weekend as kw
from aidd_agent.conformer_artifacts import META_DTYPE, FEATURE_DTYPE
from aidd_agent.chemical_companion import CHEM_META_DTYPE, ATOM_DTYPE, FEATURE_DIRECTION_DTYPE, BOND_DTYPE, TORSION_DTYPE

PREPARED=Path(__file__).resolve().parents[1]/'configs/targets/kras_sip_e116'


def test_contact_diagnostic_isolates_reference_gate_and_restores_it(tmp_path, monkeypatch):
    import diagnose_kras_contacts as diagnostic
    batch=tiny_library(tmp_path/'diagnostic-library')
    run.initialize(batch, PREPARED)
    def reject_all(state, assignments, scores):
        return np.zeros(len(assignments), dtype=bool)
    monkeypatch.setattr(kw, 'contact_mask', reject_all)
    result=diagnostic.diagnose((0, [0, 2, 3]))
    assert kw.contact_mask is reject_all
    for target, row in result.items():
        assert row['counts']['processed'] == 3
        assert row['counts']['strict_full'] == 0
        assert row['counts']['bypass_full'] >= 1
        assert row['counts']['direct_poses'] > 0
        assert row['matches']


def test_contact_diagnostic_cli_pool_and_pilot_lineage(tmp_path):
    import subprocess
    batch=tiny_library(tmp_path/'cli-library')
    pilot=tmp_path/'pilot'
    run.execute(SimpleNamespace(batch=batch, prepared=PREPARED, output=pilot,
        workers=2, chunk_size=32, sample=1000, hours=.1))
    output=tmp_path/'diagnostic'
    subprocess.run([sys.executable, str(Path(run.__file__).with_name('diagnose_kras_contacts.py')),
        '--pilot', str(pilot), '--output', str(output), '--workers', '2'],
        check=True, capture_output=True, text=True)
    report=kw.read(output/'report.json')
    assert report['status']=='complete' and report['native_controls_passed']
    for row in report['targets'].values():
        assert row['counts']['processed']==4
        assert row['counts']['bypass_full']>=row['counts']['strict_full']>=1


def tiny_library(root):
    artifact=root/'artifacts'; chemical=root/'chemical'
    for directory in (artifact,chemical): (directory/'fixture').mkdir(parents=True)
    coords=[]; feats=[]; metas=[]; atoms=[]; directions=[]; members=[]; cmetas=[]
    cids=[]; mids=[]; co=fo=ao=mo=0
    for gid,target in enumerate(['KRAS_G12D','KRAS_G12D','KRAS_G12V','KRAS_G12V']):
        state=kw.load_state(PREPARED,target); q=state['query']; chemistry=state['record']['chemistry']
        shape=q['shape_points']; origin=shape.min(axis=0)
        # The last record has no donor/acceptor chemistry and must fail both branches.
        types=q['feature_types'].copy()
        if gid==3: types[:]=5
        points=np.rint((shape-origin)*100).astype('<i2'); coords.append(points)
        f=np.zeros(len(types),dtype=FEATURE_DTYPE)
        f['xyz']=np.rint((q['feature_points']-origin)*100).astype('<i2'); f['type']=types; feats.append(f)
        m=np.zeros(1,dtype=META_DTYPE); m['global_id']=gid; m['coord_offset']=co; m['feature_offset']=fo
        m['heavy_atoms']=len(shape); m['features']=len(types); m['origin']=origin; metas.append(m)
        a=np.zeros(len(shape),dtype=ATOM_DTYPE); a['atomic_number']=chemistry['atomic_numbers']; atoms.append(a)
        d=np.zeros(len(types),dtype=FEATURE_DIRECTION_DTYPE)
        d['direction']=q['feature_directions']; d['kind']=q['feature_direction_kinds']
        local_members=[]
        for i,indices in enumerate(chemistry['feature_members']):
            d[i]['member_offset']=len(local_members); d[i]['member_count']=len(indices)
            local_members.extend(indices)
        directions.append(d); members.append(np.array(local_members,dtype='<u2'))
        cm=np.zeros(1,dtype=CHEM_META_DTYPE); cm['global_id']=gid; cm['atom_offset']=ao
        cm['feature_offset']=fo; cm['member_offset']=mo; cm['atoms']=len(shape)
        cm['features']=len(types); cm['feature_members']=len(local_members); cmetas.append(cm)
        co+=len(shape); ao+=len(shape); fo+=len(types); mo+=len(local_members)
        cids.append(f'C{gid}'); mids.append('MD' if gid<2 else f'M{gid}')
    ar=artifact/'fixture'; cr=chemical/'fixture'
    for name,arrays in [('coords.bin',coords),('feats.bin',feats),('meta.bin',metas)]:
        np.concatenate(arrays).tofile(ar/name)
    for name,arrays in [('atoms.bin',atoms),('feature-directions.bin',directions),
                        ('feature-members.bin',members),('chem-meta.bin',cmetas)]:
        np.concatenate(arrays).tofile(cr/name)
    for name,dtype in [('bonds.bin',BOND_DTYPE),('torsions.bin',TORSION_DTYPE),('torsion-members.bin','<u2')]:
        np.empty(0,dtype=dtype).tofile(cr/name)
    for folder in (ar,cr):
        np.array(cids,dtype='S16').tofile(folder/'conformer_ids.bin')
        np.array(mids,dtype='S16').tofile(folder/'molecule_ids.bin')
    for folder,kind in ((ar,'aidd-conformer-artifact-shard'),(cr,'aidd-chemical-companion-shard')):
        kw.save(folder/'manifest.json',dict(format=kind,version=1,global_id_start=0,conformers=4,coordinate_scale=100))
    for folder,kind in ((artifact,'aidd-conformer-artifact-catalog'),(chemical,'aidd-chemical-companion-catalog')):
        doc=dict(format=kind,version=1,conformers=4,shards=[dict(name='fixture',path=str(folder/'fixture'),global_id_start=0,conformers=4)])
        if folder==chemical: doc['artifact_v1_catalog_sha256']=kw.sha(artifact/'catalog.json')
        kw.save(folder/'catalog.json',doc)
    return root


def test_both_native_references_recover_and_organic_fluorine_is_not_polar():
    for target in kw.REFERENCES:
        state=kw.load_state(PREPARED,target)
        result=kw.self_control(state)
        assert result['flags']['full']
        chemistry=SimpleNamespace(**state['record']['chemistry'])
        sets=kw.chemical_sets(chemistry,state['query']['feature_types'])
        assert all(chemistry.atomic_numbers[i] in (7,8) for i in sets['acceptor'])


def test_no_cross_pose_contact_pooling_and_codon12_is_separate():
    state=kw.load_state(PREPARED,'KRAS_G12D'); q=state['query']
    sets=kw.chemical_sets(SimpleNamespace(**state['record']['chemistry']),q['feature_types'])
    matrices=np.repeat(np.eye(4)[None],2,axis=0); matrices[1,:3,3]=100
    masks=kw.direct_contact_masks(q['shape_points'],sets,state,matrices)
    assert masks['all'].tolist()==[True,False]
    altered=dict(state,record=dict(state['record'],receptor_atoms=[a for a in state['record']['receptor_atoms'] if a['residue']!=12]))
    assert not kw.direct_contact_masks(q['shape_points'],sets,altered,matrices)['all'].any()


def test_shared_process_matches_independent_evaluation(tmp_path):
    batch=tiny_library(tmp_path/'library'); run.initialize(batch,PREPARED)
    paired=run.process((0,[0,2,3]))
    reader,features,states,chemistry=run.WORK
    for target in kw.REFERENCES:
        expected=[kw.evaluate(reader.get(gid),features.get(gid),states[target],chemistry.get(gid)) for gid in [0,2,3]]
        assert paired['targets'][target]['counts']['full']==sum(r['flags']['full'] for r in expected)
        assert not expected[-1]['flags']['chemical']
        assert paired['targets'][target]['hits']==[r for r in expected if r['flags']['full']]


def test_real_pool_resume_dedup_and_corruption(tmp_path):
    batch=tiny_library(tmp_path/'library')
    args=SimpleNamespace(batch=batch,prepared=PREPARED,output=tmp_path/'run',workers=2,chunk_size=2,sample=0,hours=.1)
    first=run.execute(args)
    assert first['full_library_complete']
    assert first['targets']['KRAS_G12D']['full']>=2
    assert first['targets']['KRAS_G12D']['distinct_passing_molecules']<first['targets']['KRAS_G12D']['full']
    second=run.execute(args)
    assert second['invocation_new_conformers']==0
    assert second['targets']==first['targets']
    chunk=next((args.output/'chunks').glob('*.json.gz'))
    chunk.write_bytes(chunk.read_bytes()+b'corruption')
    with pytest.raises(ValueError,match='Corrupt completed chunk'): run.execute(args)


def test_spread_schedule_and_full_schedule():
    full=list(run.tasks(7,3,0))
    assert [gid for _,ids in full for gid in ids]==list(range(7))
    sampled=[gid for _,ids in run.tasks(100,2,5) for gid in ids]
    assert len(set(sampled))==5 and sampled[0]==0 and sampled[-1]==99


def test_reference_seal_is_portable_between_windows_and_linux(tmp_path):
    report=kw.read(PREPARED/'report.json')
    for relative,digest in report['output_hashes'].items():
        assert '\\' not in relative
        assert kw.sha(PREPARED/relative)==digest
        if relative.endswith('.json'): assert b'\r\n' not in (PREPARED/relative).read_bytes()
    output=tmp_path/'test.json'; kw.save(output,{'status':'complete'})
    assert b'\r\n' not in output.read_bytes()


def test_seal_detects_reference_tampering(tmp_path):
    import shutil
    batch=tiny_library(tmp_path/'library'); prepared=tmp_path/'prepared'
    shutil.copytree(PREPARED,prepared)
    path=prepared/'KRAS_G12D/reference.json'
    path.write_text(path.read_text()+' ')
    with pytest.raises(ValueError,match='checksum'): run.seal_inputs(batch,prepared,2,0)


def test_budget_partial_then_resume_and_pilot_lineage(tmp_path):
    batch=tiny_library(tmp_path/'library')
    args=SimpleNamespace(batch=batch,prepared=PREPARED,output=tmp_path/'pilot',workers=2,chunk_size=2,sample=2,hours=1e-12)
    partial=run.execute(args)
    assert partial['status']=='partial' and not partial['full_library_complete']
    args.hours=.1
    assert run.execute(args)['status']=='complete'
    args.pilot_report=args.output/'report.json'; args.output=tmp_path/'full'; args.sample=0
    full=run.execute(args)
    assert full['full_library_complete']
    protocol=args.pilot_report.parent/'protocol.json'
    document=kw.read(protocol); document['policy']={}; kw.save(protocol,document)
    with pytest.raises(ValueError,match='Pilot/full lineage mismatch'): run.execute(args)
