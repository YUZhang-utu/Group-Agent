"""Validate independent ligand tasks without executing AF3 inference."""
import json
from pathlib import Path
import sys
import subprocess
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools/af3_smiles_batch'))
import prepare_batch as prep
import run_batch as runner
import msa_reuse


def inputs(tmp_path):
    base=tmp_path/'base.json';table=tmp_path/'ligands.tsv'
    prep.save(base,dict(name='base',dialect='alphafold3',version=1,modelSeeds=[1],
        sequences=[{'protein':{'id':['A','B'],'sequence':'ACDEFG','unpairedMsa':'>x\nACDEFG\n'}}]))
    table.write_text('id\tsmiles\nlig1\tC[C@H](O)C(=O)O\nlig2\tC[NH3+]\n',encoding='utf-8')
    return base,table


def test_separate_ligands_and_dimer_ids(tmp_path):
    base,table=inputs(tmp_path);out=tmp_path/'batch'
    manifest=prep.prepare(base,table,out,[1,2,3],2)
    assert manifest['ligand_chain_ids']==['C','D']
    for job in manifest['jobs']:
        result=json.loads((out/job['input']).read_text())
        assert len(result['sequences'])==2
        assert result['sequences'][-1]['ligand']['smiles']==job['smiles']
        assert result['modelSeeds']==[1,2,3]
        assert result['sequences'][0]['protein']['unpairedMsa']=='>x\nACDEFG\n'


@pytest.mark.parametrize('smiles',['invalid','CC.O','CC name',''])
def test_invalid_or_disconnected_smiles(smiles):
    with pytest.raises(ValueError):prep.validate_smiles(smiles,rdkit_check=True)


def test_preparation_without_site_packages(tmp_path):
    base,table=inputs(tmp_path);out=tmp_path/'stdlib_only'
    result=subprocess.run([sys.executable,'-S',str(Path(prep.__file__)),
        '--base',str(base),'--ligands',str(table),'--output',str(out)],capture_output=True,text=True)
    assert result.returncode==0,result.stdout+result.stderr
    manifest=json.loads((out/'manifest.json').read_text())
    assert manifest['smiles_validation']=='text_only_chemical_validation_deferred_to_AF3'
    assert all(job['heavy_atoms'] is None for job in manifest['jobs'])


def test_base_ligand_requires_explicit_retention(tmp_path):
    base,table=inputs(tmp_path);value=json.loads(base.read_text())
    value['sequences'].append({'ligand':{'id':'C','ccdCodes':['NAD']}});prep.save(base,value)
    with pytest.raises(ValueError,match='Base contains a ligand'):prep.prepare(base,table,tmp_path/'bad')
    result=prep.prepare(base,table,tmp_path/'good',keep_ligands=True)
    assert result['ligand_chain_ids']==['D']


def test_unsafe_ids_and_external_paths(tmp_path):
    base,table=inputs(tmp_path);table.write_text('id\tsmiles\n../bad\tCC\n')
    with pytest.raises(ValueError,match='safe'):prep.prepare(base,table,tmp_path/'bad')
    value=json.loads(base.read_text());value['sequences'][0]['protein']['unpairedMsaPath']='msa.a3m'
    with pytest.raises(ValueError,match='inlined'):prep.validate_base(value,False)


def test_batch_verifies_model_and_hash(tmp_path,monkeypatch):
    base,table=inputs(tmp_path);out=tmp_path/'batch';prep.prepare(base,table,out)
    count=[0]
    def popen(command,**kwargs):
        count[0]+=1;folder=tmp_path/('prediction'+str(count[0]));folder.mkdir()
        (folder/'job_model.cif').write_text('data_test\n')
        prep.save(folder/'job_summary_confidences.json',{'iptm':.5})
        class Process:
            stdout=iter(['Output: '+str(folder)+'\n'])
            def wait(self):return 0
        return Process()
    monkeypatch.setattr(runner.subprocess,'Popen',popen)
    report=runner.run(out,optimize=False)
    assert report['status']=='complete' and len(report['jobs'])==2
    assert all(j['confidence']['iptm']==.5 for j in report['jobs'])
    (out/'inputs/lig1.json').write_text('{}')
    with pytest.raises(ValueError,match='checksum'):runner.run(out)


def test_no_final_model_is_failure(tmp_path,monkeypatch):
    base,table=inputs(tmp_path);out=tmp_path/'batch';prep.prepare(base,table,out)
    class Process:
        stdout=iter([])
        def wait(self):return 0
    monkeypatch.setattr(runner.subprocess,'Popen',lambda *a,**k:Process())
    with pytest.raises(RuntimeError,match='Job failed'):runner.run(out,optimize=False)


def test_one_data_pipeline_many_inference_jobs(tmp_path,monkeypatch):
    base,table=inputs(tmp_path);out=tmp_path/'batch';prep.prepare(base,table,out,[7,9])
    calls=[]
    def invoke(command,log):
        calls.append(command[-1]);data=json.loads(Path(command[2]).read_text())
        folder=tmp_path/('output'+str(len(calls)));folder.mkdir()
        if command[-1]=='--data-only':
            for e in data['sequences']:
                if 'protein' in e:e['protein'].update(unpairedMsa='>q\nACDEFG\n',pairedMsa='',templates=[])
            prep.save(folder/'lig1_data.json',data)
        else:
            assert command[-1]=='--inference-only'
            assert data['modelSeeds']==[7,9]
            assert data['sequences'][0]['protein']['pairedMsa']==''
            original=json.loads((out/'inputs'/(data['name']+'.json')).read_text())
            assert data['sequences'][-1]==original['sequences'][-1]
            (folder/'job_model.cif').write_text('data_test\n')
        return 0,folder
    monkeypatch.setattr(runner,'invoke',invoke)
    report=runner.run(out)
    assert calls==['--data-only','--inference-only','--inference-only']
    assert report['status']=='complete' and not report['shared_data']['reused_inline']


def test_inline_features_skip_shared_search(tmp_path,monkeypatch):
    base,table=inputs(tmp_path);value=json.loads(base.read_text())
    value['sequences'][0]['protein'].update(pairedMsa='',templates=[]);prep.save(base,value)
    out=tmp_path/'batch';prep.prepare(base,table,out);calls=[]
    def invoke(command,log):
        calls.append(command[-1]);folder=tmp_path/('output'+str(len(calls)));folder.mkdir()
        (folder/'job_model.cif').write_text('data_test\n');return 0,folder
    monkeypatch.setattr(runner,'invoke',invoke)
    report=runner.run(out)
    assert calls==['--inference-only','--inference-only'] and report['shared_data']['reused_inline']


def test_feature_transfer_rejects_mismatched_chain_and_sequence(tmp_path):
    base,_=inputs(tmp_path);job=json.loads(base.read_text());source=json.loads(base.read_text())
    p=source['sequences'][0]['protein'];p.update(pairedMsa='',templates=[])
    p['sequence']='ACDEFA'
    with pytest.raises(ValueError,match='sequence'):msa_reuse.transfer(job,source)
    p['sequence']='ACDEFG';p['id']='C'
    with pytest.raises(ValueError,match='chain IDs'):msa_reuse.transfer(job,source)


def test_data_pipeline_failure_is_recorded(tmp_path,monkeypatch):
    base,table=inputs(tmp_path);out=tmp_path/'batch';prep.prepare(base,table,out)
    monkeypatch.setattr(runner,'invoke',lambda *args:(1,None))
    with pytest.raises(RuntimeError,match='Shared data pipeline failed'):runner.run(out)
    report=json.loads(next(out.glob('run_*/report.json')).read_text())
    assert report['status']=='failed' and report['stage']=='shared_data'
