"""Validate independent ligand tasks without executing AF3 inference."""
import json
from pathlib import Path
import sys
import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools/af3_smiles_batch'))
import prepare_batch as prep
import run_batch as runner


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
    with pytest.raises(ValueError):prep.validate_smiles(smiles)


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
    report=runner.run(out)
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
    with pytest.raises(RuntimeError,match='Job failed'):runner.run(out)
