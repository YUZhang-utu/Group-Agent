import json

import numpy as np
import pytest

from aidd_agent.structure_diversity import resolution_bin, diverse_indices, summarize_entry
from aidd_agent.polymer_ligand_diversity import token_similarity
from aidd_agent.prompt_plan import validate_plan


def test_resolution_bins_have_exact_edges_and_unknown():
    assert [resolution_bin(v) for v in [1.5,1.5001,2,2.01,2.5,2.501,3,3.1,None]]==[
        '<=1.5','(1.5,2.0]','(1.5,2.0]','(2.0,2.5]','(2.0,2.5]','(2.5,3.0]','(2.5,3.0]','>3.0','unknown']


def test_diversity_cover_and_cap_are_distinct():
    m=np.array([[1.,.9,.1,.2],[.9,1,.2,.1],[.1,.2,1,.8],[.2,.1,.8,1]])
    assert diverse_indices(m,[0,1,2,3],8,.6)==[0,2]
    assert diverse_indices(m,[0,1,2,3],1,.6)==[0]
    assert diverse_indices(np.empty((0,0)),[])==[]


def test_polymer_ligands_are_not_mislabeled_as_no_complex():
    raw=dict(rcsb_id='1ABC',struct=dict(title='Peptide complex'),exptl=[dict(method='X-RAY DIFFRACTION')],
        rcsb_entry_info=dict(resolution_combined=[1.5]),refine=[],nonpolymer_entities=[dict(pdbx_entity_nonpoly=dict(comp_id='SO4'))],
        polymer_entities=[dict(entity_poly=dict(rcsb_entity_polymer_type='Protein',pdbx_seq_one_letter_code_can='ACDE'),
            rcsb_polymer_entity=dict(pdbx_description='Peptide'),
            rcsb_polymer_entity_container_identifiers=dict(entity_id='2',reference_sequence_identifiers=[]))])
    row=summarize_entry(raw,'Q00987')
    assert row['polymer_complex_candidate'] and not row['nonpolymer_complex_candidate']
    assert row['other_polymers'][0]['length']==4


def test_modified_residue_sequence_metric_keeps_ccd_identity():
    assert token_similarity(['ALA','A1A2K'],['ALA','ALA'])==.5
    assert token_similarity(['ALA'],['ALA'])==1
    assert token_similarity(['ALA'],[])==0


def test_chat_diversity_plan_requires_explicit_site_reference():
    plan=dict(version=1,summary='Compare MDM2 crystal ligands',clarifications=[],steps=[
        dict(id='protein',action='protein_fetch',params=dict(accession='Q00987',organism_id=9606)),
        dict(id='diversity',action='structure_diversity',params=dict(protein_step='protein',reference_pdb='5C5A'))])
    assert validate_plan(plan)==plan
    for value in [0,True,float('nan'),1.5]:
        changed=json.loads(json.dumps(plan));changed['steps'][1]['params']['similarity_threshold']=value
        with pytest.raises(ValueError):validate_plan(changed)


def test_real_mdm2_generated_report_counts_when_available():
    # Optional local artifact audit; never requires a network call for the test suite.
    from pathlib import Path
    p=Path('D:/agent/MDM2/analysis/e052/report.json')
    if not p.exists():pytest.skip('Local MDM2 census not present')
    r=json.loads(p.read_text())
    assert r['quality_site_unique_ligands']>=len(r['proposed_references'])
    assert len({x['smiles'] for x in r['proposed_references']})==len(r['proposed_references'])
    assert all(x['quality_passed'] and x['same_reference_site'] for x in r['proposed_references'])


def test_cached_public_evidence_is_not_silently_changed(tmp_path):
    import hashlib
    from aidd_agent.structure_diversity import PublicCache, save
    url='https://example.invalid/public'
    key=hashlib.sha256((url+'null').encode()).hexdigest()
    payload=b'{"verified":true}'
    (tmp_path/(key+'.json')).write_bytes(payload)
    save(tmp_path/(key+'.source.json'),dict(url=url,sha256=hashlib.sha256(payload).hexdigest()))
    cache=PublicCache(tmp_path)
    assert cache.get(url)==dict(verified=True)
    (tmp_path/(key+'.json')).write_bytes(b'{}')
    with pytest.raises(ValueError,match='Changed cached'):
        cache.get(url)
