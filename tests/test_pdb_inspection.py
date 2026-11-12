import hashlib
import json
from pathlib import Path

import pytest

from aidd_agent.pdb_inspection import inspect_file, run
from aidd_agent.prompt_plan import validate_plan
from aidd_agent.prompt_workflow import create_plan, initialize_context, run_plan
from aidd_agent.prompt_smoke import OfflineServices
from aidd_agent.structure_review import structures


CIF = '''data_1ABC
_entry.id 1ABC
loop_
_entity.id
_entity.type
_entity.pdbx_description
1 polymer 'Test protein'
2 non-polymer GDP
3 non-polymer 'Magnesium ion'
4 water Water
loop_
_entity_poly.entity_id
_entity_poly.type
_entity_poly.pdbx_seq_one_letter_code_can
_entity_poly.pdbx_strand_id
1 'polypeptide(L)' AVG A,B
loop_
_struct_asym.id
_struct_asym.entity_id
A 1
B 1
C 2
D 3
E 4
loop_
_pdbx_struct_assembly.id
_pdbx_struct_assembly.oligomeric_details
_pdbx_struct_assembly.oligomeric_count
1 monomeric 1
2 monomeric 1
loop_
_pdbx_struct_assembly_gen.assembly_id
_pdbx_struct_assembly_gen.oper_expression
_pdbx_struct_assembly_gen.asym_id_list
1 1 A,C,D
2 1 B
loop_
_struct_ref.id
_struct_ref.entity_id
_struct_ref.db_name
_struct_ref.pdbx_db_accession
1 1 UNP P01116
loop_
_struct_ref_seq.align_id
_struct_ref_seq.ref_id
_struct_ref_seq.seq_align_beg
_struct_ref_seq.seq_align_end
_struct_ref_seq.db_align_beg
_struct_ref_seq.db_align_end
1 1 1 3 11 13
loop_
_struct_ref_seq_dif.align_id
_struct_ref_seq_dif.seq_num
_struct_ref_seq_dif.mon_id
_struct_ref_seq_dif.db_mon_id
_struct_ref_seq_dif.pdbx_seq_db_seq_num
_struct_ref_seq_dif.details
1 2 VAL GLY 12 'engineered mutation'
loop_
_atom_site.group_PDB
_atom_site.label_asym_id
_atom_site.auth_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.auth_seq_id
_atom_site.label_comp_id
_atom_site.type_symbol
_atom_site.pdbx_PDB_ins_code
_atom_site.pdbx_PDB_model_num
ATOM A X 1 1 11 ALA C ? 1
ATOM A X 1 1 11 ALA N ? 1
ATOM A X 1 2 12 VAL C ? 1
HETATM B Y 1 1 11 MSE C ? 1
ATOM B Y 1 3 13 GLY C ? 1
ATOM A X 1 3 13 GLY C ? 2
HETATM C X 2 . 201 GDP C ? 1
HETATM C X 2 . 201 GDP P ? 1
HETATM D X 3 . 202 MG MG A 1
HETATM E X 4 . 301 HOH O ? 1
'''


def fixture_file(tmp_path, content=CIF):
    path = tmp_path / 'test.cif'
    path.write_text(content, encoding='utf-8')
    return path


def test_deposited_evidence_not_assembly_or_isoform_inference(tmp_path):
    result = inspect_file(fixture_file(tmp_path), '1abc')
    assert result['polymers'][0]['deposited_sequence'] == 'AVG'
    assert len(result['asymmetric_unit_chains']) == 5
    assert [r['oligomeric_count'] for r in result['assemblies']] == ['1', '1']
    assert result['assembly_generators'][0]['asym_id_list'] == 'A,C,D'
    assert result['isoform_assignment'] == 'not_independently_verified'
    assert result['sequence_differences'][0]['pdbx_seq_db_seq_num'] == '12'
    assert result['database_alignments'][0]['db_align_beg'] == '11'
    chains = result['observed_polymer_chains']
    first = next(c for c in chains if c['label_asym_id'] == 'A' and c['model'] == '1')
    assert first['auth_asym_id'] == 'X'
    assert first['observed_residue_count'] == 2
    assert first['deposited_positions_without_coordinates'] == [3]
    second_model = next(c for c in chains if c['model'] == '2')
    assert second_model['deposited_positions_without_coordinates'] == [1, 2]
    ligands = result['nonwater_nonpolymer_instances']
    assert {r['comp_id'] for r in ligands} == {'GDP', 'MG'}
    assert next(r for r in ligands if r['comp_id'] == 'MG')['insertion_code'] == 'A'
    assert result['field_status']['deposited_unobserved_residues'] == 'not_reported'
    assert result['receptor_adopted'] is False


def test_missing_annotations_and_mismatched_entry(tmp_path):
    path = fixture_file(tmp_path, 'data_1ABC\n_entry.id 1ABC\n')
    result = inspect_file(path, '1ABC')
    assert result['field_status']['assemblies'] == 'not_reported'
    assert result['database_references'] == []
    with pytest.raises(ValueError, match='entry ID'):
        inspect_file(path, '2ABC')


class Services(OfflineServices):
    @staticmethod
    def pdb_fetch(code, output):
        output = Path(output)
        output.mkdir(parents=True, exist_ok=True)
        path = fixture_file(output)
        return dict(path=str(path), pdb_id=code, url='https://files.rcsb.org/download/1ABC.cif',
                    sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def test_plan_execution_persistence_and_structure_reuse(tmp_path):
    plan = dict(version=1, summary='Inspect deposited structure', clarifications=[],
        steps=[dict(id='inspection', action='pdb_inspect', params=dict(pdb_id='1ABC'))])
    validate_plan(plan)
    ctx = initialize_context(tmp_path / 'storage')
    path = create_plan(Path(ctx['db']), ctx['user_id'], ctx['project_id'],
                      'Inspect PDB 1ABC only', local_plan=plan)
    report, _ = run_plan(Path(ctx['db']), ctx['user_id'], ctx['project_id'], path,
                         services=Services(), allow_compute=False)
    assert report['status'] == 'complete'
    result = report['steps']['inspection']['result']
    saved = Path(result['report'])
    assert saved.is_relative_to(path.parent)
    assert hashlib.sha256(saved.read_bytes()).hexdigest() == result['report_sha256']
    assert json.loads(saved.read_text())['isoform_assignment'] == 'not_independently_verified'
    job = dict(status='complete', plan=str(path), report=str(path.parent / 'execution/report.json'))
    listed = structures(job, Path(ctx['project_root']))
    assert listed[0]['path'] == result['path']
    assert listed[0]['kind'] == 'experimental'
    Path(result['path']).write_text('modified', encoding='utf-8')
    with pytest.raises(ValueError, match='provenance'):
        structures(job, Path(ctx['project_root']))


def test_bad_digest_and_invalid_id(tmp_path):
    class Corrupt(Services):
        @staticmethod
        def pdb_fetch(code, output):
            return dict(Services.pdb_fetch(code, output), sha256='0' * 64)
    with pytest.raises(ValueError, match='digest'):
        run('1ABC', tmp_path, Corrupt())
    with pytest.raises(ValueError, match='PDB ID'):
        validate_plan(dict(version=1, summary='Inspect', clarifications=[], steps=[
            dict(id='inspection', action='pdb_inspect', params=dict(pdb_id='../bad'))]))
