"""Inspect deposited mmCIF evidence without adopting a receptor or running science."""
import hashlib
import json
from pathlib import Path
import re

from Bio.PDB.MMCIF2Dict import MMCIF2Dict


def category(data, prefix):
    """Keep all deposited columns, preserving absent values as null."""
    columns = {key[len(prefix) + 1:]: value if isinstance(value, list) else [value]
               for key, value in data.items() if key.startswith(prefix + '.')}
    if not columns:
        return []
    lengths = {len(value) for value in columns.values()}
    if len(lengths) != 1:
        raise ValueError('Inconsistent mmCIF category: ' + prefix)
    return [{key: None if value[i] in ('.', '?', '') else value[i]
             for key, value in columns.items()} for i in range(next(iter(lengths)))]


def inspect_file(path, pdb_id):
    path = Path(path).resolve()
    data = MMCIF2Dict(str(path))
    entries = category(data, '_entry')
    if len(entries) != 1 or entries[0].get('id', '').upper() != pdb_id.upper():
        raise ValueError('Deposited entry ID does not match requested PDB')
    entities = {r['id']: r for r in category(data, '_entity')}
    polymers = []
    for row in category(data, '_entity_poly'):
        sequence = re.sub(r'\s+', '', row.get('pdbx_seq_one_letter_code_can') or '')
        polymers.append(dict(entity_id=row.get('entity_id'), type=row.get('type'),
            description=entities.get(row.get('entity_id'), {}).get('pdbx_description'),
            deposited_sequence=sequence or None,
            sequence_length=len(sequence) if re.fullmatch('[A-Z]+', sequence) else None,
            deposited_author_chain_ids=row.get('pdbx_strand_id'),
            sequence_provenance='Deposited canonical polymer sequence; not a fetched reference'))
    polymer_entities = {row['entity_id']: row for row in polymers}
    asym = category(data, '_struct_asym')
    asym_entities = {row['id']: row.get('entity_id') for row in asym}
    chains = {}
    nonpoly = {}
    for atom in category(data, '_atom_site'):
        label = atom.get('label_asym_id')
        entity = atom.get('label_entity_id') or asym_entities.get(label)
        model = atom.get('pdbx_PDB_model_num') or '1'
        author = atom.get('auth_asym_id')
        if entity in polymer_entities:
            key = (model, label, author, entity)
            chain = chains.setdefault(key, set())
            number = atom.get('label_seq_id')
            if number and number.isdigit():
                chain.add(int(number))
        elif entities.get(entity, {}).get('type') != 'water' and atom.get('label_comp_id') not in {'HOH', 'DOD'}:
            key = (model, label, author, atom.get('auth_seq_id'),
                   atom.get('pdbx_PDB_ins_code'), atom.get('label_comp_id'), entity)
            nonpoly.setdefault(key, set()).add(atom.get('type_symbol'))
    observed = []
    for (model, label, author, entity), numbers in sorted(chains.items(), key=lambda item: str(item[0])):
        length = polymer_entities[entity]['sequence_length']
        valid = bool(numbers) and length is not None and all(1 <= n <= length for n in numbers)
        observed.append(dict(model=model, label_asym_id=label, auth_asym_id=author,
            entity_id=entity, observed_label_seq_ids=sorted(numbers),
            observed_residue_count=len(numbers) if numbers else None,
            deposited_positions_without_coordinates=sorted(set(range(1, length + 1)) - numbers) if valid else None,
            coordinate_coverage=len(numbers) / length if valid else None,
            coverage_basis='Per-model label_seq_id versus deposited entity sequence; not full UniProt coverage'))
    components = {row['id']: row for row in category(data, '_chem_comp')}
    instances = []
    for (model, label, author, residue, insertion, comp, entity), elements in sorted(nonpoly.items(), key=lambda item: str(item[0])):
        instances.append(dict(model=model, label_asym_id=label, auth_asym_id=author,
            auth_seq_id=residue, insertion_code=insertion, comp_id=comp, entity_id=entity,
            name=components.get(comp, {}).get('name'),
            component_type=components.get(comp, {}).get('type'),
            elements=sorted(e for e in elements if e), role='Unassigned; not an adopted screening reference'))
    mappings = category(data, '_struct_ref')
    report = dict(schema_version=1, kind='pdb_inspection', pdb_id=pdb_id.upper(),
        status='inspection_completed', structure=str(path),
        structure_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        title=category(data, '_struct'), experimental_methods=category(data, '_exptl'),
        refinement=category(data, '_refine'), em_reconstruction=category(data, '_em_3d_reconstruction'),
        polymers=polymers, asymmetric_unit_chains=asym, observed_polymer_chains=observed,
        assemblies=category(data, '_pdbx_struct_assembly'),
        assembly_generators=category(data, '_pdbx_struct_assembly_gen'),
        assembly_operations=category(data, '_pdbx_struct_oper_list'),
        nonwater_nonpolymer_instances=instances,
        deposited_unobserved_residues=category(data, '_pdbx_unobs_or_zero_occ_residues'),
        database_references=mappings, database_alignments=category(data, '_struct_ref_seq'),
        sequence_differences=category(data, '_struct_ref_seq_dif'),
        isoform_assignment='not_independently_verified', receptor_adopted=False,
        limitations=[
            'Reference alignments and sequence differences are deposited annotations, not independent sequence verification.',
            'An accession or shared construct does not establish an isoform; no KRAS4A/4B assignment is inferred.',
            'Assembly annotations and operators are reported, not expanded into new coordinates or biologically validated.',
            'Unreported metadata is not evidence of absence. Coordinate coverage is relative to the deposited construct.',
            'Ligand roles, binding pockets, affinity and pose quality are not determined by this inspection.'])
    report['field_status'] = {key: 'reported' if report[key] else 'not_reported'
        for key in ('polymers', 'assemblies', 'database_references', 'database_alignments',
                    'sequence_differences', 'deposited_unobserved_residues')}
    return report


def run(pdb_id, directory, services):
    directory = Path(directory).resolve()
    directory.mkdir(parents=True, exist_ok=True)
    source = services.pdb_fetch(pdb_id.upper(), directory)
    report = inspect_file(source['path'], pdb_id)
    if source.get('sha256') and source['sha256'] != report['structure_sha256']:
        raise ValueError('Downloaded structure digest mismatch')
    report['source'] = source
    report['evidence_mode'] = services.mode
    destination = directory / 'inspection-report.json'
    destination.write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    return dict(status='inspection_completed', pdb_id=pdb_id.upper(),
        path=report['structure'], sha256=report['structure_sha256'],
        report=str(destination), report_sha256=hashlib.sha256(destination.read_bytes()).hexdigest(),
        source_url=source.get('url'), evidence_mode=services.mode,
        polymer_entity_count=len(report['polymers']),
        assembly_annotation_count=len(report['assemblies']),
        nonwater_nonpolymer_instance_count=len(report['nonwater_nonpolymer_instances']),
        isoform_assignment=report['isoform_assignment'], receptor_adopted=False)
