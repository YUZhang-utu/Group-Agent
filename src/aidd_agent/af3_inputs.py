"""Validated literal sequence and chain composition for generic AF3 prompts."""
import hashlib
import re


def sequence(value):
    if not isinstance(value, str):
        raise ValueError('Protein sequence must be text')
    cleaned = re.sub(r'\s+', '', value).upper()
    if not 1 <= len(cleaned) <= 50000 or not re.fullmatch('[ACDEFGHIKLMNPQRSTVWY]+', cleaned):
        raise ValueError('Use a nonempty protein sequence with the 20 standard amino-acid letters')
    return cleaned


def seeds(params):
    if 'seeds' in params and ('seed_count' in params or 'seed_start' in params):
        raise ValueError('Choose an explicit seeds list OR seed_count with optional seed_start')
    if 'seed_start' in params and 'seed_count' not in params:
        raise ValueError('seed_start requires seed_count')
    count, start = params.get('seed_count', 1), params.get('seed_start', 1)
    if type(count) is not int or not 1 <= count <= 100 or type(start) is not int or not 0 <= start <= 2**31-count:
        raise ValueError('seed_count must be 1..100 and seed_start must keep seeds within 0..2147483647')
    result = params.get('seeds', list(range(start, start+count)))
    if not isinstance(result, list) or not 1 <= len(result) <= 100 or any(
            type(s) is not int or not 0 <= s <= 2**31-1 for s in result) or len(set(result)) != len(result):
        raise ValueError('Invalid AF3 seeds: use 1..100 distinct nonnegative integers')
    return result


def validate(params, seen):
    modes = [key for key in ('protein_step', 'protein_steps', 'protein_sequences') if key in params]
    if len(modes) != 1:
        raise ValueError('Choose exactly one of protein_step, protein_steps or protein_sequences')
    mode = modes[0]
    values = [params[mode]] if mode == 'protein_step' else params[mode]
    if not isinstance(values, list) or not 1 <= len(values) <= 26:
        raise ValueError('Supply 1..26 protein chains')
    if mode == 'protein_sequences':
        for value in values: sequence(value)
    elif any(not isinstance(ref, str) or seen.get(ref) not in {'protein_fetch', 'protein_resolve', 'protein_from_pdb'} for ref in values):
        raise ValueError('Protein chains must refer to preceding verified protein steps')
    copies = params.get('protein_copies', 1)
    if type(copies) is not int or not 1 <= copies <= 8:
        raise ValueError('AF3 protein_copies must be an integer from 1 to 8')
    if copies != 1 and len(values) != 1:
        raise ValueError('For mixed chains repeat entries explicitly; protein_copies applies to one sequence only')
    if len(values)*copies + len(params.get('ligand_ccd', [])) + len(params.get('ligand_smiles', [])) > 26:
        raise ValueError('At most 26 total protein and ligand chains are supported')
    if len(values) > 1 and ('start' in params or 'end' in params):
        raise ValueError('Mixed-chain inputs require complete supplied constructs; global start/end is ambiguous')
    seeds(params)


def literal_check(params, prompt):
    """Only normalize whitespace/case in sequence text actually supplied by the user."""
    normalized = re.sub(r'\s+', '', prompt).upper()
    for value in params.get('protein_sequences', []):
        if sequence(value) not in normalized:
            raise ValueError('Protein sequences must be copied from the original user request')


def resolve(params, execution):
    from .library_acceptance import read
    if 'protein_sequences' in params:
        records = [dict(sequence=sequence(value), accession=None, source_url=None,
            provenance='literal_user_sequence_not_database_verified') for value in params['protein_sequences']]
    else:
        refs = params.get('protein_steps', [params.get('protein_step')])
        records = [read(execution/ref/'protein.json') for ref in refs]
    records = records * params.get('protein_copies', 1)
    sequences, evidence = [], []
    for i, record in enumerate(records):
        full = sequence(record['sequence'])
        start, end = params.get('start', 1), params.get('end', len(full))
        if end > len(full): raise ValueError('Construct exceeds supplied protein sequence')
        construct = full[start-1:end]
        sequences.append(construct)
        evidence.append(dict(chain_id=chr(65+i), accession=record.get('accession'), source_url=record.get('source_url'),
            provenance=record.get('provenance', 'verified_protein_record'), start=start, end=end,
            full_sequence_sha256=hashlib.sha256(full.encode()).hexdigest(),
            construct_sequence_sha256=hashlib.sha256(construct.encode()).hexdigest()))
    return sequences, evidence
