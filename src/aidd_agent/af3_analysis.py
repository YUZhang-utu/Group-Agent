"""Validate supplied ligand chemistry and explain recorded AF3 confidence."""
import math
from collections import Counter


def validate_smiles_list(values):
    if not isinstance(values, list) or len(values) > 8:
        raise ValueError('ligand_smiles must be a list of at most eight SMILES strings')
    for value in values:
        if not isinstance(value, str) or not 1 <= len(value) <= 6000 or any(c.isspace() for c in value):
            raise ValueError('Use a nonempty SMILES without whitespace or a molecule name')
    return values


def ligand_evidence(smiles):
    from rdkit import Chem
    from rdkit import rdBase
    validate_smiles_list([smiles])
    mol = Chem.MolFromSmiles(smiles)
    if mol is None or mol.GetNumHeavyAtoms() == 0:
        raise ValueError('Ligand SMILES failed strict RDKit sanitization')
    if len(Chem.GetMolFrags(mol)) != 1:
        raise ValueError('Supply disconnected ligand components as separate SMILES entities')
    stereo = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    return dict(input_smiles=smiles, canonical_isomeric_smiles=Chem.MolToSmiles(mol, isomericSmiles=True),
                formal_charge=Chem.GetFormalCharge(mol), heavy_atoms=mol.GetNumHeavyAtoms(),
                unassigned_tetrahedral_centers=[i for i, label in stereo if label == '?'],
                rdkit_version=rdBase.rdkitVersion,
                policy='Original SMILES passed unchanged; no protonation, tautomer or stereoisomer enumeration')


MEANINGS = {
    'iptm': 'Interface predicted TM score (0 to 1): confidence in relative placement across chains. Not affinity or proof of binding.',
    'ptm': 'Predicted TM score (0 to 1): confidence in overall complex structure; does not isolate the ligand interface.',
    'chain_pair_iptm': 'Off-diagonal entries: pairwise interface confidence. Diagonal entries: within-chain pTM. Interpret with recorded chain_ids.',
    'chain_iptm': 'Per-chain interface confidence. Not a per-chain binding energy.',
    'chain_ptm': 'Per-chain structural confidence.',
    'chain_pair_pae_min': 'Minimum predicted aligned error for chain pairs in angstroms. A minimum alone does not establish a reliable entire interface.',
    'ranking_score': 'AF3 sample-selection score: 0.8*ipTM + 0.2*pTM + 0.5*fraction_disordered - 100*has_clash. Not affinity or a probability.',
    'fraction_disordered': 'Model-estimated disordered fraction used in AF3 ranking.',
    'has_clash': 'AF3 clash flag; inspect the structure and chemistry if set.',
    'plddt': 'Local per-atom confidence (0 to 100), stored in AF3 model B-factor fields; these are not experimental temperature factors.',
    'pae': 'Predicted aligned error in angstroms between tokens. Lower values indicate greater confidence in relative placement.',
}


def assess_chain_mapping(confidence):
    """Validate summary labels without inferring matrix order from input or counts."""
    labels = confidence.get('chain_ids')
    valid_labels = (isinstance(labels, list) and bool(labels)
                    and all(isinstance(x, str) and bool(x) for x in labels))
    warnings = []
    counts = dict(Counter(labels)) if valid_labels else {}
    if not valid_labels:
        warnings.append('Summary chain_ids are absent or invalid; matrix chain mapping is unverified.')
    elif len(counts) != len(labels):
        warnings.append('Summary chain_ids contain repeated labels. They are not a unique chain list; '
                        'do not deduplicate them to assign matrix order. Token/atom granularity is unverified.')
    dimensions = {}
    for key in ('chain_iptm', 'chain_ptm', 'chain_pair_iptm', 'chain_pair_pae_min'):
        if key not in confidence:
            continue
        value = confidence[key]
        if key.startswith('chain_pair_'):
            shape_ok = (isinstance(value, list) and all(isinstance(row, list) and len(row) == len(value)
                                                       for row in value))
        else:
            shape_ok = isinstance(value, list) and all(not isinstance(x, (dict, list)) for x in value)
        if not isinstance(value, list) or not value or not shape_ok:
            warnings.append(f'{key} has an invalid chain-level shape.')
            dimensions[key] = None
        else:
            dimensions[key] = len(value)
    sizes = {size for size in dimensions.values() if size is not None}
    if len(sizes) > 1 or (valid_labels and any(size != len(labels) for size in sizes)):
        warnings.append('Chain-level score dimensions do not match each other or summary chain_ids.')
    verified = valid_labels and not warnings
    return dict(status='summary_labels_consistent' if verified else 'unverified',
                matrix_chain_ids=list(labels) if verified else None,
                label_count=len(labels) if isinstance(labels, list) else None,
                label_counts=counts, metric_dimensions=dimensions, warnings=warnings)


def compact_confidence(confidence):
    """Return a display copy; raw arrays remain in the recorded source artifacts."""
    result = dict(confidence)
    for key in ('chain_ids', 'token_chain_ids', 'atom_chain_ids'):
        labels = result.get(key)
        if (isinstance(labels, list) and all(isinstance(x, str) for x in labels)
                and (len(labels) > 16 or len(set(labels)) != len(labels))):
            result[key] = dict(representation='label_counts_not_matrix_order',
                               total_entries=len(labels), counts=dict(Counter(labels)))
    return result


def explain_confidence(confidence, chain_ids=()):
    def clean(value):
        if isinstance(value, float) and not math.isfinite(value): return None
        if isinstance(value, list): return [clean(v) for v in value]
        if isinstance(value, dict): return {k:clean(v) for k,v in value.items()}
        return value
    return dict(metrics={key:dict(value=clean(value), meaning=MEANINGS.get(key, 'Raw AF3 output field.'))
                         for key,value in compact_confidence(confidence).items()},
                chain_mapping=assess_chain_mapping(confidence),
                input_chain_ids=list(chain_ids),
                additional_metrics={key:MEANINGS[key] for key in ('plddt','pae')},
                limitations=['No affinity, activity or validated binding conclusion follows from confidence alone.',
                             'No universal ligand acceptance threshold is applied.',
                             'Only use unique, dimension-consistent summary chain_ids for matrix order. Repeated labels and input chain order are not verified matrix mappings.'],
                documentation='https://github.com/google-deepmind/alphafold3/blob/main/docs/output.md')
