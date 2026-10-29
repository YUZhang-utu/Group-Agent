"""Target-independent receptor/reference queries for the existing guided engine."""
import copy
from pathlib import Path

import numpy as np

from .gaussian_batch import write_gaussian_query, _load_query, _atomic_json
from .expanded_wee1 import fingerprint


def load_receptor(path, chains, keep_residues=()):
    from Bio.PDB import PDBParser, MMCIFParser
    from Bio.PDB.Polypeptide import is_aa
    from .anchor_extraction import protein_hbond_roles
    parser = PDBParser(QUIET=True) if Path(path).suffix.lower() == '.pdb' else MMCIFParser(QUIET=True)
    structure = parser.get_structure('receptor', str(path))
    models = list(structure)
    if len(models) != 1: raise ValueError('Supply one explicit receptor model')
    available = {chain.id for chain in models[0]}
    if not set(chains) <= available: raise ValueError('Receptor chains not found: ' + ','.join(sorted(set(chains)-available)))
    keep = {(r['chain'], r['residue'], r['name']) for r in keep_residues}
    found, atoms, points = set(), [], []
    for chain in models[0]:
        for residue in chain:
            number = str(residue.id[1]) + residue.id[2].strip()
            identity = (chain.id, number, residue.resname)
            protein = (residue.id[0] == ' ' or is_aa(residue,standard=False)) and chain.id in chains
            if not protein and identity not in keep: continue
            if identity in keep: found.add(identity)
            if residue.resname in {'HOH', 'WAT', 'DOD'}: raise ValueError('This dry-screening adapter does not retain water')
            for atom in residue:
                if atom.is_disordered():
                    atom = sorted(atom.disordered_get_list(), key=lambda a: (-(a.occupancy or 0), a.altloc))[0]
                if atom.element.upper() in {'H', 'D'}: continue
                xyz = np.asarray(atom.coord, dtype=float)
                if not np.isfinite(xyz).all(): raise ValueError('Nonfinite receptor coordinate')
                points.append(xyz)
                if protein:
                    donor, acceptor = protein_hbond_roles(residue.resname, atom.name)
                    atoms.append(dict(chain=chain.id, residue_number=number, residue=residue.resname,
                        atom_name=atom.name, element=atom.element, donor=donor, acceptor=acceptor, xyz=xyz))
    if keep-found: raise ValueError('Requested retained cofactor/ion residue was not found')
    if not atoms: raise ValueError('Selected receptor chains have no protein heavy atoms')
    return atoms, np.asarray(points)


def load_ligand(path):
    from rdkit import Chem
    path = Path(path)
    if path.suffix.lower() == '.sdf':
        mols = list(Chem.SDMolSupplier(str(path), removeHs=False, sanitize=True))
        if len(mols) != 1: raise ValueError('Reference SDF must contain exactly one posed ligand')
        mol = mols[0]
    else:
        if path.read_text(encoding='utf-8').count('@<TRIPOS>MOLECULE') != 1:
            raise ValueError('Reference MOL2 must contain exactly one posed ligand')
        mol = Chem.MolFromMol2File(str(path), removeHs=False, sanitize=True)
    if mol is None or mol.GetNumConformers() != 1 or mol.GetNumHeavyAtoms() < 3:
        raise ValueError('Reference ligand requires a valid chemical graph and one nondegenerate pose')
    if len(Chem.GetMolFrags(mol)) != 1: raise ValueError('Reference ligand must be one connected molecule')
    if not np.isfinite(mol.GetConformer().GetPositions()).all(): raise ValueError('Nonfinite ligand coordinate')
    return mol


def prepare(params, output):
    from rdkit import RDConfig
    from rdkit.Chem import ChemicalFeatures
    from .conformer_artifacts import FEATURE_TYPES
    from .chemical_companion import _feature_direction
    from .classified_features import extra_hypotheses
    from scipy.spatial import cKDTree
    output = Path(output); output.mkdir(parents=True, exist_ok=True)
    atoms, receptor = load_receptor(params['receptor'], params['receptor_chains'], params.get('keep_residues', []))
    mol = load_ligand(params['reference_ligand'])
    xyz = mol.GetConformer().GetPositions()
    shape = xyz[[a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]]
    if np.linalg.norm(np.std(shape, axis=0)) < 1e-5: raise ValueError('Degenerate reference ligand coordinates')
    distance, _ = cKDTree(receptor).query(shape)
    if distance.min() > 6: raise ValueError('Reference ligand is outside the selected receptor site; supply receptor-frame coordinates')
    factory = ChemicalFeatures.BuildFeatureFactory(str(Path(RDConfig.RDDataDir)/'BaseFeatures.fdef'))
    features = [f for f in factory.GetFeaturesForMol(mol) if f.GetFamily() in FEATURE_TYPES]
    if not features: raise ValueError('Reference ligand has no supported features')
    points = np.asarray([xyz[list(f.GetAtomIds())].mean(axis=0) for f in features])
    directions = [_feature_direction(mol, f.GetFamily(), tuple(f.GetAtomIds())) for f in features]
    query = dict(shape_points=shape, feature_points=points,
        feature_types=np.asarray([FEATURE_TYPES[f.GetFamily()] for f in features]),
        feature_directions=np.asarray([d[0] for d in directions]),
        feature_direction_kinds=np.asarray([d[1] for d in directions]))
    anchors = []
    for i, feature in enumerate(features):
        family = feature.GetFamily()
        if family not in {'Donor', 'Acceptor'}: continue
        role = 'HBD' if family == 'Donor' else 'HBA'
        partners = [a for a in atoms if a['acceptor' if role == 'HBD' else 'donor'] and np.linalg.norm(a['xyz']-points[i]) <= 3.5]
        nearest = {}
        for atom in sorted(partners, key=lambda a: np.linalg.norm(a['xyz']-points[i])):
            nearest.setdefault((atom['chain'], atom['residue_number']), atom)
        for (chain, residue), atom in nearest.items():
            anchors.append(dict(anchor_id=f'{params["name"]}:{role}:{i}:{chain}:{residue}', feature_index=i,
                feature_class=role, evidence=dict(interaction='direct_hydrogen_bond_geometry',
                    distance_angstrom=float(np.linalg.norm(atom['xyz']-points[i])),
                    protein_partner={k: atom[k] for k in ('chain','residue_number','residue','atom_name')},
                    angle_status='not_evaluated')))
    anchors += extra_hypotheses(params['name'], query, atoms)
    for anchor in anchors:
        anchor['evidence_class']='user_reference_geometric_hypothesis'
        anchor['interpretation']='Feature supported by the supplied reference pose; not proof of an experimental complex or validated binding energy'
    if not anchors: raise ValueError('No supported reference/protein contact hypotheses; inspect site and chemistry')
    path = output/'query.npz'
    write_gaussian_query(path, **query, anchored_weights=np.ones(len(features)),
        anchor_feature_indices=sorted({a['feature_index'] for a in anchors}),
        source=dict(query_id=params['name'], receptor=str(params['receptor']), reference_ligand=str(params['reference_ligand']),
            provenance='user_selected_receptor_and_posed_reference', feature_atom_indices=[list(f.GetAtomIds()) for f in features]))
    receptor_path = output/'receptor.npz'; np.savez_compressed(receptor_path, points=receptor)
    files = [Path(params['receptor']), Path(params['reference_ligand']), path, path.with_suffix('.manifest.json'), receptor_path]
    result = dict(kind='target_query', status='complete', target=params['name'],
        query=dict(query_id=params['name'], query_npz=str(path.resolve()), receptor_npz=str(receptor_path.resolve()), anchors=anchors),
        receptor_chains=params['receptor_chains'], kept_residues=params.get('keep_residues', []),
        protein_residues=[dict(chain=c, residue=r) for c,r in sorted({(a['chain'],a['residue_number']) for a in atoms})],
        source_parameters=params, sources=fingerprint(files),
        native_reference_minimum_distance=float(distance.min()),
        limitations=['User-selected receptor and posed reference; no automatic target identity or pocket-state validation',
            'Highest-occupancy atom altloc; ties sorted by altloc ID',
            'Contact hypotheses are geometric; protonation and H-bond angles are not validated',
            'Modified amino-acid side-chain chemistry may lack contact templates',
            'Reference shape and atom distances are not cavity-volume measurement'])
    _atomic_json(output/'report.json', result)
    return result


def design(source, output, contacts=(), policy=None):
    from .library_acceptance import read
    from .screening_selection import check_hashes
    from .guided_filters import GuidedFilter, rule_passes
    from .interaction_fast import match_batch
    from .joint_coarse import validate_constraints
    source = Path(source); report = read(source); check_hashes(report['sources'])
    if report.get('kind') != 'target_query' or report.get('status') != 'complete': raise ValueError('Need a completed target query')
    q = copy.deepcopy(report['query']); anchors = q['anchors']; groups = []
    if contacts:
        for contact in contacts:
            group = [a['anchor_id'] for a in anchors if a['evidence']['protein_partner']['chain'] == contact['chain']
                and a['evidence']['protein_partner']['residue_number'] == contact['residue']
                and ('residue_name' not in contact or a['evidence']['protein_partner']['residue'] == contact['residue_name'])
                and ('feature_class' not in contact or a['feature_class'] == contact['feature_class'])]
            if not group: raise ValueError('Required residue/contact has no reference-supported feature: ' + str(contact))
            groups.append(group)
    else:
        # A broad default requires any observed contact, not every observed residue.
        nearest = {}
        for anchor in sorted(anchors, key=lambda a: a['evidence']['distance_angstrom']):
            nearest.setdefault(anchor['feature_index'], anchor['anchor_id'])
        groups = [list(nearest.values())]
    order = sorted({aid for group in groups for aid in group})
    if len(order) > 20: raise ValueError('More than 20 selected anchors; specify narrower residue/contact groups')
    settings = dict(minimum_score=.5, coarse_constraints=dict(heavy_atom_ratio=[.1,10.],
        maximum_extent_distance=10., minimum_feature_coverage=.01),
        pocket=dict(minimum_heavy_atom_distance=1.2, maximum_clashing_fraction=0.))
    if policy:
        if set(policy)-set(settings): raise ValueError('Unknown target search policy field')
        settings.update(copy.deepcopy(policy))
    validate_constraints(settings['coarse_constraints'])
    score = settings['minimum_score']; pocket = settings['pocket']
    if type(score) not in (int,float) or not 0 < score <= 1: raise ValueError('minimum_score must be in (0,1]')
    if not isinstance(pocket,dict) or set(pocket) != {'minimum_heavy_atom_distance','maximum_clashing_fraction'}:
        raise ValueError('Specify both pocket distance and clashing fraction')
    if any(type(v) not in (int,float) or not np.isfinite(v) for v in pocket.values()) or not (
            0 < pocket['minimum_heavy_atom_distance'] <= 5 and 0 <= pocket['maximum_clashing_fraction'] <= 1):
        raise ValueError('Invalid pocket distance or clashing fraction')
    atoms, _ = load_receptor(report['source_parameters']['receptor'], report['receptor_chains'])
    direct = [dict(points=[a['xyz'].tolist() for a in atoms if a['chain']==c['chain'] and a['residue_number']==c['residue']],
                   distance=c.get('distance',4.5), chain=c['chain'], residue=c['residue']) for c in contacts]
    adopted = dict(**settings, mandatory_anchors=[], alternative_groups=groups, optional_anchors=[],
        receptor_npz=q['receptor_npz'], direct_contact_groups=direct,
        selection_authority='explicit_target_search_request', version='generic-receptor-reference-v1')
    _, original = _load_query(Path(q['query_npz']))
    indices = sorted({a['feature_index'] for a in anchors})
    expanded = dict(original, anchor_feature_indices=np.asarray(indices))
    gate = GuidedFilter(q, expanded, adopted)
    if not gate.pocket_mask(original['shape_points'], np.eye(4)[None])[0]:
        raise ValueError('Reference fails receptor distance/contact policy; inspect inputs rather than relaxing automatically')
    match_query = tuple(expanded[k][indices] for k in ('feature_points','feature_types','feature_directions','feature_direction_kinds','anchored_weights'))
    _, assigned, scores = match_batch(match_query, original['feature_points'][None], original['feature_types'][None],
        original['feature_directions'][None], original['feature_direction_kinds'][None], np.eye(4)[None])
    columns = [indices.index(next(a['feature_index'] for a in anchors if a['anchor_id']==aid)) for aid in order]
    mask = sum(1<<i for i,col in enumerate(columns) if assigned[0,col]>=0 and scores[0,col]>=score)
    if not rule_passes(mask,order,adopted): raise ValueError('Reference fails contact-feature self control')
    result = dict(kind='anchor_design',status='complete',readiness='ready_for_guided_funnel',query=q,design=adopted,
        anchor_order=order, reference_control=dict(pocket_passed=True,anchor_rule_passed=True),
        sources={**report['sources'],**fingerprint([source])},
        limitations=report['limitations']+['Finite rigid seed search; failure is not molecular impossibility',
            'Broad coarse defaults are exploratory; no shape/color score cutoff or affinity claim',
            'Required residue proximity is rechecked in the same pose; typed features match the reference, not validated H-bond energies'])
    output = Path(output); output.mkdir(parents=True,exist_ok=True); _atomic_json(output/'report.json',result)
    return result
