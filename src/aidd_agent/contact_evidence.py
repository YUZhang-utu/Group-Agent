"""Target-independent atom proximity ledger, separate from pharmacophore scoring.

Every supported nearby pair is retained. Chemical labels are hypotheses, never
claims of hydrogen bonds, binding energy or absence of contacts in missing atoms.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from .anchor_extraction import protein_hbond_roles
from .consensus_admission import read_structure, ligand_atoms, reference_state_check
from .expanded_wee1 import fingerprint
from .screening_selection import check_hashes

STANDARD_RESIDUES = set('ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split())


def quality(atom):
    reasons = []
    if atom['label_alt_id']: reasons.append('alternate_location')
    if float(atom['occupancy']) < .9: reasons.append('low_occupancy')
    if atom['pdbx_PDB_ins_code']: reasons.append('insertion_code')
    return reasons


def atom_record(atom, matrix):
    point = np.asarray(atom['xyz']) @ matrix[:3, :3].T + matrix[:3, 3]
    if not np.isfinite(point).all(): raise ValueError('Nonfinite atom coordinate')
    return dict(chain=atom['auth_asym_id'], author_residue=atom['auth_seq_id'],
        residue_name=atom['auth_comp_id'], atom=atom['auth_atom_id'], element=atom['type_symbol'].upper(),
        canonical_residue=atom['canonical_residue'], model=atom['pdbx_PDB_model_num'],
        altloc=atom['label_alt_id'], occupancy=float(atom['occupancy']), point=point.tolist(),
        quality_issues=quality(atom))


def ligand_roles(ccd, atoms):
    from rdkit import RDConfig
    from rdkit.Chem import ChemicalFeatures
    from .chemistry_prep import _ccd_molecule
    coordinates = [dict(atom_id=a['auth_atom_id'], x=float(a['xyz'][0]), y=float(a['xyz'][1]),
                        z=float(a['xyz'][2])) for a in atoms]
    if len({a['atom_id'] for a in coordinates}) != len(coordinates):
        raise ValueError('Ambiguous ligand atom identities')
    molecule = _ccd_molecule(ccd, coordinates)
    factory = ChemicalFeatures.BuildFeatureFactory(str(Path(RDConfig.RDDataDir)/'BaseFeatures.fdef'))
    roles = {a.GetProp('_CCDAtomName'): dict(donor=False, acceptor=False,
        aromatic=a.GetIsAromatic(), formal_charge=a.GetFormalCharge(),
        feature_families=[]) for a in molecule.GetAtoms()}
    for feature in factory.GetFeaturesForMol(molecule):
        for i in feature.GetAtomIds():
            row = roles[molecule.GetAtomWithIdx(i).GetProp('_CCDAtomName')]
            row['feature_families'].append(feature.GetFamily())
            if feature.GetFamily() in {'Donor','Acceptor'}: row[feature.GetFamily().lower()] = True
    return roles


def pair_ledger(ligand, receptor, cutoff=4.5):
    """All observed heavy-atom pairs within cutoff; never nearest-partner-only."""
    if not np.isfinite(cutoff) or not 0 < cutoff <= 8: raise ValueError('Invalid contact cutoff')
    usable = [r for r in receptor if not r['quality_issues']]
    if not usable: return []
    tree = cKDTree([r['point'] for r in usable])
    rows = []
    for atom in ligand:
        if atom['quality_issues']: continue
        for j in sorted(tree.query_ball_point(atom['point'],cutoff)):
            other = usable[j]
            distance = float(np.linalg.norm(np.array(atom['point'])-other['point']))
            classes = ['heavy_atom_proximity']
            if atom['element'] in {'C','S'} and other['element'] in {'C','S'}:
                classes.append('carbon_sulfur_proximity')
            if atom['element'] in {'F','CL','BR','I'}:
                classes.append('halogen_proximity_not_halogen_bond')
            role = atom.get('chemistry',{})
            if distance <= 3.5 and ((role.get('donor') and other.get('acceptor')) or
                                    (role.get('acceptor') and other.get('donor'))):
                classes.append('donor_acceptor_proximity')
            identity = [atom['chain'], atom['author_residue'], atom['atom'], other['chain'],
                        other['author_residue'], other['atom']]
            rows.append(dict(pair_id=hashlib.sha256(json.dumps(identity).encode()).hexdigest()[:20],
                ligand_atom=atom['atom'], protein_chain=other['chain'],
                protein_author_residue=other['author_residue'], target_residue=other['canonical_residue'],
                protein_atom=other['atom'], distance=distance, classes=classes,
                status='geometric_hypothesis', angle_status='not_evaluated',
                protonation_status='unvalidated'))
    return rows


def build(survey, output):
    """Create ledger from a verified consensus; never mutate source query packages."""
    from .gaussian_batch import _load_query
    check_hashes(survey['sources'])
    prepared = {q['query_id']:q for q in survey['prepared_complexes']}
    cache = {}; complexes = []; sources = dict(survey['sources'])
    for row in survey['cohort']['admitted']:
        if row['query_id'] not in prepared: continue
        path = Path(row['structure_path'])
        if str(path) not in cache: cache[str(path)] = read_structure(path,survey['target'])
        structure = cache[str(path)]; chain = row['admission']['target_chain']
        if chain not in structure['chains']: raise ValueError('Target chain cannot be verified')
        matrix = np.asarray(row['admission']['transform'],float)
        if matrix.shape!=(4,4) or not np.isfinite(matrix).all() or not np.allclose(matrix[3],[0,0,0,1]) or not np.allclose(matrix[:3,:3].T@matrix[:3,:3],np.eye(3),atol=1e-6) or not np.isclose(np.linalg.det(matrix[:3,:3]),1):
            raise ValueError('Invalid protein rigid transform')
        original = ligand_atoms(structure,row)
        ligand = [atom_record(a,matrix) for a in original]
        receptor = [atom_record(a,matrix) for a in structure['chains'][chain]]
        manifest,_ = _load_query(Path(prepared[row['query_id']]['query_npz']).parent/'native.npz')
        ccd = Path(manifest['source']['ccd'])
        sources.update(fingerprint([path,ccd]))
        roles = {}; chemistry_status = 'CCD_perceived_protonation_unvalidated'
        try: roles = ligand_roles(ccd,original)
        except (ValueError,RuntimeError) as exc:
            chemistry_status = 'unresolved_'+type(exc).__name__
        for atom in ligand: atom['chemistry'] = roles.get(atom['atom'],{})
        for atom in receptor:
            supported = atom['residue_name'] in STANDARD_RESIDUES
            donor,acceptor = protein_hbond_roles(atom['residue_name'],atom['atom']) if supported else (False,False)
            atom.update(donor=donor,acceptor=acceptor,
                        chemistry_status='residue_template_protonation_unvalidated' if supported else 'unsupported_residue')
        pairs = pair_ledger(ligand,receptor)
        for pair in pairs: pair['pair_id']=row['query_id']+'/'+pair['pair_id']
        complexes.append(dict(query_id=row['query_id'],pdb_id=row['pdb_id'],target_chain=chain,
            chemotype=prepared[row['query_id']].get('chemotype'),
            chemistry_status=chemistry_status, ligand_atoms=ligand,receptor_atoms=receptor,pairs=pairs,
            counts=dict(pairs=len(pairs),excluded_ligand_atoms=sum(bool(a['quality_issues']) for a in ligand),
                        excluded_receptor_atoms=sum(bool(a['quality_issues']) for a in receptor)),
            observed_partner_atoms=prepared[row['query_id']]['observed_partner_atoms']))
    reference=next(c for c in complexes if c['query_id']==survey['cohort']['reference']['query_id'])
    ref_atoms=[dict(xyz=a['point'],label_alt_id=a['altloc'],occupancy=a['occupancy'],
                   canonical_residue=a['canonical_residue'],auth_atom_id=a['atom'],
                   pdbx_PDB_ins_code='uncertain' if a['quality_issues'] else '') for a in reference['receptor_atoms']]
    for c in complexes:
        c['reference_state_check']=reference_state_check(np.asarray([a['point'] for a in c['ligand_atoms'] if not a['quality_issues']]),ref_atoms)
        if c['counts']['excluded_ligand_atoms']:
            c['reference_state_check']['status']='unknown'
            c['reference_state_check']['reason']='uncertain_ligand_geometry'
    result=dict(version='atom-contact-ledger-v1',target=survey['target']['accession'],
        coordinate_frame=survey['cohort']['reference']['coordinate_frame'], cutoff_angstrom=4.5,
        complexes=complexes, sources=sources, implementation=fingerprint([Path(__file__)]), policy=dict(quality_minimum_occupancy=.9,
        missing_geometry='unknown_not_absent', contact_classes='proximity_not_validated_interactions',
        symmetry_mates='not_enumerated', waters_metals='not_in_direct_protein_pair_ledger'))
    Path(output).write_text(json.dumps(result,indent=2),encoding='utf-8')
    return result


def attach(survey, ledger_path):
    ledger = build(survey,ledger_path)
    survey['contact_evidence'] = dict(version=ledger['version'],path=str(Path(ledger_path).resolve()),
        complexes=len(ledger['complexes']),pairs=sum(c['counts']['pairs'] for c in ledger['complexes']),
        unresolved_reference_states=[c['query_id'] for c in ledger['complexes'] if c['reference_state_check']['status']!='no_severe_overlap'],
        scope='Complete observed heavy-atom proximity within 4.5 A; not validated bonds')
    survey.setdefault('outputs',{})['contacts.json']=str(Path(ledger_path).resolve())
    survey['sources'].update(ledger['sources'])
    survey['sources'].update(fingerprint([ledger_path]))


def main():
    parser=argparse.ArgumentParser(description='Augment old consensus into NEW evidence artifacts.')
    parser.add_argument('--survey',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    survey=json.loads(args.survey.read_text());check_hashes(survey['sources'])
    args.output.mkdir(parents=True,exist_ok=False)
    attach(survey,args.output/'contacts.json')
    survey['sources'].update(fingerprint([args.survey]))
    (args.output/'report.json').write_text(json.dumps(survey,indent=2),encoding='utf-8')
    from .consensus_report import render
    render(args.output)
    print(json.dumps(survey['contact_evidence']))


if __name__=='__main__':main()
