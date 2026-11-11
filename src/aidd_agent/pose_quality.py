"""Independent, bounded chemical/geometry checks; never an affinity estimate."""
import json
from pathlib import Path

import numpy as np

from .expanded_wee1 import fingerprint


def receptor_atoms(path, chains=None, exclude=()):
    from Bio.PDB import PDBParser, MMCIFParser
    parser = PDBParser(QUIET=True) if Path(path).suffix.lower() == '.pdb' else MMCIFParser(QUIET=True)
    models = list(parser.get_structure('assembly', str(path)))
    if len(models) != 1:
        raise ValueError('QC requires one explicit assembly model')
    available = {c.id for c in models[0]}
    if chains and not set(chains) <= available:
        raise ValueError('Requested assembly chains are missing')
    atoms = []
    for chain in models[0]:
        if chains and chain.id not in chains:
            continue
        for residue in chain:
            number = str(residue.id[1]) + residue.id[2].strip()
            if (chain.id, number, residue.resname) in exclude:
                continue
            for atom in residue:
                if atom.is_disordered():
                    atom = sorted(atom.disordered_get_list(), key=lambda a: (-(a.occupancy or 0), a.altloc))[0]
                element = atom.element.upper()
                if element in {'H', 'D'}:
                    continue
                atoms.append(dict(chain=chain.id, residue=number, name=residue.resname,
                                  element=element, xyz=[float(x) for x in atom.coord]))
    if not atoms or not np.isfinite([a['xyz'] for a in atoms]).all():
        raise ValueError('Assembly has no finite heavy atom coordinates')
    return atoms


def assess(mol, atoms=(), expected_smiles=None, contacts=(), clash_ratio=.65):
    from rdkit import Chem
    if mol is None or mol.GetNumConformers() != 1:
        raise ValueError('QC requires one parsed molecular graph and pose')
    if not 0 < clash_ratio < 1:
        raise ValueError('clash_ratio must be between zero and one')
    mol = Chem.RemoveHs(Chem.Mol(mol))
    Chem.SanitizeMol(mol)
    xyz = np.asarray(mol.GetConformer().GetPositions())
    if not len(xyz) or not np.isfinite(xyz).all():
        raise ValueError('Invalid pose coordinates')
    table = Chem.GetPeriodicTable()
    numbers = [a.GetAtomicNum() for a in mol.GetAtoms()]
    radii = np.asarray([table.GetRvdw(n) for n in numbers])
    failures = []
    checks = {}
    canonical = Chem.MolToSmiles(mol, isomericSmiles=True)
    if expected_smiles is not None:
        expected = Chem.MolFromSmiles(expected_smiles)
        if expected is None:
            raise ValueError('Invalid expected SMILES')
        checks['chemical_identity'] = canonical == Chem.MolToSmiles(Chem.RemoveHs(expected), isomericSmiles=True)
        if not checks['chemical_identity']:
            failures.append('chemical_identity')
    else:
        checks['chemical_identity'] = None
    # Broad covalent-radius bounds are diagnostics, not a force-field validation.
    bad_bonds = []
    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        length = float(np.linalg.norm(xyz[i]-xyz[j]))
        scale = table.GetRcovalent(numbers[i]) + table.GetRcovalent(numbers[j])
        if not .55*scale <= length <= 1.45*scale:
            bad_bonds.append(dict(atoms=[i,j], angstrom=length))
    checks['bond_lengths'] = not bad_bonds
    if bad_bonds:
        failures.append('bond_lengths')
    geometric = Chem.Mol(mol)
    Chem.AssignAtomChiralTagsFromStructure(geometric, replaceExistingTags=True)
    Chem.AssignStereochemistry(geometric, cleanIt=True, force=True)
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    wrong_stereo = [a.GetIdx() for a in mol.GetAtoms() if a.HasProp('_CIPCode') and
        (not geometric.GetAtomWithIdx(a.GetIdx()).HasProp('_CIPCode') or
         a.GetProp('_CIPCode') != geometric.GetAtomWithIdx(a.GetIdx()).GetProp('_CIPCode'))]
    checks['specified_tetrahedral_geometry'] = not wrong_stereo
    if wrong_stereo:
        failures.append('specified_tetrahedral_geometry')
    distances = np.linalg.norm(xyz[:,None,:]-xyz[None,:,:],axis=2)
    graph_distance = Chem.GetDistanceMatrix(mol)
    internal = np.argwhere(np.triu((graph_distance > 2) & (distances < clash_ratio*(radii[:,None]+radii[None,:])),1))
    checks['internal_clashes'] = len(internal)==0
    if len(internal):
        failures.append('internal_clashes')
    collisions, contact_results = [], []
    # Chunk receptor atoms to keep memory bounded for large assemblies.
    for start in range(0,len(atoms),4096):
        subset = atoms[start:start+4096]
        rec = np.asarray([a['xyz'] for a in subset])
        rr = np.asarray([table.GetRvdw(table.GetAtomicNumber(a['element'].title())) for a in subset])
        d = np.linalg.norm(xyz[:,None,:]-rec[None,:,:],axis=2)
        for i,j in np.argwhere(d < clash_ratio*(radii[:,None]+rr[None,:])):
            collisions.append(dict(ligand_atom=int(i),chain=subset[j]['chain'],residue=subset[j]['residue'],
                                   name=subset[j]['name'],distance=float(d[i,j])))
    checks['assembly_clashes'] = not collisions if atoms else None
    if collisions:
        failures.append('assembly_clashes')
    for rule in contacts:
        selected = [a['xyz'] for a in atoms if a['chain']==rule['chain'] and a['residue']==str(rule['residue'])]
        nearest = float(np.min(np.linalg.norm(xyz[:,None,:]-np.asarray(selected)[None,:,:],axis=2))) if selected else None
        cutoff = float(rule.get('distance',4.5))
        if not np.isfinite(cutoff) or cutoff<=0:
            raise ValueError('Contact distance must be positive and finite')
        passed = nearest is not None and nearest<=cutoff
        contact_results.append(dict(rule=rule,nearest_angstrom=nearest,passed=passed))
        if not passed:
            failures.append('required_contact')
    return dict(execution_status='complete', qc_status='failed' if failures else 'passed_available_checks',
        experimental_status='not_assessed',checks=checks,failures=sorted(set(failures)),
        bad_bonds=bad_bonds,wrong_stereo_atoms=wrong_stereo,internal_clashes=internal.tolist(),
        assembly_clashes=collisions[:1000],assembly_clash_count=len(collisions),
        assembly_chains=sorted({a['chain'] for a in atoms}),contacts=contact_results,
        policy=dict(clash_ratio=clash_ratio,all_supplied_assembly_atoms=True),
        limitations=['Supplied coordinates only; biological assembly is not inferred',
                    'Broad geometry diagnostics; no strain-energy, affinity or experimental validation',
                    'Distance contacts are not hydrogen-bond or metal-coordination validation',
                    'Double-bond stereo geometry and bond-angle validation are not included'])


def execute(args,directory,root=None):
    from .target_query import load_ligand
    if args.get('set_id'):
        from .candidate_sets import export
        from .library_acceptance import read
        from .screening_selection import check_hashes
        from rdkit import Chem
        export_dir=Path(directory)/'candidate-poses';export_dir.mkdir(parents=True,exist_ok=True)
        result=export(dict(set_id=args['set_id'],count=args['count']),export_dir,root)
        manifest=read(Path(result['export_manifest']));check_hashes(manifest['outputs'])
        records=[]
        evidence=Path(result['export_manifest']).parent/'molecules.jsonl'
        for line in evidence.read_text(encoding='utf-8').splitlines():
            record=json.loads(line);folder=evidence.parent
            original=load_ligand(folder/record['original_mol2'])
            values={k:v for k,v in args.items() if k not in {'set_id','count'}}
            values.update(ligand_path=str(folder/record['posed_mol2']),expected_smiles=Chem.MolToSmiles(Chem.RemoveHs(original),isomericSmiles=True))
            records.append(dict(molecule_id=record['molecule_id'],**execute(values,directory)))
        return dict(kind='pose_qc',status='complete',set_id=args['set_id'],poses=records,
            checked=len(records),shortfall=result['shortfall'],
            qc_status='failed' if any(r['qc_status']=='failed' for r in records) else 'passed_available_checks' if records else 'not_assessed',
            experimental_status='not_assessed',export_manifest=result['export_manifest'],outputs=manifest['outputs'])
    mol = load_ligand(args['ligand_path'])
    atoms = receptor_atoms(args['receptor_path'],args.get('receptor_chains'),
        [tuple(x) for x in args.get('exclude_residues',[])]) if args.get('receptor_path') else []
    result = assess(mol,atoms,args.get('expected_smiles'),args.get('required_contacts',[]))
    result.update(kind='pose_qc',status='complete',sources=fingerprint([Path(args[k]) for k in
        ('ligand_path','receptor_path') if args.get(k)]))
    return result


def af3_ligands(model,input_path):
    """Map explicit ligand chains using the supplied SMILES; fail closed on ambiguity."""
    from Bio.PDB import MMCIFParser, PDBIO, Select
    from rdkit import Chem
    from rdkit.Chem import AllChem
    import io
    data = json.loads(Path(input_path).read_text(encoding='utf-8'))
    structure = MMCIFParser(QUIET=True).get_structure('af3',str(model))
    results = []
    for entity in data['sequences']:
        ligand = entity.get('ligand')
        if not ligand:
            continue
        chain = ligand['id']
        row = dict(chain=chain,qc_status='not_assessed')
        try:
            if not isinstance(chain,str) or not ligand.get('smiles'):
                raise ValueError('QC currently requires one chain ID and literal ligand SMILES')
            class Chain(Select):
                def accept_chain(self,item):
                    return item.id==chain
            writer=PDBIO();writer.set_structure(structure);stream=io.StringIO();writer.save(stream,Chain())
            pose=Chem.MolFromPDBBlock(stream.getvalue(),removeHs=True,sanitize=False)
            template=Chem.MolFromSmiles(ligand['smiles'])
            if pose is None or template is None or pose.GetNumHeavyAtoms()!=template.GetNumHeavyAtoms():
                raise ValueError('Ligand atom coverage does not match input chemistry')
            pose=AllChem.AssignBondOrdersFromTemplate(template,pose)
            # Verify geometry against specified stereochemistry after assigning the graph.
            matches=pose.GetSubstructMatches(template,useChirality=False,uniquify=False,maxMatches=1000)
            if not matches:
                raise ValueError('Cannot map supplied ligand chemistry to structure')
            mapped=Chem.Mol(template);mapped.RemoveAllConformers();conf=Chem.Conformer(mapped.GetNumAtoms())
            conf.Set3D(True)
            for i,j in enumerate(matches[0]):conf.SetAtomPosition(i,pose.GetConformer().GetAtomPosition(j))
            mapped.AddConformer(conf)
            atoms=receptor_atoms(model,chains=[c.id for c in next(iter(structure)) if c.id!=chain])
            row.update(assess(mapped,atoms,ligand['smiles']))
            row['atom_mapping']='first graph isomorphism; symmetric atom mappings may be ambiguous'
        except (ValueError,RuntimeError,KeyError,TypeError) as exc:
            row['error']=str(exc)
        results.append(row)
    return dict(ligands=results,experimental_status='not_assessed',
                qc_status='not_assessed' if not results or any(r['qc_status']=='not_assessed' for r in results)
                else 'failed' if any(r['qc_status']=='failed' for r in results) else 'passed_available_checks')
