"""Graph-verified cyclic peptide units for the explicitly supported naming alphabet."""
from functools import lru_cache
import re

from rdkit import Chem

QUERY = Chem.MolFromSmarts('[N]-[C;X4]-[C](=[O])')


def tokens(name):
    if not name.startswith('c--') or not name.endswith('-c'):
        raise ValueError('unsupported_name_wrapper')
    text = name[3:-2].replace('-', '')
    found = re.findall(r'd?[AWLFVPG](?:nme|NMe)?', text)
    if ''.join(found) != text or not found:
        raise ValueError('unsupported_residue_name')
    return found


def signature(mol):
    Chem.SanitizeMol(mol)
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)


@lru_cache(maxsize=64)
def reference(token):
    d = token.startswith('d')
    text = token[1:] if d else token
    mol = Chem.MolFromSequence(text[0])
    n, ca, c, _ = mol.GetSubstructMatch(QUERY)
    if d:
        mol.GetAtomWithIdx(ca).InvertChirality()
    rw = Chem.RWMol(mol)
    if len(text) > 1:
        methyl = rw.AddAtom(Chem.Atom(6))
        rw.AddBond(n, methyl, Chem.BondType.SINGLE)
    hydroxyl = [a.GetIdx() for a in rw.GetAtomWithIdx(c).GetNeighbors()
                if a.GetAtomicNum() == 8 and rw.GetBondBetweenAtoms(c, a.GetIdx()).GetBondType() == Chem.BondType.SINGLE]
    if len(hydroxyl) != 1:
        raise ValueError('invalid_reference_amino_acid')
    rw.RemoveAtom(hydroxyl[0])
    return signature(rw.GetMol())


def peptide(mol, name=None):
    """Return unit atom indices in heavy-atom order, with cyclic name correspondence."""
    mol = Chem.RemoveHs(Chem.Mol(mol))
    unit_rows = list(mol.GetSubstructMatches(QUERY))
    expected = tokens(name) if name is not None else None
    if not unit_rows:
        raise ValueError('no_peptide_units')
    if expected is not None and len(unit_rows) != len(expected):
        raise ValueError('peptide_unit_count_mismatch')
    by_n = {r[0]: r for r in unit_rows}
    if len(by_n) != len(unit_rows):
        raise ValueError('ambiguous_peptide_units')
    following = {}
    for n, ca, c, o in unit_rows:
        nxt = [a.GetIdx() for a in mol.GetAtomWithIdx(c).GetNeighbors()
               if a.GetIdx() in by_n and a.GetIdx() != n
               and mol.GetBondBetweenAtoms(c, a.GetIdx()).GetBondType() == Chem.BondType.SINGLE]
        if len(nxt) != 1:
            raise ValueError('ambiguous_peptide_link')
        following[n] = nxt[0]
    start = min(by_n); ordered = []; current = start
    while current not in ordered:
        ordered.append(current); current = following[current]
    if current != start or len(ordered) != len(by_n):
        raise ValueError('multiple_or_open_peptide_cycles')
    ring = [a for n in ordered for a in by_n[n][:3]]
    if len(set(ring)) != len(ring):
        raise ValueError('overlapping_peptide_units')
    rw = Chem.RWMol(mol)
    for n in ordered:
        rw.RemoveBond(by_n[n][2], following[n])
    # Cap each cleaved peptide endpoint with one H on the analysis copy.
    # MOL2-derived atoms may forbid implicit H, unlike SMILES-derived atoms.
    # Without explicit caps the former would become radical fragments.
    for n in ordered:
        for index in (n, by_n[n][2]):
            atom = rw.GetAtomWithIdx(index)
            atom.SetNumExplicitHs(mol.GetAtomWithIdx(index).GetTotalNumHs() + 1)
            atom.SetNoImplicit(True)
            atom.SetNumRadicalElectrons(0)
    disconnected = rw.GetMol()
    Chem.SanitizeMol(disconnected)
    atom_maps = []
    fragments = Chem.GetMolFrags(disconnected, asMols=True, fragsMolAtomMapping=atom_maps)
    if len(fragments) != len(ordered):
        raise ValueError('crosslinked_residues_require_review')
    labels = {}
    for fragment, indices in zip(fragments, atom_maps):
        ns = set(indices) & set(ordered)
        if len(ns) != 1:
            raise ValueError('ambiguous_residue_fragment')
        labels[next(iter(ns))] = signature(fragment)
    if expected is None:
        return dict(tokens=None, name_rotations=None,
                    residue_signatures=[labels[n] for n in ordered],
                    units=[list(by_n[n]) for n in ordered], ring_atoms=ring)
    wanted = [reference(t) for t in expected]
    rotations = [i for i in range(len(ordered))
                 if [labels[n] for n in ordered[i:] + ordered[:i]] == wanted]
    if not rotations:
        raise ValueError('name_structure_residue_mismatch')
    offset = rotations[0]; ordered = ordered[offset:] + ordered[:offset]
    return dict(tokens=expected, name_rotations=len(rotations),
                residue_signatures=[labels[n] for n in ordered],
                units=[list(by_n[n]) for n in ordered],
                ring_atoms=[a for n in ordered for a in by_n[n][:3]])
