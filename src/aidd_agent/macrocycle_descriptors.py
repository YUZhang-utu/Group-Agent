"""Directed peptide descriptors in proper local frames, with cyclic canonicalization."""
from functools import lru_cache
from pathlib import Path

import numpy as np
from rdkit import Chem, RDConfig
from rdkit.Chem import ChemicalFeatures

from .macrocycle_blocks import angle, amide_state, digest
from .macrocycle_peptide import peptide

VERSION = 'verified-peptide-local-frame-v1'
FAMILIES = ('Donor', 'Acceptor', 'Aromatic', 'Hydrophobe', 'LumpedHydrophobe',
            'PosIonizable', 'NegIonizable')
VARIANTS = ('backbone', 'chemistry', 'typed')


@lru_cache(maxsize=1)
def feature_factory():
    return ChemicalFeatures.BuildFeatureFactory(str(Path(RDConfig.RDDataDir) / 'BaseFeatures.fdef'))


def describe_peptide(mol, name=None, variant='typed'):
    if variant not in VARIANTS:
        raise ValueError('Unsupported descriptor variant')
    mol = Chem.RemoveHs(Chem.Mol(mol))
    mapping = peptide(mol, name)
    xyz = np.asarray(mol.GetConformer().GetPositions(), dtype=float)
    if not np.isfinite(xyz).all():
        raise ValueError('Nonfinite coordinates')
    units = mapping['units']; count = len(units)
    cuts = Chem.RWMol(mol)
    for i, (_, _, carbon, _) in enumerate(units):
        cuts.RemoveBond(carbon, units[(i+1) % count][0])
    fragments = Chem.GetMolFrags(cuts.GetMol(), sanitizeFrags=False)
    features = feature_factory().GetFeaturesForMol(mol) if variant != 'backbone' else []
    torsions = []; vectors = []; states = []; chirality = []; feature_maps = []
    for i, (n, ca, carbon, oxygen) in enumerate(units):
        previous, following = units[(i-1) % count], units[(i+1) % count]
        values = [angle(xyz[[previous[2], n, ca, carbon]]),
                  angle(xyz[[n, ca, carbon, following[0]]]),
                  angle(xyz[[ca, carbon, following[0], following[1]]])]
        torsions.append(values)
        states.append(amide_state(np.degrees(values[2])))
        atom = mol.GetAtomWithIdx(ca)
        chirality.append(atom.GetProp('_CIPCode') if atom.HasProp('_CIPCode') else 'achiral')
        vector = list(np.sin(values)) + list(np.cos(values))
        side = set(next(fragment for fragment in fragments if ca in fragment)) - {n, ca, carbon, oxygen}
        side_features = [(f.GetFamily(), list(f.GetAtomIds())) for f in features
                         if f.GetFamily() in FAMILIES and set(f.GetAtomIds()) <= side]
        feature_maps.append(side_features)
        if variant != 'backbone':
            # Fixed, versioned scaling; these weights are benchmark hypotheses.
            vector.extend([len(side)/10, sum(mol.GetAtomWithIdx(a).GetFormalCharge() for a in side),
                           float(mol.GetAtomWithIdx(n).IsInRingSize(5)),
                           sum(1 for a in mol.GetAtomWithIdx(n).GetNeighbors()
                               if a.GetIdx() in side and a.GetAtomicNum() == 6 and a.GetTotalNumHs() == 3)])
            vector.extend(sum(family == kind for family, _ in side_features)/4 for kind in FAMILIES)
        if variant == 'typed':
            e1 = xyz[carbon] - xyz[ca]; e1 /= np.linalg.norm(e1)
            e2 = xyz[n] - xyz[ca]; e2 -= e1 * np.dot(e1, e2)
            if np.linalg.norm(e2) < 1e-8:
                raise ValueError('Degenerate local backbone frame')
            e2 /= np.linalg.norm(e2); frame = np.stack((e1, e2, np.cross(e1, e2)), axis=1)
            for kind in FAMILIES:
                points = [(xyz[ids].mean(axis=0) - xyz[ca]) @ frame / 5
                          for family, ids in side_features if family == kind]
                vector.extend(np.mean(points, axis=0).tolist() if points else [0.0]*3)
                vector.extend(np.mean(np.square(points), axis=0).tolist() if points else [0.0]*3)
        vectors.append(vector)
    # Only directed rotations are allowed, never reverse traversal or reflection.
    labels = mapping['residue_signatures']
    keys = [tuple(labels[i:] + labels[:i]) for i in range(count)]
    candidates = [i for i, key in enumerate(keys) if key == min(keys)]
    offset = min(candidates, key=lambda i: (
        tuple(states[i:] + states[:i]),
        tuple(np.round(np.asarray(vectors[i:] + vectors[:i]).ravel(), 8))))
    rotate = lambda rows: rows[offset:] + rows[:offset]
    descriptor = np.asarray(rotate(vectors)).ravel()
    if not np.isfinite(descriptor).all():
        raise ValueError('Invalid peptide descriptor')
    return dict(version=VERSION, variant=variant, descriptor=descriptor.tolist(),
        hard_group=digest([VERSION, variant, count, rotate(chirality), rotate(states)]),
        units=rotate(units), ring_atoms=[a for u in rotate(units) for a in u[:3]],
        residue_signatures=rotate(labels), omega_states=rotate(states),
        torsions_radians=rotate(torsions), equivalent_rotations=len(candidates),
        typed_feature_atoms=rotate(feature_maps), name_mapping_verified=name is not None,
        metric='Euclidean distance on fixed scaled descriptor; no energy interpretation')
