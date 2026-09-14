from pathlib import Path

import pytest

pytest.importorskip("rdkit")
from rdkit import RDConfig
from rdkit.Chem import ChemicalFeatures, Descriptors3D, rdMolDescriptors

from aidd_agent.mol2 import (
    Mol2Error, RDKIT_SANITIZATION_AROMATIC, RDKIT_SANITIZATION_STRICT,
    load_rdkit_mol2,
)


def _block(atom_lines: str, bond_lines: str, atoms: int, bonds: int) -> str:
    return f"""@<TRIPOS>MOLECULE
test
{atoms} {bonds} 0 0 0
SMALL
NO_CHARGES
@<TRIPOS>ATOM
{atom_lines}@<TRIPOS>BOND
{bond_lines}"""


def test_strict_mol2_stays_on_strict_path():
    raw = _block(
        "1 C 0 0 0 C.3 1 M 0\n2 N 1 0 0 N.3 1 M 0\n",
        "1 1 2 1\n", 2, 1)
    molecule, mode = load_rdkit_mol2(raw, "strict-test")
    assert molecule.GetNumAtoms() == 2
    assert mode == RDKIT_SANITIZATION_STRICT


def test_non_kekule_aromatic_graph_uses_controlled_fallback():
    raw = _block(
        "1 N 0 0 0 N.am 1 M 0\n"
        "2 C 1 0 0 C.ar 1 M 0\n"
        "3 C 2 0 0 C.ar 1 M 0\n"
        "4 N 2 1 0 N.ar 1 M 0\n"
        "5 C 1 2 0 C.ar 1 M 0\n"
        "6 C 0 1 0 C.ar 1 M 0\n",
        "1 1 2 1\n2 2 3 ar\n3 3 4 ar\n4 4 5 ar\n"
        "5 5 6 ar\n6 6 2 ar\n", 6, 6)
    molecule, mode = load_rdkit_mol2(raw, "aromatic-test")
    assert mode == RDKIT_SANITIZATION_AROMATIC
    assert sum(bond.GetIsAromatic() for bond in molecule.GetBonds()) == 5
    factory = ChemicalFeatures.BuildFeatureFactory(
        str(Path(RDConfig.RDDataDir) / "BaseFeatures.fdef"))
    assert factory.GetFeaturesForMol(molecule)
    assert Descriptors3D.PMI1(molecule) >= 0
    assert len(rdMolDescriptors.GetUSRCAT(molecule)) == 60


def test_fallback_does_not_admit_other_sanitization_failures():
    raw = _block(
        "1 C 0 0 0 C.3 1 M 0\n"
        "2 Cl 1 0 0 Cl 1 M 0\n3 Cl -1 0 0 Cl 1 M 0\n"
        "4 Cl 0 1 0 Cl 1 M 0\n5 Cl 0 -1 0 Cl 1 M 0\n"
        "6 Cl 0 0 1 Cl 1 M 0\n",
        "1 1 2 1\n2 1 3 1\n3 1 4 1\n4 1 5 1\n5 1 6 1\n", 6, 5)
    with pytest.raises(Mol2Error, match="non-kekulization sanitization failed"):
        load_rdkit_mol2(raw, "invalid-valence-test")
