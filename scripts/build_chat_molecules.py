"""Render illustrative noncanonical cyclic peptides and a schematic solvent field."""
import json
import math
from pathlib import Path
import random

from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Geometry import Point2D

ROOT = Path(__file__).resolve().parents[1] / 'src/aidd_agent/web'
EXAMPLES = [
    ('N-methyl and para-fluorophenyl', 'CN1[C@@H](Cc2ccc(F)cc2)C(=O)NCC(=O)N[C@@H](CC2CC2)C(=O)N[C@@H](COC)C(=O)N[C@H](C)C1=O'),
    ('Thienyl and cyclopropyl side chains', 'N1[C@@H](Cc2cccs2)C(=O)N(C)[C@@H](C)C(=O)N[C@@H](CC2CC2)C(=O)NCC(=O)N[C@@H](CC(C)C)C(=O)N[C@H](C)C1=O'),
    ('Pipecolic and naphthyl building blocks', 'O=C1N2CCCC[C@H]2C(=O)N[C@@H](Cc2ccc3ccccc3c2)C(=O)N(C)CC(=O)N[C@@H](CO)C(=O)N[C@@H](C(C)C)C1'),
    ('Proline, alkyne and fluorinated side chains', 'N1[C@@H](CC#C)C(=O)N2CCC[C@H]2C(=O)N[C@@H](CC(F)(F)F)C(=O)N(C)[C@@H](C)C(=O)NCC1=O'),
    ('Beta-amino acid and pyridyl building blocks', 'N1CC[C@H](C)C(=O)N[C@@H](Cc2ccncc2)C(=O)N(C)[C@@H](C)C(=O)N[C@@H](CC2CC2)C(=O)NCC1=O'),
    ('Methoxy and biphenyl side chains', 'CN1[C@@H](COC)C(=O)N[C@@H](Cc2ccc(-c3ccccc3)cc2)C(=O)NCC(=O)N[C@H](C(C)C)C(=O)N[C@@H](C)C1=O'),
]


def draw(mol, filename, width=380, height=280, line_width=1.55):
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    opts = drawer.drawOptions()
    opts.clearBackground = False
    opts.bondLineWidth = line_width
    opts.padding = .07
    opts.setAtomPalette({-1: (.13, .27, .28)})
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    (ROOT / filename).write_text(drawer.GetDrawingText(), encoding='utf-8')


def d_coordinates(count):
    # Artistic 2D coordinates only; this outline is not an energy-minimized pose.
    perimeter = 4 + 3 * math.pi + 6
    result = []
    for i in range(count):
        distance = i * perimeter / count
        if distance < 2:
            x, y = -2 + distance, 3
        elif distance < 2 + 3 * math.pi:
            angle = math.pi / 2 - (distance - 2) / 3
            x, y = 3 * math.cos(angle), 3 * math.sin(angle)
        elif distance < 4 + 3 * math.pi:
            x, y = -(distance - 2 - 3 * math.pi), -3
        else:
            x, y = -2, -3 + distance - 4 - 3 * math.pi
        result.append(Point2D(x, y))
    return result


def solvent_field():
    rng = random.Random(20261005)
    parts = ['<svg xmlns="http://www.w3.org/2000/svg" width="1600" height="340" viewBox="0 0 1600 340">',
             '<rect width="1600" height="340" fill="#edf7fb"/>',
             '<path d="M0 176 C240 152 390 189 650 168 S1040 156 1250 175 S1480 163 1600 171 V340 H0Z" fill="#f4eee5"/>',
             '<path d="M0 176 C240 152 390 189 650 168 S1040 156 1250 175 S1480 163 1600 171" fill="none" stroke="#d5e4e3" stroke-width="1.2"/>']
    for row in range(3):
        for col in range(15):
            x, y = col * 113 + rng.uniform(-15, 22), 18 + row * 51 + rng.uniform(-9, 12)
            angle = rng.uniform(-38, 38)
            parts.append(f'<g transform="translate({x:.1f} {y:.1f}) rotate({angle:.1f})" opacity=".25" stroke="#72abc1" stroke-width="1.4" fill="none"><path d="M-17 13 L-3 3 M3 3 L17 13"/><g fill="#568ca5" stroke="none" font-family="Arial,sans-serif" font-size="10"><text x="-4" y="3">O</text><text x="-24" y="23">H</text><text x="17" y="23">H</text></g></g>')
    for row in range(3):
        for col in range(12):
            x, y = col * 143 + rng.uniform(-17, 20), 216 + row * 48 + rng.uniform(-12, 12)
            angle = rng.uniform(-25, 25)
            skeleton = 'M-30 6 L-18 -1 L-6 6 L6 -1 L18 6 L30 -1' if (col + row) % 3 else 'M-12 0 L-6 -10 L6 -10 L12 0 L6 10 L-6 10Z'
            parts.append(f'<path d="{skeleton}" transform="translate({x:.1f} {y:.1f}) rotate({angle:.1f})" fill="none" stroke="#b1a284" stroke-width="1.5" opacity=".28"/>')
    for x, y, label in [(220, 40, 'H₂O'), (1290, 133, 'H₂O'), (570, 25, 'H₂O'), (1000, 307, 'C₆H₁₄'), (300, 282, 'C₆H₁₂')]:
        parts.append(f'<text x="{x}" y="{y}" fill="#789396" opacity=".24" font-family="Georgia,serif" font-size="15">{label}</text>')
    parts.append('</svg>')
    (ROOT / 'solvent-field.svg').write_text('\n'.join(parts), encoding='utf-8')


def main():
    metadata = []
    for index, (label, smiles) in enumerate(EXAMPLES, 1):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None: raise ValueError(label)
        rdDepictor.Compute2DCoords(mol)
        draw(mol, f'peptide-{index}.svg')
        metadata.append(dict(label=label, smiles=Chem.MolToSmiles(mol), role='Illustrative chemical graph; not a screened hit or synthesis claim'))
    logo_smiles = 'CN1[C@@H](C)C(=O)NC(C)(C)C(=O)N[C@@H](C2CC2)C(=O)N(C)CC(=O)N[C@@H](C)C1=O'
    logo = Chem.MolFromSmiles(logo_smiles)
    ring = max(logo.GetRingInfo().AtomRings(), key=len)
    rdDepictor.Compute2DCoords(logo, canonOrient=False, coordMap=dict(zip(ring, d_coordinates(len(ring)))))
    draw(logo, 'peptide-d.svg', 250, 300, line_width=2.6)
    (ROOT / 'peptides.json').write_text(json.dumps(dict(examples=metadata, wordmark=dict(smiles=Chem.MolToSmiles(logo), building_blocks='N-methyl, alpha-aminoisobutyric and cyclopropyl substitutions', layout='Constrained illustrative D-shaped 2D ring; not a physical conformation'), solvent_background='Schematic water above nonpolar hydrocarbon motifs; not a density or partitioning simulation'), indent=2) + '\n', encoding='utf-8')
    solvent_field()


if __name__ == '__main__':
    main()
