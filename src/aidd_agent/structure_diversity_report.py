"""Standalone human review of crystal evidence and reference diversity."""
import csv
import html
import json
from pathlib import Path


def render(root):
    from rdkit import Chem
    from rdkit.Chem.Draw import rdMolDraw2D
    root=Path(root);r=json.loads((root/'report.json').read_text())
    entries={e['pdb_id']:e for e in json.loads((root/'entries.json').read_text())['entries']}
    e=lambda value:html.escape(str(value))
    parts=['<!doctype html><html lang="en"><meta charset="utf-8"><title>Crystal reference review</title>',
        '<style>body{font:16px system-ui;color:#182c3e;background:#f3f6f8;margin:0}main{max-width:1200px;margin:auto;padding:32px}h1{font-size:38px}h2{margin-top:36px}p{line-height:1.6}.card{background:white;border:1px solid #d9e2e9;padding:20px;border-radius:12px;margin:12px 0}table{border-collapse:collapse;width:100%;background:white}td,th{text-align:left;padding:10px;border-bottom:1px solid #dce5ea}th{background:#e6edf2}.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(330px,1fr));gap:16px}.note{border-left:4px solid #c47e19;padding:12px;background:#fff4dc}.scroll{overflow:auto;max-height:640px}.matrix td{font-size:10px;padding:4px;min-width:30px}.matrix th{font-size:10px;position:sticky;top:0}a{color:#006e8a}</style><main>',
        f'<h1>{e(r["target"])} / crystal reference review</h1>',
        f'<p>Verified target <b>{e(r["target"])}</b> · {r["discovered_entries"]} experimental PDB entries · {r["downloaded_entries"]} coordinate files downloaded.</p>',
        '<p class="note">Exploratory reference proposal. No library screening was started. Resolution is one quality indicator; ligand density and biological relevance remain review items.</p>',
        '<h2>1. Resolution census</h2><table><tr><th>Resolution / A</th><th>X-ray entries</th><th>Nonpolymer candidates</th><th>Other-polymer candidates</th><th>Organic contacts verified</th><th>Peptide contacts verified</th></tr>']
    for b in r['census']:
        parts.append('<tr>'+''.join(f'<td>{e(b[k])}</td>' for k in ('bin','xray_entries','nonpolymer_candidates','other_polymer_candidates','coordinate_checked_organic_contact_entries','coordinate_checked_peptide_contact_entries'))+'</tr>')
    parts+=['</table><p>Candidate categories overlap. Contact checks cover downloaded structures at ≤3 A; a zero outside that range means not assessed. Contacts do not establish biological binding.</p>',
        '<h2>2. Your supplied structures</h2><table><tr><th>PDB</th><th>Resolution / A</th><th>Deposited ligand class</th><th>Local coordinate chains</th></tr>']
    for local in r['local_structures']:
        entry=entries.get(local['pdb_id'],{})
        kind='; '.join(entry.get('ligand_ids',[])+[p['description'] or 'polymer' for p in entry.get('other_polymers',[])])
        parts.append(f'<tr><td><a href="https://www.rcsb.org/structure/{e(local["pdb_id"])}">{e(local["pdb_id"])}</a></td><td>{e(entry.get("resolution"))}</td><td>{e(kind)}</td><td>{e(", ".join(local["chains"]))}</td></tr>')
    parts+=['</table><p>Aligned local exports are preserved. Full original coordinates supply missing metadata and chemical-component definitions.</p>',
        '<h2>3. Diverse small-molecule proposal</h2>',
        f'<p>{r["quality_site_unique_ligands"]} unique ligands passed the resolution, occupancy, heavy-atom count and pocket-contact screen. The displayed {len(r["proposed_references"])} references are a review budget, not an exhaustive representation.</p><div class="grid">']
    for ref in r['proposed_references']:
        drawer=rdMolDraw2D.MolDraw2DSVG(360,240);drawer.DrawMolecule(Chem.MolFromSmiles(ref['smiles']));drawer.FinishDrawing()
        svg=drawer.GetDrawingText();svg=svg[svg.index('<svg'):]
        parts.append(f'<section class="card"><b>{e(ref["query_id"])}</b>{svg}<p>Resolution: {ref["resolution"]} A · R-free: {e(ref["r_free"])}<br>Heavy atoms: {ref["heavy_atoms"]} · Minimum occupancy: {ref["minimum_occupancy"]}<br>Density validation: not evaluated</p></section>')
    parts+=['</div><h2>4. Reference count versus chemical coverage</h2><table><tr><th>References</th><th>Ligands with fingerprint similarity ≥0.6</th></tr>',
        *[f'<tr><td>{x["references"]}</td><td>{x["covered"]} / {r["quality_site_unique_ligands"]}</td></tr>' for x in r['coverage_curve']],
        '</table><p>Morgan radius 2, 2048 bits, chirality-aware Tanimoto. Chemical coverage is not 3D similarity or screening recall.</p>',
        '<h2>5. Peptide / macrocycle track</h2>',
        f'<p>{r["polymer_ligands"]["unique_sequence_link_variants"]} distinct deposited residue-sequence/link variants at the reference pocket. The examples below require chemical-graph and query preparation; they are not interchangeable with CCD small-molecule references.</p>',
        '<table><tr><th>Example</th><th>Resolution / A</th><th>Residues</th><th>Deposited description</th></tr>',
        *[f'<tr><td>{e(x["query_id"])}</td><td>{x["resolution"]}</td><td>{x["length"]}</td><td>{e(x["description"])}</td></tr>' for x in r['polymer_ligands']['sequence_diversity_examples']],
        '</table><p>Peptide comparison uses modified-residue CCD-token sequence edit similarity. Cyclization links are listed separately; this is not full chemical similarity.</p>',
        '<h2>6. Full small-molecule similarity matrix</h2><div class="scroll"><table class="matrix">']
    with (root/'ligand-similarity.csv').open(encoding='utf-8') as f:
        rows=list(csv.reader(f))
    parts.append('<tr>'+''.join('<th>'+e(x)+'</th>' for x in rows[0])+'</tr>')
    for row in rows[1:]:
        parts.append('<tr><th>'+e(row[0])+'</th>'+''.join(f'<td style="background:rgba(0,140,150,{float(x)*.65:.2f})">{float(x):.2f}</td>' for x in row[1:])+'</tr>')
    parts+=['</table></div><h2>7. Unresolved mapping</h2><ul>',
        *[f'<li>{e(x["pdb_id"])}: {e(x.get("error"))}</li>' for x in r['failed_entries']],
        '</ul><p><a href="report.json">Full report JSON</a> · <a href="ligand-similarity.csv">Chemical similarity CSV</a> · <a href="polymer-ligands.json">Polymer contacts and links</a> · <a href="peptide-sequence-similarity.csv">Peptide sequence similarity CSV</a></p></main></html>']
    (root/'review.html').write_text('\n'.join(parts),encoding='utf-8')
    return root/'review.html'
