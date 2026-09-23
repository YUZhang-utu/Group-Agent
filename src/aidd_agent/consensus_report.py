"""Human-readable admission, denominator and template evidence."""
import html
import json
from pathlib import Path


def render(output):
    output=Path(output);r=json.loads((output/'report.json').read_text(encoding='utf-8'))
    e=lambda x:html.escape(str(x))
    rows=['<!doctype html><html lang="en"><meta charset="utf-8"><title>Pocket consensus review</title>',
        '<style>body{font:15px system-ui;margin:32px;color:#183045}table{border-collapse:collapse;width:100%}td,th{border:1px solid #ccd;padding:8px;text-align:left}th{background:#eef}details{margin:8px 0}pre{white-space:pre-wrap}</style>',
        '<h1>Reference pocket and independent shape templates</h1>',
        '<p>Exploratory geometric hypotheses. Frequency is not energetic necessity; crystal self-controls are not independent activity validation.</p>',
        '<h2>Reference</h2><pre>'+e(json.dumps(r['cohort']['reference'],indent=2))+'</pre>',
        '<h2>Admission decisions</h2><table><tr><th>Ligand</th><th>Status</th><th>Evidence / reasons</th></tr>']
    for decision in r['cohort']['decisions']:
        rows.append('<tr><td>'+e(decision['query_id'])+'</td><td>'+e(decision['status'])+'</td><td><details><summary>Inspect alignment and chain checks</summary><pre>'+e(json.dumps(decision,indent=2))+'</pre></details></td></tr>')
    rows+=['</table>']
    if r.get('contact_evidence'):
        rows+=['<h2>Atom contact evidence</h2><p>The feature anchors below are a sparse representation, not the complete contact map. '
               'The ledger retains observed heavy-atom pairs within 4.5 A, all compatible polar proximity alternatives, '
               'source atom identities, aligned coordinates and uncertain atoms. Proximity is not a validated bond. '
               'Missing geometry is unknown.</p><a href="contacts.json">Open complete contact ledger</a><pre>'+
               e(json.dumps(r['contact_evidence'],indent=2))+'</pre>']
    rows+=['<h2>Pocket anchors</h2><table><tr><th>ID / contact</th><th>Support</th><th>Automatic mandatory eligible</th><th>Evidence</th></tr>']
    for a in r['anchors']:
        rows.append(f'<tr><td>{e(a["anchor_id"])}<br>{e(a["target_residue"])} {e(a["protein_atom"])} / {e(a["feature_class"])}</td>'
            f'<td>{a["distinct_structures"]}/{a["eligible_structures"]} eligible PDBs; {a["distinct_chemotypes"]} chemotypes</td><td>{e(a["mandatory_proposal_eligible"])}</td>'
            '<td><details><summary>Coordinates, sources and denominator</summary><pre>'+e(json.dumps(a,indent=2))+'</pre></details></td></tr>')
    rows+=['</table><h2>Independent template panel</h2><p>Proposed: '+e(', '.join(r['proposed_template_ids']))+'</p>',
        '<pre>'+e(json.dumps(r.get('template_selection',{}),indent=2))+'</pre><table><tr><th>Template</th><th>Resolution</th><th>Structure</th></tr>']
    for q in r['templates']:
        rows.append(f'<tr><td>{e(q["query_id"])}</td><td>{e(q["resolution"])}</td><td>{e(q["smiles"])}</td></tr>')
    rows+=['</table><h2>Limitations</h2><ul>',*['<li>'+e(x)+'</li>' for x in r['limitations']],'</ul></html>']
    path=output/'review.html';path.write_text('\n'.join(rows),encoding='utf-8')
    r.setdefault('outputs',{})['review.html']=str(path.resolve())
    from .gaussian_batch import _atomic_json
    _atomic_json(output/'report.json',r)
    return path
