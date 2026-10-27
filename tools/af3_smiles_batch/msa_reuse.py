"""Reuse inline protein MSA/templates without changing ligands or protein identity."""
import copy
import json
from pathlib import Path

FIELDS=('unpairedMsa','pairedMsa','templates')


def proteins(data):
    found={}
    for entity in data['sequences']:
        if 'protein' not in entity:continue
        p=entity['protein'];ids=p['id'] if isinstance(p['id'],list) else [p['id']]
        for chain in ids:
            if chain in found:raise ValueError('Duplicate protein chain')
            found[chain]=p
    if not found:raise ValueError('No protein chains')
    return found


def identity(p):
    return (p['sequence'],json.dumps(p.get('modifications',[]),sort_keys=True))


def ready(data):
    return all(all(k in p and p[k] is not None for k in FIELDS) for p in proteins(data).values())


def transfer(job,prepared):
    source=proteins(prepared);target=proteins(job)
    if set(source)!=set(target):raise ValueError('Prepared protein chain IDs differ')
    for chain,p in target.items():
        if identity(p)!=identity(source[chain]):raise ValueError('Prepared protein sequence/modifications differ')
    if not ready(prepared):raise ValueError('Prepared input lacks inline MSA/templates; cannot skip data pipeline')
    result=copy.deepcopy(job)
    for entity in result['sequences']:
        if 'protein' not in entity:continue
        p=entity['protein'];ids=p['id'] if isinstance(p['id'],list) else [p['id']]
        values=[{k:source[c][k] for k in FIELDS} for c in ids]
        if any(v!=values[0] for v in values[1:]):
            raise ValueError('Grouped protein chains have different prepared features; use separate protein entries')
        p.update(copy.deepcopy(values[0]))
    return result


def load_data_output(output):
    paths=[] if output is None else list(Path(output).rglob('*_data.json'))
    if len(paths)!=1:raise ValueError('Expected exactly one AF3 data-pipeline JSON')
    return json.loads(paths[0].read_text(encoding='utf-8')),paths[0]
