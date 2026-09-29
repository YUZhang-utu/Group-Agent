"""Experimental PDB pocket grids, permissive clustering and sealed adoption."""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from .consensus_admission import (read_structure, ca_map, fit_protein, ligand_atoms,
                                  contact_residues, assembly_compatible)
from .expanded_wee1 import fingerprint
from .gaussian_batch import _atomic_json
from .screening_selection import check_hashes

DEFAULTS = dict(spacing=1.0, region_padding=4.0, probe_radius=.5,
                cluster_distance=.40, chemical_weight=.15, feature_radius=2.0,
                minimum_alignment_residues=12, maximum_core_rmsd=2.5,
                maximum_grid_points=150000, local_change_volume=30.0)
RADII = dict(C=1.70, N=1.55, O=1.52, S=1.80, P=1.80, SE=1.90)
SIDECHAINS = {
    'GLY':'', 'ALA':'CB', 'VAL':'CB CG1 CG2', 'LEU':'CB CG CD1 CD2',
    'ILE':'CB CG1 CG2 CD1', 'SER':'CB OG', 'THR':'CB OG1 CG2',
    'CYS':'CB SG', 'MET':'CB CG SD CE', 'PRO':'CB CG CD',
    'ASP':'CB CG OD1 OD2', 'ASN':'CB CG OD1 ND2', 'GLU':'CB CG CD OE1 OE2',
    'GLN':'CB CG CD OE1 NE2', 'LYS':'CB CG CD CE NZ',
    'ARG':'CB CG CD NE CZ NH1 NH2', 'HIS':'CB CG ND1 CD2 CE1 NE2',
    'PHE':'CB CG CD1 CD2 CE1 CE2 CZ', 'TYR':'CB CG CD1 CD2 CE1 CE2 CZ OH',
    'TRP':'CB CG CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2'}
DONORS = {'SER':{'OG'}, 'THR':{'OG1'}, 'TYR':{'OH'}, 'TRP':{'NE1'},
          'ASN':{'ND2'}, 'GLN':{'NE2'}, 'LYS':{'NZ'}, 'ARG':{'NE','NH1','NH2'}}
ACCEPTORS = {'SER':{'OG'}, 'THR':{'OG1'}, 'TYR':{'OH'}, 'ASN':{'OD1'},
             'GLN':{'OE1'}, 'ASP':{'OD1','OD2'}, 'GLU':{'OE1','OE2'}, 'MET':{'SD'}}
CHANNELS = ('donor', 'acceptor', 'hydrophobic', 'aromatic', 'positive', 'negative')


def policy(values=None):
    p = dict(DEFAULTS, **(values or {}))
    if set(p) != set(DEFAULTS) or any(type(v) not in (int, float) or not np.isfinite(v) for v in p.values()):
        raise ValueError('Invalid pocket policy')
    if not (.5 <= p['spacing'] <= 2 and 2 <= p['region_padding'] <= 8
            and 0 <= p['probe_radius'] <= 1.5 and .05 <= p['cluster_distance'] <= .8
            and 0 <= p['chemical_weight'] <= .3 and 1 <= p['feature_radius'] <= 3
            and 3 <= p['minimum_alignment_residues'] <= 100
            and .5 <= p['maximum_core_rmsd'] <= 4
            and 100 <= p['maximum_grid_points'] <= 500000
            and 0 < p['local_change_volume'] <= 200):
        raise ValueError('Pocket policy outside supported bounds')
    return p


def pocket_grid(seed, p):
    seed = np.asarray(seed, float)
    if seed.ndim != 2 or seed.shape[1] != 3 or not len(seed) or not np.isfinite(seed).all():
        raise ValueError('Invalid pocket seed coordinates')
    spacing = p['spacing']; padding = p['region_padding']
    lo = np.floor((seed.min(0)-padding)/spacing)*spacing
    hi = np.ceil((seed.max(0)+padding)/spacing)*spacing
    axes = [np.arange(a, b+spacing*.1, spacing) for a, b in zip(lo, hi)]
    if np.prod([len(a) for a in axes]) > p['maximum_grid_points'] * 8:
        raise ValueError('Pocket bounding grid too large; review binding site')
    grid = np.stack(np.meshgrid(*axes, indexing='ij'), -1).reshape(-1, 3)
    grid = grid[cKDTree(seed).query(grid)[0] <= padding]
    if len(grid) > p['maximum_grid_points']:
        raise ValueError('Pocket grid exceeds resource bound')
    return grid


def atom_features(atoms, with_elements=False):
    """Conservative residue templates; no inferred His/cysteine protonation."""
    features = {name: [] for name in CHANNELS}
    residues = {}
    for atom in atoms: residues.setdefault(atom['canonical_residue'], []).append(atom)
    for atom in atoms:
        residue = atom['auth_comp_id']; name = atom['auth_atom_id']; point = np.asarray(atom['xyz'])
        labels = []
        if (name == 'N' and residue != 'PRO') or name in DONORS.get(residue, set()): labels.append('donor')
        if name == 'O' or name in ACCEPTORS.get(residue, set()): labels.append('acceptor')
        polar_carbon={'ASP':'CG','ASN':'CG','GLU':'CD','GLN':'CD','ARG':'CZ'}.get(residue)
        if atom['type_symbol'] == 'C' and name not in {'C', 'CA', polar_carbon}: labels.append('hydrophobic')
        if residue in {'PHE','TYR','TRP','HIS'} and name not in {'N','CA','C','O','CB','OH'}: labels.append('aromatic')
        if (residue == 'LYS' and name == 'NZ') or (residue == 'ARG' and name in {'NE','NH1','NH2'}): labels.append('positive')
        if name in ACCEPTORS.get(residue, set()) and residue in {'ASP','GLU'}: labels.append('negative')
        neighbors = [np.asarray(a['xyz']) for a in residues[atom['canonical_residue']] if a is not atom
                     and .7 < np.linalg.norm(np.asarray(a['xyz'])-point) < 1.9]
        direction = point-np.mean(neighbors, axis=0) if neighbors else np.zeros(3)
        norm = np.linalg.norm(direction)
        if norm > 1e-8: direction /= norm
        for label in labels:
            # Outward heavy-atom proxy, not a reconstructed hydrogen/lone pair.
            center = point + (1.5*direction if label in {'donor','acceptor'} else 0)
            features[label].append((center, atom['type_symbol']) if with_elements else center)
    return features


def describe(grid, atoms, p):
    xyz = np.asarray([a['xyz'] for a in atoms], float)
    if not len(xyz) or not np.isfinite(xyz).all(): raise ValueError('Invalid receptor atoms')
    occupied = np.zeros(len(grid), bool)
    elements = np.array([a['type_symbol'] for a in atoms])
    for element in sorted(set(elements)):
        if element not in RADII: raise ValueError('Unsupported receptor element: '+element)
        occupied |= cKDTree(xyz[elements == element]).query(grid)[0] < RADII[element]+p['probe_radius']
    cavity = ~occupied
    fields = np.zeros((len(CHANNELS), len(grid)), bool)
    for index, (channel, features) in enumerate(atom_features(atoms, with_elements=True).items()):
        for element in sorted({element for _, element in features}):
            points = [point for point, symbol in features if symbol == element]
            # Nondirectional contact fields start at the accessible atomic surface.
            # An atom-centred 2 A ball was wholly buried inside the exclusion.
            radius = p['feature_radius']
            if channel not in {'donor', 'acceptor'}:
                radius += RADII[element] + p['probe_radius']
            fields[index] |= (cKDTree(points).query(grid)[0] <= radius) & cavity
    return cavity, fields


def iou(a, b):
    union = np.count_nonzero(a | b)
    return np.count_nonzero(a & b)/union if union else 1.0


def local_diagnostics(masks, grid, subpocket_centers, p):
    n = len(masks)
    raw = np.zeros((n,n)); legacy = raw.copy(); smooth = raw.copy(); changed = raw.copy()
    regions = [np.linalg.norm(grid-np.asarray(c), axis=1) <= 4 for c in subpocket_centers]
    for i in range(n):
        for j in range(i):
            values = [(1-iou(masks[i]&r, masks[j]&r),
                       np.count_nonzero((masks[i]^masks[j])&r)*p['spacing']**3) for r in regions]
            raw[i,j] = raw[j,i] = max((d for d,v in values), default=0)
            legacy[i,j] = legacy[j,i] = max((d for d,v in values if v >= p['local_change_volume']), default=0)
            smooth[i,j] = smooth[j,i] = max((d*v/(v+p['local_change_volume']) for d,v in values), default=0)
            changed[i,j] = changed[j,i] = max((v for d,v in values), default=0)
    return dict(raw=raw, legacy=legacy, smooth=smooth, maximum_changed_volume=changed)


def compare(masks, fields, grid, subpocket_centers, p):
    n = len(masks); overlap = np.eye(n); chemical = np.zeros((n,n)); local = np.zeros((n,n))
    local = local_diagnostics(masks, grid, subpocket_centers, p)['smooth']
    for i in range(n):
        for j in range(i):
            overlap[i,j] = overlap[j,i] = iou(masks[i], masks[j])
            # Compare chemistry only in mutually available space; geometry has its own term.
            common = masks[i] & masks[j]
            active = [k for k in range(len(CHANNELS)) if np.any((fields[i,k] | fields[j,k]) & common)]
            value = np.mean([1-iou(fields[i,k] & common, fields[j,k] & common) for k in active]) if active else 0
            chemical[i,j] = chemical[j,i] = value
    shape_distance = np.maximum(1-overlap, .5*local)
    distance = (1-p['chemical_weight'])*shape_distance + p['chemical_weight']*chemical
    np.fill_diagonal(distance, 0)
    return distance, overlap, chemical, local


def clusters(distance, rows, cutoff):
    if not rows: return []
    distance=np.asarray(distance,float)
    if distance.shape!=(len(rows),len(rows)) or not np.isfinite(distance).all() or not np.allclose(distance,distance.T):
        raise ValueError('Invalid pocket distance matrix')
    # Weighted average linkage: each PDB contributes one unit across its chains.
    # This corrects chain multiplicity, not crystallographic sampling bias.
    counts=Counter(r['pdb_id'] for r in rows)
    weights={i:1/counts[r['pdb_id']] for i,r in enumerate(rows)}
    members={i:[i] for i in range(len(rows))}; active=list(members)
    work=distance.copy(); np.fill_diagonal(work,np.inf)
    while len(active)>1:
        sub=work[np.ix_(active,active)]
        a,b=np.unravel_index(np.argmin(sub),sub.shape)
        if sub[a,b]>cutoff: break
        i,j=active[a],active[b]; total=weights[i]+weights[j]
        for k in active:
            if k not in {i,j}:
                work[i,k]=work[k,i]=(weights[i]*work[i,k]+weights[j]*work[j,k])/total
        members[i]+=members.pop(j);weights[i]=total;active.remove(j)
    groups=sorted([sorted(g) for g in members.values()],key=lambda g:min(rows[i]['id'] for i in g))
    result = []
    for group in groups:
        counts = Counter(rows[i]['pdb_id'] for i in group)
        weights = np.array([1/counts[rows[i]['pdb_id']] for i in group])
        costs = distance[np.ix_(group, group)] @ weights / weights.sum()
        minimum = float(costs.min())
        eligible = [i for i,cost in zip(group,costs) if cost <= minimum+.02]
        representative = min(eligible, key=lambda i: (not bool(rows[i].get('queries')), rows[i].get('resolution') or 99, rows[i]['id']))
        members = [rows[i]['id'] for i in group]
        result.append(dict(id='pocket-'+hashlib.sha256(json.dumps(members).encode()).hexdigest()[:12],
                           members=members, representative=rows[representative]['id'],
                           distinct_pdb_support=len(counts), small_support=len(counts)<3,
                           representative_max_distance=float(distance[representative,group].max()),
                           representative_weighted_distance=float(costs[group.index(representative)])))
    return result


def quality(atoms, residues, relevant=None):
    issues = []
    for residue in sorted(residues):
        selected = [a for a in atoms if a['canonical_residue'] == residue]
        if not selected: issues.append(f'missing_pocket_residue:{residue}'); continue
        name = selected[0]['auth_comp_id']
        if name not in SIDECHAINS: issues.append(f'unsupported_pocket_residue:{residue}:{name}'); continue
        observed = {a['auth_atom_id'] for a in selected}
        required = set(SIDECHAINS[name].split()) | {'N','CA','C','O'}
        if relevant is not None: required &= relevant.get(residue,set())
        missing = required - observed
        if missing: issues.append(f'missing_pocket_atoms:{residue}:'+','.join(sorted(missing)))
        if any((a['label_alt_id'] or float(a['occupancy']) < .9) and a['auth_atom_id'] in required for a in selected):
            issues.append(f'uncertain_pocket_occupancy:{residue}')
    return issues


def build(source, output, reference_query=None, target_chain=None, settings=None, excluded_entries=None, advice_settings=None):
    source = Path(source).resolve(); root = source.parent; out = Path(output).resolve()
    if (out/'report.json').exists(): raise FileExistsError('Use a fresh pocket output')
    out.mkdir(parents=True, exist_ok=True); p = policy(settings)
    excluded_entries = {str(k).upper(): str(v) for k,v in (excluded_entries or {}).items()}
    survey = json.loads(source.read_text()); protein = json.loads((root/'protein.json').read_text())
    instances = json.loads((root/'ligand-instances.json').read_text()); cache = {}
    def get(code):
        if code not in cache: cache[code] = read_structure(root/'structures'/f'{code}.cif', protein)
        return cache[code]
    options = []
    for row in instances:
        if row['pdb_id'] != survey['reference_site_pdb'] or not row['quality_passed']: continue
        if reference_query and row['query_id'] != reference_query: continue
        structure = get(row['pdb_id']); lig = ligand_atoms(structure, row)
        for chain, atoms in structure['chains'].items():
            if target_chain and chain != target_chain: continue
            if len(contact_residues(atoms, lig)) >= 3 and assembly_compatible(structure, chain, lig):
                options.append((row, chain))
    report = dict(kind='pocket_states', status='complete', descriptor_version=2,
                  excluded_entries=excluded_entries, policy=p, target=protein,
                  diversity_report=str(source), sources=fingerprint([source,root/'protein.json',root/'ligand-instances.json']),
                  limitations=['PDB support is not equilibrium occupancy; rare clusters are retained',
                    'Reference-ligand-expanded local region, not exhaustive cavity or access-path detection',
                    'Residue-template chemical fields with outward heavy-atom direction proxies; no affinity claim',
                    'His/Cys protonation, termini, waters, metals and transformed interfaces require separate review',
                    'Only model 1 and identity assembly 1; uncertain pocket atoms are held for review',
                    'Thresholds are exploratory defaults, not calibrated target-independent cutoffs'])
    if len(options) != 1:
        report.update(readiness='needs_reference_instance', reference_options=[dict(query_id=r['query_id'],target_chain=c) for r,c in options])
        _atomic_json(out/'report.json', report); return report
    reference, chain = options[0]; structure = get(reference['pdb_id'])
    ligand = ligand_atoms(structure, reference); seed = np.array([a['xyz'] for a in ligand])
    ra = structure['chains'][chain]; rm = ca_map(ra)
    pocket_residues = contact_residues(ra, ligand, 6.)
    grid = pocket_grid(seed, p); grid_tree=cKDTree(grid)
    relevant={r:{a['auth_atom_id'] for a in ra if a['canonical_residue']==r and
                 grid_tree.query(a['xyz'])[0] <= RADII.get(a['type_symbol'],1.8)+p['probe_radius']+.5}
              for r in pocket_residues}
    problems = quality(ra, pocket_residues, relevant)
    if problems:
        report.update(readiness='needs_reference_quality', reasons=problems)
        _atomic_json(out/'report.json', report); return report
    rows = []; masks = []; fields = []; held = []
    reference_id = reference['pdb_id']+':'+chain
    for file in sorted((root/'structures').glob('*.cif')):
        code = file.stem.upper()
        report['sources'].update(fingerprint([file]))
        if code in excluded_entries:
            held.append(dict(id=code, scope='entire_entry', reasons=['explicit_entry_exclusion'],
                             detail=excluded_entries[code])); continue
        try: s = get(code)
        except (ValueError,KeyError,OSError) as exc:
            held.append(dict(id=code,reasons=['structure_read_failed'],detail=str(exc)[:300]));continue
        if not s['chains']:
            held.append(dict(id=code,reasons=['target_identity_or_organism_unverified']));continue
        methods = s['data'].get('_exptl.method', [])
        if isinstance(methods, str): methods = [methods]
        if not methods:
            held.append(dict(id=code, reasons=['missing_experimental_method'])); continue
        for c, atoms in sorted(s['chains'].items()):
            sid = code+':'+c; reasons = quality(atoms, pocket_residues, relevant)
            names={a['canonical_residue']:a['auth_comp_id'] for a in atoms}
            reference_names={a['canonical_residue']:a['auth_comp_id'] for a in ra}
            if any(names.get(r)!=reference_names[r] for r in pocket_residues):
                reasons.append('pocket_sequence_difference_requires_review')
            if not assembly_compatible(s, c, []): reasons.append('assembly_context_unresolved')
            cm = ca_map(atoms); common = sorted(set(cm)&set(rm))
            core = [r for r in common if np.min(np.linalg.norm(seed-rm[r],axis=1)) <= 20]
            if len(core) < p['minimum_alignment_residues']: reasons.append('insufficient_alignment')
            if reasons: held.append(dict(id=sid,reasons=reasons)); continue
            # Trim displaced residues while keeping a recorded stable local core.
            for _ in range(3):
                matrix, rmsd, residuals = fit_protein([cm[r] for r in core],[rm[r] for r in core])
                keep = [r for r,d in zip(core,residuals) if d <= max(2., float(np.quantile(residuals,.8)))]
                if len(keep)<p['minimum_alignment_residues'] or len(keep)==len(core): break
                core = keep
            matrix, rmsd, _ = fit_protein([cm[r] for r in core],[rm[r] for r in core])
            if rmsd > p['maximum_core_rmsd']:
                held.append(dict(id=sid,reasons=['alignment_unresolved'],core_rmsd=rmsd)); continue
            aligned = [dict(a,xyz=np.asarray(a['xyz'])@matrix[:3,:3].T+matrix[:3,3]) for a in atoms]
            near = cKDTree(grid).query(np.array([a['xyz'] for a in aligned]))[0] < 8
            aligned = [a for a, keep in zip(aligned, near) if keep]
            try: mask, feature = describe(grid,aligned,p)
            except ValueError as exc:
                held.append(dict(id=sid,reasons=['pocket_description_failed'],detail=str(exc)));continue
            if not mask.any(): held.append(dict(id=sid,reasons=['empty_local_cavity'])); continue
            queries = []
            for q in instances:
                if q['pdb_id'] != code or not q['quality_passed']: continue
                lig = ligand_atoms(s,q)
                if len(contact_residues(atoms,lig)&pocket_residues) >= 3 and assembly_compatible(s,c,lig): queries.append(q)
            rows.append(dict(id=sid,pdb_id=code,target_chain=c,transform=matrix.tolist(),core_rmsd=rmsd,
                             structure_path=str(file.resolve()),
                             alignment_residues=core,queries=queries,resolution=min([q['resolution'] for q in queries],default=None),
                             feature_source_counts={k:len(v) for k,v in atom_features(aligned).items()},
                             field_volumes={k:float(feature[i].sum()*p['spacing']**3) for i,k in enumerate(CHANNELS)},
                             cavity_volume=float(mask.sum()*p['spacing']**3)))
            masks.append(mask); fields.append(feature)
            print(f'Pocket states: prepared {sid}; {len(rows)} structures', flush=True)
    if not rows:
        report.update(readiness='needs_usable_pockets',held_for_review=held)
        _atomic_json(out/'report.json',report); return report
    masks=np.array(masks); fields=np.array(fields)
    distance, overlap, chemical, local = compare(masks,fields,grid,[rm[r] for r in sorted(pocket_residues) if r in rm],p)
    diagnostics = local_diagnostics(masks,grid,[rm[r] for r in sorted(pocket_residues) if r in rm],p)
    groups = clusters(distance,rows,p['cluster_distance'])
    for group in groups:
        members=[i for i,r in enumerate(rows) if r['id'] in group['members']]
        group['maximum_chemical_difference']=float(chemical[np.ix_(members,members)].max())
        group['review_flags']=[]
        if group['representative_max_distance']>p['cluster_distance']:
            group['review_flags'].append('average_linkage_member_outside_representative_cutoff')
        if group['maximum_chemical_difference']>.6:
            group['review_flags'].append('inspect_within_cluster_chemical_differences')
    artifacts=out/'pocket-grids.npz'
    np.savez_compressed(artifacts,grid=grid,masks=masks,chemical_fields=fields,
                        distance=distance,overlap=overlap,chemical_distance=chemical,local_difference=local,
                        **{'local_'+k:v for k,v in diagnostics.items()})
    report.update(readiness='needs_user_adoption',reference=dict(query_id=reference['query_id'],target_chain=chain,id=reference_id),
                  pocket_residues=sorted(pocket_residues),quality_atom_scope={str(k):sorted(v) for k,v in relevant.items()},
                  structures=rows,clusters=groups,held_for_review=held,
                  artifact=str(artifacts),grid_points=len(grid),channels=CHANNELS,
                  outputs=dict(grids=str(artifacts),html=str(out/'report.html')),
                  sensitivity=[dict(distance=float(t),clusters=len(clusters(distance,rows,float(t))))
                               for t in sorted({max(.05,p['cluster_distance']-.1),p['cluster_distance'],min(.8,p['cluster_distance']+.1)})],
                  warnings=['Many clusters: inspect alignment, pocket coverage and sensitivity before changing thresholds'] if len(groups)>6 else [],
                  pairwise=dict(ids=[r['id'] for r in rows],overlap=overlap.tolist(),distance=distance.tolist(),
                                chemical_distance=chemical.tolist(),local_difference=local.tolist()),
                  code_hashes=fingerprint([Path(__file__),Path(__file__).with_name('consensus_admission.py')]))
    report['sources'].update(fingerprint([artifacts]))
    from .receptor_advice import attach, summary
    attach(report, advice_settings)
    report['code_hashes'].update(fingerprint([Path(__file__).with_name('receptor_advice.py')]))
    # Keep the saved matrix consistent with the newly versioned recommendation.
    with np.load(artifacts) as archive:
        arrays = {k: archive[k] for k in archive.files}
    arrays['legacy_distance'] = arrays['distance']
    arrays['distance'] = np.asarray(report['pairwise']['distance'])
    np.savez_compressed(artifacts, **arrays)
    report['sources'].update(fingerprint([artifacts]))
    advice_path = out/'receptor-recommendation.txt'
    advice_path.write_text(summary(report['receptor_advice']), encoding='utf-8')
    report['outputs']['recommendation'] = str(advice_path)
    _atomic_json(out/'report.json',report)
    render(report,out/'report.html',grid,masks)
    return report


def render(report, path, grid, masks):
    import html
    data=json.dumps(dict(ids=[r['id'] for r in report['structures']],grid=np.round(grid,2).tolist(),
                         masks=[np.flatnonzero(m).tolist() for m in masks])).replace('<','\\u003c')
    rows=''.join('<tr>'+''.join('<td>'+html.escape(str(c[k]))+'</td>' for k in
                             ('id','representative','distinct_pdb_support','small_support'))+'</tr>' for c in report['clusters'])
    page='''<!doctype html><meta charset="utf-8"><title>PDB pocket states</title>
<style>body{font:16px system-ui;margin:32px;max-width:1100px}td,th{padding:8px;text-align:left}canvas{border:1px solid #ccc}</style>
<h1>PDB pocket states</h1><p>Experimental support is not equilibrium occupancy. Review small clusters before adoption.</p>
<table><tr><th>State</th><th>Representative</th><th>Distinct PDBs</th><th>Small support</th></tr>ROWS</table>
<p>Compare pocket space: <select id="a"></select> <select id="b"></select> Projection <select id="axis"><option value="1">XY</option><option value="2">XZ</option></select></p>
<p>Gray: shared free space. Blue: first only. Orange: second only. Projection may hide depth; use the saved 3D grids for detailed inspection.</p>
<canvas id="plot" width="850" height="500"></canvas><pre id="summary"></pre><script>
const d=DATA; const a=document.getElementById('a'),b=document.getElementById('b'),axis=document.getElementById('axis');
for(const s of [a,b])d.ids.forEach((id,i)=>{let o=document.createElement('option');o.value=i;o.textContent=id;s.append(o)});
b.value=Math.min(1,d.ids.length-1);function draw(){let ctx=document.getElementById('plot').getContext('2d');ctx.clearRect(0,0,850,500);
const A=new Set(d.masks[+a.value]),B=new Set(d.masks[+b.value]),k=+axis.value;
let xs=d.grid.map(p=>p[0]),ys=d.grid.map(p=>p[k]);let x=xs.reduce((a,b)=>Math.min(a,b)),y=ys.reduce((a,b)=>Math.min(a,b)),scale=Math.min(790/(xs.reduce((a,b)=>Math.max(a,b))-x||1),440/(ys.reduce((a,b)=>Math.max(a,b))-y||1));
for(let i=0;i<d.grid.length;i++){if(!A.has(i)&&!B.has(i))continue;ctx.fillStyle=A.has(i)&&B.has(i)?'#b7c0c5':A.has(i)?'#1676b8':'#db843b';let p=d.grid[i];ctx.fillRect(30+(p[0]-x)*scale,470-(p[k]-y)*scale,3,3)}}
a.onchange=b.onchange=axis.onchange=draw;draw();</script>'''
    if report.get('receptor_advice'):
        from .receptor_advice import summary
        header = '<h1>Receptor selection recommendation</h1><pre style="white-space:pre-wrap">'+html.escape(summary(report['receptor_advice']))+'</pre>'
        start = page.index('<h1>'); end = page.index('<table>')
        page = page[:start]+header+'<details><summary>Detailed pocket evidence and selection IDs</summary>'+page[end:]
        page = page.replace('<th>State</th>', '<th>Selection ID</th>')+'</details>'
    Path(path).write_text(page.replace('ROWS',rows).replace('DATA',data),encoding='utf-8')


def adopt(source, output, cluster_ids=None):
    source=Path(source).resolve(); report=json.loads(source.read_text()); check_hashes(report['sources'])
    if report.get('descriptor_version', 1) < 2:
        raise ValueError('Recompute pocket fields: descriptor v1 has buried chemical channels and censored local differences')
    if report.get('readiness')!='needs_user_adoption': raise ValueError('Review a completed pocket-state report first')
    allowed={c['id'] for c in report['clusters']+report.get('receptor_options',[])}
    selected=list(cluster_ids) if cluster_ids is not None else [c['id'] for c in report['clusters']]
    if not selected or len(set(selected))!=len(selected) or not set(selected)<=allowed: raise ValueError('Invalid pocket-state selection')
    output=Path(output); output.mkdir(parents=True,exist_ok=True)
    if (output/'report.json').exists(): raise FileExistsError('Use a fresh adoption output')
    result=dict(kind='pocket_adoption',status='complete',readiness='ready_for_state_consensus',
                pocket_report=str(source),selected_cluster_ids=selected,
                sources={**report['sources'],**fingerprint([source])},
                authorization='Explicit adoption invocation; never inferred from report generation')
    _atomic_json(output/'report.json',result); return result


def state_cohort(adoption, state_id=None):
    adopted=json.loads(Path(adoption).read_text()); check_hashes(adopted['sources'])
    if adopted.get('kind')!='pocket_adoption' or adopted.get('readiness')!='ready_for_state_consensus':
        raise ValueError('Explicit pocket adoption required before consensus')
    selected=adopted['selected_cluster_ids']
    if state_id is None and len(selected)==1: state_id=selected[0]
    if state_id not in selected: raise ValueError('Choose one adopted pocket_state_id: '+', '.join(selected))
    report=json.loads(Path(adopted['pocket_report']).read_text())
    group=next(c for c in report['clusters']+report.get('receptor_options',[]) if c['id']==state_id)
    rows={r['id']:r for r in report['structures']}; rep=rows[group['representative']]
    if not rep['queries']: raise ValueError('Representative has no prepared ligand reference; review before consensus')
    reference=rep['queries'][0]; inverse=np.linalg.inv(np.array(rep['transform'])); admitted=[]
    seen=set()
    for sid in group['members']:
        row=rows[sid]
        for query in row['queries']:
            if query['query_id'] in seen: continue
            seen.add(query['query_id'])
            admitted.append(dict(query,admission=dict(target_chain=row['target_chain'],
                transform=(inverse@np.array(row['transform'])).tolist(),pocket_state_id=state_id)))
    cohort=dict(status='complete',readiness='cohort_ready',reference=dict(query_id=reference['query_id'],
                target_chain=rep['target_chain'],assembly_id='1',coordinate_frame=reference['query_id']),
                admitted=admitted,decisions=[],policy=report['policy'],pocket_state_id=state_id,
                receptor_selection_role=group.get('role','legacy_pocket_partition'),
                limitations=report['limitations']+['Consensus restricted to an explicitly adopted local cohort; a coverage neighborhood is not a physical state'])
    return report['diversity_report'],cohort,adopted


def main():
    parser=argparse.ArgumentParser(description=__doc__); sub=parser.add_subparsers(dest='command',required=True)
    build_parser=sub.add_parser('build'); build_parser.add_argument('--source',type=Path,required=True)
    build_parser.add_argument('--output',type=Path,required=True); build_parser.add_argument('--reference-query')
    build_parser.add_argument('--target-chain'); build_parser.add_argument('--cluster-distance',type=float,default=.40)
    build_parser.add_argument('--exclude-entry', action='append', default=[], help='Exclude every chain of this PDB entry')
    build_parser.add_argument('--maximum-representatives', type=int, default=32)
    build_parser.add_argument('--coverage-fraction', type=float, default=.95)
    adoption=sub.add_parser('adopt'); adoption.add_argument('--source',type=Path,required=True)
    adoption.add_argument('--output',type=Path,required=True); adoption.add_argument('--cluster-ids',nargs='+')
    consensus=sub.add_parser('consensus'); consensus.add_argument('--source',type=Path,required=True)
    consensus.add_argument('--output',type=Path,required=True); consensus.add_argument('--pocket-state-id')
    args=parser.parse_args()
    if args.command=='build':
        result=build(args.source,args.output,args.reference_query,args.target_chain,dict(cluster_distance=args.cluster_distance),
                     {code:'Explicit CLI entry exclusion' for code in args.exclude_entry},
                     dict(maximum_representatives=args.maximum_representatives,coverage_fraction=args.coverage_fraction))
    elif args.command=='adopt': result=adopt(args.source,args.output,args.cluster_ids)
    else:
        from .guided_workflow import execute
        result=execute('pocket_consensus',args.source,args.output,dict(pocket_state_id=args.pocket_state_id),{},False)
    print(json.dumps(result,indent=2))


if __name__=='__main__': main()
