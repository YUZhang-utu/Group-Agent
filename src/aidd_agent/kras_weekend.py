"""E116: native KRAS reference preparation and finite-pose screening primitives."""
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from scipy.spatial import cKDTree
from .chemistry_prep import _mmcif_dict, _column, enumerate_ligand_instances, _ccd_molecule
from .gaussian_batch import prepare_gaussian_query, _load_query, prepare_seeds, write_gaussian_query
from .full_library_screen import POSE_PARAMETERS
from .interaction_fast import match_batch
from .gaussian_overlay import score_overlay
from .anchor_extraction import protein_hbond_roles
from .classified_features import extra_hypotheses

POLICY = dict(anchor=.5, residues=3, polar=1, shape=.5, color=.3,
              minimum_distance=2., envelope=.65, region=.5, distinct_features=3,
              acidic_polar_residues=[62,69], pocket_wall_residues=[95,96])
POLICY.update(shared_regions=[95,69,62,72,68], polar_distance=3.5,
              hydrophobic_distance=4.5, g12d_polar_residues=[12,68],
              g12v_residue12='steric_compatibility', version='shared-five-regions-v1')
REFERENCES = {'KRAS_G12D': ('7RPZ', '6IC', 'ASP'), 'KRAS_G12V': ('12FF', 'A1DB0', 'VAL')}


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def save(path, value):
    path = Path(path)
    temporary = path.with_suffix(path.suffix + '.partial')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n', encoding='utf-8')
    temporary.replace(path)


def regions(reference, tree):
    contact = reference[tree.query(reference)[0] <= 4.5]
    if len(contact) < 9:
        raise ValueError('Too few reference contact atoms for three spatial regions')
    seeds = [int(np.argmax(np.sum((contact-contact.mean(axis=0))**2, axis=1)))]
    for _ in range(2):
        nearest = np.min(np.sum((contact[:, None]-contact[seeds][None, :])**2, axis=2), axis=1)
        if nearest.max() < 4:
            raise ValueError('Reference regions cannot be separated')
        seeds.append(int(np.argmax(nearest)))
    labels = np.argmin(np.sum((contact[:, None]-contact[seeds][None, :])**2, axis=2), axis=1)
    groups = [contact[labels == i] for i in range(3)]
    if any(len(g) < 3 for g in groups):
        raise ValueError('Unsupported reference contact region')
    return contact, groups


def atom_rows(path):
    data = _mmcif_dict(path)
    keys = ['group_PDB', 'auth_comp_id', 'auth_asym_id', 'auth_seq_id', 'auth_atom_id',
            'type_symbol', 'label_alt_id', 'occupancy', 'pdbx_PDB_ins_code', 'pdbx_PDB_model_num',
            'Cartn_x', 'Cartn_y', 'Cartn_z']
    n = len(data['_atom_site.group_PDB'])
    columns = {k: _column(data, '_atom_site.'+k, n) for k in keys}
    rows = []
    for i in range(n):
        a = {k: columns[k][i] for k in keys}
        if a['pdbx_PDB_model_num'] != '1' or a['type_symbol'].upper() in ('H', 'D'):
            continue
        a['xyz'] = [float(a['Cartn_'+axis]) for axis in 'xyz']
        rows.append(a)
    return rows


def prepare(sources, output):
    from Bio.SeqUtils import seq1
    sources, output = Path(sources), Path(output)
    output.mkdir(parents=True, exist_ok=False)
    canonical = ''.join(line.strip() for line in (sources/'P01116-2.fasta').read_text().splitlines()
                        if not line.startswith('>'))
    records = {}
    for target, (pdb, ccd, mutation) in REFERENCES.items():
        cif, topology = sources/(pdb+'.cif'), sources/(ccd+'.cif')
        instances = [r for r in enumerate_ligand_instances(cif, [ccd]) if r['model']=='1'
                     and not r['altloc'] and not r['insertion_code']]
        if not instances:
            raise ValueError('No complete unambiguous ligand instance: '+pdb)
        instance = sorted(instances, key=lambda r:(r['chain_id'],int(r['residue_number'])))[0]
        chain = instance['chain_id']
        query_id = f"{pdb}:{ccd}:{chain}:{instance['residue_number']}"
        folder = output/target; folder.mkdir()
        manifest_path = folder/'input.json'
        save(manifest_path, dict(query_id=query_id, anchors=[], weights={'ordinary':1.0}))
        prepare_gaussian_query(cif, topology, manifest_path, folder/'query.npz')
        _, query = _load_query(folder/'query.npz')
        mol = _ccd_molecule(topology, instance['atoms'])
        from .conformer_artifacts import _feature_factory, FEATURE_TYPES
        heavy = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum()>1]
        heavy_map = {original:i for i,original in enumerate(heavy)}
        members = [[heavy_map[j] for j in f.GetAtomIds() if j in heavy_map]
                   for f in _feature_factory().GetFeaturesForMol(mol) if f.GetFamily() in FEATURE_TYPES]
        chemical = dict(atomic_numbers=[mol.GetAtomWithIdx(i).GetAtomicNum() for i in heavy],
                        feature_members=members)
        ccd_data = _mmcif_dict(topology)
        expected = sum(e.upper() not in ('H','D') for e in ccd_data['_chem_comp_atom.type_symbol'])
        if mol.GetNumHeavyAtoms() != expected:
            raise ValueError('Incomplete native ligand: '+query_id)
        all_atoms = atom_rows(cif)
        # Use a consistent major A conformer where present; record that policy.
        chosen = {}
        for a in all_atoms:
            if a['auth_asym_id'] != chain or a['pdbx_PDB_ins_code'] or a['label_alt_id'] not in ('', 'A'):
                continue
            if float(a['occupancy']) <= 0:
                continue
            if a['group_PDB'] != 'ATOM' and a['auth_comp_id'] not in ('GDP','MG'):
                continue
            key = (a['auth_seq_id'], a['auth_comp_id'], a['auth_atom_id'])
            if key not in chosen or float(a['occupancy']) > float(chosen[key]['occupancy']):
                chosen[key] = a
        atoms = list(chosen.values())
        protein = [a for a in atoms if a['group_PDB']=='ATOM']
        resnames = {int(a['auth_seq_id']):a['auth_comp_id'] for a in protein}
        if resnames.get(12) != mutation:
            raise ValueError('Wrong mutation or residue numbering: '+pdb)
        observed_mismatches = [dict(position=i, canonical=canonical[i-1], observed=seq1(res))
                               for i,res in sorted(resnames.items()) if 1<=i<=len(canonical)
                               and seq1(res)!=canonical[i-1]]
        environment = []
        for a in protein:
            donor, acceptor = protein_hbond_roles(a['auth_comp_id'],a['auth_atom_id'])
            environment.append(dict(residue=a['auth_comp_id'],chain=chain,
                residue_number=a['auth_seq_id'],atom_name=a['auth_atom_id'],element=a['type_symbol'],
                donor=donor,acceptor=acceptor,xyz=np.array(a['xyz'])))
        anchors = []
        for i,(point,kind) in enumerate(zip(query['feature_points'], query['feature_types'])):
            if kind not in (1,2):
                continue
            for a in environment:
                distance = float(np.linalg.norm(point-a['xyz']))
                if distance<=3.5 and (a['acceptor'] if kind==1 else a['donor']):
                    anchors.append(dict(feature_index=i,target_residue=int(a['residue_number']),
                        protein_atom=a['atom_name'],feature_class='HBD' if kind==1 else 'HBA',
                        distance=distance))
        for a in extra_hypotheses(query_id,query,environment):
            partner = a['evidence']['protein_partner']
            if partner.get('chain') != chain:
                continue
            anchors.append(dict(feature_index=int(a['feature_index']),
                target_residue=int(partner['residue_number']),protein_atom=partner['atom_name'],
                feature_class=a['feature_class'],
                distance=float(a['evidence']['distance_angstrom'])))
        anchors = [a for a in anchors if a['feature_class'] not in ('HBD','HBA') or
                   any(chemical['atomic_numbers'][i] in (7,8) for i in members[a['feature_index']])]
        indices = np.array(sorted({a['feature_index'] for a in anchors}),dtype=int)
        if not len(indices):
            raise ValueError('No native contacts: '+query_id)
        query['anchor_feature_indices'] = indices
        write_gaussian_query(folder/'query.npz', **query, source=dict(query_id=query_id,
            mmcif_sha256=sha(cif),ccd_sha256=sha(topology)))
        xyz = np.array([a['xyz'] for a in atoms])
        tree = cKDTree(xyz)
        contact, groups = regions(query['shape_points'],tree)
        np.savez_compressed(folder/'geometry.npz', receptor=xyz, contact=contact,
                            **{f'region{i}':g for i,g in enumerate(groups)})
        record = dict(target=target,query_id=query_id,anchors=anchors,
            chemistry=chemical,
            receptor_atoms=[dict(residue=int(a['auth_seq_id']),resname=a['auth_comp_id'],
                atom=a['auth_atom_id'],element=a['type_symbol'],xyz=a['xyz']) for a in protein],
            observed_sequence_differences=observed_mismatches,
            receptor_scope='Native target chain plus same-chain GDP/Mg; no waters/additives',
            alternate_policy='Blank/A atom locations, highest occupancy per identity; model 1',
            nucleotide='GDP', source_hashes={p.name:sha(p) for p in (cif,topology,sources/'P01116-2.fasta')},
            policy=POLICY, limitations=['Native engineered construct, not a repaired canonical mutant.',
                'Contact transfer is geometric; hydrogen-bond angles and energetics are not validated.',
                'One chemotype/state per allele; no allele selectivity inference.'])
        save(folder/'reference.json',record)
        records[target]=record
    save(output/'report.json',dict(status='prepared',targets=records,
         output_hashes={str(p.relative_to(output)):sha(p) for p in output.rglob('*') if p.is_file()}))
    return records


def load_state(root, target):
    folder = Path(root)/target
    _, query = _load_query(folder/'query.npz')
    with np.load(folder/'geometry.npz') as z:
        geometry = {k:z[k].copy() for k in z.files}
    record=read(folder/'reference.json')
    if record['target']!=target or record['policy']!=POLICY:
        raise ValueError('Reference target/policy mismatch')
    return dict(query=query, record=record, geometry=geometry,
                tree=cKDTree(geometry['receptor']))


def contact_mask(state, assignments, scores):
    indices = state['query']['anchor_feature_indices']
    lookup = {int(index):i for i,index in enumerate(indices)}
    masks = []
    for assigned, row in zip(assignments,scores):
        matched = [a for a in state['record']['anchors'] if assigned[lookup[a['feature_index']]]>=0
                   and row[lookup[a['feature_index']]]>=POLICY['anchor']]
        residue_ids = {a['target_residue'] for a in matched}
        required=[(95,'HBA'),(69,'HBD'),(62,'HBD')]
        if state['record']['target']=='KRAS_G12D': required += [(12,'HBD'),(68,'HBA')]
        modes=all(any(a['target_residue']==res and a['feature_class']==kind for a in matched)
                  for res,kind in required)
        masks.append(len(residue_ids)>=POLICY['residues'] and modes and
                     len({a['feature_index'] for a in matched})>=POLICY['distinct_features'])
    return np.array(masks,dtype=bool)


def chemical_sets(chemistry, types):
    numbers=np.asarray(chemistry.atomic_numbers)
    if len(types)!=len(chemistry.feature_members): raise ValueError('Feature/member count mismatch')
    def members(kinds,elements):
        return np.array(sorted({int(j) for i,kind in enumerate(types) if int(kind) in kinds
            for j in chemistry.feature_members[i] if int(numbers[j]) in elements}),dtype=int)
    return dict(donor=members({1},{7,8}),acceptor=members({2},{7,8}),carbon=members({5,6},{6}))


def direct_contact_masks(points, sets, state, matrices):
    """Same-pose chemically typed proximity, not an H-bond energy/angle claim."""
    atoms=state['record']['receptor_atoms']
    moved=np.einsum('bij,nj->bni',matrices[:,:3,:3],points)+matrices[:,None,:3,3]
    def proximity(indices,residue,names,distance):
        partners=np.array([a['xyz'] for a in atoms if a['residue']==residue and a['atom'] in names])
        if not len(indices) or not len(partners): return np.zeros(len(matrices),dtype=bool)
        squared=np.sum((moved[:,indices,None,:]-partners[None,None,:,:])**2,axis=-1)
        return squared.min(axis=(1,2))<=distance**2
    polar=POLICY['polar_distance']; hydro=POLICY['hydrophobic_distance']
    masks={
        'His95':proximity(sets['acceptor'],95,{'ND1','NE2'},polar),
        'Asp69':proximity(sets['donor'],69,{'OD1','OD2'},polar),
        'Glu62':proximity(sets['donor'],62,{'OE1','OE2'},polar),
        'Met72':proximity(sets['carbon'],72,{'CB','CG','CE'},hydro),
        'Arg68_polar':proximity(sets['acceptor'],68,{'NE','NH1','NH2'},polar),
        'Arg68_carbon':proximity(sets['carbon'],68,{'CB','CG','CD'},hydro)}
    masks['Arg68']=masks['Arg68_polar'] | masks['Arg68_carbon']
    if state['record']['target']=='KRAS_G12D':
        masks['residue12']=proximity(sets['donor'],12,{'OD1','OD2'},polar) & masks['Arg68_polar']
    else:
        # No unobserved attractive Val12 interaction is imposed.
        clashes=proximity(np.arange(len(points)),12,{'N','CA','C','O','CB','CG1','CG2'},
                          POLICY['minimum_distance']-1e-9)
        masks['residue12']=~clashes
    masks['all']=np.logical_and.reduce([masks[k] for k in
                                     ('His95','Asp69','Glu62','Met72','Arg68','residue12')])
    return masks


def evaluate(candidate, features, state, chemistry, sets=None):
    query=state['query']; indices=query['anchor_feature_indices']
    sets=chemical_sets(chemistry,candidate.feature_types) if sets is None else sets
    flags=dict(chemical=all(len(v)>0 for v in sets.values()),contact=False,shape=False,color=False,geometry=False,full=False)
    def finish(n,best=None,passing=None,region_counts=None):
        return dict(gid=int(getattr(candidate,'global_id',-1)),mid=candidate.molecule_id,cid=candidate.conformer_id,
                    seeds=n,flags=flags,best_contact_shape=best,passing_pose=passing,
                    direct_region_pose_counts=region_counts or {})
    if not flags['chemical']: return finish(0)
    seeds,_=prepare_seeds(candidate,query,**POSE_PARAMETERS)
    matrices=np.array([s.transform_matrix for s in seeds]).reshape(-1,4,4)
    generated=len(matrices)
    direct=direct_contact_masks(candidate.shape_points,sets,state,matrices)
    region_counts={k:int(v.sum()) for k,v in direct.items()}
    matrices=matrices[direct['all']]
    if not len(matrices): return finish(generated,region_counts=region_counts)
    q=tuple(query[k][indices] for k in ('feature_points','feature_types','feature_directions','feature_direction_kinds'))
    n=len(matrices)
    matching_types=np.array(features.feature_types,copy=True)
    for i,kind in enumerate(matching_types):
        if kind in (1,2) and not any(chemistry.atomic_numbers[j] in (7,8) for j in chemistry.feature_members[i]):
            matching_types[i]=0
    average, assignments, scores=match_batch((*q,np.ones(len(indices))),
        np.repeat(features.feature_points[None],n,axis=0), np.repeat(matching_types[None],n,axis=0),
        np.repeat(features.feature_directions[None],n,axis=0), np.repeat(features.feature_kinds[None],n,axis=0), matrices)
    valid=contact_mask(state,assignments,scores)
    flags['contact']=bool(valid.any())
    best=None; passing=None
    for j in np.flatnonzero(valid):
        matrix=matrices[j]
        score=score_overlay(query['shape_points'],candidate.shape_points,query['feature_points'],query['feature_types'],
                            candidate.feature_points,candidate.feature_types,matrix.ravel(),cutoff=None)
        shape,color=score['shape']['tanimoto'],score['color']['tanimoto']
        flags['shape'] |= shape>=POLICY['shape']; flags['color'] |= color>=POLICY['color']
        pose=dict(shape=shape,color=color,coverage=float(average[j]),transform=matrix.ravel().tolist())
        if best is None or shape>best['shape']:
            best=pose
        if shape<POLICY['shape']:
            continue
        moved=candidate.shape_points @ matrix[:3,:3].T + matrix[:3,3]
        tree=cKDTree(moved); geometry=state['geometry']
        minimum=float(state['tree'].query(moved)[0].min())
        envelope=float(np.mean(tree.query(geometry['contact'])[0]<=1.5))
        fractions=[float(np.mean(tree.query(geometry[f'region{i}'])[0]<=1.5)) for i in range(3)]
        geometry_pass=minimum>=POLICY['minimum_distance'] and envelope>=POLICY['envelope'] and min(fractions)>=POLICY['region']
        flags['geometry'] |= geometry_pass
        pose.update(minimum_distance=minimum,envelope=envelope,regions=fractions,geometry_pass=geometry_pass)
        if geometry_pass and color>=POLICY['color']:
            flags['full']=True
            if passing is None or shape+color>passing['shape']+passing['color']:
                passing=pose
    return finish(generated,best,passing,region_counts)


def self_control(state):
    q=state['query']; rotation=np.array([[0,-1,0],[0,0,-1],[1,0,0]])
    move=lambda x:x @ rotation.T + np.array([17,-11,6])
    candidate=SimpleNamespace(global_id=-1,molecule_id='native_control',conformer_id='rigid_perturbation',
        shape_points=move(q['shape_points']),feature_points=move(q['feature_points']),feature_types=q['feature_types'])
    features=SimpleNamespace(feature_points=candidate.feature_points,feature_types=q['feature_types'],
        feature_directions=q['feature_directions'] @ rotation.T,feature_kinds=q['feature_direction_kinds'])
    chemistry=SimpleNamespace(**state['record']['chemistry'])
    result=evaluate(candidate,features,state,chemistry)
    if not result['flags']['full']:
        raise ValueError('Native search recovery failed: '+json.dumps(result))
    return result
