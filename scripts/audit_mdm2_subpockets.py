"""Read-only E053 coordinate audit; spatial probes are not screening rules."""
import argparse
from collections import defaultdict
import csv
import json
from pathlib import Path

import numpy as np
from scipy.spatial.distance import cdist

from aidd_agent.consensus_admission import read_structure, ligand_atoms, ca_map, fit_protein
from aidd_agent.screening_selection import check_hashes
from aidd_agent.expanded_wee1 import fingerprint
from aidd_agent.conformer_artifacts import FEATURE_TYPES

SITES = {'A_Trp23': (23, 'TRP'), 'B_Leu26': (26, 'LEU'), 'C_Phe19': (19, 'PHE')}
USER_GROUPS = {'A_Trp23': [57, 86, 99], 'B_Leu26': [54, 93, 99, 103],
               'C_Phe19': [61, 62, 67, 75, 72, 58]}
BACKBONE = {'N', 'CA', 'C', 'O', 'OXT'}
NONPOLAR_ELEMENTS = {'C', 'S', 'F', 'CL', 'BR', 'I'}


def observed(atom):
    return not atom['label_alt_id'] and float(atom['occupancy']) >= .9 and not atom['pdbx_PDB_ins_code']


def xyz(atoms):
    return np.array([a['xyz'] for a in atoms], dtype=float)


def move(points, matrix):
    return points @ matrix[:3, :3].T + matrix[:3, 3]


def write_csv(path, rows):
    if not rows:
        path.write_text('', encoding='utf-8')
        return
    with path.open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def spatial_assign(points, probes, cutoff, margin=.5):
    distances = np.column_stack([cdist(points, probe).min(axis=1) for probe in probes])
    order = np.argsort(distances, axis=1)
    labels = []
    for i, ranks in enumerate(order):
        if distances[i, ranks[0]] > cutoff:
            labels.append('outside')
        elif distances[i, ranks[1]] - distances[i, ranks[0]] <= margin:
            labels.append('ambiguous')
        else:
            labels.append(list(SITES)[ranks[0]])
    return labels, distances


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--survey', type=Path, required=True)
    parser.add_argument('--peptide-reference', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    survey = json.loads(args.survey.read_text())
    if survey['target']['accession'] != 'Q00987':
        raise ValueError('This audit protocol is explicitly human MDM2-specific')
    check_hashes(survey['sources'])
    args.output.mkdir(parents=True, exist_ok=False)
    prepared = {r['query_id']: r for r in survey['prepared_complexes']}
    rows = [r for r in survey['cohort']['admitted'] if r['query_id'] in prepared]
    reference = next(r for r in rows if r['query_id'] == survey['cohort']['reference']['query_id'])
    structures = {}
    def structure(row):
        path = row['structure_path']
        if path not in structures:
            structures[path] = read_structure(Path(path), survey['target'])
        return structures[path]
    fixed = structure(reference)
    peptide = read_structure(args.peptide_reference, survey['target'])
    if len(peptide['chains']) != 1:
        raise ValueError('Ambiguous MDM2 chain in peptide reference')
    peptide_chain = next(iter(peptide['chains']))
    data = peptide['data']
    refs = list(zip(data['_struct_ref.pdbx_db_accession'], data['_struct_ref.pdbx_seq_one_letter_code']))
    p53seq = next(seq for acc, seq in refs if acc == 'P04637')
    # The standard receptor mapper intentionally rejects short peptides. Verify
    # this probe through explicit database correspondence and exact sequence.
    p53seq = ''.join(p53seq.split())
    p53mapping = {}
    p53notes = []
    for i, acc in enumerate(data['_struct_ref_seq.pdbx_db_accession']):
        if acc != 'P04637':
            continue
        chain = data['_struct_ref_seq.pdbx_strand_id'][i]
        begin = int(data['_struct_ref_seq.db_align_beg'][i])
        end = int(data['_struct_ref_seq.db_align_end'][i])
        label_begin = int(data['_struct_ref_seq.seq_align_beg'][i])
        entity_index = next(j for j, c in enumerate(data['_entity_poly.pdbx_strand_id']) if c == chain)
        sequence = ''.join(data['_entity_poly.pdbx_seq_one_letter_code_can'][entity_index].split())
        if label_begin != 1 or sequence != p53seq[begin-1:end]:
            raise ValueError('Peptide reference sequence correspondence is not exact')
        for atom in peptide['atoms']:
            if atom['auth_asym_id'] == chain and atom['label_seq_id'].isdigit():
                p53mapping[(chain, atom['auth_seq_id'])] = begin + int(atom['label_seq_id']) - 1
        p53notes.append(dict(accession=acc, chain=chain, begin=begin, end=end, sequence=sequence,
                             method='Explicit mmCIF database range and exact entity sequence'))
    rm = ca_map(fixed['chains'][reference['admission']['target_chain']])
    pm = ca_map(peptide['chains'][peptide_chain])
    common = sorted(set(reference['admission']['alignment_residues']) & set(rm) & set(pm))
    if len(common) < 12:
        raise ValueError('Insufficient local backbone correspondence')
    transform, rmsd, residuals = fit_protein([pm[r] for r in common], [rm[r] for r in common])
    probes = {}
    probe_rows = []
    for site, (residue, name) in SITES.items():
        atoms = [a for a in peptide['atoms'] if p53mapping.get((a['auth_asym_id'], a['auth_seq_id'])) == residue
                 and a['auth_comp_id'] == name and a['auth_atom_id'] not in BACKBONE and observed(a)]
        expected = {'TRP': 10, 'LEU': 4, 'PHE': 7}[name]
        if len(atoms) != expected or len({a['auth_asym_id'] for a in atoms}) != 1:
            raise ValueError('Missing or ambiguous peptide side-chain probe: ' + site)
        probes[site] = move(xyz(atoms), transform)
        for atom, point in zip(atoms, probes[site]):
            probe_rows.append(dict(site=site, chain=atom['auth_asym_id'], residue=residue,
                                   atom=atom['auth_atom_id'], x=point[0], y=point[1], z=point[2]))
    mode_queries = defaultdict(set)
    for anchor in survey['anchors']:
        if anchor['feature_class'] == 'hydrophobic':
            mode_queries[anchor['target_residue']].update(o['query_id'] for o in anchor['observations'])
    contacts, pairs, spatial, polar, complexes = [], [], [], [], []
    for row in rows:
        qid = row['query_id']
        st = structure(row)
        ligand_all = ligand_atoms(st, row)
        receptor_all = st['chains'][row['admission']['target_chain']]
        ligand = [a for a in ligand_all if observed(a)]
        receptor = [a for a in receptor_all if observed(a)]
        if not ligand or not receptor:
            raise ValueError('Missing usable geometry: ' + qid)
        matrix = np.array(row['admission']['transform'])
        lp, rp = move(xyz(ligand), matrix), move(xyz(receptor), matrix)
        distances = cdist(lp, rp)
        invariant_error = float(np.abs(distances - cdist(xyz(ligand), xyz(receptor))).max())
        with np.load(prepared[qid]['query_npz'], allow_pickle=False) as query:
            shape = query['shape_points']
            donor = query['feature_points'][query['feature_types'] == FEATURE_TYPES['Donor']]
            discrepancy = cdist(move(xyz(ligand_all), matrix), shape)
            shape_error = float(max(discrepancy.min(axis=0).max(), discrepancy.min(axis=1).max()))
        if invariant_error > 1e-5 or shape_error > 1e-3:
            raise ValueError('Coordinate-frame verification failed: ' + qid)
        labels, probe_distances = spatial_assign(lp, list(probes.values()), 2.)
        occupancy = {}
        for cutoff in (1.5, 2., 2.5):
            labs, _ = spatial_assign(lp, list(probes.values()), cutoff)
            occupancy[str(cutoff)] = {site: labs.count(site) for site in SITES}
        for i, atom in enumerate(ligand):
            spatial.append(dict(query_id=qid, ligand_atom=atom['auth_atom_id'], element=atom['type_symbol'],
                x=lp[i,0], y=lp[i,1], z=lp[i,2], site_at_2A=labels[i],
                distance_A=probe_distances[i,0], distance_B=probe_distances[i,1], distance_C=probe_distances[i,2]))
        byres = defaultdict(list)
        for j, atom in enumerate(receptor):
            byres[atom['canonical_residue']].append(j)
        for residue, indices in sorted(byres.items()):
            block = distances[:, indices]
            i, col = np.unravel_index(block.argmin(), block.shape)
            j = indices[col]
            nonpolar = [(li, rj) for li, la in enumerate(ligand) for rj in indices
                        if la['type_symbol'].upper() in NONPOLAR_ELEMENTS and
                        receptor[rj]['type_symbol'].upper() in NONPOLAR_ELEMENTS and
                        receptor[rj]['auth_atom_id'] not in BACKBONE]
            nmin = min((distances[li,rj] for li,rj in nonpolar), default=None)
            contacts.append(dict(query_id=qid, pdb_id=row['pdb_id'], canonical_residue=residue,
                residue_name=receptor[j]['auth_comp_id'], minimum_heavy_distance=float(block.min()),
                ligand_atom=ligand[i]['auth_atom_id'], protein_atom=receptor[j]['auth_atom_id'],
                nonpolar_sidechain_distance=nmin, extracted_hydrophobic=qid in mode_queries[residue],
                raw_atom_count=len([a for a in receptor_all if a['canonical_residue']==residue]),
                observed_atom_count=len(indices)))
        li, rj = np.where(distances <= 4.5)
        for i, j in zip(li, rj):
            pairs.append(dict(query_id=qid, pdb_id=row['pdb_id'], ligand_atom=ligand[i]['auth_atom_id'],
                ligand_element=ligand[i]['type_symbol'], canonical_residue=receptor[j]['canonical_residue'],
                residue_name=receptor[j]['auth_comp_id'], protein_atom=receptor[j]['auth_atom_id'],
                protein_element=receptor[j]['type_symbol'], distance=float(distances[i,j]), ligand_site_at_2A=labels[i]))
        for residue, name in ((54,'O'), (72,'O'), (72,'OE1'), (93,'O')):
            atoms = [a for a in receptor if a['canonical_residue']==residue and a['auth_atom_id']==name]
            d = float(cdist(donor, move(xyz(atoms),matrix)).min()) if atoms and len(donor) else None
            polar.append(dict(query_id=qid, residue=residue, protein_atom=name, observed=bool(atoms),
                              ligand_donor_features=len(donor), nearest_donor_distance=d,
                              interpretation='Distance only; hydrogen geometry and chemistry not validated'))
        complexes.append(dict(query_id=qid, pdb_id=row['pdb_id'], ligand_atoms=len(ligand),
            excluded_ligand_atoms=len(ligand_all)-len(ligand), excluded_receptor_atoms=len(receptor_all)-len(receptor),
            transform_distance_error=invariant_error, prepared_shape_error=shape_error,
            probe_overlap_atom_counts=occupancy, observed_residues=sorted(byres)))
    residues = sorted({r for group in USER_GROUPS.values() for r in group} | {96})
    summary = []
    for residue in residues:
        records = [r for r in contacts if r['canonical_residue']==residue]
        entry = dict(residue=residue, observed_complexes=len(records),
                     extracted_hydrophobic_complexes=len(mode_queries[residue]),
                     extracted_hydrophobic_pdbs=len({q.split(':')[0] for q in mode_queries[residue]}))
        for cutoff in (4., 4.5):
            selected = [r for r in records if r['nonpolar_sidechain_distance'] is not None and r['nonpolar_sidechain_distance']<=cutoff]
            entry['nonpolar_proxy_complexes_'+str(cutoff)] = len(selected)
            entry['nonpolar_proxy_pdbs_'+str(cutoff)] = len({r['pdb_id'] for r in selected})
        summary.append(entry)
    result = dict(status='complete',scope='Exploratory coordinate audit, not adopted occupancy rules',
        complex_count=len(complexes), pdb_count=len({r['pdb_id'] for r in complexes}), user_groups=USER_GROUPS,
        peptide_alignment=dict(reference=str(args.peptide_reference), target_chain=peptide_chain, transform=transform.tolist(),
            residues=common, ca_rmsd=rmsd, maximum_ca_residual=max(residuals), p53_mapping=p53notes),
        residue_summary=summary, complexes=complexes,
        sources=fingerprint([args.survey,args.peptide_reference,Path(__file__),*survey['sources']]),
        limitations=['Raw proximity is not a chemical interaction assignment.',
            'Spatial overlap with p53 side-chain atoms is not a validated cavity boundary.',
            'Counts are restricted to the admitted prepared cohort, not universal occupancy.',
            'Unobserved geometry is unknown; no energies, induced fit, or full-library NH retention inferred.'])
    for name, values in [('residue_contacts',contacts),('atom_contacts',pairs),('spatial_atoms',spatial),
                         ('polar_distance_checks',polar),('reference_probes',probe_rows),('residue_summary',summary)]:
        write_csv(args.output/(name+'.csv'), values)
    (args.output/'report.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps({k:result[k] for k in ('status','complex_count','pdb_count','residue_summary')},indent=2))


if __name__ == '__main__':
    main()
