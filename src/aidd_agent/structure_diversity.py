"""Experimental structure census and quality-gated ligand diversity proposals."""
from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
import csv
import hashlib
import json
from pathlib import Path
import re
import threading
from datetime import datetime, timezone
from urllib.parse import urlencode

import numpy as np

from .rcsb import build_search_query, classify_nonpolymer_components
from .structure_survey import target_residues


def save(path, value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,indent=2),encoding='utf-8')


class PublicCache:
    def __init__(self, root):
        self.root=Path(root);self.root.mkdir(parents=True,exist_ok=True)
        self.locks={}

    def get(self, url, payload=None, binary=False):
        key=hashlib.sha256((url+json.dumps(payload,sort_keys=True)).encode()).hexdigest()
        with self.locks.setdefault(key,threading.Lock()):
            return self._get(url,payload,binary)

    def _get(self, url, payload, binary):
        import requests
        from requests.adapters import HTTPAdapter
        from urllib3.util.retry import Retry
        key=hashlib.sha256((url+json.dumps(payload,sort_keys=True)).encode()).hexdigest()
        path=self.root/(key+('.bin' if binary else '.json'))
        if not path.exists():
            with requests.Session() as session:
                session.mount('https://',HTTPAdapter(max_retries=Retry(total=3,backoff_factor=1,
                    status_forcelist=[429,500,502,503,504],allowed_methods=['GET','POST'])))
                response=(session.get(url,timeout=60) if payload is None else session.post(url,json=payload,timeout=60))
                response.raise_for_status()
                content=response.content if binary else json.dumps(response.json() if response.status_code!=204 else {}).encode()
                path.write_bytes(content)
                save(path.with_suffix('.source.json'),dict(url=url,payload=payload,sha256=hashlib.sha256(content).hexdigest(),
                    retrieved_utc=datetime.now(timezone.utc).isoformat()))
        metadata=json.loads(path.with_suffix('.source.json').read_text(encoding='utf-8'))
        if metadata['sha256']!=hashlib.sha256(path.read_bytes()).hexdigest():raise ValueError('Changed cached public data')
        return path.read_bytes() if binary else json.loads(path.read_text(encoding='utf-8'))


def resolution_bin(value):
    if value is None:return 'unknown'
    for edge,label in [(1.5,'<=1.5'),(2.,'(1.5,2.0]'),(2.5,'(2.0,2.5]'),(3.,'(2.5,3.0]')]:
        if value<=edge:return label
    return '>3.0'


def diverse_indices(matrix, quality_order, maximum=8, similarity=.6):
    """Deterministic quality-seeded farthest-first cover; explicitly report cap."""
    if not quality_order:return []
    chosen=[quality_order[0]]
    while len(chosen)<min(maximum,len(quality_order)):
        remaining=[i for i in quality_order if i not in chosen]
        next_id=min(remaining,key=lambda i:(float(matrix[i,chosen].max()),quality_order.index(i)))
        if matrix[next_id,chosen].max()>=similarity:break
        chosen.append(next_id)
    return chosen


def inventory(folder):
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict
    rows=[]
    for path in sorted(Path(folder).iterdir()):
        if path.suffix.lower() not in {'.pdb','.cif','.mmcif'}:continue
        raw=path.read_bytes();text=raw.decode('utf-8',errors='replace')
        match=re.search(r'(?i)([0-9][a-z0-9]{3})',path.stem)
        row=dict(file=str(path.resolve()),sha256=hashlib.sha256(raw).hexdigest(),pdb_id=match[1].upper() if match else None)
        if path.suffix.lower()=='.pdb':
            atoms=[line for line in text.splitlines() if line.startswith(('ATOM  ','HETATM'))]
            row.update(chains=sorted({line[21:22].strip() for line in atoms}),
                nonwater_hetero_residues=sorted({line[17:20].strip() for line in atoms if line.startswith('HETATM') and line[17:20].strip()!='HOH'}),
                crystallographic_metadata_present='REMARK   2 RESOLUTION.' in text)
        else:
            data=MMCIF2Dict(str(path));groups=data.get('_atom_site.group_PDB',[])
            row.update(chains=sorted(set(data.get('_atom_site.auth_asym_id',[]))),
                nonwater_hetero_residues=sorted({c for g,c in zip(groups,data.get('_atom_site.label_comp_id',[])) if g=='HETATM' and c!='HOH'}),
                crystallographic_metadata_present='_refine.ls_d_res_high' in data,
                target_mapping_metadata_present='_struct_ref.pdbx_db_accession' in data)
        rows.append(row)
    return rows


GRAPHQL='''query($ids:[String!]!){entries(entry_ids:$ids){rcsb_id struct{title}
 exptl{method} rcsb_entry_info{resolution_combined} refine{ls_R_factor_R_free}
 polymer_entities{entity_poly{rcsb_entity_polymer_type pdbx_seq_one_letter_code_can}
 rcsb_polymer_entity{pdbx_description}
 rcsb_polymer_entity_container_identifiers{entity_id reference_sequence_identifiers{database_name database_accession}}}
 nonpolymer_entities{pdbx_entity_nonpoly{comp_id}
 nonpolymer_comp{chem_comp{id name} rcsb_chem_comp_descriptor{SMILES SMILES_stereo}}}}}'''


def fetch_entries(protein, cache):
    query=build_search_query(protein['accession'],max_resolution=None,method=None)
    search=cache.get('https://search.rcsb.org/rcsbsearch/v2/query?'+urlencode({'json':json.dumps(query)}))
    ids=sorted({r['identifier'].upper() for r in search.get('result_set',[])})
    entries=[]
    for start in range(0,len(ids),25):
        raw=cache.get('https://data.rcsb.org/graphql',dict(query=GRAPHQL,variables=dict(ids=ids[start:start+25])))
        if raw.get('errors'):raise ValueError('RCSB GraphQL schema/query failure: '+json.dumps(raw['errors']))
        entries.extend(r for r in raw['data']['entries'] if r)
    if {e['rcsb_id'] for e in entries}!=set(ids):raise ValueError('Incomplete RCSB entry coverage')
    return query,entries


def summarize_entry(e, accession):
    resolutions=(e.get('rcsb_entry_info') or {}).get('resolution_combined') or []
    method=[r['method'] for r in e.get('exptl') or []]
    components=[r['pdbx_entity_nonpoly']['comp_id'] for r in e.get('nonpolymer_entities') or []]
    ligands,excluded=classify_nonpolymer_components(components)
    partners=[];targets=[]
    for p in e.get('polymer_entities') or []:
        ids=p['rcsb_polymer_entity_container_identifiers']
        refs=ids.get('reference_sequence_identifiers') or []
        if any(r['database_name']=='UniProt' and r['database_accession']==accession for r in refs):targets.append(ids['entity_id'])
        else:partners.append(dict(entity_id=ids['entity_id'],description=(p.get('rcsb_polymer_entity') or {}).get('pdbx_description'),
            polymer_type=p['entity_poly']['rcsb_entity_polymer_type'],length=len(re.sub(r'\s+','',p['entity_poly'].get('pdbx_seq_one_letter_code_can') or ''))))
    free=[r['ls_R_factor_R_free'] for r in e.get('refine') or [] if r.get('ls_R_factor_R_free') is not None]
    return dict(pdb_id=e['rcsb_id'],title=e['struct']['title'],methods=method,
        resolution=min(resolutions) if resolutions else None,resolution_bin=resolution_bin(min(resolutions) if resolutions else None),
        r_free=max(free) if free else None,ligand_ids=ligands,excluded_components=excluded,
        target_entities=targets,other_polymers=partners,
        nonpolymer_complex_candidate=bool(ligands),polymer_complex_candidate=bool(partners))


def inspect_structure(entry, protein, output, cache):
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict
    from scipy.spatial import cKDTree
    from rdkit import Chem
    code=entry['pdb_id'];path=output/'structures'/f'{code}.cif';path.parent.mkdir(exist_ok=True)
    original=cache.get(f'https://files.rcsb.org/download/{code}.cif',binary=True)
    if path.exists() and path.read_bytes()!=original:raise ValueError('Downloaded coordinate file changed')
    if not path.exists():path.write_bytes(original)
    data=MMCIF2Dict(str(path));mapping,notes=target_residues(data,protein)
    n=len(data['_atom_site.group_PDB'])
    def col(name,default=''):return data.get('_atom_site.'+name,[default]*n)
    groups,chains,resnums,comp,labels,elements,models,alt,occ=[col(k,d) for k,d in [
        ('group_PDB',''),('auth_asym_id',''),('auth_seq_id',''),('auth_comp_id',''),('auth_atom_id',''),
        ('type_symbol',''),('pdbx_PDB_model_num','1'),('label_alt_id','.'),('occupancy','1')]]
    xyz=np.column_stack([np.asarray(col(k),float) for k in ('Cartn_x','Cartn_y','Cartn_z')])
    target=[i for i in range(n) if groups[i]=='ATOM' and models[i]=='1' and elements[i] not in {'H','D'} and (chains[i],resnums[i]) in mapping]
    if not target:return [],dict(pdb_id=code,error='No mapped target heavy atoms',mapping=notes)
    tree=cKDTree(xyz[target]);instances={}
    for i in range(n):
        if groups[i]=='HETATM' and models[i]=='1' and comp[i] in entry['ligand_ids'] and elements[i] not in {'H','D'}:
            instances.setdefault((comp[i],chains[i],resnums[i]),[]).append(i)
    rows=[]
    for (ccd,chain,residue),indices in instances.items():
        chemical=cache.get(f'https://data.rcsb.org/rest/v1/core/chemcomp/{ccd}')
        desc=chemical.get('rcsb_chem_comp_descriptor') or {}
        smiles=desc.get('SMILES_stereo') or desc.get('SMILES')
        mol=Chem.MolFromSmiles(smiles) if smiles else None
        if mol is None or not any(a.GetAtomicNum()==6 for a in mol.GetAtoms()):continue
        contacts=set();contact_chains=set()
        for nearby in tree.query_ball_point(xyz[indices],4.5):
            for j in nearby:
                index=target[j];contacts.add(mapping[chains[index],resnums[index]]);contact_chains.add(chains[index])
        expected=mol.GetNumHeavyAtoms();observed=len({labels[i] for i in indices})
        occupancy=min(float(occ[i]) for i in indices)
        alternates=sorted({alt[i] for i in indices if alt[i] not in {'.','?',''}})
        reasons=[]
        if entry['resolution'] is None or entry['resolution']>2.5:reasons.append('resolution_above_2.5_or_unknown')
        if observed!=expected:reasons.append('incomplete_or_unexpected_heavy_atom_count')
        if occupancy<.9:reasons.append('occupancy_below_0.9')
        if alternates:reasons.append('alternate_ligand_conformers')
        if len(contacts)<3:reasons.append('fewer_than_three_target_contact_residues')
        rows.append(dict(query_id=f'{code}:{ccd}:{chain}:{residue}',pdb_id=code,ccd_id=ccd,chain=chain,residue=residue,
            resolution=entry['resolution'],r_free=entry['r_free'],name=chemical['chem_comp'].get('name'),
            smiles=Chem.MolToSmiles(mol,isomericSmiles=True),heavy_atoms=expected,observed_heavy_atoms=observed,
            minimum_occupancy=occupancy,alternate_locations=alternates,target_contacts=sorted(contacts),target_chains=sorted(contact_chains),
            quality_passed=not reasons,quality_rejections=reasons,structure_path=str(path.resolve()),
            density_validation='not_evaluated'))
    return rows,None


def run(protein, output, local_folder=None, reference_pdb='5C5A', maximum=8, similarity=.6):
    from rdkit import Chem, DataStructs, rdBase
    from rdkit.Chem import rdFingerprintGenerator
    output=Path(output);output.mkdir(parents=True,exist_ok=True);cache=PublicCache(output/'raw')
    local=inventory(local_folder) if local_folder else [];save(output/'local-inventory.json',local)
    save(output/'protein.json',protein)
    query,raw=fetch_entries(protein,cache);entries=[summarize_entry(e,protein['accession']) for e in raw]
    save(output/'entries.json',dict(query=query,entries=entries))
    bins=[]
    for label in ['<=1.5','(1.5,2.0]','(2.0,2.5]','(2.5,3.0]','>3.0','unknown']:
        rows=[e for e in entries if 'X-RAY DIFFRACTION' in e['methods'] and e['resolution_bin']==label]
        bins.append(dict(bin=label,xray_entries=len(rows),nonpolymer_candidates=sum(e['nonpolymer_complex_candidate'] for e in rows),
            other_polymer_candidates=sum(e['polymer_complex_candidate'] for e in rows)))
    save(output/'census.json',dict(total_entries=len(entries),methods=dict(Counter(m for e in entries for m in e['methods'])),bins=bins,
        interpretation='Metadata complex candidates before target-site contact verification; categories overlap'))
    print(f'Discovered {len(entries)} entries. Census: {json.dumps(bins)}',flush=True)
    local_ids={r['pdb_id'] for r in local}
    chosen=[e for e in entries if ('X-RAY DIFFRACTION' in e['methods'] and e['resolution'] is not None and e['resolution']<=3) or e['pdb_id'] in local_ids]
    instances=[];failures=[]
    def inspect(e):
        try:return inspect_structure(e,protein,output,cache)
        except Exception as exc:return [],dict(pdb_id=e['pdb_id'],error=type(exc).__name__,detail=str(exc)[:500])
    with ThreadPoolExecutor(max_workers=4) as pool:
        for index,(rows,error) in enumerate(pool.map(inspect,chosen),1):
            instances.extend(rows)
            if error:failures.append(error)
            print(f'Coordinates checked {index}/{len(chosen)}; ligand instances {len(instances)}; failures {len(failures)}',flush=True)
    site=set().union(*(set(r['target_contacts']) for r in instances if r['pdb_id']==reference_pdb and r['quality_passed']))
    for r in instances:
        r['site_overlap_residues']=sorted(site.intersection(r['target_contacts']))
        r['same_reference_site']=len(r['site_overlap_residues'])>=3
    save(output/'ligand-instances.json',instances);save(output/'failures.json',failures)
    unique={}
    for r in sorted(instances,key=lambda r:(r['resolution'] or 999,r['query_id'])):
        if r['quality_passed'] and r['same_reference_site']:unique.setdefault(r['smiles'],[]).append(r)
    molecules=list(unique);generator=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048,includeChirality=True)
    fps=[generator.GetFingerprint(Chem.MolFromSmiles(s)) for s in molecules]
    matrix=np.array([DataStructs.BulkTanimotoSimilarity(fp,fps) for fp in fps]) if fps else np.empty((0,0))
    names=[unique[s][0]['query_id'] for s in molecules]
    with (output/'ligand-similarity.csv').open('w',newline='',encoding='utf-8') as stream:
        writer=csv.writer(stream);writer.writerow(['query_id',*names])
        for name,row in zip(names,matrix):writer.writerow([name,*map(float,row)])
    selected=diverse_indices(matrix,list(range(len(molecules))),maximum,similarity)
    uncovered=[names[i] for i in range(len(names)) if not selected or matrix[i,selected].max()<similarity]
    refs=[dict(unique[molecules[i]][0],crystal_occurrences=[r['query_id'] for r in unique[molecules[i]]]) for i in selected]
    report=dict(kind='structure_diversity',status='complete',classification='exploratory_reference_proposal',target=protein['accession'],
        generated_utc=datetime.now(timezone.utc).isoformat(),
        quality_policy=dict(download_xray_resolution_max=3.,proposal_resolution_max=2.5,minimum_occupancy=.9,
            complete_ligand_heavy_atom_count=True,no_ligand_altloc=True,contact_distance_angstrom=4.5,
            minimum_target_contact_residues=3,minimum_shared_reference_contact_residues=3),
        discovered_entries=len(entries),downloaded_entries=sum((output/'structures'/f'{e["pdb_id"]}.cif').exists() for e in chosen),
        successfully_checked_entries=len(chosen)-len(failures),failed_entries=failures,
        local_structures=local,local_ids_not_in_target_search=sorted(local_ids-{e['pdb_id'] for e in entries}),
        census=bins,reference_site_pdb=reference_pdb,reference_site_residues=sorted(site),
        organic_ligand_instances=len(instances),quality_site_unique_ligands=len(unique),
        proposed_references=refs,uncovered_at_similarity_threshold=uncovered,
        selection=dict(maximum_references=maximum,similarity_threshold=similarity,method='quality-seeded farthest-first',
                       fingerprint='Morgan radius 2, 2048 bits, chirality',rdkit=rdBase.rdkitVersion),
        limitations=['Nonpolymeric organic ligand diversity only; polymer ligands reported separately',
            'Resolution/occupancy/completeness screen is not density or biological validation',
            'No new library scan; multi-reference selection is a proposal for human review'])
    report['code_hashes']={name:hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
        for name in ('structure_diversity.py','polymer_ligand_diversity.py','structure_diversity_report.py','structure_survey.py','rcsb.py')}
    save(output/'report.json',report)
    from .polymer_ligand_diversity import analyze
    peptides=analyze(output)
    report['polymer_ligands']=dict(unique_sequence_link_variants=peptides['unique_sequence_link_variants'],
        sequence_diversity_examples=peptides['sequence_diversity_examples'],limitations=peptides['limitations'])
    report['coverage_curve']=[dict(references=len(panel),covered=int(sum(matrix[i,panel].max()>=similarity for i in range(len(names)))))
        for limit in sorted({4,8,12,16,len(names)}) if limit>0
        for panel in [diverse_indices(matrix,list(range(len(names))),limit,similarity)] if panel]
    report['full_cover_reference_ids']=[names[i] for i in diverse_indices(matrix,list(range(len(names))),len(names),similarity)]
    for b in bins:
        in_bin={e['pdb_id'] for e in entries if e['resolution_bin']==b['bin'] and 'X-RAY DIFFRACTION' in e['methods']}
        b['coordinate_checked_organic_contact_entries']=len({r['pdb_id'] for r in instances if r['pdb_id'] in in_bin and len(r['target_contacts'])>=3})
        b['coordinate_checked_peptide_contact_entries']=len({r['pdb_id'] for r in peptides['instances'] if r['pdb_id'] in in_bin and r['length']<=50 and len(r['target_contacts'])>=3})
        b['coordinate_scope']='Downloaded <=3 A structures only; higher-resolution-number bins not assessed'
    report['outputs']={name:str((output/name).resolve()) for name in ('report.md','ligand-similarity.csv','ligand-instances.json','polymer-ligands.json','peptide-sequence-similarity.csv')}
    save(output/'report.json',report)
    lines=['# Crystal census and reference proposal','',f'Target: {protein["accession"]}. Experimental entries: {len(entries)}.',
        '', '| Resolution (A) | X-ray entries | Nonpolymer candidates | Other-polymer candidates |','|---|---:|---:|---:|']
    lines += [f'| {b["bin"]} | {b["xray_entries"]} | {b["nonpolymer_candidates"]} | {b["other_polymer_candidates"]} |' for b in bins]
    lines+=['','Metadata candidate categories overlap. Target contacts are checked separately.',
        '',f'Quality/site-qualified unique organic ligands: {len(unique)}. Proposals: {len(refs)}. Uncovered: {len(uncovered)}.',
        '', '| Reference | Resolution | Ligand heavy atoms | Target contact residues |','|---|---:|---:|---|']
    lines += [f'| {r["query_id"]} | {r["resolution"]} | {r["heavy_atoms"]} | {r["target_contacts"]} |' for r in refs]
    lines+=['','## Limits','',*['- '+x for x in report['limitations']]]
    lines+=['','## Polymer ligands','',f"Sequence/link variants at the same pocket (<=2.5 A, <=50 residues): {peptides['unique_sequence_link_variants']}.",
        'These are sequence diagnostics, not chemical similarity or prepared 3D search references.','',
        *[f"- {r['query_id']}: {r['resolution']} A; {r['length']} residues; {r['description']}" for r in peptides['sequence_diversity_examples']],
        '', '## Coverage sensitivity','',json.dumps(report['coverage_curve']),
        '',f"The greedy full-cover proposal uses {len(report['full_cover_reference_ids'])} representatives at similarity {similarity}. This is not a minimum-set proof or screening recall."]
    (output/'report.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    from .structure_diversity_report import render
    report['outputs']['review.html']=str((output/'review.html').resolve())
    save(output/'report.json',report)
    render(output)
    files=[p for p in output.rglob('*') if p.is_file() and p.name not in {'artifact-manifest.json','run.log'}]
    save(output/'artifact-manifest.json',{str(p.relative_to(output)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files})
    return report


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True);parser.add_argument('--local-folder',type=Path)
    parser.add_argument('--accession',default='Q00987');parser.add_argument('--reference-pdb',default='5C5A')
    args=parser.parse_args()
    from .protein_data import fetch_protein
    cache=PublicCache(args.output/'raw')
    protein=fetch_protein(args.accession,9606,cache.get)
    report=run(protein,args.output,args.local_folder,args.reference_pdb)
    print(json.dumps({k:report[k] for k in ('status','discovered_entries','quality_site_unique_ligands','uncovered_at_similarity_threshold')},indent=2))


if __name__=='__main__':main()
