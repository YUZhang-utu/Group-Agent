"""Target-grounded crystal interaction evidence, recurrence and literature fallback."""
from __future__ import annotations

from collections import defaultdict
from itertools import islice
import json
from pathlib import Path
import re
from urllib.parse import urlencode

import numpy as np

from . import library_acceptance as ev
from .chemistry_prep import _mmcif_dict, _column, enumerate_ligand_instances, fetch_ccd
from .gaussian_batch import _atomic_json, prepare_gaussian_query, _load_query
from .expanded_wee1 import fingerprint


def target_residues(data, protein):
    """Map author residues through a uniquely aligned target-associated entity.

    Require an explicit mmCIF UniProt association and high identity/coverage;
    do not equate author residue numbers across structures.
    """
    from Bio.Align import PairwiseAligner
    refs = data.get('_struct_ref.id', [])
    if isinstance(refs,str): refs=[refs]
    accessions = _column(data,'_struct_ref.pdbx_db_accession',len(refs))
    entities = _column(data,'_struct_ref.entity_id',len(refs))
    databases = _column(data,'_struct_ref.db_name',len(refs))
    allowed = {e for e,a,d in zip(entities,accessions,databases)
               if a == protein['accession'] and d.upper() in ('UNP','UNIPROT')}
    poly = data.get('_entity_poly.entity_id',[])
    if isinstance(poly,str):poly=[poly]
    sequences = _column(data,'_entity_poly.pdbx_seq_one_letter_code_can',len(poly))
    maps = {}; notes=[]
    for entity, sequence in zip(poly,sequences):
        if entity not in allowed:continue
        sequence = re.sub(r'\s+','',sequence)
        if not re.fullmatch('[A-Z]+',sequence):continue
        aligner=PairwiseAligner(mode='local',match_score=2,mismatch_score=-3,
                               open_gap_score=-5,extend_gap_score=-.5)
        alignments=list(islice(iter(aligner.align(protein['sequence'],sequence)),2))
        if len(alignments)!=1:
            notes.append(dict(entity=entity,status='ambiguous_alignment'));continue
        alignment=alignments[0]; mapping={}; aligned=0
        for (a,b),(c,d) in zip(*alignment.aligned):
            for target,query in zip(range(a,b),range(c,d)):
                aligned+=1
                if protein['sequence'][target]==sequence[query]:mapping[query+1]=target+1
        identity=len(mapping)/max(aligned,1);coverage=aligned/max(len(sequence),1)
        if aligned<20 or identity<.95 or coverage<.8:
            notes.append(dict(entity=entity,status='insufficient_alignment',identity=identity,coverage=coverage));continue
        maps[entity]=mapping
        notes.append(dict(entity=entity,status='mapped',identity=identity,coverage=coverage,
                          method='UniProt-associated entity; unique sequence alignment; matching residues only'))
    labels=data.get('_atom_site.group_PDB',[])
    if isinstance(labels,str):labels=[labels]
    n=len(labels)
    cols={k:_column(data,'_atom_site.'+k,n) for k in
          ('label_entity_id','label_seq_id','auth_asym_id','auth_seq_id','pdbx_PDB_ins_code')}
    result={}
    for i,group in enumerate(labels):
        if group!='ATOM' or cols['pdbx_PDB_ins_code'][i]:continue
        entity=cols['label_entity_id'][i]
        seq=cols['label_seq_id'][i]
        if entity in maps and seq.isdigit() and int(seq) in maps[entity]:
            result[(cols['auth_asym_id'][i],cols['auth_seq_id'][i])]=maps[entity][int(seq)]
    return result,notes


def literature(protein, fetch):
    genes=[x.get('geneName',{}).get('value','') for x in protein.get('genes',[])]
    names=[protein['accession'],*[g for g in genes if g]]
    terms=' OR '.join('"'+re.sub(r'[^A-Za-z0-9_. -]','',x)+'"' for x in names)
    query=f'({terms}) AND (ligand OR inhibitor OR binding OR pharmacophore)'
    url='https://www.ebi.ac.uk/europepmc/webservices/rest/search?'+urlencode(
        dict(query=query,format='json',resultType='core',pageSize=20))
    raw=fetch(url)
    rows=[]
    for item in raw.get('resultList',{}).get('result',[]):
        if not item.get('id') or not item.get('source'):continue
        rows.append(dict(evidence_id='paper:'+item['source']+':'+item['id'],
            title=item.get('title',''),doi=item.get('doi'),year=item.get('pubYear'),
            abstract=item.get('abstractText',''),
            url=f"https://europepmc.org/article/{item['source']}/{item['id']}",
            interpretation='Search result; target relevance and experimental claims require review'))
    return dict(query=query,url=url,total_hits=raw.get('hitCount'),papers=rows)


def recurrence(queries):
    groups=defaultdict(lambda:dict(pdbs=set(),ligands=set(),anchors=[]))
    for q in queries:
        for a in q['anchors']:
            residue=a.get('target_residue')
            if residue is None:continue
            p=a['evidence'].get('protein_partner',{})
            key=(residue,a['feature_class'],p.get('atom_name',''))
            group=groups[key];group['pdbs'].add(q['pdb_id']);group['ligands'].add(q['ccd_id'])
            group['anchors'].append(a['anchor_id'])
    return [dict(evidence_id=f'contact:{residue}:{kind}:{atom}',target_residue=residue,
                 feature_class=kind,protein_atom=atom,distinct_structures=len(g['pdbs']),
                 pdb_ids=sorted(g['pdbs']),distinct_ligands=len(g['ligands']),anchor_ids=sorted(set(g['anchors'])),
                 eligible_structures=len({q['pdb_id'] for q in queries if residue in q['mapped_residues']}),
                 interpretation='Observed recurrence in analyzed structures; not proof of energetic necessity')
            for (residue,kind,atom),g in sorted(groups.items())]


def analyze_entry(entry, protein, output, services):
    from .anchor_extraction import extract_and_write_query_manifest, _protein_atoms
    from .screening_selection import anchor_evidence
    from .classified_features import crystal_environment, extra_hypotheses
    record=entry.get('download_record') or services.pdb_fetch(entry['pdb_id'],output)
    mmcif=Path(record['path'])
    mapping,notes=target_residues(_mmcif_dict(mmcif),protein)
    queries=[]; failures=[]
    instances=enumerate_ligand_instances(mmcif,entry['ligand_ids'])
    for instance in instances:
        query_id=f"{entry['pdb_id']}:{instance['ccd_id']}:{instance['chain_id']}:{instance['residue_number']}"
        if instance['altloc'] or instance['insertion_code'] or instance['model']!='1':
            failures.append(dict(query_id=query_id,reason='Alternate/insertion/model ambiguity requires explicit structure preparation'));continue
        if sum(a['element'].upper()!='H' for a in instance['atoms'])<6:continue
        folder=output/('query-'+str(len(queries))+'-'+str(len(failures)));folder.mkdir(parents=True,exist_ok=True)
        try:
            if not mapping:raise ValueError('No verified target residue mapping')
            ccd=fetch_ccd(instance['ccd_id'],output/'ccd')
            manifest=folder/'query_manifest.json'
            extract_and_write_query_manifest(mmcif,ccd,instance['ccd_id'],query_id,manifest)
            query_path=folder/'gaussian-query.npz'
            prepare_gaussian_query(mmcif,ccd,manifest,query_path)
            anchors,paths,query=anchor_evidence(query_path)
            atoms,waters,metals=crystal_environment(mmcif,query_id)
            anchors+=extra_hypotheses(query_id,query,atoms,waters,metals)
            retained=[]
            for a in anchors:
                p=a['evidence'].get('protein_partner',{})
                key=(p.get('chain'),p.get('residue_number'))
                if key not in mapping:continue
                a['target_residue']=mapping[key]
                a['diagnostic_counts']=[]
                retained.append(a)
            if not retained:raise ValueError('No mapped target contacts; ligand may belong to another chain/site')
            receptor=folder/'receptor.npz'
            np.savez_compressed(receptor,points=np.array([a['xyz'] for a in atoms]),
                                elements=np.array([a['element'] for a in atoms]))
            files=[*paths,receptor,manifest]
            queries.append(dict(query_id=query_id,pdb_id=entry['pdb_id'],ccd_id=instance['ccd_id'],
                resolution_angstrom=entry.get('resolution_angstrom'),query_npz=str(query_path.resolve()),
                receptor_npz=str(receptor.resolve()),anchors=retained,mapped_residues=sorted(set(mapping.values())),
                mapping=notes,sources=fingerprint(files),
                evidence_id='crystal:'+query_id,
                limitations=['Geometric contacts, not energetic validation','Pocket completeness and protonation unvalidated']))
        except Exception as exc:
            failures.append(dict(query_id=query_id,reason=str(exc) if isinstance(exc,ValueError) else type(exc).__name__))
    return queries,failures


def survey(protein, output, services, max_structures=12, max_resolution=3., pdb_ids=None):
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    report=dict(kind='structure_survey',status='running',target={k:protein[k] for k in
        ('accession','organism','sequence_sha256','source_url')},queries=[],failures=[],sources={})
    _atomic_json(output/'protein.json',protein)
    try:
        try:
            query,entries=services.pdb_search(protein['accession'],max_resolution=max_resolution)
            report['pdb_search_status']='complete'
        except Exception as exc:
            query={};entries=[]
            report['pdb_search_status']='failed'
            report['failures'].append(dict(stage='pdb_search',reason=type(exc).__name__))
        for code in pdb_ids or []:
            code=code.upper()
            if any(e['pdb_id']==code for e in entries):continue
            try:
                from .rcsb import classify_nonpolymer_components
                record=services.pdb_fetch(code,output/code)
                data=_mmcif_dict(Path(record['path']))
                components=data.get('_pdbx_entity_nonpoly.comp_id',[])
                if isinstance(components,str):components=[components]
                ligands,excluded=classify_nonpolymer_components(components)
                entries.insert(0,dict(pdb_id=code,ligand_ids=ligands,download_record=record,
                                      excluded_nonpolymer_components=excluded,user_supplied=True))
            except Exception as exc:
                report['failures'].append(dict(pdb_id=code,reason=type(exc).__name__))
        report.update(pdb_query=query,structures=entries,discovered_structures=len(entries),
                      analyzed_structure_limit=max_structures)
        eligible=[e for e in entries if e.get('ligand_ids')]
        report['unanalyzed_structures']=[e['pdb_id'] for e in eligible[max_structures:]]
        for entry in eligible[:max_structures]:
            try:
                rows,failures=analyze_entry(entry,protein,output/entry['pdb_id'],services)
                report['queries']+=rows;report['failures']+=failures
                for row in rows:report['sources'].update(row['sources'])
            except Exception as exc:
                report['failures'].append(dict(pdb_id=entry['pdb_id'],reason=type(exc).__name__))
            _atomic_json(output/'report.json',report)
        report['recurrence']=recurrence(report['queries'])
        if not report['queries']:
            try:report['literature']=literature(protein,services.fetch_json)
            except Exception as exc:report['literature']=dict(status='failed',error=type(exc).__name__,papers=[])
            report['needs_input']=['Provide a ligand-bound PDB ID associated with this target, or prepared receptor and reference ligand coordinates with chain/site identity. Literature alone cannot define 3D anchor positions.']
        report.update(status='complete',readiness='proposal_ready' if report['queries'] else 'needs_structure_input')
        report['sources'].update(fingerprint([output/'protein.json']))
        _atomic_json(output/'report.json',report)
        return report
    except Exception as exc:
        report.update(status='failed',error=type(exc).__name__)
        _atomic_json(output/'report.json',report)
        raise
