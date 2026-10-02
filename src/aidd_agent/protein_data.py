"""Public protein evidence; sequence identity is supplied by UniProt, never LLM."""
import hashlib
import json
import re
from urllib.parse import urlencode
from urllib.request import Request, urlopen


def get_json(url):
    with urlopen(Request(url, headers={"User-Agent": "aidd-agent/0.2", "Accept": "application/json"}), timeout=60) as response:
        raw = response.read(8*1024*1024+1)
    if len(raw) > 8*1024*1024: raise ValueError("Protein response too large")
    return json.loads(raw)


def resolve_pdb_target(pdb_id, entity_id=None, fetch=get_json):
    """Resolve target identity from RCSB evidence, never from a guessed PDB label."""
    if not re.fullmatch(r'[0-9][A-Za-z0-9]{3}',pdb_id): raise ValueError('Invalid PDB ID')
    code=pdb_id.upper();base='https://data.rcsb.org/rest/v1/core/'
    entry=fetch(base+'entry/'+code);entities=entry['rcsb_entry_container_identifiers']['polymer_entity_ids']
    if entity_id is not None and entity_id not in entities: raise ValueError('Requested entity absent from PDB entry')
    candidates=[]
    for eid in entities:
        if entity_id is not None and eid!=entity_id:continue
        raw=fetch(base+f'polymer_entity/{code}/{eid}')
        identifiers=raw.get('rcsb_polymer_entity_container_identifiers',{})
        refs=identifiers.get('reference_sequence_identifiers') or []
        for ref in refs:
            if ref.get('database_name')=='UniProt':
                candidates.append(dict(entity_id=eid,accession=ref['database_accession'],
                    chains=identifiers.get('auth_asym_ids',[]),description=raw.get('rcsb_polymer_entity',{}).get('pdbx_description'),
                    source_url=base+f'polymer_entity/{code}/{eid}',source_payload_sha256=hashlib.sha256(json.dumps(raw,sort_keys=True).encode()).hexdigest()))
    accessions={c['accession'] for c in candidates}
    if len(accessions)!=1:
        return dict(status='needs_target_identity',pdb_id=code,candidates=candidates,
                    reason='Choose the target polymer entity; the entry has ambiguous or missing UniProt mappings')
    protein=fetch_protein(next(iter(accessions)),fetch=fetch)
    return dict(status='complete',pdb_id=code,candidates=candidates,protein=protein)


def fetch_protein(accession, organism_id=None, fetch=get_json):
    if not re.fullmatch(r"[A-Z0-9]{6}(?:[A-Z0-9]{4})?", accession): raise ValueError("Invalid accession")
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.json"
    raw = fetch(url)
    if raw.get("primaryAccession") != accession: raise ValueError("UniProt accession mismatch or obsolete accession")
    organism = raw.get("organism", {})
    if organism_id is not None and organism.get("taxonId") != organism_id: raise ValueError("Protein organism mismatch")
    sequence = raw.get("sequence", {}).get("value", "")
    if not sequence or not re.fullmatch(r"[A-Z]+", sequence) or len(sequence) != raw["sequence"].get("length"):
        raise ValueError("Invalid UniProt sequence or length")
    return dict(accession=accession, organism=organism, genes=raw.get("genes", []),
                protein_description=raw.get("proteinDescription", {}), sequence=sequence, length=len(sequence),
                sequence_sha256=hashlib.sha256(sequence.encode()).hexdigest(), source_url=url,
                source_payload_sha256=hashlib.sha256(json.dumps(raw, sort_keys=True).encode()).hexdigest(),
                reviewed=raw.get("entryType") == "UniProtKB reviewed (Swiss-Prot)", source="UniProt")


def resolve_protein(gene, organism_id, fetch=get_json):
    if not re.fullmatch(r"[A-Za-z][A-Za-z0-9_.-]{0,30}", gene) or type(organism_id) is not int or organism_id <= 0:
        raise ValueError("Invalid gene/taxonomy query")
    query = f"(gene_exact:{gene}) AND (organism_id:{organism_id}) AND (reviewed:true)"
    url = "https://rest.uniprot.org/uniprotkb/search?" + urlencode({"query": query, "format": "json", "size": 25})
    result = fetch(url)
    records = result.get("results")
    if not isinstance(records, list): raise ValueError("Malformed UniProt search response")
    matches = [r.get("primaryAccession") for r in records]
    if len(matches) != 1:
        raise ValueError(f"Target unresolved: {len(matches)} reviewed matches in first page; specify accession. Candidates: {matches}")
    protein = fetch_protein(matches[0], organism_id, fetch)
    names = [g.get("geneName", {}).get("value", "") for g in protein["genes"]]
    names += [s.get("value", "") for g in protein["genes"] for s in g.get("synonyms", [])]
    if gene.casefold() not in [n.casefold() for n in names]: raise ValueError("Resolved gene annotation mismatch")
    protein["resolution_query"] = dict(gene=gene, organism_id=organism_id, source_url=url)
    return protein
