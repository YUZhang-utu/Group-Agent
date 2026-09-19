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
