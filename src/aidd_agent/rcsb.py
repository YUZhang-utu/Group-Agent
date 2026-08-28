from __future__ import annotations

import json
from typing import Any, Callable
from urllib.request import Request, urlopen

SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
GRAPHQL_URL = "https://data.rcsb.org/graphql"


def build_search_query(uniprot_id: str, *, max_resolution: float | None = 3.0,
                       method: str | None = "X-RAY DIFFRACTION") -> dict[str, Any]:
    nodes: list[dict[str, Any]] = [{
        "type": "terminal", "service": "text",
        "parameters": {"attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
                       "operator": "exact_match", "value": uniprot_id},
    }]
    if method:
        nodes.append({"type": "terminal", "service": "text", "parameters": {
            "attribute": "exptl.method", "operator": "exact_match", "value": method}})
    if max_resolution is not None:
        nodes.append({"type": "terminal", "service": "text", "parameters": {
            "attribute": "rcsb_entry_info.resolution_combined", "operator": "less_or_equal",
            "value": max_resolution}})
    return {"query": {"type": "group", "logical_operator": "and", "nodes": nodes},
            "return_type": "entry", "request_options": {"return_all_hits": True}}


def _post_json(url: str, payload: dict[str, Any]) -> dict[str, Any]:
    request = Request(url, data=json.dumps(payload).encode(),
                      headers={"Content-Type": "application/json", "User-Agent": "aidd-macrocycle-agent/0.1"})
    with urlopen(request, timeout=30) as response:
        return json.load(response)


def search_structures(uniprot_id: str, *, max_resolution: float | None = 3.0,
                      method: str | None = "X-RAY DIFFRACTION",
                      post_json: Callable[[str, dict[str, Any]], dict[str, Any]] = _post_json
                      ) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    query = build_search_query(uniprot_id, max_resolution=max_resolution, method=method)
    result = post_json(SEARCH_URL, query)
    ids = [row["identifier"].upper() for row in result.get("result_set", [])]
    if not ids:
        return query, []
    graphql = """query($ids:[String!]!){entries(entry_ids:$ids){rcsb_id struct{title}
      exptl{method} rcsb_entry_info{resolution_combined}
      rcsb_accession_info{deposit_date}
      nonpolymer_entities{pdbx_entity_nonpoly{comp_id}}}}"""
    detail = post_json(GRAPHQL_URL, {"query": graphql, "variables": {"ids": ids}})
    candidates = []
    for entry in detail.get("data", {}).get("entries", []) or []:
        if not entry or not entry.get("rcsb_id"):
            continue
        resolutions = (entry.get("rcsb_entry_info") or {}).get("resolution_combined") or []
        experiments = entry.get("exptl") or []
        nonpolymer_entities = entry.get("nonpolymer_entities") or []
        candidates.append({
            "pdb_id": entry["rcsb_id"].upper(),
            "title": (entry.get("struct") or {}).get("title", ""),
            "experimental_method": next(
                (item.get("method") for item in experiments if item and item.get("method")), None
            ),
            "resolution_angstrom": min(resolutions) if resolutions else None,
            "deposition_date": (entry.get("rcsb_accession_info") or {}).get("deposit_date"),
            "ligand_ids": sorted({(x.get("pdbx_entity_nonpoly") or {}).get("comp_id")
                                  for x in nonpolymer_entities if x
                                  if (x.get("pdbx_entity_nonpoly") or {}).get("comp_id")}),
        })
    candidates.sort(key=lambda x: (x["resolution_angstrom"] is None,
                                   x["resolution_angstrom"] or 999.0, x["pdb_id"]))
    return query, candidates
