from __future__ import annotations

import hashlib
import json
from pathlib import Path
import re
from typing import Callable
from urllib.request import Request, urlopen

from .project_context import ensure_within

PDB_ID = re.compile(r"^[0-9][A-Za-z0-9]{3}$")


def _download(url: str) -> tuple[bytes, str]:
    request = Request(url, headers={"User-Agent": "aidd-macrocycle-agent/0.2"})
    with urlopen(request, timeout=60) as response:
        return response.read(), response.headers.get_content_type()


def acquire_rcsb_mmcif(pdb_id: str, project_root: Path, *,
                       fetch: Callable[[str], tuple[bytes, str]] = _download) -> dict:
    code = pdb_id.upper()
    if not PDB_ID.fullmatch(code):
        raise ValueError("PDB ID must contain exactly four alphanumeric characters and start with a digit")
    root = project_root.resolve()
    destination = ensure_within(root / "inputs" / "structures" / f"{code}.cif", root)
    url = f"https://files.rcsb.org/download/{code}.cif"
    payload, content_type = fetch(url)
    prefix = payload[:4096].decode("utf-8", errors="ignore")
    if len(payload) < 40 or f"data_{code}".lower() not in prefix.lower():
        raise ValueError("RCSB response is not the requested mmCIF structure")
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_bytes(payload)
    digest = hashlib.sha256(payload).hexdigest()
    metadata = {
        "schema_version": 1, "source": "RCSB PDB", "pdb_id": code,
        "url": url, "content_type": content_type, "sha256": digest,
        "size_bytes": len(payload), "path": str(destination),
    }
    metadata_path = destination.with_suffix(".metadata.json")
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    return metadata

