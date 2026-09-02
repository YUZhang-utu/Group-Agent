from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def _sha256(path: Path) -> str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda:f.read(1024*1024),b""):h.update(chunk)
    return h.hexdigest()


def build_catalog(root: Path, library_id: str) -> dict:
    shards=[]
    for manifest_path in sorted(root.glob("*/manifest.json")):
        manifest=json.loads(manifest_path.read_text(encoding="utf-8"))
        if manifest.get("format")!="aidd-conformer-artifact-shard": continue
        if manifest["library_id"]!=library_id: raise ValueError("Cross-library shard")
        for name,record in manifest["files"].items():
            path=manifest_path.parent/name
            if path.stat().st_size!=record["bytes"] or _sha256(path)!=record["sha256"]:
                raise ValueError(f"Shard checksum mismatch: {path}")
        shards.append({"name":manifest_path.parent.name,"path":str(manifest_path.parent.resolve()),
                       "manifest_sha256":_sha256(manifest_path),"global_id_start":manifest["global_id_start"],
                       "conformers":manifest["conformers"],"heavy_atoms":manifest["heavy_atoms"],
                       "features":manifest["features"]})
    if not shards: raise ValueError("No artifact shards")
    shards.sort(key=lambda x:x["global_id_start"])
    expected=shards[0]["global_id_start"]
    if expected!=0: raise ValueError("Catalog global IDs must start at zero")
    for shard in shards:
        if shard["global_id_start"]!=expected: raise ValueError("Shard global IDs overlap or have a gap")
        expected+=shard["conformers"]
    catalog={"format":"aidd-conformer-artifact-catalog","version":1,"library_id":library_id,
             "conformers":expected,"heavy_atoms":sum(x["heavy_atoms"] for x in shards),
             "features":sum(x["features"] for x in shards),"shards":shards}
    path=root/"catalog.json";path.write_text(json.dumps(catalog,indent=2),encoding="utf-8")
    return catalog


def main():
    p=argparse.ArgumentParser();p.add_argument("--root",type=Path,required=True);p.add_argument("--library",required=True)
    a=p.parse_args();print(json.dumps(build_catalog(a.root,a.library),indent=2));return 0

if __name__=="__main__":raise SystemExit(main())
