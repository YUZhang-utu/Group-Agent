import hashlib
import json

import pytest

from aidd_agent.artifact_catalog import build_catalog


def _shard(root,name,start,count):
    path=root/name;path.mkdir();data=path/"x.bin";data.write_bytes(name.encode())
    digest=hashlib.sha256(data.read_bytes()).hexdigest()
    manifest={"format":"aidd-conformer-artifact-shard","library_id":"LIB-X",
              "global_id_start":start,"conformers":count,"heavy_atoms":count*2,"features":count,
              "files":{"x.bin":{"bytes":data.stat().st_size,"sha256":digest}}}
    (path/"manifest.json").write_text(json.dumps(manifest))


def test_catalog_orders_contiguous_shards(tmp_path):
    _shard(tmp_path,"b",2,3);_shard(tmp_path,"a",0,2)
    result=build_catalog(tmp_path,"LIB-X")
    assert result["conformers"]==5
    assert [x["name"] for x in result["shards"]]==["a","b"]


def test_catalog_rejects_global_id_gap(tmp_path):
    _shard(tmp_path,"a",0,2);_shard(tmp_path,"b",3,1)
    with pytest.raises(ValueError,match="gap"):
        build_catalog(tmp_path,"LIB-X")
