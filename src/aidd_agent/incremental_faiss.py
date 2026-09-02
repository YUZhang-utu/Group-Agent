from __future__ import annotations

import argparse
import json
from pathlib import Path
import time

import numpy as np


def build_incremental(catalog_path: Path, output_dir: Path, *, nlist: int=512,
                      m: int=20, nbits: int=8, seed: int=20260902) -> dict:
    import faiss
    catalog=json.loads(catalog_path.read_text(encoding="utf-8"));shards=catalog["shards"]
    output_dir.mkdir(parents=True,exist_ok=True);rng=np.random.default_rng(seed)
    first_path=Path(shards[0]["path"])/"usrcat.f32.bin"
    first=np.memmap(first_path,dtype="<f4",mode="r").reshape(-1,60)
    train_size=min(len(first),max(30000,len(first)//10));ids=rng.choice(len(first),train_size,replace=False)
    sample=np.asarray(first[ids],dtype=np.float32);mean=sample.mean(0);std=sample.std(0);std[std<1e-6]=1
    training=np.ascontiguousarray((sample-mean)/std,dtype=np.float32)
    index=faiss.IndexIVFPQ(faiss.IndexFlatL2(60),60,nlist,m,nbits)
    t0=time.perf_counter();index.train(training);train_seconds=time.perf_counter()-t0
    additions=[]
    for position,shard in enumerate(shards):
        vector_path=Path(shard["path"])/"usrcat.f32.bin"
        raw=np.memmap(vector_path,dtype="<f4",mode="r").reshape(-1,60)
        transformed=np.ascontiguousarray((np.asarray(raw)-mean)/std,dtype=np.float32)
        global_ids=np.arange(shard["global_id_start"],shard["global_id_start"]+len(raw),dtype=np.int64)
        t0=time.perf_counter();index.add_with_ids(transformed,global_ids);seconds=time.perf_counter()-t0
        stage_path=output_dir/f"index_after_{position+1}_shards.faiss";faiss.write_index(index,str(stage_path))
        additions.append({"shard":shard["name"],"vectors":len(raw),"seconds":seconds,
                          "ntotal":index.ntotal,"index_path":str(stage_path.resolve()),
                          "index_bytes":stage_path.stat().st_size})
        # Reload between additions to validate persisted incremental operation.
        del index;index=faiss.read_index(str(stage_path))
    np.savez(output_dir/"transform.npz",mean=mean,std=std)
    manifest={"format":"aidd-incremental-faiss","version":1,"catalog":str(catalog_path.resolve()),
              "parameters":{"dimension":60,"nlist":nlist,"m":m,"nbits":nbits,"seed":seed,
                            "train_size":train_size,"direct_map":"disabled"},
              "train_seconds":train_seconds,"additions":additions,"final_ntotal":index.ntotal}
    (output_dir/"manifest.json").write_text(json.dumps(manifest,indent=2),encoding="utf-8");return manifest


def main():
    p=argparse.ArgumentParser();p.add_argument("--catalog",type=Path,required=True);p.add_argument("--output-dir",type=Path,required=True)
    p.add_argument("--nlist",type=int,default=512);p.add_argument("--m",type=int,default=20);a=p.parse_args()
    print(json.dumps(build_incremental(a.catalog,a.output_dir,nlist=a.nlist,m=a.m),indent=2));return 0

if __name__=="__main__":raise SystemExit(main())
