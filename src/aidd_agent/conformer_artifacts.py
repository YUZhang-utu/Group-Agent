from __future__ import annotations

import argparse
import hashlib
import json
import multiprocessing as mp
from pathlib import Path
import sqlite3
import time

import numpy as np

from .mol2 import iter_mol2_records


SCHEMA_VERSION = 1
COORD_SCALE = 100.0
FEATURE_TYPES = {"Donor": 1, "Acceptor": 2, "PosIonizable": 3,
                 "NegIonizable": 4, "Hydrophobe": 5, "Aromatic": 6}
META_DTYPE = np.dtype([
    ("global_id", "<i8"), ("coord_offset", "<u8"), ("feature_offset", "<u8"),
    ("heavy_atoms", "<u2"), ("features", "<u2"), ("origin", "<f4", (3,)),
    ("bbox_extent", "<f4", (3,)), ("pmi", "<f4", (3,)),
    ("feature_counts", "u1", (6,)), ("reserved", "u1", (2,)),
])
FEATURE_DTYPE = np.dtype([("xyz", "<i2", (3,)), ("type", "u1"), ("parent", "<u2")])
_WORKER_FACTORY = None


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""): h.update(chunk)
    return h.hexdigest()


def _feature_factory():
    from rdkit import Chem, RDConfig
    from rdkit.Chem import ChemicalFeatures
    return ChemicalFeatures.BuildFeatureFactory(str(Path(RDConfig.RDDataDir) / "BaseFeatures.fdef"))


def _worker_init():
    global _WORKER_FACTORY
    _WORKER_FACTORY = _feature_factory()


def _heavy_coordinates(mol) -> tuple[np.ndarray, list[int]]:
    conf = mol.GetConformer()
    indices = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
    xyz = np.asarray([[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y,
                       conf.GetAtomPosition(i).z] for i in indices], dtype=np.float32)
    return xyz, indices


def _feature_template(mol, factory) -> tuple[tuple[str, tuple[int, ...]], ...]:
    return tuple((feature.GetFamily(), tuple(feature.GetAtomIds()))
                 for feature in factory.GetFeaturesForMol(mol)
                 if feature.GetFamily() in FEATURE_TYPES)


def _features_from_template(mol, template) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    conf = mol.GetConformer()
    points, types, counts = [], [], np.zeros(6, dtype=np.uint8)
    for family, atom_ids in template:
        type_id = FEATURE_TYPES.get(family)
        if type_id is None: continue
        xyz = np.asarray([[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y,
                           conf.GetAtomPosition(i).z] for i in atom_ids], dtype=np.float32)
        points.append(tuple(xyz.mean(axis=0)))
        types.append(type_id)
        if counts[type_id - 1] < 255: counts[type_id - 1] += 1
    return np.asarray(points, dtype=np.float32).reshape(-1, 3), np.asarray(types, dtype=np.uint8), counts


def _quantize(points: np.ndarray, origin: np.ndarray) -> np.ndarray:
    values = np.rint((points - origin) * COORD_SCALE)
    if values.size and (values.min() < -32768 or values.max() > 32767):
        raise ValueError("Coordinate range exceeds int16 at 0.01 angstrom")
    return values.astype("<i2")


def _process_group(task):
    from rdkit import Chem
    from rdkit.Chem import Descriptors3D, rdMolDescriptors
    molecule_id, records = task
    template = None; output = []
    for conformer_id, raw_text in records:
        mol = Chem.MolFromMol2Block(raw_text, sanitize=True, removeHs=False)
        if mol is None: raise ValueError(f"RDKit rejected {conformer_id}")
        if template is None: template = _feature_template(mol, _WORKER_FACTORY)
        xyz, _ = _heavy_coordinates(mol); origin = xyz.min(axis=0)
        feature_xyz, feature_types, counts = _features_from_template(mol, template)
        feature_records = np.empty(len(feature_types), dtype=FEATURE_DTYPE)
        feature_records["xyz"] = _quantize(feature_xyz, origin)
        feature_records["type"] = feature_types
        feature_records["parent"] = np.arange(len(feature_types), dtype=np.uint16)
        pmi = np.asarray([Descriptors3D.PMI1(mol), Descriptors3D.PMI2(mol),
                          Descriptors3D.PMI3(mol)], dtype=np.float32)
        output.append((conformer_id, molecule_id, _quantize(xyz, origin), feature_records,
                       origin, xyz.max(axis=0)-origin, pmi, counts,
                       np.asarray(rdMolDescriptors.GetUSRCAT(mol), dtype="<f4")))
    return output


def _grouped_tasks(source_path: Path, lookup: dict[int, sqlite3.Row]):
    current_molecule = None; records = []
    for record in iter_mol2_records(source_path):
        row = lookup.get(record.record_index)
        if row is None: continue
        if current_molecule is not None and row["molecule_id"] != current_molecule:
            yield current_molecule, tuple(records); records = []
        current_molecule = row["molecule_id"]
        records.append((row["conformer_id"], record.raw_text))
    if records: yield current_molecule, tuple(records)


def build_shard_parallel(db_path: Path, library_id: str, source_path: Path,
                         output_root: Path, global_id_start: int,
                         workers: int = 24) -> dict:
    from rdkit import RDLogger, rdBase
    RDLogger.DisableLog("rdApp.warning")
    source_path=source_path.resolve(); shard_name=source_path.stem
    final_dir=output_root/shard_name; partial_dir=output_root/f".{shard_name}.partial"
    if final_dir.exists(): return json.loads((final_dir/"manifest.json").read_text())
    partial_dir.mkdir(parents=True,exist_ok=True)
    with sqlite3.connect(db_path) as connection:
        connection.row_factory=sqlite3.Row
        rows=connection.execute("""SELECT c.id conformer_id,c.molecule_id,c.source_record_index
          FROM conformer c JOIN molecule m ON m.id=c.molecule_id
          WHERE m.library_id=? AND c.source_path=? ORDER BY c.source_record_index""",
          (library_id,str(source_path))).fetchall()
    lookup={r["source_record_index"]:r for r in rows}; started=time.perf_counter()
    coord_offset=feature_offset=valid=0
    context=mp.get_context("spawn")
    paths={name:partial_dir/name for name in ("coords.bin","feats.bin","meta.bin",
           "conformer_ids.bin","molecule_ids.bin","usrcat.f32.bin")}
    streams={name:path.open("wb") for name,path in paths.items()}
    try:
        with context.Pool(workers,initializer=_worker_init) as pool:
            for molecule_results in pool.imap(_process_group,_grouped_tasks(source_path,lookup),chunksize=8):
                for conformer_id,molecule_id,qxyz,features,origin,extent,pmi,counts,usr in molecule_results:
                    qxyz.tofile(streams["coords.bin"]);features.tofile(streams["feats.bin"])
                    meta=np.asarray([(global_id_start+valid,coord_offset,feature_offset,len(qxyz),len(features),
                                     origin,extent,pmi,counts,(0,0))],dtype=META_DTYPE)
                    meta.tofile(streams["meta.bin"])
                    np.asarray([conformer_id.encode()],dtype="S16").tofile(streams["conformer_ids.bin"])
                    np.asarray([molecule_id.encode()],dtype="S16").tofile(streams["molecule_ids.bin"])
                    usr.tofile(streams["usrcat.f32.bin"])
                    coord_offset+=len(qxyz);feature_offset+=len(features);valid+=1
    finally:
        for stream in streams.values(): stream.close()
    if valid!=len(rows): raise ValueError(f"Expected {len(rows)} records, built {valid}")
    files={p.name:{"bytes":p.stat().st_size,"sha256":_sha256(p)} for p in paths.values()}
    manifest={"format":"aidd-conformer-artifact-shard","version":SCHEMA_VERSION,
      "library_id":library_id,"source_path":str(source_path),"source_sha256":_sha256(source_path),
      "global_id_start":global_id_start,"conformers":valid,"heavy_atoms":coord_offset,
      "features":feature_offset,"coordinate_scale":COORD_SCALE,"coordinates":"heavy_atoms_only",
      "feature_types":FEATURE_TYPES,"workers":workers,"rdkit_version":rdBase.rdkitVersion,
      "seconds":time.perf_counter()-started,"files":files}
    (partial_dir/"manifest.json").write_text(json.dumps(manifest,indent=2),encoding="utf-8")
    partial_dir.replace(final_dir);return manifest


def build_shard(db_path: Path, library_id: str, source_path: Path,
                output_root: Path, global_id_start: int) -> dict:
    from rdkit import Chem, RDLogger, rdBase
    from rdkit.Chem import Descriptors3D, rdMolDescriptors
    RDLogger.DisableLog("rdApp.warning")
    source_path = source_path.resolve()
    shard_name = source_path.stem
    final_dir, partial_dir = output_root / shard_name, output_root / f".{shard_name}.partial"
    if final_dir.exists():
        return json.loads((final_dir / "manifest.json").read_text(encoding="utf-8"))
    partial_dir.mkdir(parents=True, exist_ok=True)
    with sqlite3.connect(db_path) as connection:
        connection.row_factory = sqlite3.Row
        rows = connection.execute(
            """SELECT c.id conformer_id,c.molecule_id,c.source_record_index
               FROM conformer c JOIN molecule m ON m.id=c.molecule_id
               WHERE m.library_id=? AND c.source_path=? ORDER BY c.source_record_index""",
            (library_id, str(source_path))).fetchall()
    lookup = {r["source_record_index"]: r for r in rows}
    factory = _feature_factory(); meta_rows=[]; conf_ids=[]; mol_ids=[]; usrcat=[]
    coord_offset=feature_offset=valid=0; started=time.perf_counter()
    previous_molecule_id = None; template = None
    with (partial_dir/"coords.bin").open("wb") as coord_file, \
         (partial_dir/"feats.bin").open("wb") as feat_file:
        for record in iter_mol2_records(source_path):
            row=lookup.get(record.record_index)
            if row is None: continue
            mol=Chem.MolFromMol2Block(record.raw_text,sanitize=True,removeHs=False)
            if mol is None: raise ValueError(f"RDKit rejected {row['conformer_id']}")
            if row["molecule_id"] != previous_molecule_id:
                template = _feature_template(mol, factory)
                previous_molecule_id = row["molecule_id"]
            xyz,_=_heavy_coordinates(mol); origin=xyz.min(axis=0)
            qxyz=_quantize(xyz,origin); qxyz.tofile(coord_file)
            feature_xyz,feature_types,counts=_features_from_template(mol,template)
            qfeatures=_quantize(feature_xyz,origin)
            feature_records=np.empty(len(feature_types),dtype=[("xyz","<i2",(3,)),("type","u1"),("parent","<u2")])
            feature_records["xyz"]=qfeatures; feature_records["type"]=feature_types
            feature_records["parent"]=np.arange(len(feature_types),dtype=np.uint16)
            feature_records.tofile(feat_file)
            pmi=np.asarray([Descriptors3D.PMI1(mol),Descriptors3D.PMI2(mol),Descriptors3D.PMI3(mol)],dtype=np.float32)
            meta_rows.append((global_id_start+valid,coord_offset,feature_offset,len(xyz),len(feature_types),
                              origin,xyz.max(axis=0)-origin,pmi,counts,(0,0)))
            conf_ids.append(row["conformer_id"].encode());mol_ids.append(row["molecule_id"].encode())
            usrcat.append(rdMolDescriptors.GetUSRCAT(mol));coord_offset+=len(xyz);feature_offset+=len(feature_types);valid+=1
    if valid != len(rows): raise ValueError(f"Expected {len(rows)} records, built {valid}")
    np.asarray(meta_rows,dtype=META_DTYPE).tofile(partial_dir/"meta.bin")
    np.asarray(conf_ids,dtype="S16").tofile(partial_dir/"conformer_ids.bin")
    np.asarray(mol_ids,dtype="S16").tofile(partial_dir/"molecule_ids.bin")
    np.asarray(usrcat,dtype="<f4").tofile(partial_dir/"usrcat.f32.bin")
    files={}
    for p in partial_dir.iterdir(): files[p.name]={"bytes":p.stat().st_size,"sha256":_sha256(p)}
    manifest={"format":"aidd-conformer-artifact-shard","version":SCHEMA_VERSION,
              "library_id":library_id,"source_path":str(source_path),"source_sha256":_sha256(source_path),
              "global_id_start":global_id_start,"conformers":valid,"coordinate_scale":COORD_SCALE,
              "coordinates": "heavy_atoms_only", "feature_types":FEATURE_TYPES,
              "rdkit_version":rdBase.rdkitVersion,"seconds":time.perf_counter()-started,"files":files}
    (partial_dir/"manifest.json").write_text(json.dumps(manifest,indent=2),encoding="utf-8")
    partial_dir.replace(final_dir); return manifest


def main():
    p=argparse.ArgumentParser();p.add_argument("--db",type=Path,required=True);p.add_argument("--library",required=True)
    p.add_argument("--source",type=Path,required=True);p.add_argument("--output-root",type=Path,required=True)
    p.add_argument("--global-id-start",type=int,required=True);p.add_argument("--workers",type=int,default=1);a=p.parse_args()
    function=build_shard_parallel if a.workers>1 else build_shard
    print(json.dumps(function(a.db,a.library,a.source,a.output_root,a.global_id_start,**({"workers":a.workers} if a.workers>1 else {})),indent=2));return 0

if __name__=="__main__":raise SystemExit(main())
