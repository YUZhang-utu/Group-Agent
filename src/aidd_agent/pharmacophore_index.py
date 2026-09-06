from __future__ import annotations

from collections import OrderedDict
from dataclasses import asdict, dataclass
import hashlib
import json
import math
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import numpy as np

from .conformer_artifacts import COORD_SCALE, FEATURE_DTYPE, META_DTYPE


INDEX_FORMAT = "aidd-local-pharmacophore-index"
INDEX_VERSION = 1
QUERY_FORMAT = "aidd-local-pharmacophore-query"
RESULT_FORMAT = "aidd-local-pharmacophore-results"

FEATURE_ALIASES = {
    "HBD": 1, "DONOR": 1,
    "HBA": 2, "ACCEPTOR": 2,
    "POSIONIZABLE": 3, "POSITIVE": 3, "CATION": 3,
    "NEGIONIZABLE": 4, "NEGATIVE": 4, "ANION": 4,
    "HYDROPHOBE": 5, "HYDROPHOBIC": 5,
    "AROMATIC": 6, "AROMATICRING": 6,
}


@dataclass(frozen=True)
class SearchProfile:
    name: str
    tolerance_angstrom: float
    required_pairs: int

    def __post_init__(self) -> None:
        if not self.name or self.tolerance_angstrom < 0 or self.required_pairs <= 0:
            raise ValueError("invalid pharmacophore search profile")


DEFAULT_PROFILES = (
    SearchProfile("loose", 2.0, 1),
    SearchProfile("balanced", 1.0, 2),
    SearchProfile("strict", 0.5, 3),
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _canonical_json_hash(value: object) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(payload).hexdigest()


def encode_pair_key(first_type: int, second_type: int, distance_bin: int) -> int:
    first, second = sorted((int(first_type), int(second_type)))
    if not (1 <= first <= 255 and 1 <= second <= 255):
        raise ValueError("feature type IDs must fit uint8 and be positive")
    if not (0 <= distance_bin <= 65535):
        raise ValueError("distance bin must fit uint16")
    return (first << 24) | (second << 16) | distance_bin


def decode_pair_key(key: int) -> tuple[int, int, int]:
    return (int(key) >> 24) & 255, (int(key) >> 16) & 255, int(key) & 65535


def _distance_bin(distance: float, bin_width: float) -> int:
    return int(math.floor(float(distance) / bin_width))


def _feature_points(meta: np.void, features: np.ndarray, scale: float) -> np.ndarray:
    return np.asarray(meta["origin"], dtype=np.float64) + (
        np.asarray(features["xyz"], dtype=np.float64) / scale)


def _conformer_pair_keys(points: np.ndarray, types: np.ndarray, bin_width: float,
                         max_distance: float) -> np.ndarray:
    """Return unique uint32 keys using one vectorized upper-triangle pass."""
    if len(points) < 2:
        return np.empty(0, dtype=np.uint32)
    left, right = np.triu_indices(len(points), 1)
    distances = np.linalg.norm(points[left] - points[right], axis=1)
    keep = distances <= max_distance
    if not np.any(keep):
        return np.empty(0, dtype=np.uint32)
    left, right, distances = left[keep], right[keep], distances[keep]
    first = np.minimum(types[left], types[right]).astype(np.uint32)
    second = np.maximum(types[left], types[right]).astype(np.uint32)
    bins = np.floor(distances / bin_width).astype(np.uint32)
    if len(bins) and int(bins.max()) > 65535:
        raise ValueError("distance bin must fit uint16")
    encoded = (first << 24) | (second << 16) | bins
    return np.unique(encoded)


def build_shard_pharmacophore_index(shard_directory: Path, output_directory: Path,
                                     bin_width: float = 0.5,
                                     max_distance: float = 20.0,
                                     flush_records: int = 1_000_000) -> dict:
    """Build one immutable postings index from precomputed local artifacts."""
    if bin_width <= 0 or max_distance <= 0 or flush_records <= 0:
        raise ValueError("bin width, maximum distance and flush size must be positive")
    shard_directory = shard_directory.resolve()
    source_manifest_path = shard_directory / "manifest.json"
    source_manifest = json.loads(source_manifest_path.read_text(encoding="utf-8"))
    source_hash = _sha256(source_manifest_path)
    final = output_directory.resolve()
    if final.exists():
        existing = json.loads((final / "manifest.json").read_text(encoding="utf-8"))
        expected = (source_hash, float(bin_width), float(max_distance))
        observed = (existing.get("source_manifest_sha256"), existing.get("bin_width_angstrom"),
                    existing.get("max_distance_angstrom"))
        if observed != expected:
            raise ValueError("existing pharmacophore index has different inputs or parameters")
        return existing
    partial = final.with_name(f".{final.name}.partial")
    if partial.exists():
        raise ValueError(f"partial pharmacophore index already exists: {partial}")
    postings_dir = partial / "postings"
    postings_dir.mkdir(parents=True)

    meta = np.memmap(shard_directory / "meta.bin", dtype=META_DTYPE, mode="r")
    feats = np.memmap(shard_directory / "feats.bin", dtype=FEATURE_DTYPE, mode="r")
    scale = float(source_manifest.get("coordinate_scale", COORD_SCALE))
    buffered_keys: list[np.ndarray] = []
    buffered_ids: list[np.ndarray] = []
    counts: dict[int, int] = {}
    streams: OrderedDict[int, object] = OrderedDict()
    max_open_streams = 256
    buffered = 0

    def posting_stream(key: int):
        stream = streams.pop(key, None)
        if stream is None:
            if len(streams) >= max_open_streams:
                _, oldest = streams.popitem(last=False)
                oldest.close()
            stream = (postings_dir / f"{key:08x}.i64").open("ab")
        streams[key] = stream
        return stream

    def flush() -> None:
        nonlocal buffered
        keys = np.concatenate(buffered_keys)
        gids = np.concatenate(buffered_ids)
        order = np.argsort(keys, kind="stable")
        keys, gids = keys[order], gids[order]
        unique_keys, starts, key_counts = np.unique(
            keys, return_index=True, return_counts=True)
        for key, start, key_count in zip(unique_keys, starts, key_counts):
            values = np.asarray(gids[start:start + key_count], dtype="<i8")
            integer_key = int(key)
            values.tofile(posting_stream(integer_key))
            counts[integer_key] = counts.get(integer_key, 0) + len(values)
        buffered_keys.clear(); buffered_ids.clear()
        buffered = 0

    previous_gid = -1
    row = None
    try:
        for row in meta:
            gid = int(row["global_id"])
            if gid <= previous_gid:
                raise ValueError("artifact global IDs must be strictly increasing within a shard")
            previous_gid = gid
            offset, count = int(row["feature_offset"]), int(row["features"])
            record = feats[offset:offset + count]
            points = _feature_points(row, record, scale)
            keys = _conformer_pair_keys(points, record["type"], bin_width, max_distance)
            if len(keys):
                buffered_keys.append(keys)
                buffered_ids.append(np.full(len(keys), gid, dtype="<i8"))
                buffered += len(keys)
            if buffered >= flush_records:
                flush()
        if buffered:
            flush()
    finally:
        for stream in streams.values():
            stream.close()
        streams.clear()
    # Windows will not atomically rename the directory while memmaps are open.
    del row, meta, feats

    keys = []
    for key in sorted(counts):
        path = postings_dir / f"{key:08x}.i64"
        values = np.memmap(path, dtype="<i8", mode="r")
        if len(values) > 1 and np.any(values[1:] <= values[:-1]):
            raise ValueError(f"postings are not strictly increasing for key {key}")
        del values
        first, second, distance_bin = decode_pair_key(key)
        keys.append({
            "key": key, "first_type": first, "second_type": second,
            "distance_bin": distance_bin, "count": counts[key],
            "file": f"postings/{key:08x}.i64", "sha256": _sha256(path),
        })
    manifest = {
        "format": INDEX_FORMAT, "version": INDEX_VERSION,
        "library_id": source_manifest["library_id"],
        "source_shard": shard_directory.name,
        "source_manifest": str(source_manifest_path),
        "source_manifest_sha256": source_hash,
        "global_id_start": int(source_manifest["global_id_start"]),
        "conformers": int(source_manifest["conformers"]),
        "bin_width_angstrom": float(bin_width),
        "max_distance_angstrom": float(max_distance),
        "coordinate_scale": scale,
        "feature_types": source_manifest.get("feature_types", {}),
        "keys": keys,
    }
    (partial / "manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8")
    partial.replace(final)
    return manifest


def build_pharmacophore_index(artifact_catalog_path: Path, output_root: Path,
                               bin_width: float = 0.5,
                               max_distance: float = 20.0) -> dict:
    """Index every artifact shard independently; existing matching shards are reused."""
    catalog_path = artifact_catalog_path.resolve()
    catalog = json.loads(catalog_path.read_text(encoding="utf-8"))
    if catalog.get("format") != "aidd-conformer-artifact-catalog":
        raise ValueError("not an AIDD conformer artifact catalog")
    output_root.mkdir(parents=True, exist_ok=True)
    shard_records = []
    for shard in catalog["shards"]:
        target = output_root / "shards" / shard["name"]
        built = build_shard_pharmacophore_index(
            Path(shard["path"]), target, bin_width=bin_width,
            max_distance=max_distance)
        manifest_path = target / "manifest.json"
        shard_records.append({
            "name": shard["name"], "path": str(target.resolve()),
            "manifest_sha256": _sha256(manifest_path),
            "global_id_start": built["global_id_start"],
            "conformers": built["conformers"],
        })
    result = {
        "format": INDEX_FORMAT + "-catalog", "version": INDEX_VERSION,
        "library_id": catalog["library_id"],
        "artifact_catalog": str(catalog_path),
        "artifact_catalog_sha256": _sha256(catalog_path),
        "bin_width_angstrom": float(bin_width),
        "max_distance_angstrom": float(max_distance),
        "conformers": int(catalog["conformers"]), "shards": shard_records,
    }
    (output_root / "catalog.json").write_text(
        json.dumps(result, indent=2), encoding="utf-8")
    return result


def _feature_type_id(name: str) -> int | None:
    return FEATURE_ALIASES.get(name.replace("_", "").replace("-", "").upper())


def compile_pharmacophore_query(query_manifest_path: Path, output_path: Path,
                                 profiles: Sequence[SearchProfile] = DEFAULT_PROFILES) -> dict:
    """Compile a small reusable co-crystal query; no library artifact is read."""
    source_path = query_manifest_path.resolve()
    source = json.loads(source_path.read_text(encoding="utf-8"))
    anchors = []
    for anchor in source.get("anchors", []):
        type_id = _feature_type_id(str(anchor["feature_type"]))
        if type_id is None:
            continue
        point = np.asarray(anchor["atom_center"], dtype=np.float64)
        if point.shape != (3,) or not np.isfinite(point).all():
            raise ValueError(f"invalid atom center for anchor {anchor.get('anchor_id')}")
        anchors.append((str(anchor["anchor_id"]), type_id, point))
    pairs = []
    for left in range(len(anchors)):
        for right in range(left + 1, len(anchors)):
            first, second = anchors[left], anchors[right]
            pairs.append({
                "pair_id": f"{first[0]}::{second[0]}",
                "anchor_ids": [first[0], second[0]],
                "feature_types": sorted([first[1], second[1]]),
                "distance_angstrom": float(np.linalg.norm(first[2] - second[2])),
            })
    compiled = {
        "format": QUERY_FORMAT, "version": 1,
        "query_id": source.get("query_id"),
        "source_manifest": str(source_path),
        "source_manifest_sha256": _sha256(source_path),
        "supported_anchors": len(anchors), "pairs": pairs,
        "profiles": [asdict(profile) for profile in profiles],
    }
    # Cache identity is content/configuration based and independent of where the
    # query manifest happens to be mounted on a workstation.
    compiled["query_hash"] = _canonical_json_hash({
        "format": compiled["format"], "version": compiled["version"],
        "query_id": compiled["query_id"],
        "source_manifest_sha256": compiled["source_manifest_sha256"],
        "supported_anchors": compiled["supported_anchors"],
        "pairs": compiled["pairs"], "profiles": compiled["profiles"],
    })
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(compiled, indent=2), encoding="utf-8")
    return compiled


def _keys_for_pair(pair: Mapping, tolerance: float, bin_width: float,
                   max_distance: float) -> list[int]:
    distance = float(pair["distance_angstrom"])
    low, high = max(0.0, distance - tolerance), min(max_distance, distance + tolerance)
    if low > high:
        return []
    first, second = pair["feature_types"]
    return [encode_pair_key(first, second, value)
            for value in range(_distance_bin(low, bin_width),
                               _distance_bin(high, bin_width) + 1)]


def _shard_hits(shard_path: Path, pairs: Sequence[Mapping], tolerance: float,
                bin_width: float, max_distance: float) -> tuple[np.ndarray, np.ndarray]:
    manifest = json.loads((shard_path / "manifest.json").read_text(encoding="utf-8"))
    key_files = {int(row["key"]): row["file"] for row in manifest["keys"]}
    ids, pair_ids = [], []
    for pair_index, pair in enumerate(pairs):
        pair_hits = []
        for key in _keys_for_pair(pair, tolerance, bin_width, max_distance):
            relative = key_files.get(key)
            if relative:
                pair_hits.append(np.memmap(shard_path / relative, dtype="<i8", mode="r"))
        if pair_hits:
            unique = np.unique(np.concatenate(pair_hits))
            ids.append(unique)
            pair_ids.append(np.full(len(unique), pair_index, dtype=np.int32))
    if not ids:
        return np.empty(0, dtype=np.int64), np.empty(0, dtype=np.int32)
    return np.concatenate(ids), np.concatenate(pair_ids)


def _profile_matches(index_catalog: Mapping, pairs: Sequence[Mapping],
                     profile: SearchProfile) -> tuple[np.ndarray, np.ndarray]:
    all_ids, all_pairs = [], []
    for shard in index_catalog["shards"]:
        ids, pair_ids = _shard_hits(
            Path(shard["path"]), pairs, profile.tolerance_angstrom,
            float(index_catalog["bin_width_angstrom"]),
            float(index_catalog["max_distance_angstrom"]))
        if len(ids):
            all_ids.append(ids); all_pairs.append(pair_ids)
    if not all_ids:
        return np.empty(0, dtype=np.int64), np.empty(0, dtype=np.uint16)
    ids, pair_ids = np.concatenate(all_ids), np.concatenate(all_pairs)
    order = np.lexsort((pair_ids, ids)); ids, pair_ids = ids[order], pair_ids[order]
    keep = np.ones(len(ids), dtype=bool)
    keep[1:] = (ids[1:] != ids[:-1]) | (pair_ids[1:] != pair_ids[:-1])
    ids = ids[keep]
    unique_ids, counts = np.unique(ids, return_counts=True)
    required = min(profile.required_pairs, len(pairs))
    accepted = counts >= required
    return unique_ids[accepted], counts[accepted].astype(np.uint16)


def search_pharmacophore_index(index_catalog_path: Path, query_plan_path: Path,
                                output_path: Path,
                                external_l1_ids: Iterable[int] = ()) -> dict:
    """Return nested evidence tiers plus an anchor-safe union with external L1."""
    index_path, query_path = index_catalog_path.resolve(), query_plan_path.resolve()
    index = json.loads(index_path.read_text(encoding="utf-8"))
    query = json.loads(query_path.read_text(encoding="utf-8"))
    if index.get("format") != INDEX_FORMAT + "-catalog":
        raise ValueError("not a pharmacophore index catalog")
    if query.get("format") != QUERY_FORMAT:
        raise ValueError("not a compiled pharmacophore query")
    pairs = query["pairs"]
    profiles = tuple(SearchProfile(**row) for row in query["profiles"])
    profile_results: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for profile in profiles:
        profile_results[profile.name] = _profile_matches(index, pairs, profile)

    # Enforce the scientific invariant even if a caller supplies custom profiles.
    if all(name in profile_results for name in ("loose", "balanced", "strict")):
        loose = set(profile_results["loose"][0].tolist())
        balanced = set(profile_results["balanced"][0].tolist())
        strict = set(profile_results["strict"][0].tolist())
        if not strict <= balanced <= loose:
            raise ValueError("search profiles are not nested")

    external = np.unique(np.fromiter(external_l1_ids, dtype=np.int64))
    pharma = np.unique(np.concatenate([
        ids for ids, _ in profile_results.values()
    ])) if profile_results else np.empty(0, dtype=np.int64)
    union = np.union1d(external, pharma)
    source_flags = np.zeros(len(union), dtype=np.uint8)
    source_flags[np.isin(union, pharma)] |= 1
    source_flags[np.isin(union, external)] |= 2
    tier = np.zeros(len(union), dtype=np.uint8)
    tier_names = ("loose", "balanced", "strict")
    for value, name in enumerate(tier_names, start=1):
        if name in profile_results:
            tier[np.isin(union, profile_results[name][0])] = value
    matched = np.zeros((len(union), len(profiles)), dtype=np.uint16)
    positions = {int(gid): idx for idx, gid in enumerate(union)}
    for column, profile in enumerate(profiles):
        ids, counts = profile_results[profile.name]
        for gid, count in zip(ids, counts):
            matched[positions[int(gid)], column] = count

    output_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output_path, global_ids=union, source_flags=source_flags,
                        highest_tier=tier, matched_pair_counts=matched)
    result_path = output_path if output_path.suffix == ".npz" else output_path.with_suffix(
        output_path.suffix + ".npz")
    manifest = {
        "format": RESULT_FORMAT, "version": 1,
        "library_id": index["library_id"], "query_id": query.get("query_id"),
        "index_catalog": str(index_path), "index_catalog_sha256": _sha256(index_path),
        "query_plan": str(query_path), "query_plan_sha256": _sha256(query_path),
        "result_file": str(result_path.resolve()), "result_sha256": _sha256(result_path),
        "profiles": query["profiles"],
        "profile_counts": {name: int(len(rows[0])) for name, rows in profile_results.items()},
        "external_l1_count": int(len(external)),
        "pharmacophore_count": int(len(pharma)), "union_count": int(len(union)),
        "source_flags": {"pharmacophore": 1, "external_l1": 2},
        "highest_tier": {"external_only": 0, "loose": 1, "balanced": 2, "strict": 3},
    }
    manifest_path = result_path.with_suffix(".manifest.json")
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return manifest


def load_external_l1_ids(path: Path | None) -> np.ndarray:
    """Load a baseline candidate set without coupling this index to FAISS format."""
    if path is None:
        return np.empty(0, dtype=np.int64)
    suffix = path.suffix.lower()
    if suffix == ".npy":
        values = np.load(path)
    elif suffix == ".npz":
        archive = np.load(path)
        name = next((candidate for candidate in ("global_ids", "ids")
                     if candidate in archive.files), None)
        if name is None:
            raise ValueError("NPZ L1 candidates require a global_ids or ids array")
        values = archive[name]
    elif suffix == ".json":
        payload = json.loads(path.read_text(encoding="utf-8"))
        if isinstance(payload, dict):
            payload = payload.get("global_ids", payload.get("ids"))
        if payload is None:
            raise ValueError("JSON L1 candidates require global_ids or ids")
        values = payload
    else:
        values = np.loadtxt(path, dtype=np.int64, ndmin=1)
    result = np.asarray(values, dtype=np.int64).reshape(-1)
    if np.any(result < 0):
        raise ValueError("external L1 candidate IDs must be non-negative")
    return np.unique(result)
