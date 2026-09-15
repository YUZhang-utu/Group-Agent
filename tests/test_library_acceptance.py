import json
from pathlib import Path
import sqlite3
from types import SimpleNamespace

import numpy as np
import pytest

from aidd_agent import library_acceptance as ev
from aidd_agent.registry import initialize


def dump(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    ev.write(path, value)


def manifest(path, base):
    payload = dict(base, files={p.name: dict(bytes=p.stat().st_size, sha256=ev.sha(p))
                               for p in path.iterdir() if p.is_file() and p.name != 'manifest.json'})
    dump(path/'manifest.json', payload)
    return ev.sha(path/'manifest.json')


@pytest.fixture
def library(tmp_path):
    batch = tmp_path/'batch'; batch.mkdir()
    initialize(batch/'registry.sqlite3')
    rng = np.random.default_rng(71)
    raw = rng.normal(size=(12, 60)).astype(np.float32)
    # Equal-distance boundary and same molecule across a shard boundary.
    raw[8] = raw[7]
    molecules = np.asarray([f'MOL-{i//2:012d}' for i in range(12)], dtype='S16')
    molecules[6] = molecules[0]
    conformers = np.asarray([f'CNF-{i:012d}' for i in range(12)], dtype='S16')
    rows = {k: [] for k in ('artifacts', 'chemical', 'pharmacophore')}
    sources = []
    with sqlite3.connect(batch/'registry.sqlite3') as db:
        db.execute("INSERT INTO library VALUES ('L','test','now')")
        db.execute('CREATE TABLE batch_source(path TEXT PRIMARY KEY,sha TEXT,records INTEGER,inserted INTEGER,duplicates INTEGER)')
        for mid in np.unique(molecules):
            db.execute('INSERT INTO molecule VALUES (?,?,?,?)', (mid.decode(), 'L', mid.decode(), 'now'))
        for si, start in enumerate((0, 6)):
            source = str(tmp_path/f's{si}.mol2'); sources.append(dict(path=source))
            db.execute('INSERT INTO batch_source VALUES (?,?,?,?,?)', (source, f'hash{si}', 6, 6, 0))
            for i in range(start, start+6):
                db.execute('INSERT INTO conformer VALUES (?,?,?,?,?,?,?,?,?,?,?,?)',
                           (conformers[i].decode(), molecules[i].decode(), i, f'record{i}', 2, 1,
                            f'content{i}', 'topology', source, i-start, '[]', 'now'))
            for kind in rows:
                p = batch/kind/f's{si}'; p.mkdir(parents=True)
                base = dict(library_id='L', conformers=6, global_id_start=start,
                            source_path=source, source_sha256=f'hash{si}')
                if kind != 'pharmacophore':
                    molecules[start:start+6].tofile(p/'molecule_ids.bin')
                    conformers[start:start+6].tofile(p/'conformer_ids.bin')
                    dtype = ev.META_DTYPE if kind == 'artifacts' else ev.CHEM_META_DTYPE
                    meta = np.zeros(6, dtype=dtype); meta['global_id'] = np.arange(start, start+6)
                    meta.tofile(p/('meta.bin' if kind == 'artifacts' else 'chem-meta.bin'))
                if kind == 'artifacts':
                    raw[start:start+6].tofile(p/'usrcat.f32.bin')
                else:
                    key = 'artifact_v1_manifest_sha256' if kind == 'chemical' else 'source_manifest_sha256'
                    base[key] = rows['artifacts'][-1]['manifest_sha256']
                h = manifest(p, base)
                rows[kind].append(dict(name=f's{si}', path=str(p), manifest_sha256=h,
                                       conformers=6, global_id_start=start))
    for kind in rows:
        cat = dict(library_id='L', conformers=12, shards=rows[kind])
        if kind != 'artifacts':
            key = 'artifact_v1_catalog_sha256' if kind == 'chemical' else 'artifact_catalog_sha256'
            cat[key] = ev.sha(batch/'artifacts'/'catalog.json')
        dump(batch/kind/'catalog.json', cat)
    dump(batch/'batch.json', dict(files=sources))
    dump(batch/'COMPLETE.json', dict(conformers=12, shards=2, batch_sha256=ev.sha(batch/'batch.json')))
    f = batch/'faiss'; trained = f/'trained'; trained.mkdir(parents=True)
    (f/'index.faiss').write_bytes(b'fake index for integrity tests only')
    np.savez(f/'transform.npz', mean=np.zeros(60, dtype=np.float32), std=np.ones(60, dtype=np.float32))
    (trained/'transform.npz').write_bytes((f/'transform.npz').read_bytes())
    (trained/'template.faiss').write_bytes(b'fake template')
    manifest(trained, dict(catalog_sha256=ev.sha(batch/'artifacts'/'catalog.json')))
    dump(f/'manifest.json', dict(final_ntotal=12, catalog_sha256=ev.sha(batch/'artifacts'/'catalog.json'),
                                  index_sha256=ev.sha(f/'index.faiss')))
    return batch, raw, molecules


def test_acceptance_and_no_library_mutation(library):
    batch, _, _ = library
    before = {str(p): ev.sha(p) for p in batch.rglob('*') if p.is_file()}
    catalog, report = ev.accept_library(batch)
    assert report['status'] == 'passed' and report['library_conformers'] == 12
    assert report['library_molecules'] == 6
    assert before == {str(p): ev.sha(p) for p in batch.rglob('*') if p.is_file()}


def test_missing_complete_is_not_registration_success(tmp_path):
    with pytest.raises(ValueError, match='registration alone'):
        ev.accept_library(tmp_path)


def test_corrupt_payload_rejected(library):
    batch, _, _ = library
    path = batch/'artifacts/s0/usrcat.f32.bin'
    data = bytearray(path.read_bytes()); data[0] ^= 1; path.write_bytes(data)
    with pytest.raises(ValueError, match='Checksum mismatch'):
        ev.accept_library(batch)
    _, report = ev.accept_library(batch, full=False)
    assert report['status'] == 'metadata_passed_hashes_not_checked'


def test_registry_missing_source_and_bad_range_rejected(library):
    batch, _, _ = library
    with sqlite3.connect(batch/'registry.sqlite3') as db:
        db.execute('DELETE FROM batch_source WHERE path=(SELECT path FROM batch_source LIMIT 1)')
    with pytest.raises(ValueError, match='frozen sources'):
        ev.accept_library(batch)


def test_registry_identity_corruption_is_detected(library):
    batch, _, _ = library
    with sqlite3.connect(batch/'registry.sqlite3') as db:
        db.execute("UPDATE conformer SET id='CNF-CORRUPTED' WHERE id=(SELECT id FROM conformer LIMIT 1)")
    with pytest.raises(ValueError, match='identity mismatch'):
        ev.accept_library(batch)


def test_immutable_old_manifest_paths_can_differ_from_registered_source(library):
    batch, _, _ = library
    # Move only the logical source paths in the frozen inventory and registry;
    # manifests remain immutable, as when E032 restores old shards.
    inventory = ev.read(batch/'batch.json')
    with sqlite3.connect(batch/'registry.sqlite3') as db:
        for row in inventory['files']:
            old = row['path']; new = str(Path(old).parent/'relocated'/Path(old).name)
            db.execute('UPDATE batch_source SET path=? WHERE path=?', (new, old))
            db.execute('UPDATE conformer SET source_path=? WHERE source_path=?', (new, old))
            row['path'] = new
    dump(batch/'batch.json', inventory)
    complete = ev.read(batch/'COMPLETE.json'); complete['batch_sha256'] = ev.sha(batch/'batch.json')
    dump(batch/'COMPLETE.json', complete)
    assert ev.accept_library(batch)[1]['status'] == 'passed'


def test_failed_gate_does_not_recommend_a_fast_configuration():
    row = dict(requested_candidate_budget=1000, nprobe=64, recall={'top1000_strict': .7},
               warm_search_seconds=[.001], warm_search_plus_single_fetch_rerank_seconds=.002,
               retained_molecules=50, molecule_reduction_fraction=.99)
    summary, selected = ev.summarize([row], .95, 1000)
    assert selected is None and not summary[0]['gate_passed']


def test_streaming_exact_matches_dense_excludes_whole_molecule_and_ties(library):
    batch, raw, mols = library
    corpus = ev.Corpus(ev.read(batch/'artifacts/catalog.json'))
    q = raw[[0, 7]]
    exclusions = [bytes(mols[0]), b'']
    for chunk in (1, 4, 100):
        result, excluded, _ = ev.exact_truth(corpus, q, exclusions, np.zeros(60), np.ones(60), 5, chunk)
        assert excluded == [3, 0]
        for qi in range(2):
            delta = raw-q[qi]; dist = np.einsum('ij,ij->i', delta, delta)
            if exclusions[qi]:
                dist[mols == exclusions[qi]] = np.inf
            expected = np.lexsort((np.arange(12), dist))[:5]
            assert np.array_equal(result[qi][1], expected)
            assert np.allclose(result[qi][0], dist[expected])
    found, vectors = corpus.fetch([11, 0, 6, 0], vectors=True)
    assert np.array_equal(found, mols[[11, 0, 6, 0]])
    assert np.array_equal(vectors, raw[[11, 0, 6, 0]])


def test_topk_stable_edge_and_caps_are_molecule_level():
    distances = np.asarray([1., 1., 0., 1., np.inf])
    ids = np.asarray([9, 1, 3, 0, 4])
    assert ev.topk(distances, ids, 3)[1].tolist() == [3, 0, 1]
    chosen, mols = ev.unique_representatives(np.asarray([8, 2, 4, 5]), np.asarray([b'A', b'A', b'B', b'C']), 2)
    assert chosen.tolist() == [8, 4] and mols.tolist() == [b'A', b'B']


class ExactIndex:
    def __init__(self, vectors):
        self.vectors = vectors; self.ntotal = len(vectors); self.d = 60; self.nlist = 2; self.nprobe = 1

    def search(self, queries, k):
        distances = ((queries[:, None]-self.vectors[None])**2).sum(axis=2)
        ids = np.argsort(distances, axis=1, kind='stable')[:, :k]
        return np.take_along_axis(distances, ids, axis=1), ids


def test_search_denominator_exclusion_and_candidate_mapping(library):
    batch, raw, mols = library
    corpus = ev.Corpus(ev.read(batch/'artifacts/catalog.json'))
    mean, std = np.zeros(60, dtype=np.float32), np.ones(60, dtype=np.float32)
    truth, excluded, _ = ev.exact_truth(corpus, raw[[0]], [bytes(mols[0])], mean, std, 1000, 4)
    result, arrays = ev.evaluate_setting(ExactIndex(raw), corpus, raw[0], bytes(mols[0]), excluded[0], truth[0],
                                        mean, std, 2, 1000, 2, [100, 1000], [2, 100], 6)
    assert result['retained_conformers'] == 9
    assert result['recall']['top1000_strict'] == 1
    assert bytes(mols[0]) not in arrays['molecule_ids']
    assert result['cap_scenarios'][0]['retained_molecules'] == 2
    assert result['cap_scenarios'][1]['exact_truth_molecule_coverage'] == 1


def test_full_report_pipeline_without_faiss_dependency(library, tmp_path):
    batch, raw, _ = library
    fake = SimpleNamespace(__version__='test-exact', omp_set_num_threads=lambda _: None,
                           read_index=lambda _: ExactIndex(raw))
    output = tmp_path/'assessment'; output.mkdir()
    args = SimpleNamespace(metadata_only=False, threads=1, queries=None, seed=71, query_count=2,
                           nprobes=[1, 2], budgets=[3, 12], truth_ks=[2, 1000], caps=[2, 12], repeats=1,
                           chunk_size=4, recall_gate=.95)
    ev._run_locked(args, batch, output, fake)
    report = ev.read(output/'report.json')
    assert report['calibration_status'] == 'gate_passed_on_panel'
    assert report['provisional_setting']['candidate_budget'] == 12
    assert len(report['results']) == 8
    assert (output/'metrics.csv').is_file() and (output/'report.md').is_file()
    assert ev.read(output/'EVALUATION_COMPLETE.json')['report_sha256'] == ev.sha(output/'report.json')


def test_real_faiss_ivf_search(library):
    faiss = pytest.importorskip('faiss')
    batch, raw, mols = library
    faiss.omp_set_num_threads(1)
    index = faiss.IndexIVFFlat(faiss.IndexFlatL2(60), 60, 1)
    index.train(raw); index.add_with_ids(raw, np.arange(len(raw), dtype=np.int64))
    corpus = ev.Corpus(ev.read(batch/'artifacts/catalog.json'))
    mean, std = np.zeros(60, dtype=np.float32), np.ones(60, dtype=np.float32)
    truth, exclusions, _ = ev.exact_truth(corpus, raw[[0]], [bytes(mols[0])], mean, std, 1000, 4)
    result, _ = ev.evaluate_setting(index, corpus, raw[0], bytes(mols[0]), exclusions[0], truth[0],
                                    mean, std, 1, 1000, 1, [1000], [10], 6)
    assert result['recall']['top1000_strict'] == 1
