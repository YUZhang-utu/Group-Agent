"""Isolated Python 3.9 EquiScore worker; only consumes coordinator-created jobs."""
import argparse
from contextlib import closing
from functools import lru_cache
import hashlib
import io
import importlib.metadata
import json
from pathlib import Path
import pickle
import sqlite3
import sys
import time


WEIGHT_SHA = 'd4367bb73686b2363e238abb778fab55e2924458ec1ced561072bd82f711695d'


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''): h.update(block)
    return h.hexdigest()


def load_model(profile):
    repo = Path(profile['repository']).resolve()
    weights = repo / 'workdir/official_weight/save_model_screen.pt'
    if sha(weights) != WEIGHT_SHA:
        raise ValueError('Expected the pinned official screening checkpoint')
    sys.path.insert(0, str(repo))
    import torch
    import numpy as np
    import random
    # Match upstream import order: dataset and utils have a circular dependency.
    import utils.utils
    from utils.parsing import parse_train_args
    from model.equiscore import EquiScore
    argv = sys.argv
    try:
        sys.argv = ['equiscore_worker', '--test']
        args = parse_train_args()
    finally:
        sys.argv = argv
    if not torch.cuda.is_available():
        raise RuntimeError('CUDA is unavailable in the isolated EquiScore environment')
    args.local_rank = 0
    args.N_atom_features = 39 if args.FP else 28
    random.seed(42); np.random.seed(42); torch.manual_seed(42)
    torch.cuda.manual_seed_all(42)
    model = EquiScore(args)
    checkpoint = torch.load(str(weights), map_location='cpu')
    model.load_state_dict(checkpoint['model'], strict=True)
    model = model.to('cuda:0').eval()
    return model, args


@lru_cache(maxsize=4)
def mol2_records(filename):
    records, lines = [], []
    with open(filename, encoding='utf-8') as stream:
        for line in stream:
            if line.startswith('@<TRIPOS>MOLECULE') and lines:
                records.append(''.join(lines)); lines = []
            if lines or line.startswith('@<TRIPOS>MOLECULE'): lines.append(line)
    if lines: records.append(''.join(lines))
    return records


def ligand_from_pose(row, root):
    import numpy as np
    from rdkit import Chem
    path = Path(row['pose_file'])
    if not path.is_absolute(): path = root / path
    path = path.resolve()
    path.relative_to(root.resolve())
    index = int(row['pose_index'])
    if index < 0: raise ValueError('Pose index must be zero-based and nonnegative')
    raw = mol2_records(str(path))[index]
    if raw.splitlines()[1].strip() != row['pose_name']:
        raise ValueError('Pose name/index mismatch')
    ligand = Chem.MolFromMol2Block(raw, sanitize=True, removeHs=False)
    if ligand is None:
        raise ValueError('Pose cannot be read without chemical repair')
    original = next((r for r in mol2_records(str(path.parent.parent / 'ligands.mol2'))
                     if r.splitlines()[1].strip() == row['alias']), None)
    if original is None: raise ValueError('Original docking-input ligand not found')
    source = Chem.MolFromMol2Block(original, sanitize=True, removeHs=False)
    if source is None: raise ValueError('Original ligand cannot be read without repair')
    for mol in (ligand, source):
        Chem.AssignStereochemistryFrom3D(mol)
    if Chem.MolToSmiles(Chem.RemoveHs(ligand), isomericSmiles=True) != Chem.MolToSmiles(Chem.RemoveHs(source), isomericSmiles=True):
        raise ValueError('Pose chemistry/stereochemistry differs from the docking input')
    heavy_positions = np.array([ligand.GetConformer().GetAtomPosition(a.GetIdx())
                               for a in ligand.GetAtoms() if a.GetAtomicNum() > 1])
    ligand = Chem.RemoveHs(ligand)
    if not np.allclose(heavy_positions, ligand.GetConformer().GetPositions(), atol=1e-7, rtol=0):
        raise ValueError('Heavy-atom coordinates changed during H removal')
    return ligand


def receptor_structure(pdb, mol2):
    import numpy as np
    from Bio.PDB import PDBParser
    from scipy.spatial import cKDTree
    structure = PDBParser(QUIET=True).get_structure('receptor', pdb)
    pdb_atoms = [a for a in structure.get_atoms() if a.element.upper() not in {'H', 'D'}]
    rows, active = [], False
    for line in mol2_records(str(mol2))[0].splitlines():
        if line == '@<TRIPOS>ATOM': active = True; continue
        if line.startswith('@<TRIPOS>'): active = False
        if active and line.strip():
            f = line.split(); element = f[5].split('.')[0].upper()
            if element not in {'H', 'D'}: rows.append((element, [float(v) for v in f[2:5]]))
    if len(pdb_atoms) != len(rows):
        raise ValueError('PDB and docked receptor have different heavy-atom counts')
    for element in {e for e, _ in rows}:
        a = np.array([v for e, v in rows if e == element])
        b = np.array([v.coord for v in pdb_atoms if v.element.upper() == element])
        if len(a) != len(b): raise ValueError('Receptor element inventory differs')
        distances, indices = cKDTree(b).query(a)
        if max(distances) > 0.03 or len(set(indices)) != len(indices):
            raise ValueError('Receptor PDB is not in the saved PLANTS coordinate frame')
    return structure


def environment():
    import torch
    return dict(python=sys.version, packages={dist.metadata['Name']: dist.version
                for dist in importlib.metadata.distributions()}, gpu=torch.cuda.get_device_name(0),
                cuda=torch.version.cuda)


def verify_pose_inputs(row, root, hashes, verified):
    path = Path(row['pose_file'])
    path = (path if path.is_absolute() else root / path).resolve()
    for item in (path, path.parent.parent / 'ligands.mol2'):
        relative = item.relative_to(root).as_posix()
        if relative not in hashes: raise ValueError('Pose/input is not sealed in the docking report')
        if relative not in verified:
            if sha(item) != hashes[relative]: raise ValueError('Saved docking pose/input changed')
            verified.add(relative)


def make_pair(ligand, structure, filename):
    import numpy as np
    from Bio.PDB import PDBIO, Select
    from rdkit import Chem
    from scipy.spatial import cKDTree
    from utils.ifp_construct import get_nonBond_pair
    tree = cKDTree(ligand.GetConformer().GetPositions())
    class Pocket(Select):
        def accept_residue(self, residue):
            points = np.array([a.coord for a in residue.get_atoms()])
            return bool(len(points) and min(tree.query(points)[0]) < 8.0)
    buffer = io.StringIO(); writer = PDBIO(); writer.set_structure(structure)
    writer.save(buffer, Pocket())
    pocket = Chem.MolFromPDBBlock(buffer.getvalue(), sanitize=True, removeHs=True)
    if pocket is None or pocket.GetNumAtoms() == 0:
        raise ValueError('Pocket cannot be parsed; no atoms/metals are silently removed')
    if ligand.GetNumAtoms() + pocket.GetNumAtoms() > 2000:
        raise ValueError('Complex exceeds the declared 2000-atom graph budget')
    # Precompute strictly. Upstream otherwise catches fingerprint errors silently.
    atom_pairs, interaction_types = get_nonBond_pair(ligand, pocket)
    with open(filename, 'wb') as stream:
        pickle.dump((ligand, pocket, atom_pairs, interaction_types), stream, pickle.HIGHEST_PROTOCOL)


def run(job_path):
    import csv
    job_path = Path(job_path).resolve()
    job = json.loads(job_path.read_text())
    for name, digest in job['input_hashes'].items():
        if sha(name) != digest: raise ValueError('Worker input changed: ' + name)
    model, args = load_model(job['profile'])
    import torch
    from dataset.dataset import ESDataset
    structures = {key: receptor_structure(r['pdb'], r['mol2']) for key, r in job['receptors'].items()}
    output = job_path.parent
    receipt = output / 'environment.json'
    current_env = environment()
    if receipt.exists() and json.loads(receipt.read_text()) != current_env:
        raise ValueError('Environment changed; use a fresh output and repeat the pilot')
    receipt.write_text(json.dumps(current_env, indent=2) + '\n', encoding='utf-8')
    if job.get('pilot_environment') is not None and job['pilot_environment'] != current_env:
        raise ValueError('Full scoring environment differs from the successful pilot')
    verified = set()
    selected = set(job['selected_cids']) if job['selected_cids'] is not None else None
    started = time.monotonic()
    with closing(sqlite3.connect(output / 'predictions.sqlite')) as db:
        db.execute('CREATE TABLE IF NOT EXISTS scores(cid TEXT,receptor TEXT,status TEXT,score REAL,error TEXT,PRIMARY KEY(cid,receptor))')
        done = set(db.execute('SELECT cid,receptor FROM scores'))
        with open(job['scores'], encoding='utf-8', newline='') as stream:
            for row in csv.DictReader(stream):
                if selected is not None and row['cid'] not in selected: continue
                key = row['cid'], row['receptor']
                if key in done: continue
                value, error, status = None, '', 'unsupported_input'
                pair_file = output / 'current-pair.pkl'
                try:
                    if row['status'] != 'ok': raise ValueError('Source docking pair is not scored')
                    verify_pose_inputs(row, Path(job['source_root']), job['source_hashes'], verified)
                    ligand = ligand_from_pose(row, Path(job['source_root']))
                    make_pair(ligand, structures[row['receptor']], pair_file)
                    dataset = ESDataset([str(pair_file)], args, str(output), False)
                    graph = dataset[0]
                    if graph is None: raise ValueError('Upstream graph construction rejected this pose')
                    g, full_g, _ = dataset.collate([graph])
                    status = 'model_failed'
                    with torch.no_grad():
                        pred = model(g.to('cuda:0'), full_g.to('cuda:0'))
                        if tuple(pred.shape) != (1, 2): raise ValueError('Unexpected screening head shape')
                        value = float(torch.softmax(pred, dim=-1)[0, 1].cpu())
                    if not 0 <= value <= 1: raise ValueError('Nonfinite/out-of-range EquiScore probability')
                    status = 'ok'
                except Exception as exc:
                    value = None; error = type(exc).__name__ + ': ' + str(exc)
                    if isinstance(exc, RuntimeError) and 'out of memory' in str(exc).lower(): torch.cuda.empty_cache()
                finally:
                    if pair_file.exists(): pair_file.unlink()
                db.execute('INSERT INTO scores VALUES(?,?,?,?,?)', (*key, status, value, error))
                db.commit(); done.add(key)
                if len(done) % 100 == 0:
                    print('Scored or recorded failure for %d pairs; %.1f seconds' % (len(done), time.monotonic()-started), flush=True)
    print('EquiScore worker completed; %.1f seconds' % (time.monotonic()-started), flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--job')
    parser.add_argument('--doctor', help='Profile JSON to check without scientific inference')
    args = parser.parse_args()
    if args.doctor:
        profile = json.loads(Path(args.doctor).read_text())
        load_model(profile)
        import torch, dgl, rdkit, prolif
        import dgl.function as fn
        # Exercise CUDA kernels, not just driver discovery or checkpoint loading.
        graph = dgl.graph(([0, 1], [1, 0])).to('cuda:0')
        graph.ndata['x'] = torch.ones((2, 8), device='cuda:0')
        graph.update_all(fn.copy_u('x', 'm'), fn.sum('m', 'y'))
        assert bool(torch.isfinite(graph.ndata['y']).all())
        assert bool(torch.isfinite(torch.ones((8, 8), device='cuda:0') @ torch.ones((8, 8), device='cuda:0')).all())
        torch.cuda.synchronize()
        print(json.dumps(dict(status='ready', checkpoint_loaded_strictly=True,
            cuda_graph_smoke_passed=True,
            torch=torch.__version__, dgl=dgl.__version__, rdkit=rdkit.__version__,
            prolif=prolif.__version__, gpu=torch.cuda.get_device_name(0))))
    elif args.job: run(args.job)
    else: parser.error('Choose --doctor or --job')


if __name__ == '__main__': main()
