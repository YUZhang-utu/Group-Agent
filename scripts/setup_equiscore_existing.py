"""Reuse working workstation Torch/DGL through a non-mutating venv overlay."""
import argparse
import importlib.metadata
import json
import os
from pathlib import Path
import subprocess
import sys
import venv


def package_constraints():
    values = {}
    for dist in importlib.metadata.distributions():
        name = dist.metadata['Name']
        if name:
            values[name.lower().replace('_', '-')] = dist.version
    return [f'{key}=={value}' for key, value in sorted(values.items())]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prefix', default='/mnt/local/hand/yuzhang/aidd/tools/equiscore-dl4s')
    parser.add_argument('--source-profile', default='/mnt/local/hand/yuzhang/aidd/tools/equiscore/profile.json')
    args = parser.parse_args()
    if os.name != 'posix': raise RuntimeError('This setup targets the Linux workstation')
    root = Path(__file__).resolve().parents[1]
    worker = root / 'scripts/equiscore_worker.py'
    env = os.environ.copy(); env['DGLBACKEND'] = 'pytorch'
    # This gate runs before creating files or downloading any packages.
    print('Checking the existing Torch and DGL CUDA kernels first.', flush=True)
    subprocess.run([sys.executable, str(worker), '--gpu-check'], check=True, env=env)
    profile = json.loads(Path(args.source_profile).read_text())
    prefix = Path(args.prefix).absolute()
    if prefix == Path('/'): raise ValueError('Use a dedicated setup prefix')
    marker = prefix / 'base-python.json'
    base = dict(python=os.path.abspath(sys.executable), version=sys.version)
    if prefix.exists() and (not marker.is_file() or json.loads(marker.read_text()) != base):
        raise ValueError('Existing unmanaged or different-base prefix; choose a fresh --prefix')
    prefix.mkdir(parents=True, exist_ok=True)
    marker.write_text(json.dumps(base, indent=2) + '\n')
    constraints = prefix / 'base-constraints.txt'
    constraints.write_text('\n'.join(package_constraints()) + '\n')
    venv_path = prefix / 'env'
    python = venv_path / 'bin/python'
    if not python.exists():
        venv.EnvBuilder(system_site_packages=True, with_pip=True).create(venv_path)
    # Pin every visible base distribution. Pip can add missing dependencies only;
    # conflicts fail instead of replacing Torch, DGL, NumPy or RDKit.
    subprocess.run([str(python), '-m', 'pip', '--isolated', 'install', '--constraint', str(constraints),
                    'prolif==1.1.0', 'lmdb', 'prefetch-generator', 'biopython', 'scikit-learn'], check=True, env=env)
    profile['python'] = str(python)  # Do not resolve the venv's interpreter symlink.
    target = prefix / 'profile.json'
    if target.exists() and json.loads(target.read_text()) != profile:
        raise ValueError('Existing profile differs; choose a fresh setup prefix')
    target.write_text(json.dumps(profile, indent=2) + '\n')
    with (prefix / 'doctor.log').open('w') as log:
        result = subprocess.run([str(python), str(worker), '--doctor', str(target)], stdout=log, stderr=subprocess.STDOUT, env=env)
    print((prefix / 'doctor.log').read_text(), flush=True)
    if result.returncode:
        raise RuntimeError('Model/chemistry preflight failed; inspect ' + str(prefix / 'doctor.log'))
    print('Ready for a new pilot. EQUISCORE_PROFILE=' + str(target))


if __name__ == '__main__': main()
