"""Package the reviewed MEDCHEM update without data, credentials or old deliveries."""
import hashlib
import json
from pathlib import Path
import subprocess
import zipfile

ROOT = Path(__file__).resolve().parents[1]
NEW = [
    'src/aidd_agent/block_results.py', 'src/aidd_agent/web/med.png',
    'src/aidd_agent/web/peptide-1.svg', 'src/aidd_agent/web/peptide-2.svg',
    'src/aidd_agent/web/peptide-3.svg', 'src/aidd_agent/web/peptides.json',
    'tests/test_block_results.py', 'scripts/build_chat_molecules.py',
    'scripts/qa_medchem_server.py', 'scripts/qa_medchem_ui.cjs',
    'scripts/build_e098_delivery.py', 'experiments/E098-medchem-chat.md',
    'to_human/20261005_MEDCHEM_CHAT_HANDOFF.md',
    'src/aidd_agent/web/peptide-4.svg', 'src/aidd_agent/web/peptide-5.svg',
    'src/aidd_agent/web/peptide-6.svg', 'src/aidd_agent/web/peptide-d.svg',
    'src/aidd_agent/web/solvent-field.svg', 'to_human/20261005_MEDCHEM_VISUAL_REFRESH.md',
]


def digest(raw):
    try: raw = raw.decode('utf-8-sig').replace('\r\n', '\n').encode('utf-8')
    except UnicodeError: pass
    return hashlib.sha256(raw).hexdigest()


INSTALLER = '''"""Apply only hash-matched source updates; never touch scientific data."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import time
import uuid

def digest(raw):
    try: raw = raw.decode('utf-8-sig').replace('\\r\\n', '\\n').encode('utf-8')
    except UnicodeError: pass
    return hashlib.sha256(raw).hexdigest()

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--repo', required=True, type=Path)
p.add_argument('--apply', action='store_true', help='Apply after stopping the Chat server; default is check only')
args = p.parse_args()
root = args.repo.resolve()
package = Path(__file__).resolve().parent
if not (root / 'src/aidd_agent').is_dir() or not (root / 'pyproject.toml').is_file():
    raise SystemExit('Choose the existing Group-Agent repository, not the data workspace')
manifest = json.loads((package / 'manifest.json').read_text(encoding='utf-8'))
changes = []
for row in manifest['files']:
    rel = Path(row['path'])
    target = (root / rel).resolve()
    payload = (package / 'files' / rel).resolve()
    if not target.is_relative_to(root) or not payload.is_relative_to(package / 'files'):
        raise SystemExit('Package path escapes its root')
    raw = payload.read_bytes()
    if hashlib.sha256(raw).hexdigest() != row['sha256']:
        raise SystemExit('Package checksum mismatch: ' + str(rel))
    current = digest(target.read_bytes()) if target.exists() else None
    if current == row['new_normalized_sha256']: continue
    if current not in row.get('accepted_previous_sha256', [row['old_normalized_sha256']]):
        raise SystemExit('Existing local edit or different version; no files changed: ' + str(rel))
    changes.append((rel, target, raw))
print(str(len(changes)) + ' source files ready; project data and credentials are excluded.')
if not args.apply:
    print('Check only. Stop Chat and rerun with --apply to install.')
    raise SystemExit(0)
backup = root / '.medchem-updates' / (time.strftime('%Y%m%d-%H%M%S') + '-' + uuid.uuid4().hex[:8])
for rel, target, raw in changes:
    if target.exists():
        saved = backup / rel
        saved.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(target, saved)
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_name(target.name + '.e098-update')
    temporary.write_bytes(raw)
    temporary.replace(target)
print('Update applied. Restart the same Chat workspace. Backup: ' + str(backup))
'''


def main():
    changed = subprocess.check_output(['git', 'diff', '--name-only', 'HEAD'], cwd=ROOT).decode().splitlines()
    names = sorted(set(changed + NEW))
    manifest = dict(base_commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT).decode().strip(), files=[])
    baseline = {}
    for name in ('MEDCHEM_E098_UPDATE.zip', 'MEDCHEM_E099_UI_UPDATE.zip'):
        previous = ROOT / 'to_human' / name
        if previous.is_file():
            with zipfile.ZipFile(previous) as archive:
                for row in json.loads(archive.read('manifest.json'))['files']:
                    baseline.setdefault(row['path'], []).append(row['new_normalized_sha256'])
    output = ROOT / 'to_human/MEDCHEM_E100_TYPOGRAPHY_UPDATE.zip'
    with zipfile.ZipFile(output, 'w', zipfile.ZIP_DEFLATED) as archive:
        for name in names:
            raw = (ROOT / name).read_bytes()
            prior = subprocess.run(['git', 'show', 'HEAD:' + name], cwd=ROOT, capture_output=True)
            manifest['files'].append(dict(path=name, sha256=hashlib.sha256(raw).hexdigest(),
                                         new_normalized_sha256=digest(raw), old_normalized_sha256=digest(prior.stdout) if prior.returncode == 0 else None))
            accepted = [manifest['files'][-1]['old_normalized_sha256']]
            if name in baseline: accepted.extend(baseline[name])
            manifest['files'][-1]['accepted_previous_sha256'] = list(dict.fromkeys(accepted))
            archive.writestr('files/' + name, raw)
        archive.writestr('manifest.json', json.dumps(manifest, indent=2))
        archive.writestr('install_e098.py', INSTALLER)
        archive.writestr('install_medchem.py', INSTALLER)
        for name in ('desktop.png', 'masthead.png', 'mobile.png', 'results.png'):
            preview = ROOT / 'data/e100-ui' / name
            if preview.is_file(): archive.write(preview, 'previews/' + name)
        archive.writestr('README.txt', 'MEDCHEM Agent source update. No data or credentials included.\n'
                          'Stop Chat after confirming no scientific tasks are active.\n'
                          'Check: python install_medchem.py --repo /path/to/Group-Agent\n'
                          'Apply: python install_medchem.py --repo /path/to/Group-Agent --apply\n'
                          'See files/to_human/20261005_MEDCHEM_VISUAL_REFRESH.md for this revision.\n')
    print(json.dumps(dict(package=str(output), files=len(names), sha256=sha256_file(output))))


def sha256_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


if __name__ == '__main__':
    main()
