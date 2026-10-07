"""Verified complete-block reuse; independent arm ranking and stopping remain local."""
import hashlib
import json
from pathlib import Path
import shutil

import numpy as np

from .expanded_wee1 import fingerprint
from .final_work_blocks import read
from .joint_spatial_profiles import save
from .prompt_workflow import file_lock
from .screening_selection import check_hashes


def refine_cached(batch, output, design, ids, settings, refine):
    config = design['_block_cache']
    query = {k: v for k, v in design.items() if k != '_block_cache'}
    identity = dict(version=1, query=query, inputs=config['identity'],
                    ids_sha256=hashlib.sha256(np.asarray(ids, dtype='<i8').tobytes()).hexdigest(),
                    count=len(ids), chunk=settings['refine_chunk'])
    key = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()
    root = Path(config['root'])
    root.mkdir(parents=True, exist_ok=True)
    cache = root / key
    with file_lock(root / (key + '.lock')):
        cache.mkdir(exist_ok=True)
        seal = cache / 'complete.json'
        reused = seal.exists()
        if reused:
            record = read(seal)
            if record['identity'] != identity:
                raise ValueError('Block cache identity mismatch')
            check_hashes(record['files'])
        else:
            refine(batch, cache, query, ids, settings['workers'], settings['refine_chunk'])
            # Validate every scheduled chunk before publishing a reusable seal.
            from .budget_screen import verified_chunk
            for ti in range(len(query['templates'])):
                for start in range(0, len(ids), settings['refine_chunk']):
                    target = cache/'chunks'/f'{ti:02d}'/f'{start:010d}.npz'
                    if not verified_chunk(target, ids[start:start+settings['refine_chunk']]):
                        raise ValueError('Incomplete block cache')
            files = sorted(p for p in (cache/'chunks').rglob('*') if p.is_file())
            save(seal, dict(identity=identity, files=fingerprint(files)))
        # Local receipts must hash local copies, not just their cached originals.
        for source in sorted((cache/'chunks').glob('*/*.receipt.json')):
            record = read(source)
            check_hashes(record['files'])
            target = output/'chunks'/source.parent.name/source.name
            target.parent.mkdir(parents=True, exist_ok=True)
            copied = []
            for original in record['files']:
                original = Path(original)
                if original.parent.resolve() != source.parent.resolve():
                    raise ValueError('Block cache artifact outside chunk directory')
                destination = target.parent/original.name
                shutil.copyfile(original, destination)
                copied.append(destination)
            record['files'] = fingerprint(copied)
            save(target, record)
        save(output/'cache-reuse.json', dict(key=key, reused=reused,
             timing_scope='copied_worker_times_are_source_compute_not_current_wall_time'))
    return reused
