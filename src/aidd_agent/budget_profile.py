"""Read completed budget receipts without changing or rescoring an active run."""
import argparse
from collections import Counter
from itertools import islice
import json
from pathlib import Path


def summarize(run,limit=200):
    run=Path(run)
    if limit<1:raise ValueError('Positive receipt limit required')
    protocol=json.loads((run/'protocol.json').read_text())
    seconds=Counter();counts=Counter();receipts=0
    for path in islice((run/'chunks').glob('*/*.receipt.json'),limit):
        row=json.loads(path.read_text());seconds.update(row.get('worker_seconds',{}))
        counts.update(row.get('counts',{}));receipts+=1
    total=sum(seconds.values())
    return dict(receipts_sampled=receipts,configuration={k:protocol.get(k) for k in
        ('workers','chunk_conformers','assignment_backend','budget','ranking_policy')},
        worker_seconds=dict(seconds),worker_fraction={k:v/total for k,v in seconds.items()} if total else {},
        counts=dict(counts),scope='First available completed receipts; not a random sample, full run or wall-time benchmark')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run',type=Path,required=True);parser.add_argument('--receipts',type=int,default=200)
    args=parser.parse_args();print(json.dumps(summarize(args.run,args.receipts),indent=2))


if __name__=='__main__':main()
