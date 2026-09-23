"""Rebuild fresh consensus artifacts from local evidence with networking disabled."""
import argparse
import json
from pathlib import Path
import shutil
from unittest.mock import patch

from aidd_agent.consensus_model import build
from aidd_agent.expanded_wee1 import fingerprint


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,required=True)
    parser.add_argument('--cached-public-data',type=Path,required=True)
    parser.add_argument('--reference-query',required=True)
    parser.add_argument('--target-chain',required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    shutil.copytree(args.cached_public_data,args.output/'raw')
    with patch('requests.sessions.Session.request',side_effect=RuntimeError('Offline validation: missing cached public data')):
        result=build(args.source,args.output,args.reference_query,args.target_chain)
    receipt=dict(scope='Offline consensus rebuild; no live LLM or library execution',
        readiness=result['readiness'],prepared=len(result.get('prepared_complexes',[])),
        anchors=len(result.get('anchors',[])),templates=len(result.get('templates',[])),
        failures=result.get('failures',[]),contact_evidence=result.get('contact_evidence'),
        implementation=fingerprint([Path(__file__),*Path('src/aidd_agent').glob('*.py')]))
    (args.output/'rebuild_receipt.json').write_text(json.dumps(receipt,indent=2),encoding='utf-8')
    print(json.dumps(receipt,indent=2)[:2500])
    if receipt['failures'] or receipt['readiness']!='proposal_ready':raise SystemExit(1)


if __name__=='__main__':main()
