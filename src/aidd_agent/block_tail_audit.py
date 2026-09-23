"""Summarize sampled block tails without selecting or discarding library blocks."""
import argparse
import csv
import json
import math


def summarize(rows,k=5,threshold=None):
    if type(k) is not int or k<=0:raise ValueError('Positive top-k required')
    if threshold is not None and not math.isfinite(threshold):raise ValueError('Nonfinite threshold')
    blocks={}
    for row in rows:
        block,mid=row['block_id'],row['molecule_id'];score=float(row['score'])
        if not block or not mid or not math.isfinite(score):raise ValueError('Invalid sampled block row')
        members=blocks.setdefault(block,{})
        members[mid]=max(members.get(mid,-math.inf),score)
    result=[]
    for block,members in sorted(blocks.items()):
        ranked=sorted(members.items(),key=lambda item:(-item[1],item[0]))
        scores=[score for _,score in ranked];n=len(scores)
        result.append(dict(block_id=block,sampled_molecules=n,maximum=scores[0],
            top_k_requested=k,top_k_used=min(k,n),top_k_mean=sum(scores[:k])/min(k,n),
            top_molecule_ids=[mid for mid,_ in ranked[:k]],
            high_score_fraction=sum(s>=threshold for s in scores)/n if threshold is not None else None))
    return dict(blocks=result,threshold=threshold,
        comparable_sample_sizes=len({r['sampled_molecules'] for r in result})<=1,
        scope='Sampled unique-molecule tails only; unequal samples and search budgets are not directly comparable; no block rejection or recall guarantee')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--csv',required=True);p.add_argument('--top-k',type=int,default=5)
    p.add_argument('--threshold',type=float);args=p.parse_args()
    with open(args.csv,newline='',encoding='utf-8') as stream:
        result=summarize(csv.DictReader(stream),args.top_k,args.threshold)
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
