"""Inspect stored omega angles and propose chirality-preserving definite-class candidates."""
import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import itertools
import json
import math
from pathlib import Path
import sqlite3

VERSION='verified-peptide-local-frame-v1'
WIDTHS={'backbone':6,'chemistry':17,'typed':59}


def digest(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':')).encode()).hexdigest()


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def rotate(values,offset):return values[offset:]+values[:offset]


def recover_chirality(groups,variant,maximum_trials=2000000):
    """Verify finite R/S/achiral preimages; unsupported labels remain unresolved."""
    patterns=defaultdict(set)
    for g,states in groups.items():patterns[tuple(states)].add(g)
    recovered={};trials=0
    for number,(states,targets) in enumerate(sorted(patterns.items()),1):
        needed=3**len(states)
        if trials+needed>maximum_trials:continue
        for labels in itertools.product(('R','S','achiral'),repeat=len(states)):
            trials+=1
            value=digest([VERSION,variant,len(states),list(labels),list(states)])
            if value in targets:recovered[value]=list(labels)
        if number%32==0:print(f'Chirality verification: {number}/{len(patterns)} patterns',flush=True)
    return recovered,trials


def omega_angles(row,variant):
    states=row['provenance']['omega_states'];vector=row['descriptor'];width=WIDTHS[variant]
    if len(vector)!=len(states)*width:raise ValueError('Descriptor dimension mismatch')
    angles=[]
    for i,state in enumerate(states):
        sine=float(vector[i*width+2]);cosine=float(vector[i*width+5])
        if not math.isfinite(sine+cosine) or abs(sine*sine+cosine*cosine-1)>1e-6:
            raise ValueError('Invalid stored omega sin/cos')
        value=abs(math.degrees(math.atan2(sine,cosine)))
        observed='cis' if value<=30 else ('trans' if value>=150 else 'boundary')
        # Allow only floating-point roundtrip error at an exact classification edge.
        if observed!=state and min(abs(value-30),abs(value-150))>1e-7:
            raise ValueError('Stored omega state disagrees with angle')
        angles.append(value)
    return angles


def candidates(states,angles,chirality,variant,definite_groups):
    if chirality is None:return [],None
    uncertain=[i for i,s in enumerate(states) if s=='boundary']
    deviation=max((min(angles[i],180-angles[i]) for i in uncertain),default=0)
    if any(abs(angles[i]-90)<1e-7 for i in uncertain):return [],deviation
    assigned=[('cis' if a<90 else 'trans') if s=='boundary' else s for s,a in zip(states,angles)]
    result=[]
    for offset in range(len(states)):
        group=digest([VERSION,variant,len(states),rotate(chirality,offset),rotate(assigned,offset)])
        if group in definite_groups:result.append(dict(hard_group=group,rotation=offset))
    return result,deviation


def run(build,diagnostic,output,thresholds=(35,45,60,75)):
    build=Path(build).resolve();diagnostic=Path(diagnostic).resolve();output=Path(output).resolve()
    if any(output==p or p in output.parents for p in (build,diagnostic)):
        raise ValueError('Use a fresh output outside input trees')
    if not thresholds or any(not 30<t<90 for t in thresholds):raise ValueError('Angular thresholds must lie strictly between 30 and 90 degrees')
    receipt=json.loads((build/'report.json').read_text());audit=json.loads((diagnostic/'report.json').read_text())
    variant=receipt['schema']['variant']
    if receipt['status']!='complete' or receipt['schema']['version']!=VERSION or variant not in WIDTHS:
        raise ValueError('Unsupported or incomplete descriptor build')
    if Path(audit['descriptors']).resolve()!=build/'descriptors.sqlite':raise ValueError('Diagnostic belongs to another descriptor database')
    model=Path(audit['model'])
    if sha(model/'report.json')!=audit['provenance']['model_report_sha256']:raise ValueError('Model report changed')
    for name,value in audit['output_hashes'].items():
        if sha(diagnostic/name)!=value:raise ValueError('Diagnostic artifact changed')
    with (diagnostic/'hard_groups.csv').open(newline='') as stream:rows=list(csv.DictReader(stream))
    groups={r['hard_group']:r['omega_pattern'].split(';') for r in rows}
    definite={g for g,s in groups.items() if 'boundary' not in s}
    boundary={r['hard_group']:int(r['conformers']) for r in rows if 'boundary' in groups[r['hard_group']]}
    if len(groups)!=audit['hard_groups'] or sum(boundary.values())!=audit['boundary_attribution'].get('with_boundary',{}).get('conformers',0):
        raise ValueError('Diagnostic group counts mismatch')
    chirality,trials=recover_chirality(groups,variant)
    print(f'Recovered {len(chirality)}/{len(groups)} chirality signatures using {trials} hash checks',flush=True)
    output.mkdir(parents=True,exist_ok=False)
    counts=Counter();hist=Counter();sensitivity=Counter();positions=Counter()
    with sqlite3.connect((build/'descriptors.sqlite').as_uri()+'?mode=ro',uri=True) as source, \
         sqlite3.connect((model/'blocks.sqlite').resolve().as_uri()+'?mode=ro',uri=True) as blocks, \
         sqlite3.connect(output/'candidates.sqlite') as dest:
        dest.execute('CREATE TABLE candidate(cid TEXT PRIMARY KEY,original_group TEXT,maximum_deviation REAL,candidates TEXT,reason TEXT)')
        for group,expected in sorted(boundary.items()):
            checked=0
            for (cid,) in blocks.execute('SELECT cid FROM point WHERE group_id=?',(group,)):
                value=source.execute('SELECT payload FROM descriptor WHERE cid=?',(cid,)).fetchone()
                if value is None:raise ValueError('Missing source descriptor')
                row=json.loads(value[0]);states=row['provenance']['omega_states']
                if row['hard_group']!=group or states!=groups[group]:raise ValueError('Group/state changed')
                angles=omega_angles(row,variant);matches,deviation=candidates(states,angles,chirality.get(group),variant,definite)
                if group not in chirality:reason='chirality_unresolved'
                elif any(s=='boundary' and abs(a-90)<1e-7 for s,a in zip(states,angles)):reason='angular_midpoint'
                elif not matches:reason='no_existing_compatible_definite_class'
                else:reason='angular_candidates_require_property_and_geometry_review'
                counts[reason]+=1;positions[states.count('boundary')]+=1
                for s,a in zip(states,angles):
                    if s=='boundary':hist[min(175,int(a//5)*5)]+=1
                for t in thresholds:
                    if matches and deviation<=t:sensitivity[t]+=1
                dest.execute('INSERT INTO candidate VALUES(?,?,?,?,?)',(cid,group,deviation,json.dumps(matches),reason))
                counts['checked']+=1;checked+=1
                if counts['checked']%10000==0:
                    dest.commit();print(f"Checked {counts['checked']} boundary conformers",flush=True)
            if checked!=expected:raise ValueError('Boundary class population mismatch')
        dest.commit()
        if dest.execute('SELECT count(*) FROM candidate').fetchone()[0]!=sum(boundary.values()):raise ValueError('Boundary coverage mismatch')
    total=audit['conformers'];base=total-counts['checked']
    result=dict(status='complete',readiness='candidate_review_only_not_final_membership',variant=variant,
        total_conformers=total,definite_conformers=base,boundary_conformers=counts['checked'],
        recovered_chirality_groups=len(chirality),unresolved_chirality_groups=len(groups)-len(chirality),
        counts=dict(counts),boundary_position_counts=dict(positions),boundary_angle_histogram_5deg=dict(sorted(hist.items())),
        sensitivity=[dict(maximum_deviation_from_cis_or_trans=t,angular_candidate_conformers=sensitivity[t],
            minimum_remaining_special_after_angular_gate=counts['checked']-sensitivity[t],
            upper_bound_regular_coverage=(base+sensitivity[t])/total) for t in thresholds],
        sibling_property_builds={v:(build.parent/v/'report.json').exists() for v in ('chemistry','typed')},
        sources={str(build/'report.json'):sha(build/'report.json'),str(diagnostic/'report.json'):sha(diagnostic/'report.json'),
                 str(diagnostic/'hard_groups.csv'):sha(diagnostic/'hard_groups.csv')},
        output_hashes={'candidates.sqlite':sha(output/'candidates.sqlite')},implementation_sha256=sha(__file__),
        limitations=['Coverage is an angular candidate upper bound, not validated final block membership or search recall',
            'Only original definite cis/trans classes are candidate targets; no union through uncertain nodes',
            'Chirality is verified by original class hash using R/S/achiral labels; unresolved classes stay special',
            'Directed cyclic rotations only; exact midpoint states are not assigned',
            'Large input hashes rely on the existing full validator; source coordinates are not reparsed',
            'Chemical property similarity and side-chain steric compatibility must be evaluated before adoption'])
    (output/'report.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('build','diagnostic','output'):p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();print(json.dumps(run(a.build,a.diagnostic,a.output),indent=2))


if __name__=='__main__':main()
