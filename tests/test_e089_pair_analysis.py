import csv
import importlib.util
from pathlib import Path
import json


def test_target_control_subtraction_and_query_weighting(tmp_path):
    spec=importlib.util.spec_from_file_location('pair_analysis',Path(__file__).parents[1]/'scripts/analyze_e089_pairs.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    rows=[]
    def add(kind,q,ref,g,value):
        rows.append(dict(kind=kind,band='control' if kind=='within_definite_control' else 'up_to_35',
            query=q,reference=ref,target_group=g,rotation=0,maximum_deviation=34,same_molecule=False,
            **dict.fromkeys(module.METRICS,value)))
    add('within_definite_control','r1','r2','g1',10)
    add('within_definite_control','s1','s2','g2',100)
    for index,value in enumerate([11,12,13]):add('boundary_vs_reference','q1',f'r{index}','g1',value)
    add('boundary_vs_reference','q2','s1','g2',98)
    source=tmp_path/'pairs.csv'
    with source.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    module.run(source,tmp_path/'out')
    result=json.loads((tmp_path/'out/report.json').read_text())
    for stats in result['bands']['up_to_35']['equal_query_matched_deltas'].values():
        assert stats['median']==0 and stats['n']==2
    assert result['bands']['up_to_35']['queries_with_any_observed_envelope_candidate']==1
