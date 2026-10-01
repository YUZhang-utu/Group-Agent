"""Regression coverage for explainable replay of rejected property splits."""
import importlib.util
from pathlib import Path
import numpy as np

spec=importlib.util.spec_from_file_location('diagnose_e094',Path(__file__).parents[1]/'scripts/diagnose_e094_properties.py')
diagnostic=importlib.util.module_from_spec(spec)
spec.loader.exec_module(diagnostic)


def test_distinguishes_population_and_constant_features():
    assert diagnostic.root_probe(np.zeros((8,69)),5,.35)['reason']=='population_floor'
    assert diagnostic.root_probe(np.zeros((20,69)),5,.35)['reason']=='no_distinct_projection_cut_meeting_floor'


def test_single_channel_dilution_is_reported():
    x=np.zeros((20,69),dtype=np.float32)
    x[10:,0]=1
    r=diagnostic.root_probe(x,5,.35)
    assert r['reason']=='below_contrast_cutoff'
    assert r['largest_channel_difference']==1
    assert np.isclose(r['centroid_contrast'],1/np.sqrt(69))


def test_matches_existing_partition_decision():
    from aidd_agent.work_block_properties import partition
    for amplitude in (.1,2):
        x=np.zeros((20,69),dtype=np.float32)
        x[10:]=amplitude
        leaves,_=partition(x,minimum=5,maximum_children=2,minimum_contrast=.35)
        assert (len(leaves)==2)==(diagnostic.root_probe(x,5,.35)['reason']=='root_split_would_pass')
