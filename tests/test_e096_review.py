import csv
import importlib.util
import json
from pathlib import Path

spec=importlib.util.spec_from_file_location('review_e096',Path(__file__).resolve().parents[1]/'scripts/review_e096_results.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


def fixture(root):
    b=root/'blocks';b.mkdir();p=root/'profiles';p.mkdir()
    population=24530851;base,remainder=divmod(population,384)
    counts=[base+(i<remainder) for i in range(384)]
    with (b/'blocks.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['block_id','parent_block','conformers'])
        for i,n in enumerate(counts):w.writerow([f'b{i}',f'p{i}',n])
    (b/'transform.json').write_text(json.dumps(dict(mean=[0]*317,std=[1]*317,weights=[1]*317)))
    (b/'split-decisions.json').write_text('[]')
    (p/'report.json').write_text(json.dumps(dict(status='complete',coverage=1,width=317,total_regular=population)))
    families=['Donor','Acceptor','Aromatic','Hydrophobe','LumpedHydrophobe','PosIonizable','NegIonizable','SideHeavyAtoms']
    report=dict(status='complete',membership_gate='passed_against_complete_cached_profiles',regular_conformers=population,special_conformers=132885,
        regular_work_blocks=384,parent_regular_blocks=384,smallest_regular_block=min(counts),largest_regular_block=max(counts),blocks_below_minimum=0,
        check_grouping='source_molecule_id_crc32_mod5',active_channels_by_group=dict(broad_chemistry=33,local_sterics=36,backbone_geometry=24,typed_spatial=224),
        mean_residue_feature_counts=dict.fromkeys(families,0),decision_counts={},output_hashes={n:module.sha(b/n) for n in ['blocks.csv','transform.json','split-decisions.json']},sources={str((p/'report.json').resolve()):module.sha(p/'report.json')})
    (b/'report.json').write_text(json.dumps(report))


def test_metadata_only_review_detects_corruption(tmp_path):
    fixture(tmp_path)
    r=module.review(tmp_path)
    assert not r['errors'] and r['regular_blocks']==384
    assert r['warnings'] and not list(tmp_path.rglob('*.sqlite'))
    with (tmp_path/'blocks/blocks.csv').open('a') as f:f.write('duplicate,p0,4\n')
    r=module.review(tmp_path)
    assert 'Hash mismatch: blocks.csv' in r['errors']
    assert not r['checks']['bounded_fragmentation']


def test_malformed_transform_cannot_pass(tmp_path):
    import pytest
    fixture(tmp_path)
    (tmp_path/'blocks/transform.json').write_text(json.dumps(dict(mean=[0],std=[1],weights=[1])))
    with pytest.raises(ValueError,match='317-channel'):module.review(tmp_path)
