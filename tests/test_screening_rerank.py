import json
import numpy as np

from aidd_agent.screening_rerank import (
    AnchorDefinition, CandidateScores, PoseScores, QueryManifest,
    ranking_diagnostic, rerank, retain_l1_candidates,
)


def pose(objective: str, atom: float, projected: float, anchored: float) -> PoseScores:
    return PoseScores(objective, tuple(np.eye(4).ravel()), .5, .5, .1, .1,
                      anchored, atom, projected, (anchored,), (1.0,))


def row(conf_id: int, atom: float, projected: float, anchored: float) -> CandidateScores:
    return CandidateScores("M", conf_id, (pose("atomcentered_joint", atom, projected, anchored),
                                           pose("projected_joint", atom / 2, projected * 2, anchored)))


def test_l1_retention_has_no_anchor_or_feature_input():
    assert retain_l1_candidates([7, 8, -1], [.2, .1, 0], 2).tolist() == [8, 7]


def test_objective_specific_reranking_preserves_set_and_uses_named_pose():
    rows = [row(1, .4, .1, .2), row(2, .1, .8, .9)]
    for score in ("shape", "unweighted", "atomcentered", "projected", "max_color", "anchored"):
        assert {r.conf_id for r in rerank(rows, "atomcentered_joint", score)} == {1, 2}
    assert rerank(rows, "atomcentered_joint", "atomcentered")[0].conf_id == 1
    assert rerank(rows, "projected_joint", "projected")[0].conf_id == 2


def test_query_manifest_labels_raw_anchor_overlaps(tmp_path):
    anchor = AnchorDefinition("A0:HBA:12", (12,), "HBA", (1., 2., 3.), ((2., 2., 3.),),
                              {"interaction": "hydrogen_bond", "distance_angstrom": 2.8,
                               "angle_degrees": 165., "burial": 0.7}, 1.5)
    manifest = QueryManifest("8BJU:QT9:A:601", "projected-color-v1",
                             {"anchor": 1.5, "ordinary": 1., "solvent_exposed": .5}, (anchor,))
    path = manifest.write(tmp_path)
    assert json.loads(path.read_text())["anchors"][0]["anchor_id"] == "A0:HBA:12"
    labelled = row(1, .4, .1, .7).labelled_anchor_overlaps(manifest, "atomcentered_joint")
    assert labelled["A0:HBA:12"] == {"cross_overlap_raw": .7, "query_self_overlap_raw": 1.}


def test_anchor_sensitivity_diagnostic():
    result = ranking_diagnostic([1, 2, 3, 4], [1, 2, 4, 3], 2)
    assert result["top_2_overlap"] == 1.0
    assert 0 < result["spearman_rho"] < 1
