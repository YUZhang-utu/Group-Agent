"""E115 grouping, kernel attribution and feasibility regression tests."""
import time
import unittest
import sys
import tempfile
from pathlib import Path
from unittest.mock import patch
from types import SimpleNamespace
import numpy as np
from aidd_agent import contact_region_review as review
from aidd_agent.gaussian_overlay import score_overlay


class ReviewTests(unittest.TestCase):
    def fixture(self):
        points = np.array([[0., 0, 0], [1.5, 0, 0], [3., 0, 0]])
        data = dict(feature_points=points, feature_types=np.full(3, 5),
                    feature_directions=np.zeros((3, 3)), feature_direction_kinds=np.zeros(3, dtype=int))
        anchors = [dict(feature_index=i, target_residue=93, protein_atom="CG1",
                        feature_class="hydrophobic") for i in range(3)]
        return dict(aligned=data, anchors=anchors)

    def test_complete_link_no_chain_merging_and_shared_label_uniqueness(self):
        state = self.fixture()
        self.assertEqual([g["reference_features"] for g in review.reference_groups(state)], [[0, 1], [2]])
        state["anchors"].append(dict(state["anchors"][0], target_residue=61))
        groups = review.reference_groups(state)
        self.assertEqual(sorted(i for g in groups for i in g["reference_features"]), [0, 1, 2])
        self.assertEqual(len(groups[0]["labels"]), 2)

    def test_groups_cannot_reuse_candidate_and_ignore_extra_features(self):
        state = self.fixture()
        cf = SimpleNamespace(feature_points=np.array([[1.5, 0., 0.]]), feature_types=np.array([5]),
                             feature_directions=np.zeros((1, 3)), feature_kinds=np.array([0]))
        groups = review.reference_groups(state)
        before = review.grouped_score(state, cf, np.eye(4), groups)
        self.assertEqual(sum(r["candidate_feature"] >= 0 for r in before["groups"]), 1)
        cf.feature_points = np.vstack([cf.feature_points, [100., 0, 0]])
        cf.feature_types = np.array([5, 1])
        cf.feature_directions = np.zeros((2, 3)); cf.feature_kinds = np.zeros(2, dtype=int)
        self.assertEqual(review.grouped_score(state, cf, np.eye(4), groups)["mean"], before["mean"])

    def test_direction_prevents_invalid_merge(self):
        state = self.fixture()
        state["aligned"]["feature_direction_kinds"][:] = 1
        state["aligned"]["feature_directions"][:] = [[1, 0, 0], [-1, 0, 0], [1, 0, 0]]
        self.assertEqual(len(review.reference_groups(state)), 3)

    def test_color_decomposition_replays_kernel(self):
        rng = np.random.default_rng(8)
        a, b = rng.normal(size=(7, 3)), rng.normal(size=(9, 3))
        at, bt = np.arange(7) % 3, np.arange(9) % 3
        result = review.color_attribution(a, at, b, bt)
        expected = score_overlay(a, b, a, at, b, bt, np.eye(4).ravel(), cutoff=None)
        self.assertAlmostEqual(result["tanimoto"], expected["color"]["tanimoto"])
        self.assertAlmostEqual(sum(r["cross"] for r in result["candidate_features"]), result["cross"])

    def test_feasible_original_and_policy_positive_survive_better_invalid_pose(self):
        def evaluate(matrix):
            x = matrix[0, 3]
            return dict(shape=float(x), original_policy_pass=abs(x) < 1e-12,
                        flags=dict(contacts=True, steric=x <= .1, envelope=True, regions=True))
        original = np.eye(4); invalid = np.eye(4); invalid[0, 3] = 2
        result = review.refine_feasible(np.zeros((1, 3)), [original, invalid],
                                       evaluate, time.monotonic()+10, max_trials=24)
        self.assertTrue(review.feasible(result["selected"]["score"]))
        self.assertTrue(result["original_policy_preserved"])
        self.assertEqual(len([r for r in result["archive"] if r["origin"] == "endpoint"]), 2)
        self.assertLessEqual(result["trial_evaluations"], 24)

    def test_expired_budget_still_retains_endpoints(self):
        score = dict(shape=.5, original_policy_pass=False,
                     flags=dict(contacts=False, steric=True, envelope=True, regions=True))
        result = review.refine_feasible(np.zeros((1, 3)), [np.eye(4)], lambda m: score, 0)
        self.assertIsNone(result["selected"])
        self.assertTrue(result["budget_exhausted"])
        self.assertEqual(result["trial_evaluations"], 0)
        self.assertEqual(len(result["archive"]), 1)

    def test_polar_sasa_deduplicates_shared_atom_and_detects_occlusion(self):
        chemistry = SimpleNamespace(atomic_numbers=np.array([8]), feature_members=([0], [0]))
        result = review.polar_exposure(np.zeros((1, 3)), chemistry, np.array([1, 2]),
                                       np.array([[2., 0, 0]]), ["C"])
        self.assertEqual(len(result["atoms"]), 1)
        self.assertGreater(result["isolated_polar_sasa"], result["bound_polar_sasa"])

    def test_primary_export_frame_and_serialization(self):
        sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
        import review_contact_regions as runner
        from scipy.spatial import cKDTree
        points = np.array([[0., 0, 0]])
        atom = dict(label_alt_id="", occupancy="1", pdbx_PDB_ins_code="",
                    xyz=[3., 0, 0], auth_atom_id="O", auth_comp_id="LEU",
                    auth_seq_id="54", auth_asym_id="B", canonical_residue=54, type_symbol="O")
        state = dict(record=dict(native_mmcif="fixture.cif", chain="B", transform=np.eye(4)),
                     tree=cKDTree([[3., 0, 0]]),
                     aligned=dict(feature_points=points, feature_types=np.array([1])))
        chemistry = SimpleNamespace(atomic_numbers=np.array([7]), feature_members=([0],),
            molecule_id="fixture", global_id=1, conformer_id="c1", formal_charges=[0], bonds=[])
        candidate = SimpleNamespace(shape_points=points, feature_points=points,
                                    feature_types=np.array([1]), molecule_id="fixture")
        with tempfile.TemporaryDirectory() as directory, patch(
                "aidd_agent.consensus_admission.read_structure",
                return_value=dict(chains={"B": [atom]})):
            out = Path(directory)
            runner.export_primary(out, candidate, chemistry, state, dict(target={}),
                                  dict(selected=np.eye(4)))
            result = runner.base.read(out / "primary_review.json")
            self.assertEqual(result["poses"]["selected"]["color"]["tanimoto"], 1.0)
            self.assertEqual(len(result["poses"]["selected"][
                "complementary_polar_neighbors_distance_only"]), 1)
            self.assertTrue((out / "primary_selected.sdf").exists())
            self.assertEqual(runner.base.read(out / "primary_receptor_mapping.json")[0]["native_chain"], "B")


if __name__ == "__main__":
    unittest.main()
