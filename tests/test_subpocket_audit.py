"""Geometry checks for the standalone exploratory audit, not screening scoring."""
import importlib.util
from pathlib import Path

import numpy as np

spec = importlib.util.spec_from_file_location('subpocket_audit',
    Path(__file__).resolve().parents[1] / 'scripts' / 'audit_mdm2_subpockets.py')
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)


def test_shared_boundary_is_ambiguous_not_two_occupied_sites():
    probes = [np.array([[0., 0., 0.]]), np.array([[2., 0., 0.]]), np.array([[0., 8., 0.]])]
    labels, _ = audit.spatial_assign(np.array([[1., 0., 0.], [0., 0., 0.], [9., 9., 9.]]), probes, 2.)
    assert labels == ['ambiguous', 'A_Trp23', 'outside']


def test_spatial_labels_are_invariant_under_joint_rigid_transform():
    probes = [np.array([[0., 0., 0.]]), np.array([[4., 0., 0.]]), np.array([[0., 5., 0.]])]
    points = np.array([[.1, 0., 0.], [4., 0., .2], [0., 5., .3]])
    matrix = np.array([[0., -1., 0., 17.], [1., 0., 0., -9.], [0., 0., 1., 3.], [0., 0., 0., 1.]])
    before, distances = audit.spatial_assign(points, probes, 2.)
    after, transformed = audit.spatial_assign(audit.move(points, matrix),
        [audit.move(p, matrix) for p in probes], 2.)
    assert before == after == list(audit.SITES)
    np.testing.assert_allclose(distances, transformed, atol=1e-12)


def test_low_occupancy_and_alternate_atoms_are_not_silent_absence():
    good = dict(label_alt_id='', occupancy='1.0', pdbx_PDB_ins_code='')
    assert audit.observed(good)
    assert not audit.observed(dict(good, occupancy='0.5'))
    assert not audit.observed(dict(good, label_alt_id='A'))
    assert not audit.observed(dict(good, pdbx_PDB_ins_code='A'))
