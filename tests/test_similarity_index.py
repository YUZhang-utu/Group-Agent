from __future__ import annotations

import numpy as np
import pytest

from aidd_agent.similarity import query_usrcat_index, save_usrcat_index


def test_usrcat_index_round_trip(tmp_path):
    path = tmp_path / "shape.npz"
    save_usrcat_index([("near", np.zeros(60)), ("far", np.ones(60))], path)
    hits = query_usrcat_index(np.zeros(60), path)
    assert [hit.molecule_id for hit in hits] == ["near", "far"]
    assert hits[0].score == 1.0


def test_usrcat_index_validates_shape(tmp_path):
    with pytest.raises(ValueError, match="60"):
        save_usrcat_index([("bad", [1, 2])], tmp_path / "bad.npz")
