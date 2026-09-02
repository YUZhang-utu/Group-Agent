import numpy as np
import pytest

from aidd_agent.conformer_artifacts import META_DTYPE, _quantize


def test_coordinate_quantization_has_centangstrom_precision():
    origin=np.asarray([10.,20.,30.],dtype=np.float32)
    points=np.asarray([[10.004,20.006,31.234]],dtype=np.float32)
    encoded=_quantize(points,origin)
    decoded=encoded.astype(np.float32)/100+origin
    assert np.max(np.abs(decoded-points)) <= 0.0051


def test_coordinate_quantization_rejects_overflow():
    with pytest.raises(ValueError,match="int16"):
        _quantize(np.asarray([[400.,0.,0.]],dtype=np.float32),np.zeros(3,dtype=np.float32))


def test_meta_schema_contains_incremental_global_id_and_offsets():
    assert META_DTYPE["global_id"].itemsize == 8
    assert META_DTYPE["coord_offset"].itemsize == 8
    assert META_DTYPE["feature_counts"].shape == (6,)
