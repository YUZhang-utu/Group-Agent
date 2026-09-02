import numpy as np

from aidd_agent.faiss_benchmark import candidate_recall, synthetic_usrcat


def test_synthetic_usrcat_is_deterministic_and_scaled():
    first = synthetic_usrcat(20, 7)
    second = synthetic_usrcat(20, 7)
    assert first.shape == (20, 60)
    assert first.dtype == np.float32
    assert np.array_equal(first, second)
    assert first[:, 0].std() > first[:, 2].std()


def test_candidate_recall():
    truth = np.asarray([[1, 2], [3, 4]])
    candidates = np.asarray([[2, 8, 9], [3, 4, 5]])
    assert candidate_recall(truth, candidates) == 0.75
