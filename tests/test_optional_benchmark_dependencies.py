from aidd_agent.faiss_benchmark import candidate_recall, synthetic_usrcat


def test_faiss_helpers_import_without_psutil():
    # Importing and using pure benchmark helpers must not require RSS monitoring.
    vectors = synthetic_usrcat(4, 7)
    assert vectors.shape == (4, 60)
    assert candidate_recall(vectors[:, :1].astype(int), vectors[:, :1].astype(int)) == 1.0
