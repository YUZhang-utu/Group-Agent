import pytest

from aidd_agent.block_exploration import prepare, replay


def panel():
    return [dict(block_id=str(i//20), molecule_id=str(i), score=str(i), reference_rank=str(100-i))
            for i in range(100)]


def test_equal_budget_determinism_and_full_recovery():
    a = replay(panel(), initial=2, batch=3, top_k=2, top=10)
    assert a == replay(list(reversed(panel())), initial=2, batch=3, top_k=2, top=10)
    assert a['adaptive']['evaluated_molecules'] == a['uniform']['evaluated_molecules'] == 25
    assert sum(a['block_allocations'].values()) == 25
    full = replay(panel(), initial=2, fraction=1, top=10)
    assert full['adaptive']['reference_top_recall'] == 1


def test_invalid_and_unscored_rows():
    with pytest.raises(ValueError, match='unique'):
        prepare(panel()+panel()[:1])
    with pytest.raises(ValueError, match='exceed'):
        replay(panel(), initial=10)
    rows = panel()
    rows[0].update(score='', reference_rank='', status='no_surviving_pose')
    assert replay(rows, initial=2)['reference_top_size'] == 99
    rows[0]['status'] = 'not_evaluated'
    with pytest.raises(ValueError, match='Blank'):
        prepare(rows)


def test_external_rank_is_reference_not_contact_order():
    rows = panel()
    for i, row in enumerate(rows):
        row['reference_rank'] = str(i+1)
    members, reference, objective = prepare(rows)
    assert reference[0] == '0'
    assert objective == 'supplied_full_run_rank'
