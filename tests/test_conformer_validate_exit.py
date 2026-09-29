import sys

import pytest

from aidd_agent import conformer_block_validate as validator


@pytest.mark.parametrize('allow,status,structural,errors,expected',[
    (False,'complete','passed',[],2),
    (True,'complete','passed',[],0),
    (True,'failed','failed',['bad membership'],2),
    (True,'complete','passed',['corrupt'],2)])
def test_offline_review_continuation_keeps_gate(monkeypatch,allow,status,structural,errors,expected):
    report=dict(status=status,structural_gate=structural,errors=errors,release_gate='review_required')
    monkeypatch.setattr(validator,'validate',lambda *a:report)
    monkeypatch.setattr(sys,'argv',['validate','--build','source','--output','result']+
                        (['--allow-review-for-offline'] if allow else []))
    with pytest.raises(SystemExit) as exc:validator.main()
    assert exc.value.code==expected
    assert report['release_gate']=='review_required'
