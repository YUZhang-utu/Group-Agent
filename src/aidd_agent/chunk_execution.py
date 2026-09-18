"""Bounded scheduling; results carry their input index for deterministic assembly."""
from concurrent.futures import FIRST_COMPLETED, wait


def bounded_results(executor, function, tasks, limit):
    if limit < 1:
        raise ValueError("in-flight limit must be positive")
    iterator = iter(enumerate(tasks))
    pending = {}
    try:
        while True:
            while len(pending) < limit:
                try:
                    index, task = next(iterator)
                except StopIteration:
                    break
                pending[executor.submit(function, task)] = index
            if not pending:
                break
            done, _ = wait(pending, return_when=FIRST_COMPLETED)
            for future in done:
                index = pending.pop(future)
                yield index, future.result()
    finally:
        for future in pending:
            future.cancel()
