"""Shared exception identities for imported and `python -m` workflow execution."""


class Blocked(RuntimeError):
    """A known prerequisite prevents execution; no success is implied."""
