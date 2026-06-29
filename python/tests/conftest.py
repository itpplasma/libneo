"""Pytest configuration for libneo Python tests.

Extends libneo.__path__ to include the worktree source directory so that
modules added in the current branch (but not yet reflected in the editable
install's cached manifest) are importable during development.
"""

from pathlib import Path
import libneo


def pytest_configure(config):
    src = str(Path(__file__).resolve().parents[1] / "libneo")
    if src not in libneo.__path__:
        libneo.__path__.insert(0, src)
