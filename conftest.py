"""Configuration for pytest."""

from __future__ import annotations

import pytest


@pytest.fixture(autouse=True)
def _add_standard_imports(doctest_namespace):
    """Add pygeoogc namespace for doctest."""
    import pygeoogc as ogc

    doctest_namespace["ogc"] = ogc
