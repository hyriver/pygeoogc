"""Configuration for pytest."""

import pytest


@pytest.fixture(autouse=True)
def add_standard_imports(doctest_namespace):
    """Add pygeoogc namespace for doctest."""
    import pygeoogc as ogc

    doctest_namespace["ogc"] = ogc
