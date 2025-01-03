"""Top-level package for PyGeoOGC."""

from __future__ import annotations

import os
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

from pygeoogc import exceptions
from pygeoogc.cache_keys import create_request_key
from pygeoogc.print_versions import show_versions
from pygeoogc.pygeoogc import WFS, WMS, ArcGISRESTful, ServiceURL
from pygeoogc.utils import RetrySession, match_crs, streaming_download, traverse_json, validate_crs

cert_path = os.getenv("HYRIVER_SSL_CERT")
if cert_path is not None:
    from pyproj.network import set_ca_bundle_path

    if not Path(cert_path).exists():
        raise FileNotFoundError(cert_path)
    set_ca_bundle_path(cert_path)

try:
    __version__ = version("pygeoogc")
except PackageNotFoundError:
    __version__ = "999"

__all__ = [
    "WFS",
    "WMS",
    "ArcGISRESTful",
    "RetrySession",
    "ServiceURL",
    "__version__",
    "create_request_key",
    "exceptions",
    "match_crs",
    "show_versions",
    "streaming_download",
    "traverse_json",
    "validate_crs",
]
