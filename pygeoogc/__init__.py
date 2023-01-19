"""Top-level package for PyGeoOGC."""
from importlib.metadata import PackageNotFoundError, version

from pygeoogc.exceptions import (
    InputTypeError,
    InputValueError,
    MissingInputError,
    ServiceError,
    ServiceUnavailableError,
    ZeroMatchedError,
)
from pygeoogc.print_versions import show_versions
from pygeoogc.pygeoogc import WFS, WMS, ArcGISRESTful, ServiceURL
from pygeoogc.utils import RetrySession, match_crs, streaming_download, traverse_json, validate_crs

try:
    __version__ = version("pygeoogc")
except PackageNotFoundError:
    __version__ = "999"

__all__ = [
    "ArcGISRESTful",
    "WFS",
    "WMS",
    "ServiceURL",
    "RetrySession",
    "InputTypeError",
    "InputValueError",
    "MissingInputError",
    "ServiceError",
    "ServiceUnavailableError",
    "ZeroMatchedError",
    "traverse_json",
    "streaming_download",
    "match_crs",
    "validate_crs",
    "show_versions",
]
