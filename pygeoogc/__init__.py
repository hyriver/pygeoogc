"""Top-level package for PyGeoOGC."""
from .exceptions import (
    InvalidInputType,
    InvalidInputValue,
    MissingInputs,
    ServiceError,
    ServiceUnavailable,
    ZeroMatched,
)
from .print_versions import show_versions
from .pygeoogc import WFS, WMS, ArcGISRESTful, ServiceURL
from .utils import RetrySession

try:
    import importlib.metadata
except ImportError:
    import importlib_metadata

    __version__ = importlib_metadata.version("pygeoogc")
else:
    __version__ = importlib.metadata.version("pygeoogc")

__all__ = [
    "ArcGISRESTful",
    "WFS",
    "WMS",
    "ServiceURL",
    "RetrySession",
    "InvalidInputType",
    "InvalidInputValue",
    "MissingInputs",
    "ServiceError",
    "ServiceUnavailable",
    "ZeroMatched",
    "show_versions",
]
