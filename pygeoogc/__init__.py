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
    import importlib.metadata as metadata
except ImportError:
    import importlib_metadata as metadata  # type: ignore[no-redef]

try:
    __version__ = metadata.version("pygeoogc")
except Exception:
    __version__ = "999"

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
