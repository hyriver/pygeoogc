"""Top-level package for PyGeoOGC."""
import importlib.metadata

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
