"""Top-level package for PyGeoOGC."""
from importlib.metadata import PackageNotFoundError, version

from .exceptions import (
    InputTypeError,
    InputValueError,
    MissingInputError,
    ServiceError,
    ServiceUnavailableError,
    ZeroMatchedError,
)
from .print_versions import show_versions
from .pygeoogc import WFS, WMS, ArcGISRESTful, ServiceURL
from .utils import RetrySession

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
    "show_versions",
]
