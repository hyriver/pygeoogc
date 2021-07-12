"""Top-level package for PyGeoOGC."""
from pkg_resources import DistributionNotFound, get_distribution

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
    __version__ = get_distribution(__name__).version
except DistributionNotFound:
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
