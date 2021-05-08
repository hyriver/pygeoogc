"""Top-level package for PyGeoOGC."""
import sys

from pkg_resources import DistributionNotFound, get_distribution

from .exceptions import (
    InvalidInputType,
    InvalidInputValue,
    MissingInputs,
    ServerError,
    ThreadingException,
    ZeroMatched,
)
from .print_versions import show_versions
from .pygeoogc import WFS, WMS, ArcGISRESTful, ServiceURL
from .utils import MatchCRS, RetrySession

try:
    __version__ = get_distribution(__name__).version
except DistributionNotFound:
    # package is not installed
    pass
