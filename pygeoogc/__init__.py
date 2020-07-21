"""Top-level package for PyGeoOGC."""
from pkg_resources import DistributionNotFound, get_distribution

from .exceptions import InvalidInputType, InvalidInputValue, MissingInputs, ServerError, ZeroMatched
from .pygeoogc import WFS, ArcGISRESTful, MatchCRS, RetrySession, ServiceURL, wms_bybox

try:
    __version__ = get_distribution(__name__).version
except DistributionNotFound:
    # package is not installed
    pass
