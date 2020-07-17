"""Top-level package for PyGeoOGC."""


from .exceptions import InvalidInputType, InvalidInputValue, MissingInputs, ServerError, ZeroMatched
from .pygeoogc import WFS, ArcGISRESTful, MatchCRS, RetrySession, ServiceURL, wms_bybox

__author__ = """Taher Chegini"""
__email__ = "cheginit@gmail.com"
__version__ = "0.1.2"
