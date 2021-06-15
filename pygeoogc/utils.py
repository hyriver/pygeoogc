"""Some utilities for PyGeoOGC."""
import math
import socket
from concurrent import futures
from dataclasses import dataclass
from pathlib import Path
from sqlite3 import OperationalError
from typing import Any, Callable, Dict, Iterable, List, Mapping, Optional, Tuple, Union
from unittest.mock import _patch, patch

import defusedxml.ElementTree as etree
import pyproj
import simplejson as json
from requests import Response
from requests.adapters import HTTPAdapter
from requests.exceptions import RequestException
from requests_cache import CachedSession
from requests_cache.backends.sqlite import DbCache
from shapely import ops
from shapely.geometry import (
    LineString,
    MultiLineString,
    MultiPoint,
    MultiPolygon,
    Point,
    Polygon,
    box,
)
from urllib3 import Retry

from .exceptions import InvalidInputType, ThreadingException, ZeroMatched

DEF_CRS = "epsg:4326"
BOX_ORD = "(west, south, east, north)"
EXPIRE = 24 * 60 * 60


class RetrySession:
    """Configures the passed-in session to retry on failed requests.

    The fails can be due to connection errors, specific HTTP response
    codes and 30X redirections. The code is based on:
    https://github.com/bustawin/retry-requests

    Parameters
    ----------
    retries : int, optional
        The number of maximum retries before raising an exception, defaults to 5.
    backoff_factor : float, optional
        A factor used to compute the waiting time between retries, defaults to 0.5.
    status_to_retry : tuple, optional
        A tuple of status codes that trigger the reply behaviour, defaults to (500, 502, 504).
    prefixes : tuple, optional
        The prefixes to consider, defaults to ("http://", "https://")
    cache_name : str, optional
        Path to a folder for caching the session, default to None which uses
        system's temp directory.
    """

    def __init__(
        self,
        retries: int = 3,
        backoff_factor: float = 0.3,
        status_to_retry: Tuple[int, ...] = (500, 502, 504),
        prefixes: Tuple[str, ...] = ("https://",),
        cache_name: Optional[Union[str, Path]] = None,
    ) -> None:
        self.cache_name = (
            Path("cache", "http_cache.sqlite") if cache_name is None else Path(cache_name)
        )
        try:
            backend = DbCache(self.cache_name, fast_save=True, timeout=2)
            self.session = CachedSession(expire_after=EXPIRE, backend=backend)
            self.session.remove_expired_responses()
        except OperationalError:
            self._backup_db()

        self.retries = retries
        retry_args = {
            "total": retries,
            "read": retries,
            "connect": retries,
            "backoff_factor": backoff_factor,
            "status_forcelist": status_to_retry,
        }
        if hasattr(Retry.DEFAULT, "allowed_methods"):
            retry_args.update({"allowed_methods": False})
        else:
            retry_args.update({"method_whitelist": False})

        r = Retry(**retry_args)
        adapter = HTTPAdapter(max_retries=r)
        for prefix in prefixes:
            self.session.mount(prefix, adapter)
        self.session.hooks = {"response": [lambda r, *args, **kwargs: r.raise_for_status()]}

    def get(
        self,
        url: str,
        payload: Optional[Mapping[str, Any]] = None,
        headers: Optional[Mapping[str, Any]] = None,
    ) -> Response:
        """Retrieve data from a url by GET and return the Response."""
        try:
            return self.session.get(url, params=payload, headers=headers)
        except OperationalError:
            self._backup_db()
            return self.session.get(url, params=payload, headers=headers)
        except (ConnectionError, RequestException):
            raise ConnectionError(f"Connection failed after {self.retries} retries.")

    def post(
        self,
        url: str,
        payload: Optional[Mapping[str, Any]] = None,
        headers: Optional[Mapping[str, Any]] = None,
    ) -> Response:
        """Retrieve data from a url by POST and return the Response."""
        try:
            return self.session.post(url, data=payload, headers=headers)
        except OperationalError:
            self._backup_db()
            return self.session.post(url, data=payload, headers=headers)
        except (ConnectionError, RequestException):
            raise ConnectionError(f"Connection failed after {self.retries} retries.")

    def _backup_db(self) -> None:
        """Use a backup database if the current database is locked."""
        self.cache_name = Path(f"{self.cache_name.parent}/{self.cache_name.stem}_back_up.sqlite")
        backend = DbCache(self.cache_name, fast_save=True, timeout=2)
        self.session = CachedSession(expire_after=EXPIRE, backend=backend)
        self.session.remove_expired_responses()

    @staticmethod
    def onlyipv4() -> _patch:
        """Disable IPv6 and only use IPv4."""
        orig_getaddrinfo = socket.getaddrinfo

        def getaddrinfo_ipv4(
            host: str,
            port: str,
            family: socket.AddressFamily = socket.AF_INET,
            ptype: int = 0,
            proto: int = 0,
            flags: int = 0,
        ) -> List[
            Tuple[
                socket.AddressFamily,
                socket.SocketKind,
                int,
                str,
                Union[Tuple[str, int], Tuple[str, int, int, int]],
            ]
        ]:
            return orig_getaddrinfo(
                host=host,
                port=port,
                family=family,
                type=ptype,
                proto=proto,
                flags=flags,
            )

        return patch("socket.getaddrinfo", side_effect=getaddrinfo_ipv4)


def threading(
    func: Callable,
    iter_list: Iterable,
    param_list: Optional[List[Any]] = None,
    max_workers: int = 8,
) -> List[Any]:
    """Run a function in parallel with threading.

    Notes
    -----
    This function is suitable for IO intensive functions.

    Parameters
    ----------
    func : function
        The function to be ran in threads
    iter_list : list
        The iterable for the function
    param_list : list, optional
        List of other parameters, defaults to an empty list
    max_workers : int, optional
        Maximum number of threads, defaults to 8

    Returns
    -------
    list
        A list of function returns for each iterable. The list is not ordered.
    """
    data = []
    param_list = [] if param_list is None else param_list
    with futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_itr = {executor.submit(func, itr, *param_list): itr for itr in iter_list}
        for future in futures.as_completed(future_to_itr):
            itr = future_to_itr[future]
            try:
                data.append(future.result())
            except Exception as exc:  # noqa: B902
                raise ThreadingException(itr, exc)
    return data


def traverse_json(
    obj: Union[Dict[str, Any], List[Dict[str, Any]]], path: Union[str, List[str]]
) -> List[Any]:
    """Extract an element from a JSON file along a specified path.

    This function is based on `bcmullins <https://bcmullins.github.io/parsing-json-python/>`__.

    Parameters
    ----------
    obj : dict
        The input json dictionary
    path : list
        The path to the requested element

    Returns
    -------
    list
        The items founds in the JSON
    """

    def extract(
        obj: Optional[Union[List[Any], Dict[str, Any]]],
        path: Union[str, List[str]],
        ind: int,
        arr: List[Any],
    ) -> List[Any]:
        key = path[ind]
        if ind + 1 < len(path):
            if isinstance(obj, dict):
                if key in obj:
                    extract(obj.get(key), path, ind + 1, arr)
                else:
                    arr.append(None)
            elif isinstance(obj, list):
                if not obj:
                    arr.append(None)
                else:
                    for item in obj:
                        extract(item, path, ind, arr)
            else:
                arr.append(None)
        if ind + 1 == len(path):
            if isinstance(obj, list):
                if not obj:
                    arr.append(None)
                else:
                    for item in obj:
                        arr.append(item.get(key))
            elif isinstance(obj, dict):
                arr.append(obj.get(key))
            else:
                arr.append(None)
        return arr

    if isinstance(obj, dict):
        return extract(obj, path, 0, [])

    if isinstance(obj, list):
        outer_arr = []
        for item in obj:
            outer_arr.append(extract(item, path, 0, []))
        return outer_arr
    return []


@dataclass
class ESRIGeomQuery:
    """Generate input geometry query for ArcGIS RESTful services.

    Parameters
    ----------
    geometry : tuple or Polygon or Point or LineString
        The input geometry which can be a point (x, y), a list of points [(x, y), ...],
        bbox (xmin, ymin, xmax, ymax), or a Shapely's Polygon.
    wkid : int
        The Well-known ID (WKID) of the geometry's spatial reference e.g., for EPSG:4326,
        4326 should be passed. Check
        `ArcGIS <https://developers.arcgis.com/rest/services-reference/geographic-coordinate-systems.htm>`__
        for reference.
    """

    geometry: Union[
        Tuple[float, float],
        List[Tuple[float, float]],
        Tuple[float, float, float, float],
        Polygon,
    ]
    wkid: int

    def point(self) -> Dict[str, Union[str, bytes]]:
        """Query for a point."""
        if not (isinstance(self.geometry, tuple) and len(self.geometry) == 2):
            raise InvalidInputType("geometry (point)", "tuple", "(x, y)")

        geo_type = "esriGeometryPoint"
        geo_json = dict(zip(("x", "y"), self.geometry))
        return self.get_payload(geo_type, geo_json)

    def multipoint(self) -> Dict[str, Union[str, bytes]]:
        """Query for a multi-point."""
        if not (isinstance(self.geometry, list) and all(len(g) == 2 for g in self.geometry)):
            raise InvalidInputType("geometry (multi-point)", "list of tuples", "[(x, y), ...]")

        geo_type = "esriGeometryMultipoint"
        geo_json = {"points": [[x, y] for x, y in self.geometry]}
        return self.get_payload(geo_type, geo_json)

    def bbox(self) -> Dict[str, Union[str, bytes]]:
        """Query for a bbox."""
        if not (isinstance(self.geometry, (tuple, list)) and len(self.geometry) == 4):
            raise InvalidInputType("geometry (bbox)", "tuple or list", BOX_ORD)

        geo_type = "esriGeometryEnvelope"
        geo_json = dict(zip(("xmin", "ymin", "xmax", "ymax"), self.geometry))
        return self.get_payload(geo_type, geo_json)

    def polygon(self) -> Dict[str, Union[str, bytes]]:
        """Query for a polygon."""
        if not isinstance(self.geometry, Polygon):
            raise InvalidInputType("geometry", "Shapely's Polygon")

        geo_type = "esriGeometryPolygon"
        geo_json = {"rings": [[[x, y] for x, y in zip(*self.geometry.exterior.coords.xy)]]}
        return self.get_payload(geo_type, geo_json)

    def polyline(self) -> Dict[str, Union[str, bytes]]:
        """Query for a polyline."""
        if not isinstance(self.geometry, LineString):
            raise InvalidInputType("geometry", "Shapely's LineString")

        geo_type = "esriGeometryPolyline"
        geo_json = {"paths": [[[x, y] for x, y in zip(*self.geometry.coords.xy)]]}
        return self.get_payload(geo_type, geo_json)

    def get_payload(self, geo_type: str, geo_json: Dict[str, Any]) -> Dict[str, Union[str, bytes]]:
        """Generate a request payload based on ESRI template.

        Parameters
        ----------
        geo_type : str
            Type of the input geometry
        geo_json : dict
            Geometry in GeoJson format.

        Returns
        -------
        dict
            An ESRI geometry payload.
        """
        esri_json = json.dumps({**geo_json, "spatialRelference": {"wkid": str(self.wkid)}})
        return {
            "geometryType": geo_type,
            "geometry": esri_json,
            "inSR": str(self.wkid),
        }


class MatchCRS:
    """Reproject a geometry to another CRS.

    Parameters
    ----------
    in_crs : str
        Spatial reference of the input geometry
    out_crs : str
        Target spatial reference
    """

    def __init__(self, in_crs: str, out_crs: str):
        self.project = pyproj.Transformer.from_crs(in_crs, out_crs, always_xy=True).transform

    def geometry(
        self, geom: Union[Polygon, LineString, MultiLineString, MultiPolygon, Point, MultiPoint]
    ) -> Union[Polygon, LineString, MultiLineString, MultiPolygon, Point, MultiPoint]:
        """Reproject a geometry to the specified output CRS.

        Parameters
        ----------
        geom : LineString, MultiLineString, Polygon, MultiPolygon, Point, or MultiPoint
            Input geometry.

        Returns
        -------
        LineString, MultiLineString, Polygon, MultiPolygon, Point, or MultiPoint
            Input geometry in the specified CRS.

        Examples
        --------
        >>> from pygeoogc import MatchCRS
        >>> from shapely.geometry import Point
        >>> point = Point(-7766049.665, 5691929.739)
        >>> MatchCRS("epsg:3857", "epsg:4326").geometry(point).xy
        (array('d', [-69.7636111130079]), array('d', [45.44549114818127]))
        """
        if not isinstance(
            geom, (Polygon, LineString, MultiLineString, MultiPolygon, Point, MultiPoint)
        ):
            types = "LineString, MultiLineString, Polygon, MultiPolygon, Point, or MultiPoint"
            raise InvalidInputType("geom", types)

        return ops.transform(self.project, geom)

    def bounds(self, geom: Tuple[float, float, float, float]) -> Tuple[float, float, float, float]:
        """Reproject a bounding box to the specified output CRS.

        Parameters
        ----------
        geom : tuple
            Input bounding box (xmin, ymin, xmax, ymax).

        Returns
        -------
        tuple
            Input bounding box in the specified CRS.

        Examples
        --------
        >>> from pygeoogc import MatchCRS
        >>> bbox = (-7766049.665, 5691929.739, -7763049.665, 5696929.739)
        >>> MatchCRS("epsg:3857", "epsg:4326").bounds(bbox)
        (-69.7636111130079, 45.44549114818127, -69.73666165448431, 45.47699468552394)
        """
        if not (isinstance(geom, tuple) and len(geom) == 4):
            raise InvalidInputType("geom", "tuple", BOX_ORD)

        return ops.transform(self.project, box(*geom)).bounds

    def coords(self, geom: List[Tuple[float, float]]) -> List[Tuple[Any, ...]]:
        """Reproject a list of coordinates to the specified output CRS.

        Parameters
        ----------
        geom : list of tuple
            Input coords [(x1, y1), ...].

        Returns
        -------
        tuple
            Input list of coords in the specified CRS.

        Examples
        --------
        >>> from pygeoogc import MatchCRS
        >>> coords = [(-7766049.665, 5691929.739)]
        >>> MatchCRS("epsg:3857", "epsg:4326").coords(coords)
        [(-69.7636111130079, 45.44549114818127)]
        """
        if not (isinstance(geom, list) and all(len(c) == 2 for c in geom)):
            raise InvalidInputType("geom", "list of tuples", "[(x1, y1), ...]")

        return list(zip(*self.project(*zip(*geom))))


def check_bbox(bbox: Tuple[float, float, float, float]) -> None:
    """Check if an input inbox is a tuple of length 4."""
    if not (isinstance(bbox, tuple) and len(bbox) == 4):
        raise InvalidInputType("bbox", "tuple", BOX_ORD)


def bbox_resolution(
    bbox: Tuple[float, float, float, float], resolution: float, bbox_crs: str = DEF_CRS
) -> Tuple[int, int]:
    """Image size of a bounding box WGS84 for a given resolution in meters.

    Parameters
    ----------
    bbox : tuple
        A bounding box in WGS84 (west, south, east, north)
    resolution : float
        The resolution in meters
    bbox_crs : str, optional
        The spatial reference of the input bbox, default to EPSG:4326.

    Returns
    -------
    tuple
        The width and height of the image
    """
    check_bbox(bbox)

    bbox = MatchCRS(bbox_crs, DEF_CRS).bounds(bbox)
    west, south, east, north = bbox
    geod = pyproj.Geod(ellps="WGS84")

    linex = LineString([Point(west, south), Point(east, south)])
    delx = geod.geometry_length(linex)

    liney = LineString([Point(west, south), Point(west, north)])
    dely = geod.geometry_length(liney)

    return int(delx / resolution), int(dely / resolution)


def bbox_decompose(
    bbox: Tuple[float, float, float, float],
    resolution: float,
    box_crs: str = DEF_CRS,
    max_px: int = 8000000,
) -> List[Tuple[Tuple[float, float, float, float], str, int, int]]:
    """Split the bounding box vertically for WMS requests.

    Parameters
    ----------
    bbox : tuple
        A bounding box; (west, south, east, north)
    resolution : float
        The target resolution for a WMS request in meters.
    box_crs : str, optional
        The spatial reference of the input bbox, default to EPSG:4326.
    max_px : int, opitonal
        The maximum allowable number of pixels (width x height) for a WMS requests,
        defaults to 8 million based on some trial-and-error.

    Returns
    -------
    tuple
        The first element is a list of bboxes and the second one is width of the last bbox
    """
    check_bbox(bbox)
    _bbox = MatchCRS(box_crs, DEF_CRS).bounds(bbox)
    width, height = bbox_resolution(_bbox, resolution, DEF_CRS)

    n_px = width * height
    if n_px < max_px:
        return [(bbox, "0_0", width, height)]

    geod = pyproj.Geod(ellps="WGS84")
    west, south, east, north = _bbox

    def directional_split(
        az: float, origin: float, dest: float, lvl: float, xy: bool, px: int
    ) -> Tuple[List[Tuple[float, Any]], List[int]]:
        divs = [0]
        mul = 1.0
        coords = []

        def get_args(dst: float, dx: float) -> Tuple[float, float, float, float, int]:
            return (dst, lvl, az, dx, 0) if xy else (lvl, dst, az, dx, 1)

        while divs[-1] < 1:
            dim = int(math.sqrt(max_px) * mul)
            step = (dim - 1) * resolution

            _dest = origin
            while _dest < dest:
                args = get_args(_dest, step)
                coords.append((_dest, geod.fwd(*args[:-1])[args[-1]]))

                args = get_args(coords[-1][-1], resolution)
                _dest = geod.fwd(*args[:-1])[args[-1]]

            coords[-1] = (coords[-1][0], dest)

            n_dim = len(coords)
            divs = [dim for _ in range(n_dim)]
            divs[-1] = px - (n_dim - 1) * dim
            mul -= 0.1
        return coords, divs

    az_x = geod.inv(west, south, east, south)[0]
    lons, widths = directional_split(az_x, west, east, south, True, width)

    az_y = geod.inv(west, south, west, north)[0]
    lats, heights = directional_split(az_y, south, north, west, False, height)

    bboxs = []
    for i, ((bottom, top), h) in enumerate(zip(lats, heights)):
        for j, ((left, right), w) in enumerate(zip(lons, widths)):
            bx_crs = MatchCRS(DEF_CRS, box_crs).bounds((left, bottom, right, top))
            bboxs.append((bx_crs, f"{i}_{j}", w, h))
    return bboxs


def check_response(resp: Response) -> None:
    """Check if a ``requests.Resonse`` returned an error message."""
    if resp.headers["Content-Type"] == "application/xml":
        root = etree.fromstring(resp.text)
        raise ZeroMatched(root[0][0].text.strip())
