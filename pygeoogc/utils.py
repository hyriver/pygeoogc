"""Some utilities for PyGeoOGC."""
import math
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Tuple, TypeVar, Union

import async_retriever as ar
import defusedxml.ElementTree as etree
import pyproj
import shapely.geometry as sgeom
import ujson as json
from requests.adapters import HTTPAdapter
from requests.exceptions import RequestException
from requests_cache import CachedSession, Response
from requests_cache.backends.sqlite import SQLiteCache
from shapely import ops
from urllib3 import Retry

from .exceptions import InvalidInputType, ServiceError

DEF_CRS = "epsg:4326"
BOX_ORD = "(west, south, east, north)"
EXPIRE = -1
G = TypeVar(
    "G",
    sgeom.Point,
    sgeom.MultiPoint,
    sgeom.Polygon,
    sgeom.MultiPolygon,
    sgeom.LineString,
    sgeom.MultiLineString,
    Tuple[float, float, float, float],
    List[Tuple[float, float]],
)


def check_response(resp: str) -> str:
    """Extract error message from a response, if any."""
    try:
        root = etree.fromstring(resp)
    except etree.ParseError:
        return resp
    else:
        try:
            return str(root[-1][0].text).strip()
        except IndexError:
            return str(root[-1].text).strip()


class RetrySession:
    """Configures the passed-in session to retry on failed requests.

    The fails can be due to connection errors, specific HTTP response
    codes and 30X redirections. The code is was originally based on:
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
    expire_after : int, optional
        Expiration time for the cache in seconds, defaults to -1 (never expire).
    """

    def __init__(
        self,
        retries: int = 3,
        backoff_factor: float = 0.3,
        status_to_retry: Tuple[int, ...] = (500, 502, 504),
        prefixes: Tuple[str, ...] = ("https://",),
        cache_name: Optional[Union[str, Path]] = None,
        expire_after: int = EXPIRE,
    ) -> None:
        self.cache_name = (
            Path("cache", "http_cache.sqlite") if cache_name is None else Path(cache_name)
        )
        backend = SQLiteCache(self.cache_name, fast_save=True, timeout=1)
        self.session = CachedSession(
            expire_after=int(os.getenv("HYRIVER_CACHE_EXPIRE", expire_after)), backend=backend
        )

        adapter = HTTPAdapter(
            max_retries=Retry(
                total=retries,
                read=retries,
                connect=retries,
                backoff_factor=backoff_factor,
                status_forcelist=status_to_retry,
                allowed_methods=False,
            )
        )
        for prefix in prefixes:
            self.session.mount(prefix, adapter)

    def get(
        self,
        url: str,
        payload: Optional[Mapping[str, Any]] = None,
        headers: Optional[Mapping[str, Any]] = None,
    ) -> Response:
        """Retrieve data from a url by GET and return the Response."""
        resp = self.session.get(url, params=payload, headers=headers)
        try:
            resp.raise_for_status()
        except RequestException as ex:
            raise ServiceError(check_response(resp.text)) from ex
        else:
            return resp

    def post(
        self,
        url: str,
        payload: Optional[Mapping[str, Any]] = None,
        headers: Optional[Mapping[str, Any]] = None,
    ) -> Response:
        """Retrieve data from a url by POST and return the Response."""
        resp = self.session.post(url, data=payload, headers=headers)
        try:
            resp.raise_for_status()
        except RequestException as ex:
            raise ServiceError(check_response(resp.text)) from ex
        else:
            return resp


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

    Examples
    --------
    >>> from pygeoogc.utils import traverse_json
    >>> data = [{
    ...     "employees": [
    ...         {"name": "Alice", "role": "dev", "nbr": 1},
    ...         {"name": "Bob", "role": "dev", "nbr": 2}],
    ...     "firm": {"name": "Charlie's Waffle Emporium", "location": "CA"},
    ... },]
    >>> traverse_json(data, ["employees", "name"])
    [['Alice', 'Bob']]
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

    outer_arr = []
    for item in obj:
        outer_arr.append(extract(item, path, 0, []))
    return outer_arr


@dataclass
class ESRIGeomQuery:
    """Generate input geometry query for ArcGIS RESTful services.

    Parameters
    ----------
    geometry : tuple or sgeom.Polygon or sgeom.Point or sgeom.LineString
        The input geometry which can be a point (x, y), a list of points [(x, y), ...],
        bbox (xmin, ymin, xmax, ymax), or a Shapely's sgeom.Polygon.
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
        sgeom.Polygon,
    ]
    wkid: int

    def _get_payload(self, geo_type: str, geo_json: Dict[str, Any]) -> Mapping[str, str]:
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

    def point(self) -> Mapping[str, str]:
        """Query for a point."""
        if not (isinstance(self.geometry, tuple) and len(self.geometry) == 2):
            raise InvalidInputType("geometry", "tuple", "(x, y)")

        geo_type = "esriGeometryPoint"
        geo_json = dict(zip(("x", "y"), self.geometry))
        return self._get_payload(geo_type, geo_json)

    def multipoint(self) -> Mapping[str, str]:
        """Query for a multi-point."""
        if not (isinstance(self.geometry, list) and all(len(g) == 2 for g in self.geometry)):
            raise InvalidInputType("geometry", "list of tuples", "[(x, y), ...]")

        geo_type = "esriGeometryMultipoint"
        geo_json = {"points": [[x, y] for x, y in self.geometry]}
        return self._get_payload(geo_type, geo_json)

    def bbox(self) -> Mapping[str, str]:
        """Query for a bbox."""
        if not (isinstance(self.geometry, (tuple, list)) and len(self.geometry) == 4):
            raise InvalidInputType("geometry", "tuple", BOX_ORD)

        geo_type = "esriGeometryEnvelope"
        geo_json = dict(zip(("xmin", "ymin", "xmax", "ymax"), self.geometry))
        return self._get_payload(geo_type, geo_json)

    def polygon(self) -> Mapping[str, str]:
        """Query for a polygon."""
        if not isinstance(self.geometry, sgeom.Polygon):
            raise InvalidInputType("geometry", "Polygon")

        geo_type = "esriGeometryPolygon"
        geo_json = {"rings": [[[x, y] for x, y in zip(*self.geometry.exterior.coords.xy)]]}
        return self._get_payload(geo_type, geo_json)

    def polyline(self) -> Mapping[str, str]:
        """Query for a polyline."""
        if not isinstance(self.geometry, sgeom.LineString):
            raise InvalidInputType("geometry", "LineString")

        geo_type = "esriGeometryPolyline"
        geo_json = {"paths": [[[x, y] for x, y in zip(*self.geometry.coords.xy)]]}
        return self._get_payload(geo_type, geo_json)


def match_crs(geom: G, in_crs: str, out_crs: str) -> G:
    """Reproject a geometry to another CRS.

    Parameters
    ----------
    geom : list or tuple or geometry
        Input geometry which could be a list of coordinates such as ``[(x1, y1), ...]``,
        a bounding box like so ``(xmin, ymin, xmax, ymax)``, or any valid ``shapely``'s
        geometry such as ``Polygon``, ``MultiPolygon``, etc..
    in_crs : str
        Spatial reference of the input geometry
    out_crs : str
        Target spatial reference

    Returns
    -------
    same type as the input geometry
        Transformed geometry in the target CRS.

    Examples
    --------
    >>> from pygeoogc.utils import match_crs
    >>> from shapely.geometry import Point
    >>> point = Point(-7766049.665, 5691929.739)
    >>> match_crs(point, "epsg:3857", "epsg:4326").xy
    (array('d', [-69.7636111130079]), array('d', [45.44549114818127]))
    >>> bbox = (-7766049.665, 5691929.739, -7763049.665, 5696929.739)
    >>> match_crs(bbox, "epsg:3857", "epsg:4326")
    (-69.7636111130079, 45.44549114818127, -69.73666165448431, 45.47699468552394)
    >>> coords = [(-7766049.665, 5691929.739)]
    >>> match_crs(coords, "epsg:3857", "epsg:4326")
    [(-69.7636111130079, 45.44549114818127)]
    """
    project = pyproj.Transformer.from_crs(in_crs, out_crs, always_xy=True).transform

    if isinstance(
        geom,
        (
            sgeom.Polygon,
            sgeom.LineString,
            sgeom.MultiLineString,
            sgeom.MultiPolygon,
            sgeom.Point,
            sgeom.MultiPoint,
        ),
    ):
        return ops.transform(project, geom)  # type: ignore

    if isinstance(geom, tuple) and len(geom) == 4:
        return ops.transform(project, sgeom.box(*geom)).bounds  # type: ignore

    if isinstance(geom, list) and all(len(c) == 2 for c in geom):
        xx, yy = zip(*geom)
        return list(zip(*project(xx, yy)))

    gtypes = (
        "a list of coordinates such as [(x1, y1), ...],"
        + "a bounding box like so (xmin, ymin, xmax, ymax), or any valid shapely's geometry."
    )
    raise InvalidInputType("geom", gtypes)


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

    bbox = match_crs(bbox, bbox_crs, DEF_CRS)
    west, south, east, north = bbox
    geod = pyproj.Geod(ellps="WGS84")

    linex = sgeom.LineString([sgeom.Point(west, south), sgeom.Point(east, south)])
    delx = geod.geometry_length(linex)

    liney = sgeom.LineString([sgeom.Point(west, south), sgeom.Point(west, north)])
    dely = geod.geometry_length(liney)

    return int(delx / resolution), int(dely / resolution)


def bbox_decompose(
    bbox: Tuple[float, float, float, float],
    resolution: float,
    box_crs: str = DEF_CRS,
    max_px: int = 8000000,
) -> List[Tuple[Tuple[float, float, float, float], str, int, int]]:
    r"""Split the bounding box vertically for WMS requests.

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
    list of tuples
        Each tuple includes the following elements:

        * Tuple of length 4 that represents a bounding box (west, south, east, north) of a cell,
        * A label that represents cell ID starting from bottom-left to top-right, for example a
          2x2 decomposition has the following labels::

          |---------|---------|
          |         |         |
          |   0_1   |   1_1   |
          |         |         |
          |---------|---------|
          |         |         |
          |   0_0   |   1_0   |
          |         |         |
          |---------|---------|

        * Raster width of a cell,
        * Raster height of a cell.

    """
    check_bbox(bbox)
    _bbox = match_crs(bbox, box_crs, DEF_CRS)
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
            bx_crs = match_crs((left, bottom, right, top), DEF_CRS, box_crs)
            bboxs.append((bx_crs, f"{i}_{j}", w, h))
    return bboxs


def validate_crs(val: Union[str, int, pyproj.CRS]) -> str:
    """Validate a CRS.

    Parameters
    ----------
    val : str or int
        Input CRS.

    Returns
    -------
    str
        Validated CRS as a string.
    """
    try:
        crs: str = pyproj.CRS(val).to_string()
    except pyproj.exceptions.CRSError as ex:
        raise InvalidInputType("crs", "a valid CRS") from ex
    return crs


def valid_wms_crs(url: str) -> List[str]:
    """Get valid CRSs from a WMS service version 1.3.0."""
    ns = "http://www.opengis.net/wms"

    def get_path(tag_list: List[str]) -> str:
        return f"/{{{ns}}}".join([""] + tag_list)[1:]

    kwds = {"params": {"service": "wms", "request": "GetCapabilities"}}
    root = etree.fromstring(ar.retrieve_text([url], [kwds], ssl=False)[0])
    return [t.text.lower() for t in root.findall(get_path(["Capability", "Layer", "CRS"]))]
