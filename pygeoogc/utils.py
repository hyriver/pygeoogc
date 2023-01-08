"""Some utilities for PyGeoOGC."""
from __future__ import annotations

import itertools
import math
import os
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Generator, List, Mapping, Tuple, TypeVar, Union

import async_retriever as ar
import cytoolz as tlz
import defusedxml.ElementTree as ETree
import joblib
import pyproj
import requests
import ujson
import urllib3
from pyproj.exceptions import CRSError as ProjCRSError
from requests.adapters import HTTPAdapter
from requests.exceptions import RequestException
from requests_cache import CachedSession, Response
from requests_cache.backends.sqlite import SQLiteCache
from shapely import ops
from shapely.geometry import LineString, MultiLineString, MultiPoint, MultiPolygon, Point, Polygon
from shapely.geometry import box as shapely_box
from urllib3.exceptions import InsecureRequestWarning

from pygeoogc import cache_keys
from pygeoogc.exceptions import InputTypeError, InputValueError, ServiceError

CRSTYPE = Union[int, str, pyproj.CRS]
BOX_ORD = "(west, south, east, north)"
G = TypeVar(
    "G",
    Point,
    MultiPoint,
    Polygon,
    MultiPolygon,
    LineString,
    MultiLineString,
    Tuple[float, float, float, float],
    List[Tuple[float, float]],
)
MAX_CONN = 10
CHUNK_SIZE = int(100 * 1024 * 1024)  # 100 MB

__all__ = ["RetrySession", "traverse_json", "streaming_download", "match_crs", "validate_crs"]

warnings.filterwarnings("ignore", message=".*too short worker timeout.*")


def check_response(resp: str) -> str:
    """Extract error message from a response, if any."""
    try:
        root = ETree.fromstring(resp)
    except ETree.ParseError:
        return resp
    else:
        try:
            return str(root[-1][0].text).strip()
        except IndexError:
            return str(root[-1].text).strip()


class RetrySession:
    """Configures the passed-in session to retry on failed requests.

    Notes
    -----
    The fails can be due to connection errors, specific HTTP response
    codes and 30X redirections. The code was originally based on:
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
    disable : bool, optional
        If ``True`` temporarily disable caching request/responses, defaults to ``False``.
    ssl : bool, optional
        If ``True`` verify SSL certificates, defaults to ``True``.
    """

    def __init__(
        self,
        retries: int = 3,
        backoff_factor: float = 0.3,
        status_to_retry: tuple[int, ...] = (500, 502, 504),
        prefixes: tuple[str, ...] = ("https://",),
        cache_name: str | Path | None = None,
        expire_after: int = -1,
        disable: bool = False,
        ssl: bool = True,
    ) -> None:
        self.disable = os.getenv("HYRIVER_CACHE_DISABLE", f"{disable}").lower() == "true"
        if self.disable:
            self.session = requests.Session()
        else:
            self.cache_name = os.getenv(
                "HYRIVER_CACHE_NAME_HTTP", cache_name or Path("cache", "http_cache.sqlite")
            )
            backend = SQLiteCache(self.cache_name, fast_save=True, timeout=1)
            self.session = CachedSession(
                expire_after=int(os.getenv("HYRIVER_CACHE_EXPIRE", expire_after)), backend=backend
            )

        if not ssl:
            urllib3.disable_warnings(InsecureRequestWarning)
            self.session.verify = False

        adapter = HTTPAdapter(
            max_retries=urllib3.Retry(
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
        payload: Mapping[str, Any] | None = None,
        params: Mapping[str, Any] | None = None,
        headers: Mapping[str, Any] | None = None,
        stream: bool | None = None,
    ) -> Response:
        """Retrieve data from a url by GET and return the Response."""
        if stream and not self.disable:
            msg = ". ".join(
                (
                    "Streaming is not supported with caching enabled",
                    "Reinstantiate the class like so: RetrySession(disable=True).",
                )
            )
            raise ValueError(msg)
        params = params or payload
        resp = self.session.get(url, params=params, headers=headers, stream=stream)
        try:
            resp.raise_for_status()
        except RequestException as ex:
            raise ServiceError(check_response(resp.text)) from ex
        else:
            return resp

    def post(
        self,
        url: str,
        payload: Mapping[str, Any] | None = None,
        data: Mapping[str, Any] | None = None,
        json: Mapping[str, Any] | None = None,
        headers: Mapping[str, Any] | None = None,
        stream: bool | None = None,
    ) -> Response:
        """Retrieve data from a url by POST and return the Response."""
        if stream and not self.disable:
            msg = ". ".join(
                (
                    "Streaming is not supported with caching enabled",
                    "Reinstantiate the class like so: RetrySession(disable=True).",
                )
            )
            raise ValueError(msg)
        data = data or payload
        resp = self.session.post(url, data=data, json=json, headers=headers, stream=stream)
        try:
            resp.raise_for_status()
        except RequestException as ex:
            raise ServiceError(check_response(resp.text)) from ex
        else:
            return resp

    def head(
        self,
        url: str,
        params: Mapping[str, Any] | None = None,
        data: Mapping[str, Any] | None = None,
        json: Mapping[str, Any] | None = None,
        headers: Mapping[str, Any] | None = None,
    ) -> Response:
        """Retrieve data from a url by POST and return the Response."""
        resp = self.session.head(url, data=data, params=params, json=json, headers=headers)
        try:
            resp.raise_for_status()
        except RequestException as ex:
            raise ServiceError(check_response(resp.text)) from ex
        else:
            return resp

    def close(self) -> None:
        """Close the session."""
        self.session.close()


def _prepare_requests_args(
    urls: list[str] | str,
    kwds: list[dict[str, dict[Any, Any]]] | dict[str, dict[Any, Any]] | None,
    method: str,
    fnames: str | Path | list[str | Path] | None,
    file_prefix: str,
    file_extention: str,
) -> tuple[
    tuple[str, ...], tuple[dict[str, None | dict[Any, Any]], ...], Generator[Path, None, None]
]:
    """Get url and kwds for streaming download."""
    url_list = tuple(urls) if isinstance(urls, (list, tuple)) else (urls,)

    if kwds is None:
        if method == "GET":
            kwd_list = ({"params": None},) * len(url_list)
        else:
            kwd_list = ({"data": None},) * len(url_list)
    else:
        kwd_list = tuple(kwds) if isinstance(kwds, (list, tuple)) else (kwds,)  # type: ignore
    key_list = {k for keys in kwd_list for k in keys}
    valid_keys = ("params", "data", "json", "headers")
    if any(k not in valid_keys for k in key_list):
        raise InputValueError("kwds", valid_keys)

    if len(url_list) != len(kwd_list):
        raise InputTypeError("urls/kwds", "list of same length")

    fex = file_extention.replace(".", "")
    if fnames is None:
        cache_dir = os.getenv("HYRIVER_CACHE_NAME", str(Path("cache", "tmp")))
        cache_dir = Path(cache_dir).parent
        files = (
            Path(cache_dir, f"{file_prefix}{cache_keys.create_key(method, u, **p)}.{fex}")
            for u, p in zip(url_list, kwd_list)
        )
    else:
        f_list = tuple(fnames) if isinstance(fnames, (list, tuple)) else (fnames,)
        if len(url_list) != len(f_list):
            raise InputTypeError("urls/fnames", "lists of same length")
        files = (Path(f) for f in f_list)

    return url_list, kwd_list, files  # type: ignore[return-value]


def _download(
    session_func: Callable[[str], Response],
    url: str,
    kwd: Mapping[str, None | dict[Any, Any]],
    fname: Path,
    chunk_size: int,
) -> Path:
    """Download a single file."""
    resp = session_func(url, **kwd)
    fsize = int(resp.headers.get("Content-Length", -1))
    if not fname.exists() or fname.stat().st_size != fsize:
        fname.parent.mkdir(exist_ok=True, parents=True)
        with fname.open("wb") as f:
            f.writelines(resp.iter_content(chunk_size))
    return fname


def streaming_download(
    urls: list[str] | str,
    kwds: list[dict[str, dict[Any, Any]]] | dict[str, dict[Any, Any]] | None = None,
    fnames: str | Path | list[str | Path] | None = None,
    file_prefix: str = "",
    file_extention: str = "",
    method: str = "GET",
    ssl: bool = True,
    chunk_size: int = CHUNK_SIZE,
    n_jobs: int = MAX_CONN,
) -> Path | list[Path]:
    """Download and store files in parallel from a list of URLs/Keywords.

    Notes
    -----
    This function uses ``joblib`` with ``loky`` backend.

    Parameters
    ----------
    urls : tuple or list
        A list of URLs to download.
    kwds : tuple or list, optional
        A list of keywords associated with each URL, e.g.,
        ({"params": ..., "headers": ...}, ...). Defaults to ``None``.
    fnames : tuple or list, optional
        A list of filenames associated with each URL, e.g.,
        ("file1.zip", ...). Defaults to ``None``. If not provided,
        random unique filenames will be generated based on
        URL and keyword pairs.
    file_prefix : str, optional
        Prefix to add to filenames when storing the files, defaults
        to ``None``, i.e., no prefix. This argument will be only be
        used if ``fnames`` is not passed.
    file_extention : str, optional
        Extension to use for storing the files, defaults to ``None``,
        i.e., no extension if ``fnames`` is not provided otherwise. This
        argument will be only be used if ``fnames`` is not passed.
    method : str, optional
        HTTP method to use, i.e, ``GET`` or ``POST``, by default "GET".
    ssl : bool, optional
        Whether to use SSL verification, defaults to ``True``.
    chunk_size : int, optional
        Chunk size to use when downloading, defaults to 100 * 1024 * 1024
        i.e., 100 MB.
    n_jobs: int, optional
        The maximum number of concurrent downloads, defaults to 10.

    Returns
    -------
    list
        A list of ``pathlib.Path`` objects associated with URLs in the
        same order.
    """
    method = method.upper()
    valid_methods = ("GET", "POST")
    if method not in valid_methods:
        raise InputValueError("method", valid_methods)

    url_list, kwd_list, files = _prepare_requests_args(
        urls, kwds, method, fnames, file_prefix, file_extention
    )

    session = RetrySession(disable=True, ssl=ssl)
    if method == "GET":
        func = tlz.partial(session.get, stream=True)
    else:
        func = tlz.partial(session.post, stream=True)

    n_jobs = min(n_jobs, len(url_list))
    fpaths: list[Path] = joblib.Parallel(n_jobs=n_jobs)(
        joblib.delayed(_download)(func, u, k, f, chunk_size)
        for u, k, f in zip(url_list, kwd_list, files)
    )
    if n_jobs == 1:
        return fpaths[0]
    return fpaths


def traverse_json(
    items: dict[str, Any] | list[dict[str, Any]], ipath: str | list[str]
) -> list[Any]:
    """Extract an element from a JSON file along a specified ipath.

    This function is based on `bcmullins <https://bcmullins.github.io/parsing-json-python/>`__.

    Parameters
    ----------
    items : dict
        The input json dictionary
    ipath : list
        The ipath to the requested element

    Returns
    -------
    list
        The sub_items founds in the JSON

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
        sub_items: list[Any] | dict[str, Any] | None,
        path: str | list[str],
        ind: int,
        arr: list[Any],
    ) -> list[Any]:
        key = path[ind]
        if ind + 1 < len(path):
            if isinstance(sub_items, dict):
                if key in sub_items:
                    extract(sub_items.get(key), path, ind + 1, arr)
                else:
                    arr.append(None)
            elif isinstance(sub_items, list):
                if not sub_items:
                    arr.append(None)
                else:
                    for i in sub_items:
                        extract(i, path, ind, arr)
            else:
                arr.append(None)
        if ind + 1 == len(path):
            if isinstance(sub_items, list):
                if not sub_items:
                    arr.append(None)
                else:
                    for i in sub_items:
                        arr.append(i.get(key))
            elif isinstance(sub_items, dict):
                arr.append(sub_items.get(key))
            else:
                arr.append(None)
        return arr

    if isinstance(items, dict):
        return extract(items, ipath, 0, [])

    outer_arr = [extract(item, ipath, 0, []) for item in items]
    return outer_arr


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

    geometry: (
        tuple[float, float]
        | list[tuple[float, float]]
        | tuple[float, float, float, float]
        | Polygon
        | LineString
    )
    wkid: int

    def _get_payload(self, geo_type: str, geo_json: dict[str, Any]) -> Mapping[str, str]:
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
        esri_json = ujson.dumps({**geo_json, "spatialRelference": {"wkid": str(self.wkid)}})
        return {
            "geometryType": geo_type,
            "geometry": esri_json,
            "inSR": str(self.wkid),
        }

    def point(self) -> Mapping[str, str]:
        """Query for a point."""
        if isinstance(self.geometry, tuple) and len(self.geometry) == 2:
            geo_type = "esriGeometryPoint"
            geo_json = dict(zip(("x", "y"), self.geometry))  # type: ignore
            return self._get_payload(geo_type, geo_json)

        raise InputTypeError("geometry", "tuple", "(x, y)")

    def multipoint(self) -> Mapping[str, str]:
        """Query for a multi-point."""
        if isinstance(self.geometry, list) and all(len(g) == 2 for g in self.geometry):
            geo_type = "esriGeometryMultipoint"
            geo_json = {"points": [[x, y] for x, y in self.geometry]}  # type: ignore
            return self._get_payload(geo_type, geo_json)

        raise InputTypeError("geometry", "list of tuples", "[(x, y), ...]")

    def bbox(self) -> Mapping[str, str]:
        """Query for a bbox."""
        if isinstance(self.geometry, (tuple, list)) and len(self.geometry) == 4:
            geo_type = "esriGeometryEnvelope"
            geo_json = dict(zip(("xmin", "ymin", "xmax", "ymax"), self.geometry))  # type: ignore
            return self._get_payload(geo_type, geo_json)

        raise InputTypeError("geometry", "tuple", BOX_ORD)

    def polygon(self) -> Mapping[str, str]:
        """Query for a polygon."""
        if isinstance(self.geometry, Polygon):
            geo_type = "esriGeometryPolygon"
            geo_json = {
                "rings": [
                    [[x, y] for x, y in zip(*self.geometry.exterior.coords.xy)]  # type: ignore
                ]
            }
            return self._get_payload(geo_type, geo_json)

        raise InputTypeError("geometry", "Polygon")

    def polyline(self) -> Mapping[str, str]:
        """Query for a polyline."""
        if isinstance(self.geometry, LineString):
            geo_type = "esriGeometryPolyline"
            geo_json = {
                "paths": [[[x, y] for x, y in zip(*self.geometry.coords.xy)]]  # type: ignore
            }
            return self._get_payload(geo_type, geo_json)

        raise InputTypeError("geometry", "LineString")


def match_crs(geom: G, in_crs: CRSTYPE, out_crs: CRSTYPE) -> G:
    """Reproject a geometry to another CRS.

    Parameters
    ----------
    geom : list or tuple or geometry
        Input geometry which could be a list of coordinates such as ``[(x1, y1), ...]``,
        a bounding box like so ``(xmin, ymin, xmax, ymax)``, or any valid ``shapely``'s
        geometry such as ``Polygon``, ``MultiPolygon``, etc..
    in_crs : str, int, or pyproj.CRS
        Spatial reference of the input geometry
    out_crs : str, int, or pyproj.CRS
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
    if pyproj.CRS(in_crs) == pyproj.CRS(out_crs):
        return geom

    project = pyproj.Transformer.from_crs(in_crs, out_crs, always_xy=True).transform

    if isinstance(
        geom,
        (
            Polygon,
            LineString,
            MultiLineString,
            MultiPolygon,
            Point,
            MultiPoint,
        ),
    ):
        return ops.transform(project, geom)

    if isinstance(geom, tuple) and len(geom) == 4:
        return ops.transform(project, shapely_box(*geom)).bounds

    if isinstance(geom, list) and all(len(c) == 2 for c in geom):
        xx, yy = zip(*geom)
        return list(zip(*project(xx, yy)))  # type: ignore

    gtypes = (
        "a list of coordinates such as [(x1, y1), ...],"
        + "a bounding box like so (xmin, ymin, xmax, ymax), or any valid shapely's geometry."
    )
    raise InputTypeError("geom", gtypes)


def check_bbox(bbox: tuple[float, float, float, float]) -> None:
    """Check if an input inbox is a tuple of length 4."""
    if not (isinstance(bbox, tuple) and len(bbox) == 4):
        raise InputTypeError("bbox", "tuple", BOX_ORD)


def bbox_decompose(
    bbox: tuple[float, float, float, float],
    resolution: float,
    box_crs: CRSTYPE = 4326,
    max_px: int = 8000000,
) -> list[tuple[tuple[float, float, float, float], str, int, int]]:
    r"""Split the bounding box vertically for WMS requests.

    Parameters
    ----------
    bbox : tuple
        A bounding box; (west, south, east, north)
    resolution : float
        The target resolution for a WMS request in meters.
    box_crs : str, int, or pyproj.CRS, optional
        The spatial reference of the input bbox, default to ``epsg:4326``.
    max_px : int, opitonal
        The maximum allowable number of pixels (width x height) for a WMS requests,
        defaults to 8 million based on some trial-and-error.

    Returns
    -------
    list of tuples
        Each tuple includes the following elements:

        * Tuple of px_tot 4 that represents a bounding box (west, south, east, north) of a cell,
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

    geod = pyproj.Geod(ellps="GRS80")

    west, south, east, north = bbox

    xmin, ymin, xmax, ymax = match_crs(bbox, box_crs, 4326)

    x_dist = geod.geometry_length(LineString([(xmin, ymin), (xmax, ymin)]))
    y_dist = geod.geometry_length(LineString([(xmin, ymin), (xmin, ymax)]))
    width = int(math.ceil(x_dist / resolution))
    height = int(math.ceil(y_dist / resolution))

    if width * height <= max_px:
        bboxs = [(bbox, "0_0", width, height)]

    n_px = int(math.sqrt(max_px))

    def _split_directional(low: float, high: float, px_tot: int) -> tuple[list[int], list[float]]:
        npt = [n_px for _ in range(int(px_tot / n_px))] + [px_tot % n_px]
        xd = abs(high - low)
        dx = [xd * n / sum(npt) for n in npt]
        xs = [low + d for d in itertools.accumulate(dx)]
        xs.insert(0, low)
        return npt, xs

    nw, xs = _split_directional(west, east, width)
    nh, ys = _split_directional(south, north, height)

    bboxs = []
    for j, h in enumerate(nh):
        for i, w in enumerate(nw):
            bx_crs = (xs[i], ys[j], xs[i + 1], ys[j + 1])
            bboxs.append((bx_crs, f"{i}_{j}", w, h))
    return bboxs


def validate_crs(crs: CRSTYPE) -> str:
    """Validate a CRS.

    Parameters
    ----------
    crs : str, int, or pyproj.CRS
        Input CRS.

    Returns
    -------
    str
        Validated CRS as a string.
    """
    try:
        return pyproj.CRS(crs).to_string()  # type: ignore
    except ProjCRSError as ex:
        raise InputTypeError("crs", "a valid CRS") from ex


def valid_wms_crs(url: str) -> list[str]:
    """Get valid CRSs from a WMS service version 1.3.0."""
    ns = "http://www.opengis.net/wms"

    def get_path(tag_list: list[str]) -> str:
        return f"/{{{ns}}}".join([""] + tag_list)[1:]

    kwds = {"params": {"service": "wms", "request": "GetCapabilities"}}
    root = ETree.fromstring(ar.retrieve_text([url], [kwds], ssl=False)[0])
    return [t.text.lower() for t in root.findall(get_path(["Capability", "Layer", "CRS"]))]
