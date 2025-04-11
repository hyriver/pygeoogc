"""Some utilities for PyGeoOGC."""

from __future__ import annotations

import contextlib
import itertools
import json
import math
import os
import warnings
from dataclasses import dataclass
from pathlib import Path
from threading import Lock
from typing import (
    TYPE_CHECKING,
    Any,
    Callable,
    Literal,
    TypeVar,
    cast,
    overload,
)

import cytoolz.curried as tlz
import defusedxml.ElementTree as ETree
import joblib
import pyproj
import requests
import shapely
import urllib3
from pyproj.exceptions import CRSError as ProjCRSError
from requests.adapters import HTTPAdapter
from requests.exceptions import RequestException
from requests_cache import CachedSession, Response
from requests_cache.backends.sqlite import SQLiteCache
from shapely import LineString, MultiLineString, MultiPoint, MultiPolygon, Point, Polygon, ops
from urllib3.exceptions import InsecureRequestWarning

import async_retriever as ar
from pygeoogc import cache_keys
from pygeoogc.exceptions import (
    InputTypeError,
    InputValueError,
    ServiceError,
    ServiceUnavailableError,
)

if TYPE_CHECKING:
    from collections.abc import Generator, Mapping, Sequence

    from pyproj import CRS
    from typing_extensions import Self

    CRSType = int | str | CRS
    GeomType = TypeVar(
        "GeomType",
        Point,
        MultiPoint,
        Polygon,
        MultiPolygon,
        LineString,
        MultiLineString,
        "tuple[float, float, float, float]",
        "list[tuple[float, float]]",
    )
BOX_ORD = "(west, south, east, north)"
MAX_CONN = 10
CHUNK_SIZE = 100 * 1024 * 1024  # 100 MB
EXPIRE_AFTER = 60 * 60 * 24 * 7  # 1 week
__all__ = ["RetrySession", "match_crs", "streaming_download", "traverse_json", "validate_crs"]

warnings.filterwarnings("ignore", message=".*too short worker timeout.*")


def check_response(resp: str) -> str:
    """Extract error message from a response, if any."""
    try:
        root = ETree.fromstring(resp)
    except ETree.ParseError:
        return resp
    else:
        try:
            return str(root[-1][0].text).strip()  # pyright: ignore[reportIndexIssue]
        except IndexError:
            try:
                return str(root[-1].text).strip()  # pyright: ignore[reportIndexIssue]
            except IndexError:
                return str(root.text).strip()  # pyright: ignore[reportIndexIssue,reportAttributeAccessIssue]


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

    _lock = Lock()

    def __init__(
        self,
        retries: int = 3,
        backoff_factor: float = 0.3,
        status_to_retry: tuple[int, ...] = (500, 502, 504),
        prefixes: tuple[str, ...] = ("https://",),
        cache_name: str | Path | None = None,
        expire_after: int = EXPIRE_AFTER,
        disable: bool = False,
        ssl: bool = True,
    ) -> None:
        if cache_name is None:
            self.cache_name = Path(
                os.getenv("HYRIVER_CACHE_NAME_HTTP", Path("cache", "http_cache.sqlite"))
            )
        else:
            self.cache_name = Path(cache_name)
        self.cache_name.parent.mkdir(exist_ok=True, parents=True)

        self.expire_after = expire_after
        if self.expire_after == EXPIRE_AFTER:
            self.expire_after = int(os.getenv("HYRIVER_CACHE_EXPIRE", EXPIRE_AFTER))

        self.ssl_cert = os.getenv("HYRIVER_SSL_CERT")

        self.disable = disable

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
                allowed_methods=None,
            )
        )
        for prefix in prefixes:
            self.session.mount(prefix, adapter)

    def __del__(self) -> None:
        """Ensure resources are cleaned up."""
        # Suppress errors during garbage collection
        with contextlib.suppress(Exception):
            self.close()

    @property
    def disable(self) -> bool:
        """Disable caching request/responses."""
        return self._disable

    @disable.setter
    def disable(self, value: bool) -> None:
        with self._lock:
            self._disable = value
            if not self._disable:
                self._disable = os.getenv("HYRIVER_CACHE_DISABLE", "false").lower() == "true"

            if self._disable:
                self.session = requests.Session()
            else:
                backend = SQLiteCache(self.cache_name, fast_save=True, timeout=1)
                self.session = CachedSession(expire_after=self.expire_after, backend=backend)
            self.session.cert = self.ssl_cert

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

    def __enter__(self) -> Self:
        return self

    def __exit__(self, *args: object) -> None:
        self.close()


def _prepare_kwd_list_and_files(
    kwds: list[dict[str, dict[Any, Any]]] | dict[str, dict[Any, Any]] | None,
    method: Literal["GET", "POST", "get", "post"],
    url_list: tuple[str, ...],
    fnames: str | Path | Sequence[str | Path] | None,
    root_dir: str | Path | None,
    file_prefix: str,
    file_extention: str,
) -> tuple[tuple[dict[str, None | dict[Any, Any]], ...], Generator[Path, None, None]]:
    if kwds is None:
        if method == "GET":
            kwd_list = ({"params": None},) * len(url_list)
        else:
            kwd_list = ({"data": None},) * len(url_list)
    else:
        kwd_list = (kwds,) if isinstance(kwds, dict) else tuple(kwds)

    f_ext = file_extention.replace(".", "")
    f_ext = f".{f_ext}" if f_ext else ""

    if fnames is None:
        if root_dir is None:
            root_dir = os.getenv("HYRIVER_CACHE_NAME", str(Path("cache", "tmp")))
            root_dir = Path(root_dir).parent
        else:
            root_dir = Path(root_dir)
        root_dir.mkdir(exist_ok=True, parents=True)
        files = (
            Path(root_dir, f"{file_prefix}{cache_keys.create_request_key(method, u, **p)}{f_ext}")
            for u, p in zip(url_list, kwd_list)
        )
    else:
        f_list = (fnames,) if isinstance(fnames, (str, Path)) else tuple(fnames)
        if len(url_list) != len(f_list):
            raise InputTypeError("urls/fnames", "lists of same length")
        files = (Path(f) for f in f_list)

    return kwd_list, files


def _prepare_requests_args(
    urls: list[str] | str,
    kwds: list[dict[str, dict[Any, Any]]] | dict[str, dict[Any, Any]] | None,
    method: Literal["GET", "POST", "get", "post"],
    fnames: str | Path | Sequence[str | Path] | None,
    root_dir: str | Path | None,
    file_prefix: str,
    file_extention: str,
) -> tuple[
    tuple[str, ...], tuple[dict[str, None | dict[Any, Any]], ...], Generator[Path, None, None]
]:
    """Get url and kwds for streaming download."""
    url_list = (urls,) if isinstance(urls, str) else tuple(urls)
    kwd_list, files = _prepare_kwd_list_and_files(
        kwds, method, url_list, fnames, root_dir, file_prefix, file_extention
    )
    key_list = set(itertools.chain.from_iterable(k.keys() for k in kwd_list))
    valid_keys = ("params", "data", "json", "headers")
    if any(k not in valid_keys for k in key_list):
        raise InputValueError("kwds", valid_keys)

    if len(url_list) != len(kwd_list):
        raise InputTypeError("urls/kwds", "list of same length")

    url_list = cast("tuple[str, ...]", url_list)
    kwd_list = cast("tuple[dict[str, None | dict[Any, Any]], ...]", kwd_list)
    return url_list, kwd_list, files


def _download(
    session_func: Callable[[str], Response],
    url: str,
    kwd: Mapping[str, None | dict[Any, Any]],
    fname: Path,
    chunk_size: int,
) -> Path | None:
    """Download a single file."""
    try:
        resp = session_func(url, **kwd)
    except ServiceError as ex:
        msg = f"An error occurred when downloading {url}:\n{ex!s}"
        warnings.warn(msg, RuntimeWarning, stacklevel=2)
        return None
    fsize = int(resp.headers.get("Content-Length", -1))
    if not fname.exists() or fname.stat().st_size != fsize:
        fname.parent.mkdir(exist_ok=True, parents=True)
        with fname.open("wb") as f:
            f.writelines(resp.iter_content(chunk_size))
    return fname


@overload
def streaming_download(
    urls: str,
    kwds: dict[str, dict[Any, Any]] | None = None,
    fnames: str | Path | None = None,
    root_dir: str | Path | None = None,
    file_prefix: str = "",
    file_extention: str = "",
    method: Literal["GET", "POST", "get", "post"] = "GET",
    ssl: bool = True,
    chunk_size: int = CHUNK_SIZE,
    n_jobs: int = MAX_CONN,
) -> Path | None: ...


@overload
def streaming_download(
    urls: list[str],
    kwds: list[dict[str, dict[Any, Any]]] | None = None,
    fnames: Sequence[str | Path] | None = None,
    root_dir: str | Path | None = None,
    file_prefix: str = "",
    file_extention: str = "",
    method: Literal["GET", "POST", "get", "post"] = "GET",
    ssl: bool = True,
    chunk_size: int = CHUNK_SIZE,
    n_jobs: int = MAX_CONN,
) -> list[Path | None]: ...


def streaming_download(
    urls: list[str] | str,
    kwds: list[dict[str, dict[Any, Any]]] | dict[str, dict[Any, Any]] | None = None,
    fnames: str | Path | Sequence[str | Path] | None = None,
    root_dir: str | Path | None = None,
    file_prefix: str = "",
    file_extention: str = "",
    method: Literal["GET", "POST", "get", "post"] = "GET",
    ssl: bool = True,
    chunk_size: int = CHUNK_SIZE,
    n_jobs: int = MAX_CONN,
) -> Path | None | list[Path | None]:
    """Download and store files in parallel from a list of URLs/Keywords.

    Notes
    -----
    This function runs asynchronously in parallel using ``n_jobs`` threads.

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
    root_dir : str or Path, optional
        Root directory to store the files, defaults to ``None`` which
        uses HyRiver's cache directory. Note that you should either
        provide ``root_dir`` or ``fnames``. If both are provided,
        ``root_dir`` will be ignored.
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
        urls, kwds, method, fnames, root_dir, file_prefix, file_extention
    )

    session = RetrySession(disable=True, ssl=ssl)
    if method == "GET":
        func = tlz.partial(session.get, stream=True)
    else:
        func = tlz.partial(session.post, stream=True)

    n_jobs = min(n_jobs, len(url_list))
    fpaths = joblib.Parallel(n_jobs=n_jobs, prefer="threads")(
        joblib.delayed(_download)(func, u, k, f, chunk_size)
        for u, k, f in zip(url_list, kwd_list, files)
    )
    fpaths = cast("list[Path | None]", fpaths)
    session.close()
    if isinstance(urls, str):
        return fpaths[0]
    return fpaths


def traverse_json(json_data: dict[str, Any] | list[dict[str, Any]], ipath: list[str]) -> list[Any]:
    """Extract an element from a JSON-like object along a specified ipath.

    This function is based on
    `bcmullins <https://bcmullins.github.io/parsing-json-python/>`__.

    Parameters
    ----------
    json_data : dict or list of dicts
        The input json dictionary.
    ipath : list
        The ipath to the requested element.

    Returns
    -------
    list
        The sub-items founds in the JSON.

    Examples
    --------
    >>> data = [
    ...     {"employees": [
    ...         {"name": "Alice", "role": "dev", "nbr": 1},
    ...         {"name": "Bob", "role": "dev", "nbr": 2},
    ...         ],},
    ...     {"firm": {"name": "Charlie's Waffle Emporium", "location": "CA"}},
    ... ]
    >>> traverse_json(data, ["employees", "name"])
    [['Alice', 'Bob'], [None]]
    """

    def extract(
        sub_items: list[Any] | dict[str, Any] | None,
        path: list[str],
        ind: int,
        arr: list[Any],
    ) -> list[Any]:
        key = path[ind]
        if isinstance(sub_items, dict):
            if key in sub_items:
                if ind + 1 == len(path):
                    arr.append(sub_items[key])
                else:
                    extract(sub_items[key], path, ind + 1, arr)
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
        return arr

    if isinstance(json_data, dict):
        return extract(json_data, ipath, 0, [])

    outer_arr = [extract(item, ipath, 0, []) for item in json_data]
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

    def _get_payload(
        self,
        geo_type: Literal[
            "esriGeometryPoint",
            "esriGeometryMultipoint",
            "esriGeometryEnvelope",
            "esriGeometryPolygon",
            "esriGeometryPolyline",
        ],
        geo_json: dict[str, Any],
    ) -> Mapping[str, str]:
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
        esri_json = json.dumps(geo_json | {"spatialRelference": {"wkid": str(self.wkid)}})
        return {
            "geometryType": geo_type,
            "geometry": esri_json,
            "inSR": str(self.wkid),
        }

    def point(self) -> Mapping[str, str]:
        """Query for a point."""
        try:
            point = Point(*self.geometry)
        except (TypeError, AttributeError, ValueError) as ex:
            raise InputTypeError("geometry", "tuple", "(x, y)") from ex
        point = shapely.set_precision(point, 1e-6)
        geo_type = "esriGeometryPoint"
        geo_json = {"x": point.x, "y": point.y}
        return self._get_payload(geo_type, geo_json)

    def multipoint(self) -> Mapping[str, str]:
        """Query for a multi-point."""
        try:
            mp = MultiPoint(self.geometry)
        except (TypeError, AttributeError, ValueError) as ex:
            raise InputTypeError("geometry", "list of tuples", "[(x, y), ...]") from ex
        mp = shapely.set_precision(mp, 1e-6)
        geo_type = "esriGeometryMultipoint"
        geo_json = {"points": [[p.x, p.y] for p in mp.geoms]}
        return self._get_payload(geo_type, geo_json)

    def bbox(self) -> Mapping[str, str]:
        """Query for a bbox."""
        try:
            bbox = shapely.box(*self.geometry)  # pyright: ignore[reportArgumentType]
        except (TypeError, AttributeError, ValueError) as ex:
            raise InputTypeError("geometry", "tuple", BOX_ORD) from ex
        bbox = shapely.set_precision(bbox, 1e-6)
        geo_type = "esriGeometryEnvelope"
        geo_json = dict(zip(("xmin", "ymin", "xmax", "ymax"), bbox.bounds))
        return self._get_payload(geo_type, geo_json)  # pyright: ignore[reportArgumentType]

    def polygon(self) -> Mapping[str, str]:
        """Query for a polygon."""
        if isinstance(self.geometry, Polygon):
            geom = shapely.set_precision(self.geometry, 1e-6)
            geo_type = "esriGeometryPolygon"
            geo_json = {"rings": [[[x, y] for x, y in zip(*geom.exterior.coords.xy)]]}
            return self._get_payload(geo_type, geo_json)

        raise InputTypeError("geometry", "Polygon")

    def polyline(self) -> Mapping[str, str]:
        """Query for a polyline."""
        if isinstance(self.geometry, LineString):
            geom = shapely.set_precision(self.geometry, 1e-6)
            geo_type = "esriGeometryPolyline"
            geo_json = {"paths": [[[x, y] for x, y in zip(*geom.coords.xy)]]}
            return self._get_payload(geo_type, geo_json)

        raise InputTypeError("geometry", "LineString")


def match_crs(geom: GeomType, in_crs: CRSType, out_crs: CRSType) -> GeomType:
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
    >>> from shapely import Point
    >>> point = Point(-7766049.665, 5691929.739)
    >>> match_crs(point, 3857, 4326).xy
    (array('d', [-69.7636111130079]), array('d', [45.44549114818127]))
    >>> bbox = (-7766049.665, 5691929.739, -7763049.665, 5696929.739)
    >>> match_crs(bbox, 3857, 4326)
    (-69.7636111130079, 45.44549114818127, -69.73666165448431, 45.47699468552394)
    >>> coords = [(-7766049.665, 5691929.739)]
    >>> match_crs(coords, 3857, 4326)
    [(-69.7636111130079, 45.44549114818127)]
    """
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
        if pyproj.CRS(in_crs) == pyproj.CRS(out_crs):
            return geom
        return ops.transform(project, geom)

    with contextlib.suppress(TypeError, AttributeError, ValueError):
        if len(geom) > 4:
            raise TypeError
        if pyproj.CRS(in_crs) == pyproj.CRS(out_crs):
            bbox = shapely.box(*geom)  # pyright: ignore[reportArgumentType]
        else:
            bbox = ops.transform(project, shapely.box(*geom))  # pyright: ignore[reportArgumentType]
        bbox = cast("Polygon", bbox)
        return tuple(float(p) for p in bbox.bounds)  # pyright: ignore[reportReturnType]

    with contextlib.suppress(TypeError, AttributeError, ValueError):
        if pyproj.CRS(in_crs) == pyproj.CRS(out_crs):
            point = Point(geom)
        else:
            point = cast("Point", ops.transform(project, Point(geom)))
        return [(float(point.x), float(point.y))]  # pyright: ignore[reportReturnType]

    with contextlib.suppress(TypeError, AttributeError, ValueError):
        if pyproj.CRS(in_crs) == pyproj.CRS(out_crs):
            mp = MultiPoint(geom)
        else:
            mp = ops.transform(project, MultiPoint(geom))
        return [(float(p.x), float(p.y)) for p in mp.geoms]  # pyright: ignore[reportReturnType]

    gtypes = " ".join(
        (
            "a list of coordinates such as [(x1, y1), ...],",
            "a bounding box like so (xmin, ymin, xmax, ymax),",
            "or any valid shapely's geometry.",
        )
    )
    raise InputTypeError("geom", gtypes)


def esri_query(
    geom: GeomType,
    geo_crs: CRSType,
    out_crs: CRSType,
) -> Mapping[str, str]:
    """Generate geometry queries based on ESRI template."""
    geom = match_crs(geom, geo_crs, out_crs)
    wkid = cast("int", pyproj.CRS(out_crs).to_epsg())

    with contextlib.suppress(TypeError, AttributeError, ValueError):
        geom = Point(geom)  # pyright: ignore[reportAssignmentType]

    with contextlib.suppress(TypeError, AttributeError, ValueError):
        geom = MultiPoint(geom)  # pyright: ignore[reportAssignmentType]

    if isinstance(geom, Point):
        return ESRIGeomQuery((geom.x, geom.y), wkid).point()

    if isinstance(geom, MultiPoint):
        return ESRIGeomQuery([(g.x, g.y) for g in geom.geoms], wkid).multipoint()

    if isinstance(geom, Polygon):
        return ESRIGeomQuery(geom, wkid).polygon()

    if isinstance(geom, LineString):
        return ESRIGeomQuery(geom, wkid).polyline()

    with contextlib.suppress(InputTypeError):
        return ESRIGeomQuery(geom, wkid).bbox()

    raise InputTypeError("geom", "LineString, Polygon, Point, MultiPoint, tuple")


def check_bbox(bbox: tuple[float, float, float, float]) -> None:
    """Check if an input inbox is a tuple of length 4."""
    try:
        _ = shapely.box(*bbox)
    except (TypeError, AttributeError, ValueError) as ex:
        raise InputTypeError("bbox", "tuple", BOX_ORD) from ex


def bbox_decompose(
    bbox: tuple[float, float, float, float],
    resolution: float,
    box_crs: CRSType = 4326,
    max_px: int = 8_000_000,
) -> tuple[list[tuple[float, float, float, float]], int, int]:
    """Split the bounding box vertically for WMS requests.

    Parameters
    ----------
    bbox : tuple
        A bounding box; (west, south, east, north)
    resolution : float
        The target resolution for a WMS request in meters.
    box_crs : str, int, or pyproj.CRS, optional
        The spatial reference of the input bbox, default to ``epsg:4326``.
    max_px : int, optional
        The maximum allowable number of pixels (width x height) for a WMS requests,
        defaults to 8 million based on some trial-and-error.

    Returns
    -------
    boxes : list of tuple
        List of sub-bboxes in the form (west, south, east, north).
    sub_width : int
        Width of each sub-bbox in degrees.
    sub_height : int
        Height of each sub-bbox in degrees.
    """
    check_bbox(bbox)

    geod = pyproj.Geod(ellps="GRS80")

    west, south, east, north = bbox

    xmin, ymin, xmax, ymax = match_crs(bbox, box_crs, 4326)

    x_dist = geod.geometry_length(LineString([(xmin, ymin), (xmax, ymin)]))
    y_dist = geod.geometry_length(LineString([(xmin, ymin), (xmin, ymax)]))

    width = math.ceil(x_dist / resolution)
    height = math.ceil(y_dist / resolution)
    if width * height <= max_px:
        return [bbox], width, height

    # Divisions in each direction maintaining aspect ratio
    aspect_ratio = width / height
    n_boxes = math.ceil((width * height) / max_px)
    nx = math.ceil(math.sqrt(n_boxes * aspect_ratio))
    ny = math.ceil(n_boxes / nx)
    dx = (east - west) / nx
    dy = (north - south) / ny

    # Calculate buffer sizes in degrees
    sub_width = math.ceil(width / nx)
    sub_height = math.ceil(height / ny)

    bboxs = []
    for i in range(nx):
        box_west = west + (i * dx)
        box_east = min(west + ((i + 1) * dx), east)
        for j in range(ny):
            box_south = south + (j * dy)
            box_north = min(south + ((j + 1) * dy), north)
            bboxs.append((box_west, box_south, box_east, box_north))
    return bboxs, sub_width, sub_height


def validate_crs(crs: CRSType) -> str:
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
        return pyproj.CRS(crs).to_string()
    except ProjCRSError as ex:
        raise InputTypeError("crs", "a valid CRS") from ex


def valid_wms_crs(url: str) -> list[str]:
    """Get valid CRSs from a WMS service version 1.3.0."""
    ns = "http://www.opengis.net/wms"

    def get_path(tag_list: list[str]) -> str:
        return f"/{{{ns}}}".join(["", *tag_list])[1:]

    kwds = {"params": {"service": "wms", "request": "GetCapabilities"}}
    try:
        root = ETree.fromstring(ar.retrieve_text([url], [kwds], ssl=False)[0])
    except ETree.ParseError as ex:
        raise ServiceUnavailableError(url) from ex
    return [
        t.text.lower()  # pyright: ignore[reportOptionalMemberAccess]
        for t in root.findall(get_path(["Capability", "Layer", "CRS"]))
    ]
