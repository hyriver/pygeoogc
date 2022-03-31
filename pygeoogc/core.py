"""Base classes and function for REST, WMS, and WMF services."""
import contextlib
import logging
import os
import sys
import urllib.parse as urlparse
import uuid
from pathlib import Path
from typing import Any, Dict, Iterator, List, Mapping, Optional, Sequence, Tuple, Union

import async_retriever as ar
import cytoolz as tlz
import pyproj
from owslib.wfs import WebFeatureService
from owslib.wms import WebMapService
from pydantic import BaseModel, validator
from shapely.geometry import LineString, MultiPoint, Point, Polygon

from . import utils
from .exceptions import (
    InvalidInputType,
    InvalidInputValue,
    MissingInputs,
    ServiceError,
    ServiceUnavailable,
    ZeroMatched,
)

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
handler = logging.StreamHandler(sys.stdout)
handler.setFormatter(logging.Formatter(""))
logger.handlers = [handler]
logger.propagate = False
DEF_CRS = "epsg:4326"
EXPIRE = -1


def validate_version(val: str, valid_versions: List[str]) -> str:
    """Validate version from a list of valid versions.

    Parameters
    ----------
    val : str
        Input version value.
    valid_versions : list of str
        List of valid versions.

    Returns
    -------
    str
        Validated version value.
    """
    if val not in valid_versions:
        raise InvalidInputValue("version", valid_versions)
    return val


class RESTValidator(BaseModel):
    """Validate ArcGISRESTful inputs.

    Parameters
    ----------
    base_url : str, optional
        The ArcGIS RESTful service url. The URL must either include a layer number
        after the last ``/`` in the url or the target layer must be passed as an argument.
    layer : int, optional
        Target layer number, defaults to None. If None layer number must be included as after
        the last ``/`` in ``base_url``.
    outformat : str, optional
        One of the output formats offered by the selected layer. If not correct
        a list of available formats is shown, defaults to ``geojson``.
    outfields : str or list
        The output fields to be requested. Setting ``*`` as outfields requests
        all the available fields which is the default setting.
    crs : str, optional
        The spatial reference of the output data, defaults to EPSG:4326
    max_workers : int, optional
        Max number of simultaneous requests, default to 2. Note
        that some services might face issues when several requests are sent
        simultaneously and will return the requests partially. It's recommended
        to avoid using too many workers unless you are certain the web service
        can handle it.
    verbose : bool, optional
        If True, prints information about the requests and responses,
        defaults to False.
    disable_retry : bool, optional
        If ``True`` in case there are any failed queries, no retrying attempts
        is done and object IDs of the failed requests is saved to a text file
        which its path can be accessed via ``self.failed_path``.
    """

    base_url: str
    layer: Optional[int] = None
    outformat: str = "geojson"
    outfields: Union[List[str], str] = "*"
    crs: str = DEF_CRS
    max_workers: int = 1
    verbose: bool = False
    disable_retry: bool = False

    @staticmethod
    def _split_url(url: str) -> Tuple[str, Optional[int]]:
        """Check if layer is included in url, if so separate and return them."""
        url_split = urlparse.urlsplit(url.rstrip("/"))
        url_path, lyr = url_split.path.rsplit("/", 1)
        url_split = url_split._replace(path=url_path)
        try:
            return urlparse.urlunsplit(url_split), int(lyr)
        except ValueError:
            return url, None

    @validator("base_url")
    def _layer_from_url(cls, v: str) -> Tuple[str, Optional[int]]:
        return cls._split_url(urlparse.unquote(v))

    @validator("layer")
    def _integer_layer(cls, v: Optional[int], values: Dict[str, Any]) -> int:
        if v is None:
            try:
                return int(values["base_url"][1])
            except TypeError as ex:
                msg = "Either layer must be passed as an argument or be included in ``base_url``"
                raise MissingInputs(msg) from ex
        return v

    @validator("outfields")
    def _outfields_to_list(cls, v: Union[List[str], str]) -> List[str]:
        return v if isinstance(v, list) else [v]

    @validator("crs")
    def _valid_crs(cls, v: str) -> str:
        return utils.validate_crs(v)

    @validator("max_workers")
    def _positive_integer_threads(cls, v: int) -> int:
        if v <= 0:
            raise InvalidInputType("max_workers", "positive integer")
        return v


class ArcGISRESTfulBase:
    """Access to an ArcGIS REST service.

    Parameters
    ----------
    base_url : str, optional
        The ArcGIS RESTful service url. The URL must either include a layer number
        after the last ``/`` in the url or the target layer must be passed as an
        argument.
    layer : int, optional
        Target layer number, defaults to None. If None layer number must be
        included as after the last ``/`` in ``base_url``.
    outformat : str, optional
        One of the output formats offered by the selected layer. If not correct
        a list of available formats is shown, defaults to ``geojson``.
        It defaults to ``esriSpatialRelIntersects``.
    outfields : str or list
        The output fields to be requested. Setting ``*`` as outfields requests
        all the available fields which is the default setting.
    crs : str, optional
        The spatial reference of the output data, defaults to EPSG:4326
    max_workers : int, optional
        Max number of simultaneous requests, default to 2. Note
        that some services might face issues when several requests are sent
        simultaneously and will return the requests partially. It's recommended
        to avoid using too many workers unless you are certain the web service
        can handle it.
    verbose : bool, optional
        If True, prints information about the requests and responses,
        defaults to False.
    disable_retry : bool, optional
        If ``True`` in case there are any failed queries, no retrying attempts
        is done and object IDs of the failed requests is saved to a text file
        which its path can be accessed via ``self.failed_path``.
    """

    def __init__(
        self,
        base_url: str,
        layer: Optional[int] = None,
        outformat: str = "geojson",
        outfields: Union[List[str], str] = "*",
        crs: str = DEF_CRS,
        max_workers: int = 1,
        verbose: bool = False,
        disable_retry: bool = False,
    ) -> None:
        validated = RESTValidator(
            base_url=base_url,
            layer=layer,
            outformat=outformat,
            outfields=outfields,
            crs=crs,
            max_workers=max_workers,
            verbose=verbose,
            disable_retry=disable_retry,
        )
        self.base_url = validated.base_url[0]
        self.layer = validated.layer
        self.outformat = validated.outformat
        self.outfields = validated.outfields
        self.crs = validated.crs
        self.max_workers = validated.max_workers
        self.verbose = validated.verbose
        self.disable_retry = validated.disable_retry

        self.out_sr = pyproj.CRS(self.crs).to_epsg()
        self.n_features = 0
        self.url = f"{self.base_url}/{self.layer}"
        self.query_url = f"{self.url}/query"
        self.valid_layers: Dict[str, str] = {}
        self.query_formats: List[str] = []
        self.extent: Optional[Tuple[float, float, float, float]] = None
        self.units: Optional[str] = None
        self.max_nrecords: int = 1000
        self.valid_fields: List[str] = []
        self.field_types: Dict[str, str] = {}
        self.feature_types: Optional[Dict[int, str]] = None
        self.return_m: bool = False
        self.return_geom: bool = True
        self.n_missing: int = 0
        self.total_n_features: int = 0
        self.failed_path: Union[str, Path] = ""
        self.request_id: Optional[str] = None

        self.initialize_service()

    def _set_service_properties(self) -> None:
        rjson = self.get_response(self.base_url, [{"f": "json"}])[0]
        try:
            self.valid_layers = {f"{lyr['id']}": lyr["name"] for lyr in rjson["layers"]}
            self.query_formats = rjson["supportedQueryFormats"].replace(" ", "").lower().split(",")

            extent = rjson["extent"] if "extent" in rjson else rjson["fullExtent"]
            bounds = (extent["xmin"], extent["ymin"], extent["xmax"], extent["ymax"])
        except (ValueError, KeyError) as ex:
            raise ServiceError(self.base_url) from ex

        crs_iter = iter(["latestWkid", "wkid", "wkt"])
        while True:
            try:
                crs = utils.validate_crs(extent["spatialReference"].get(next(crs_iter)))
                self.extent = utils.match_crs(bounds, crs, DEF_CRS)
                break
            except InvalidInputType:
                continue
            except StopIteration as ex:
                raise ServiceError(self.base_url) from ex

        with contextlib.suppress(KeyError):
            self.units = rjson["units"].replace("esri", "").lower()
            self.max_nrecords = int(rjson["maxRecordCount"])

    def _set_layer_properties(self) -> None:
        """Set properties of the target layer."""
        rjson = self.get_response(self.url, [{"f": "json"}])[0]
        self.valid_fields = list(
            set(
                utils.traverse_json(rjson, ["fields", "name"])
                + utils.traverse_json(rjson, ["fields", "alias"])
                + ["*"]
            )
        )
        self.field_types = {
            f["name"]: f["type"].replace("esriFieldType", "").lower() for f in rjson["fields"]
        }

        with contextlib.suppress(KeyError):
            self.feature_types = dict(
                zip((tlz.pluck("id", rjson["types"])), tlz.pluck("name", rjson["types"]))
            )

    def initialize_service(self) -> None:
        """Initialize the RESTFul service."""
        self._set_service_properties()
        if f"{self.layer}" not in self.valid_layers:
            valids = [f'"{i}" for {n}' for i, n in self.valid_layers.items()]
            raise InvalidInputValue("layer", valids)

        self._set_layer_properties()

        if self.outformat.lower() not in self.query_formats:
            raise InvalidInputValue("outformat", self.query_formats)

        if any(f not in self.valid_fields for f in self.outfields):
            raise InvalidInputValue("outfields", self.valid_fields)

    def partition_oids(self, oids: Union[List[int], int]) -> Iterator[Tuple[str, ...]]:
        """Partition feature IDs based on ``self.max_nrecords``."""
        oid_list = [oids] if isinstance(oids, int) else set(oids)
        if len(oid_list) == 0:
            raise ZeroMatched

        self.n_features = len(oid_list)
        if not self.disable_retry and self.verbose:
            logger.info(f"Found {self.n_features:,} features in the requested region.")
        return tlz.partition_all(self.max_nrecords, [str(i) for i in oid_list])  # type: ignore

    def _cleanup_resp(
        self, resp: List[Dict[str, Any]], payloads: Sequence[Dict[str, str]]
    ) -> List[Dict[str, Any]]:
        """Remove failed responses."""
        fails = [i for i, r in enumerate(resp) if "error" in r]

        if len(fails) > 0:
            err = resp[fails[0]]["error"]["message"]
            resp = [r for i, r in enumerate(resp) if i not in fails]
            if len(resp) == 0 and self.disable_retry:
                raise ServiceError(err)

            if "objectIds" in payloads[fails[0]]:
                oids = list(tlz.concat(payloads[i]["objectIds"].split(",") for i in fails))
                self.n_missing = len(oids)

                self.failed_path = Path("cache", f"failed_ids_{self.request_id}.txt")
                with open(self.failed_path, "w") as f:
                    f.write("\n".join(oids))

                if not self.disable_retry:
                    self.return_m = bool(payloads[0]["ReturnM"])
                    self.return_geom = bool(payloads[0]["returnGeometry"])
                    self.total_n_features = self.n_features
                    logger.warning(f"Found {len(fails)} failed requests. Retrying ...")
                    resp.extend(self.retry_failed_requests())

                    logger.warning(
                        " ".join(
                            [
                                f"Total of {self.n_missing} out of {self.total_n_features}",
                                "requested features are not available in the dataset.",
                                "Returning the successfully retrieved features.",
                                "The failed object IDs have been saved in the",
                                f"file {self.failed_path}. The service returned the",
                                "following error message for the failed requests:\n",
                                err,
                            ]
                        )
                    )

        return resp

    def get_features(
        self,
        featureids: Iterator[Tuple[str, ...]],
        return_m: bool = False,
        return_geom: bool = True,
    ) -> List[Dict[str, Any]]:
        """Get features based on the feature IDs.

        Parameters
        ----------
        featureids : list
            List of feature IDs.
        return_m : bool, optional
            Whether to activate the Return M (measure) in the request,
            defaults to ``False``.
        return_geom : bool, optional
            Whether to return the geometry of the feature, defaults to ``True``.

        Returns
        -------
        dict
            (Geo)json response from the web service.
        """
        payloads = [
            {
                "objectIds": ",".join(ids),
                "returnGeometry": f"{return_geom}".lower(),
                "outSR": f"{self.out_sr}",
                "outfields": ",".join(self.outfields),
                "ReturnM": f"{return_m}".lower(),
                "f": self.outformat,
            }
            for ids in featureids
        ]
        if self.request_id is None:
            self.request_id = uuid.uuid4().hex

        resp = self.get_response(self.query_url, payloads, "POST")

        resp = self._cleanup_resp(resp, payloads)
        if len(resp) == 0:
            raise ZeroMatched
        return resp

    def esri_query(
        self,
        geom: Union[
            LineString,
            Polygon,
            Point,
            MultiPoint,
            Tuple[float, float, float, float],
        ],
        geo_crs: Union[str, pyproj.CRS] = DEF_CRS,
    ) -> Mapping[str, str]:
        """Generate geometry queries based on ESRI template."""
        geom = utils.match_crs(geom, geo_crs, self.crs)

        if isinstance(geom, tuple) and len(geom) == 4:
            return utils.ESRIGeomQuery(geom, self.out_sr).bbox()

        if isinstance(geom, Point):
            return utils.ESRIGeomQuery((geom.x, geom.y), self.out_sr).point()

        if isinstance(geom, MultiPoint):
            return utils.ESRIGeomQuery([(g.x, g.y) for g in geom.geoms], self.out_sr).multipoint()

        if isinstance(geom, Polygon):
            return utils.ESRIGeomQuery(geom, self.out_sr).polygon()

        if isinstance(geom, LineString):
            return utils.ESRIGeomQuery(geom, self.out_sr).polyline()

        raise InvalidInputType("geom", "LineString, Polygon, Point, MultiPoint, tuple")

    def _retry(
        self, return_m: bool, return_geo: bool, partition_fac: float
    ) -> List[Dict[str, Any]]:
        """Retry failed requests."""
        with open(self.failed_path) as f:
            oids = [int(i) for i in f.read().splitlines()]

        max_nrecords = self.max_nrecords
        self.max_nrecords = int(max(partition_fac * max_nrecords, 1))
        features = self.get_features(self.partition_oids(oids), return_m, return_geo)
        self.max_nrecords = max_nrecords

        return features

    def retry_failed_requests(self) -> List[Dict[str, Any]]:
        """Retry failed requests."""
        retry = self.disable_retry
        self.disable_retry = True
        caching = os.getenv("HYRIVER_CACHE_DISABLE", "false").lower() == "true"
        os.environ["HYRIVER_CACHE_DISABLE"] = "true"

        fac_min = 1.0 / self.max_nrecords
        n_retry = 4
        partition_fac = [(0.5 - fac_min) / (n_retry - 1) * i + fac_min for i in range(n_retry)]
        features = []
        for f in partition_fac[::-1]:
            try:
                features.append(self._retry(self.return_m, self.return_geom, f))
            except ServiceError:
                continue

        self.disable_retry = retry
        os.environ["HYRIVER_CACHE_DISABLE"] = f"{caching}".lower()
        return list(tlz.concat(features))

    def get_response(
        self, url: str, payloads: Sequence[Dict[str, str]], method: str = "GET"
    ) -> List[Dict[str, Any]]:
        """Send payload and get the response."""
        req_key = "params" if method == "GET" else "data"
        try:
            return ar.retrieve_json(
                [url] * len(payloads),
                [{req_key: p} for p in payloads],
                request_method=method,
                max_workers=self.max_workers,
            )
        except ValueError:
            raise ZeroMatched

    def __repr__(self) -> str:
        """Print the service configuration."""
        extent = ", ".join(f"{c:.3f}" for c in self.extent) if self.extent else None
        ftypes = ", ".join(self.feature_types.values()) if self.feature_types else None
        return "\n".join(
            [
                "Service configurations:",
                f"    URL: {self.base_url}",
                f"    Layer: {self.valid_layers[str(self.layer)]} ({self.layer})",
                f"    Max record count: {self.max_nrecords}",
                f"    Query format: {self.outformat}",
                f"    Units: {self.units}",
                f"    Extent: ({extent})",
                f"    Feature types: {ftypes}",
            ]
        )


class WMSBase(BaseModel):
    """Base class for accessing a WMS service.

    Parameters
    ----------
    url : str
        The base url for the WMS service e.g., https://www.mrlc.gov/geoserver/mrlc_download/wms
    layers : str or list
        A layer or a list of layers from the service to be downloaded. You can pass an empty
        string to get a list of available layers.
    outformat : str
        The data format to request for data from the service. You can pass an empty
        string to get a list of available output formats.
    version : str, optional
        The WMS service version which should be either 1.1.1 or 1.3.0, defaults to 1.3.0.
    crs : str, optional
        The spatial reference system to be used for requesting the data, defaults to
        ``epsg:4326``.
    """

    url: str
    layers: Union[str, List[str]]
    outformat: str
    version: str = "1.3.0"
    crs: str = DEF_CRS

    @validator("crs")
    def _valid_crs(cls, v: str) -> str:
        return utils.validate_crs(v)

    @validator("version")
    def _version(cls, v: str) -> str:
        return validate_version(v, ["1.1.1", "1.3.0"])

    def validate_wms(self) -> None:
        """Validate input arguments with the WMS service."""
        try:
            wms = WebMapService(self.url, version=self.version)
        except AttributeError as ex:
            raise ServiceUnavailable(self.url) from ex

        layers = [self.layers] if isinstance(self.layers, str) else self.layers
        valid_layers = {wms[lyr].name: wms[lyr].title for lyr in list(wms.contents)}
        if any(lyr not in valid_layers.keys() for lyr in layers):
            raise InvalidInputValue("layers", (f"{n} for {t}" for n, t in valid_layers.items()))

        valid_outformats = wms.getOperationByName("GetMap").formatOptions
        if self.outformat not in valid_outformats:
            raise InvalidInputValue("outformat", valid_outformats)

        valid_crss = {lyr: [s.lower() for s in wms[lyr].crsOptions] for lyr in layers}
        if any(self.crs.lower() not in valid_crss[lyr] for lyr in layers):
            _valid_crss = (f"{lyr}: {', '.join(cs)}\n" for lyr, cs in valid_crss.items())
            raise InvalidInputValue("CRS", _valid_crss)

    def get_validlayers(self) -> Dict[str, str]:
        """Get the layers supported by the WMS service."""
        try:
            wms = WebMapService(self.url, version=self.version)
        except AttributeError as ex:
            raise ServiceUnavailable(self.url) from ex

        return {wms[lyr].name: wms[lyr].title for lyr in list(wms.contents)}

    def __repr__(self) -> str:
        """Print the services properties."""
        layers = self.layers if isinstance(self.layers, list) else [self.layers]
        return (
            "Connected to the WMS service with the following properties:\n"
            + f"URL: {self.url}\n"
            + f"Version: {self.version}\n"
            + f"Layers: {', '.join(lyr for lyr in layers)}\n"
            + f"Output Format: {self.outformat}\n"
            + f"Output CRS: {self.crs}"
        )


class WFSBase(BaseModel):
    """Base class for WFS service.

    Parameters
    ----------
    url : str
        The base url for the WFS service, for examples:
        https://hazards.fema.gov/nfhl/services/public/NFHL/MapServer/WFSServer
    layer : str
        The layer from the service to be downloaded, defaults to None which throws
        an error and includes all the available layers offered by the service.
    outformat : str
        The data format to request for data from the service, defaults to None which
         throws an error and includes all the available format offered by the service.
    version : str, optional
        The WFS service version which should be either ``1.0.0``, ``1.1.0``, or
        ``2.0.0``. Defaults to ``2.0.0``.
    crs: str, optional
        The spatial reference system to be used for requesting the data, defaults to
        ``epsg:4326``.
    read_method : str, optional
        Method for reading the retrieved data, defaults to ``json``. Valid options are
        ``json``, ``binary``, and ``text``.
    max_nrecords : int, optional
        The maximum number of records in a single request to be retrieved from the service,
        defaults to 1000. If the number of requested records is greater than this value,
        the query will be split into multiple requests.
    """

    url: str
    layer: Optional[str] = None
    outformat: Optional[str] = None
    version: str = "2.0.0"
    crs: str = DEF_CRS
    read_method: str = "json"
    max_nrecords: int = 1000

    @validator("read_method")
    def _read_method(cls, v: str) -> str:
        valid_methods = ["json", "binary", "text"]
        if v not in valid_methods:
            raise InvalidInputValue("read_method", valid_methods)
        return v

    @validator("crs")
    def _valid_crs(cls, v: str) -> str:
        return utils.validate_crs(v)

    @validator("version")
    def _version(cls, v: str) -> str:
        return validate_version(v, ["1.0.0", "1.1.0", "2.0.0"])

    def validate_wfs(self) -> None:
        """Validate input arguments with the WFS service."""
        try:
            wfs = WebFeatureService(self.url, version=self.version)
        except AttributeError as ex:
            raise ServiceUnavailable(self.url) from ex

        valid_layers = list(wfs.contents)
        if self.layer is None:
            raise MissingInputs(
                "The layer argument is missing."
                + " The following layers are available:\n"
                + ", ".join(valid_layers)
            )
        if self.layer not in valid_layers:
            raise InvalidInputValue("layers", valid_layers)

        wfs_features = wfs.getOperationByName("GetFeature")
        if self.version == "1.0.0":
            valid_outformats = [f.rsplit("}", 1)[-1] for f in wfs_features.formatOptions]
        else:
            valid_outformats = wfs_features.parameters["outputFormat"]["values"]

        valid_outformats = [v.lower() for v in valid_outformats]
        if self.outformat is None:
            raise MissingInputs(
                "The outformat argument is missing."
                + " The following output formats are available:\n"
                + ", ".join(valid_outformats)
            )

        if self.outformat not in valid_outformats:
            raise InvalidInputValue("outformat", valid_outformats)

        valid_crss = [f"{s.authority.lower()}:{s.code}" for s in wfs[self.layer].crsOptions]
        if self.crs.lower() not in valid_crss:
            raise InvalidInputValue("crs", valid_crss)

    def get_validnames(self) -> List[str]:
        """Get valid column names for a layer."""
        max_features = "count" if self.version == "2.0.0" else "maxFeatures"

        payload = {
            "service": "wfs",
            "version": self.version,
            "outputFormat": "application/json",
            "request": "GetFeature",
            "typeName": self.layer,
            max_features: 1,
        }
        rjson = ar.retrieve_json([self.url], [{"params": payload}])[0]
        valid_fields = list(
            set(
                utils.traverse_json(rjson, ["fields", "name"])
                + utils.traverse_json(rjson, ["fields", "alias"])
                + ["*"]
            )
        )

        if None in valid_fields:
            valid_fields = list(utils.traverse_json(rjson, ["features", "properties"])[0].keys())

        return valid_fields

    def __repr__(self) -> str:
        """Print the services properties."""
        return "\n".join(
            [
                "Connected to the WFS service with the following properties:",
                f"URL: {self.url}",
                f"Version: {self.version}",
                f"Layer: {self.layer}",
                f"Output Format: {self.outformat}",
                f"Output CRS: {self.crs}",
                f"Read Method: {self.read_method}",
            ]
        )
