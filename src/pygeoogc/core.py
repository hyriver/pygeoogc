"""Base classes and function for REST, WMS, and WMF services."""

from __future__ import annotations

import contextlib
import math
import os
import uuid
import warnings
from pathlib import Path
from typing import TYPE_CHECKING, Any, Literal, cast

import cytoolz.curried as tlz
import pyproj
from owslib.wfs import WebFeatureService
from owslib.wms import WebMapService
from yarl import URL

import async_retriever as ar
from pygeoogc import utils
from pygeoogc.exceptions import (
    InputTypeError,
    InputValueError,
    MissingInputError,
    ServiceError,
    ServiceUnavailableError,
    ZeroMatchedError,
)

if TYPE_CHECKING:
    from collections.abc import Iterator, Sequence

    from pyproj import CRS

    CRSType = int | str | CRS

__all__ = ["ArcGISRESTfulBase", "WFSBase", "WMSBase"]


def _extract_layer(url: str, layer: int | None) -> tuple[str, int]:
    """Check if layer is included in url, if so separate and return them."""
    url_obj = URL(url.removesuffix("/"))
    if layer is None:
        try:
            return url_obj.parent.human_repr(), int(url_obj.parts[-1])
        except ValueError as ex:
            msg = "Either layer must be passed as an argument or be included in ``base_url``"
            raise MissingInputError(msg) from ex
    return url_obj.human_repr(), layer


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
    crs : str, int, or pyproj.CRS, optional
        The spatial reference of the output data, defaults to ``epsg:4326``
    verbose : bool, optional
        If True, prints information about the requests and responses,
        defaults to False.
    disable_retry : bool, optional
        If ``True`` in case there are any failed queries, no retrying attempts
        is done and object IDs of the failed requests is saved to a text file
        which its ipath can be accessed via ``self.failed_path``.
    """

    def __init__(
        self,
        base_url: str,
        layer: int | None = None,
        outformat: str = "geojson",
        outfields: list[str] | str = "*",
        crs: CRSType = 4326,
        verbose: bool = False,
        disable_retry: bool = False,
    ) -> None:
        self.base_url, self.layer = _extract_layer(base_url, layer)
        self.outformat = outformat
        self.outfields = outfields if isinstance(outfields, (list, tuple)) else [outfields]
        self.crs = utils.validate_crs(crs)
        self.verbose = verbose
        self.disable_retry = disable_retry

        self.out_sr = cast("int", pyproj.CRS(self.crs).to_epsg())
        self.n_features = 0
        self.url = f"{self.base_url}/{self.layer}"
        self.query_url = f"{self.url}/query"
        self.valid_layers: dict[str, str] = {}
        self.query_formats: list[str] = []
        self.extent: tuple[float, float, float, float] | None = None
        self.units: str | None = None
        self.max_nrecords = 1000
        self.valid_fields: list[str] = []
        self.field_types: dict[str, str] = {}
        self.feature_types: dict[int, str] | None = None
        self.return_m = False
        self.return_geom = True
        self.n_missing = 0
        self.total_n_features = 0
        self.failed_path = Path("cache", "failed_ids.txt")
        self.request_id: str | None = None

        self.initialize_service()

    def _set_service_properties(self) -> None:
        rjson = self.get_response(self.base_url, [{"f": "json"}])[0]
        try:
            self.valid_layers = {str(lyr["id"]): lyr["name"] for lyr in rjson["layers"]}
            self.query_formats = rjson["supportedQueryFormats"].replace(" ", "").lower().split(",")

            extent = rjson["extent"] if "extent" in rjson else rjson["fullExtent"]
            bounds = (extent["xmin"], extent["ymin"], extent["xmax"], extent["ymax"])
            if any(math.isnan(float(b)) for b in bounds):
                if rjson.get("initialExtent"):
                    extent = rjson["initialExtent"]
                    bounds = (extent["xmin"], extent["ymin"], extent["xmax"], extent["ymax"])
                else:
                    raise ServiceError(self.base_url)
        except (ValueError, KeyError) as ex:
            raise ServiceError(self.base_url) from ex

        crs_iter = iter(["latestWkid", "wkid", "wkt"])
        while True:
            try:
                crs = utils.validate_crs(extent["spatialReference"].get(next(crs_iter)))
                self.extent = utils.match_crs(bounds, crs, 4326)
                break
            except InputTypeError:
                continue
            except StopIteration as ex:
                raise ServiceError(self.base_url) from ex

        with contextlib.suppress(KeyError):
            self.units = rjson["units"].replace("esri", "").lower()
            self.max_nrecords = int(rjson["maxRecordCount"])

    def _set_layer_properties(self) -> None:
        """Set properties of the target layer."""
        rjson = self.get_response(self.url, [{"f": "json"}])[0]
        if "fields" not in rjson:
            return
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
        if "alias" in rjson["fields"][0]:
            self.field_types.update(
                {
                    f["alias"]: f["type"].replace("esriFieldType", "").lower()
                    for f in rjson["fields"]
                }
            )

        with contextlib.suppress(KeyError):
            self.feature_types = dict(
                zip((tlz.pluck("id", rjson["types"])), tlz.pluck("name", rjson["types"]))
            )

    def initialize_service(self) -> None:
        """Initialize the RESTFul service."""
        self._set_service_properties()
        if str(self.layer) not in self.valid_layers:
            valids = [f'"{i}" for {n}' for i, n in self.valid_layers.items()]
            raise InputValueError("layer", valids)

        self._set_layer_properties()

        if self.outformat.lower() not in self.query_formats:
            raise InputValueError("outformat", self.query_formats)

        if self.valid_fields and any(f not in self.valid_fields for f in self.outfields):
            raise InputValueError("outfields", self.valid_fields)

    def partition_oids(self, oids: list[int] | int) -> Iterator[tuple[str, ...]]:
        """Partition feature IDs based on ``self.max_nrecords``."""
        oid_list = [oids] if isinstance(oids, int) else set(oids)
        if len(oid_list) == 0:
            raise ZeroMatchedError

        self.n_features = len(oid_list)
        return tlz.partition_all(self.max_nrecords, [str(i) for i in oid_list])

    def _cleanup_resp(
        self, resp: list[dict[str, Any]], payloads: list[dict[str, str]]
    ) -> list[dict[str, Any]]:
        """Remove failed responses."""
        fails = [i for i, r in enumerate(resp) if "error" in r]

        if fails:
            err = resp[fails[0]]["error"]["message"]
            resp = [r for i, r in enumerate(resp) if i not in fails]
            if not resp and self.disable_retry:
                raise ServiceError(err)

            if "objectIds" in payloads[fails[0]]:
                oids = list(tlz.concat(payloads[i]["objectIds"].split(",") for i in fails))
                self.n_missing = len(oids)

                self.failed_path = Path("cache", f"failed_ids_{self.request_id}.txt")
                with self.failed_path.open("w") as f:
                    f.write("\n".join(oids))

                if not self.disable_retry:
                    self.return_m = bool(payloads[0]["ReturnM"])
                    self.return_geom = bool(payloads[0]["returnGeometry"])
                    self.total_n_features = self.n_features

                    if len(fails) > 1:
                        msg = f"Found {len(fails)} failed requests. Retrying ..."
                        warnings.warn(msg, UserWarning, stacklevel=2)
                    elif len(fails) == 1:
                        msg = "Found 1 failed request. Retrying ..."
                        warnings.warn(msg, UserWarning, stacklevel=2)

                    resp.extend(self.retry_failed_requests())

                    if self.n_missing > 0:
                        msg = " ".join(
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
                        warnings.warn(msg, UserWarning, stacklevel=2)

        return resp

    def get_features(
        self,
        featureids: Iterator[tuple[str, ...]],
        return_m: bool = False,
        return_geom: bool = True,
    ) -> list[dict[str, Any]]:
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
                "returnGeometry": str(return_geom).lower(),
                "outSR": str(self.out_sr),
                "outfields": ",".join(self.outfields),
                "ReturnM": str(return_m).lower(),
                "f": self.outformat,
            }
            for ids in featureids
        ]
        if self.request_id is None:
            self.request_id = uuid.uuid4().hex

        resp = self.get_response(self.query_url, payloads, "POST")

        resp = self._cleanup_resp(resp, payloads)
        if not resp:
            raise ZeroMatchedError
        return resp

    def _retry(
        self, return_m: bool, return_geo: bool, partition_fac: float
    ) -> list[dict[str, Any]]:
        """Retry failed requests."""
        with self.failed_path.open() as f:
            oids = [int(i) for i in f.read().splitlines()]

        max_nrecords = self.max_nrecords
        self.max_nrecords = int(max(partition_fac * max_nrecords, 1))
        features = self.get_features(self.partition_oids(oids), return_m, return_geo)
        self.max_nrecords = max_nrecords

        return features

    def retry_failed_requests(self) -> list[dict[str, Any]]:
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
            with contextlib.suppress(ServiceError):
                features.append(self._retry(self.return_m, self.return_geom, f))

        features = list(tlz.concat(features))
        self.n_missing -= len(features)
        self.disable_retry = retry
        os.environ["HYRIVER_CACHE_DISABLE"] = str(caching).lower()
        return features

    def get_response(
        self,
        url: str,
        payloads: list[dict[str, str]],
        method: Literal["get", "GET", "post", "POST"] = "GET",
    ) -> list[dict[str, Any]]:
        """Send payload and get the response."""
        req_key = "params" if method in ("GET", "get") else "data"
        try:
            resp = ar.retrieve_json(
                [url] * len(payloads),
                [{req_key: p} for p in payloads],
                request_method=method,
            )
            resp = cast("list[dict[str, Any]]", resp)
        except ValueError as ex:
            raise ZeroMatchedError from ex
        return resp

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


class WMSBase:
    """Base class for accessing a WMS service.

    Parameters
    ----------
    url : str
        The base url for the WMS service e.g., https://www.mrlc.gov/geoserver/mrlc_download/wms
    layers : str or list, optional
        A layer or a list of layers from the service to be downloaded. You can pass an empty
        string to get a list of available layers.
    outformat : str, optional
        The data format to request for data from the service. You can pass an empty
        string to get a list of available output formats.
    version : str, optional
        The WMS service version which should be either 1.1.1 or 1.3.0, defaults to 1.3.0.
    crs : str, int, or pyproj.CRS, optional
        The spatial reference system to be used for requesting the data, defaults to
        ``epsg:4326``.
    validation : bool, optional
        Validate the input arguments from the WMS service, defaults to True. Set this
        to False if you are sure all the WMS settings such as layer and crs are correct
        to avoid sending extra requests.
    """

    def __init__(
        self,
        url: str,
        layers: str | int | Sequence[str | int] = "",
        outformat: str = "",
        version: str = "1.3.0",
        crs: CRSType = 4326,
        validation: bool = True,
    ) -> None:
        self.url = url
        self.outformat = outformat
        self.version = version
        self.crs = crs
        self.outformat = self.outformat.lower()
        if isinstance(layers, (str, int)):
            self.layers = [str(layers)]
        else:
            self.layers = [str(lyr) for lyr in layers]
        self.crs_str = utils.validate_crs(self.crs).lower()
        if self.version not in ("1.1.1", "1.3.0"):
            raise InputValueError("version", ("1.1.1", "1.3.0"))

        if validation:
            self.get_service_options()
            self.validate_wms()

    def get_service_options(self) -> None:
        """Validate input arguments with the WMS service."""
        try:
            wms = WebMapService(self.url, version=self.version)
        except AttributeError as ex:
            raise ServiceUnavailableError(self.url) from ex

        self.available_layer = {wms[lyr].name: wms[lyr].title for lyr in list(wms.contents)}
        self.available_outformat = wms.getOperationByName("GetMap").formatOptions
        self.available_outformat = [f.lower() for f in self.available_outformat]
        self.available_crs = {
            lyr: [s.lower() for s in wms[lyr].crsOptions] for lyr in self.available_layer
        }

    def validate_wms(self) -> None:
        """Validate input arguments with the WMS service."""
        if any(lyr not in self.available_layer for lyr in self.layers):
            raise InputValueError(
                "layers", (f"{n} for {t}" for n, t in self.available_layer.items())
            )

        if self.outformat not in self.available_outformat:
            raise InputValueError("outformat", self.available_outformat)

        if any(self.crs_str not in self.available_crs[lyr] for lyr in self.layers):
            _valid_crss = (
                f"{lyr}: {', '.join(cs)}\n"
                for lyr, cs in self.available_crs.items()
                if lyr in self.layers
            )
            raise InputValueError("CRS", _valid_crss)

    def get_validlayers(self) -> dict[str, str]:
        """Get the layers supported by the WMS service."""
        try:
            wms = WebMapService(self.url, version=self.version)
        except AttributeError as ex:
            raise ServiceUnavailableError(self.url) from ex

        return {wms[lyr].name: wms[lyr].title for lyr in list(wms.contents)}

    def __repr__(self) -> str:
        """Print the services properties."""
        return "\n".join(
            (
                "Connected to the WMS service with the following properties:",
                f"URL: {self.url}",
                f"Version: {self.version}",
                f"Layers: {', '.join(self.layers)}",
                f"Output Format: {self.outformat}",
                f"Output CRS: {self.crs_str}",
            )
        )


class WFSBase:
    """Base class for WFS service.

    Parameters
    ----------
    url : str
        The base url for the WFS service, for examples:
        https://hazards.fema.gov/arcgis/rest/services/public/NFHL/MapServer/WFSServer
    layer : str
        The layer from the service to be downloaded, defaults to None which throws
        an error and includes all the available layers offered by the service.
    outformat : str
        The data format to request for data from the service, defaults to None which
         throws an error and includes all the available format offered by the service.
    version : str, optional
        The WFS service version which should be either ``1.0.0``, ``1.1.0``, or
        ``2.0.0``. Defaults to ``2.0.0``.
    crs : str, int, or pyproj.CRS, optional
        The spatial reference system to be used for requesting the data, defaults to
        ``epsg:4326``.
    read_method : str, optional
        Method for reading the retrieved data, defaults to ``json``. Valid options are
        ``json``, ``binary``, and ``text``.
    max_nrecords : int, optional
        The maximum number of records in a single request to be retrieved from the service,
        defaults to 1000. If the number of requested records is greater than this value,
        the query will be split into multiple requests.
    validation : bool, optional
        Validate the input arguments from the WFS service, defaults to True. Set this
        to False if you are sure all the WFS settings such as layer and crs are correct
        to avoid sending extra requests.
    """

    def __init__(
        self,
        url: str,
        layer: str | None = None,
        outformat: str | None = None,
        version: str = "2.0.0",
        crs: CRSType = 4326,
        read_method: str = "json",
        max_nrecords: int = 1000,
        validation: bool = True,
    ) -> None:
        self.url = url
        self.layer = layer
        self.outformat = outformat
        self.version = version
        self.crs = crs
        self.read_method = read_method
        self.max_nrecords = max_nrecords
        self.crs_str = utils.validate_crs(self.crs)

        if self.version == "2.0.0":
            self.nfeat_key = "numberMatched="
            self.count_key = "count"
        elif self.version == "1.1.0":
            self.nfeat_key = "numberOfFeatures="
            self.count_key = "maxFeatures"
        else:
            raise InputValueError("version", ("1.1.0", "2.0.0"))

        if self.read_method == "text":
            self.retrieve = ar.retrieve_text
        elif self.read_method == "bytes":
            self.retrieve = ar.retrieve_binary
        elif self.read_method == "json":
            self.retrieve = ar.retrieve_json
        else:
            raise InputValueError("read_method", ("json", "binary", "text"))

        self.get_service_options()
        if validation:
            self.validate_wfs()

    def get_service_options(self) -> None:
        """Validate input arguments with the WFS service."""
        try:
            wfs = WebFeatureService(self.url, version=self.version)
        except AttributeError as ex:
            raise ServiceUnavailableError(self.url) from ex

        if wfs is None:
            raise ServiceUnavailableError(self.url)

        self.available_layer = list(wfs.contents)

        wfs_features = wfs.getOperationByName("GetFeature")
        self.available_outformat = wfs_features.parameters["outputFormat"]["values"]
        self.available_outformat = [f.lower() for f in self.available_outformat]

        self.available_crs = {
            lyr: [f"{s.authority.lower()}:{s.code}" for s in wfs[lyr].crsOptions]
            for lyr in self.available_layer
        }

        if self.layer:
            if self.layer not in self.available_layer:
                raise InputValueError("layers", self.available_layer)
            layers = [self.layer]
        else:
            layers = self.available_layer
        self.schema = dict.fromkeys(layers)
        with contextlib.suppress(KeyError):
            for lyr in layers:
                self.schema[lyr] = wfs.get_schema(lyr)

    def validate_wfs(self) -> None:
        """Validate input arguments with the WFS service."""
        if self.layer is None:
            raise MissingInputError(
                "\n".join(
                    (
                        "The layer argument is missing. The following layers are available:",
                        ", ".join(self.available_layer),
                    )
                )
            )
        if self.layer not in self.available_layer:
            raise InputValueError("layers", self.available_layer)

        if self.outformat not in self.available_outformat:
            raise InputValueError("outformat", self.available_outformat)

        if self.crs_str.lower() not in self.available_crs[self.layer]:
            raise InputValueError("crs", self.available_crs[self.layer])

    def sort_params(
        self, sort_attr: str | None, nfeatures: int, start_index: int
    ) -> dict[str, str]:
        """Get sort parameters for a WFS request."""
        schema = self.schema[self.layer]
        valid_attrs = schema["properties"] if schema else []
        if sort_attr is None:
            if schema is None:
                msg = "sort_attr is None and no id column found in the schema."
                msg += " Please set the sort_attr manually."
                raise MissingInputError(msg)

            sort_attr = next((a for a in valid_attrs if "id" in a.lower()), None)
            if sort_attr is None:
                msg = "sort_attr is None and no id column found in the schema."
                msg += " Please set the sort_attr manually. Available columns are:\n"
                msg += ", ".join(list(valid_attrs))
                raise MissingInputError(msg)

        if valid_attrs and sort_attr not in valid_attrs:
            raise InputValueError("sort_attr", list(valid_attrs))

        if nfeatures <= self.max_nrecords:
            return {
                "startIndex": str(start_index),
                self.count_key: str(nfeatures),
                "sortBy": sort_attr,
            }
        return {
            "startIndex": str(start_index),
            self.count_key: str(self.max_nrecords),
            "sortBy": sort_attr,
        }

    def __repr__(self) -> str:
        """Print the services properties."""
        return "\n".join(
            [
                "Connected to the WFS service with the following properties:",
                f"URL: {self.url}",
                f"Version: {self.version}",
                f"Layer: {self.layer}",
                f"Output Format: {self.outformat}",
                f"Output CRS: {self.crs_str}",
                f"Read Method: {self.read_method}",
            ]
        )
