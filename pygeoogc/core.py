"""Base classes and function for REST, WMS, and WMF services."""
from dataclasses import dataclass
from itertools import zip_longest
from typing import Dict, List, Optional, Tuple, Union

import pyproj
from orjson import JSONDecodeError
from owslib.wfs import WebFeatureService
from owslib.wms import WebMapService

from . import utils
from .exceptions import InvalidInputType, InvalidInputValue, MissingInputs, ServerError, ZeroMatched
from .utils import RetrySession

DEF_CRS = "epsg:4326"


class ArcGISRESTfulBase:
    """Access to an ArcGIS REST service.

    Parameters
    ----------
    base_url : str, optional
        The ArcGIS RESTful service url.
    outformat : str, optional
        One of the output formats offered by the selected layer. If not correct
        a list of available formats is shown, defaults to ``geojson``.
    spatial_relation : str, optional
        The spatial relationship to be applied on the input geometry
        while performing the query. If not correct a list of available options is shown.
        It defaults to ``esriSpatialRelIntersects``.
    outfields : str or list
        The output fields to be requested. Setting ``*`` as outfields requests
        all the available fields which is the default behaviour.
    crs : str, optional
        The spatial reference of the output data, defaults to EPSG:4326
    n_threads : int, optional
        Number of simultaneous download, default to 1 i.e., no threading. Note
        that some services might face issues when several requests are sent
        simultaniously and will return the requests partially. It's recommended
        to avoid performing threading unless you are certain the web service can handle it.
    """

    def __init__(
        self,
        base_url: str,
        outformat: str = "geojson",
        outfields: Union[List[str], str] = "*",
        spatial_relation: str = "esriSpatialRelIntersects",
        crs: str = DEF_CRS,
        n_threads: int = 1,
    ) -> None:

        self.session = RetrySession()
        self.base_url = base_url[:-1] if base_url[-1] == "/" else base_url
        self.test_url()

        self._layer = ""
        self._outformat = outformat
        self._spatial_relation = spatial_relation
        self._outfields = outfields if isinstance(outfields, list) else [outfields]
        self._n_threads = n_threads
        self.nfeatures = 0
        self.crs = crs
        self.out_sr = pyproj.CRS(self.crs).to_epsg()
        self.valid_layers = self.get_validlayers()

        self._zeromatched = "No feature ID were found within the requested region."

    @property
    def layer(self) -> str:
        return self._layer

    @layer.setter
    def layer(self, value: int) -> None:
        try:
            existing_lyr = int(self.base_url.split("/")[-1])
            self.base_url = self.base_url.replace(f"/{existing_lyr}", "")
        except ValueError:
            pass

        self.valid_layers = self.get_validlayers()
        if value not in self.valid_layers:
            valids = [f'"{i}" for {n}' for i, n in self.valid_layers.items()]
            raise InvalidInputValue("layer", valids)

        self._layer = f"{value}"
        try:
            existing_lyr = int(self.base_url.split("/")[-1])
            self.base_url = self.base_url.replace(f"/{existing_lyr}", f"/{value}")
        except ValueError:
            self.base_url = f"{self.base_url}/{value}"

        self.test_url()

    @property
    def outformat(self) -> str:
        return self._outformat

    @outformat.setter
    def outformat(self, value: str) -> None:
        if self.base_url is not None and value.lower() not in self.query_formats:
            raise InvalidInputValue("outformat", self.query_formats)

        self._outformat = value

    @property
    def spatial_relation(self) -> str:
        return self._spatial_relation

    @spatial_relation.setter
    def spatial_relation(self, value: str) -> None:
        valid_spatialrels = [
            "esriSpatialRelIntersects",
            "esriSpatialRelContains",
            "esriSpatialRelCrosses",
            "esriSpatialRelEnvelopeIntersects",
            "esriSpatialRelIndexIntersects",
            "esriSpatialRelOverlaps",
            "esriSpatialRelTouches",
            "esriSpatialRelWithin",
            "esriSpatialRelRelation",
        ]
        if value not in valid_spatialrels:
            raise InvalidInputValue("spatial_rel", valid_spatialrels)
        self._spatial_relation = value

    @property
    def outfields(self) -> List[str]:
        return self._outfields

    @outfields.setter
    def outfields(self, value: Union[List[str], str]) -> None:
        if not isinstance(value, (list, str)):
            raise InvalidInputType("outfields", "str or list")

        self._outfields = value if isinstance(value, list) else [value]

    def get_validlayers(self) -> Dict[str, str]:
        try:
            existing_lyr = int(self.base_url.split("/")[-1])
            url = self.base_url.replace(f"/{existing_lyr}", "")
        except ValueError:
            url = self.base_url

        resp = self.session.get(url, {"f": "json"})
        utils.check_response(resp)

        layers = {"No layer": ""}
        try:
            layers = {lyr["id"]: lyr["name"] for lyr in resp.json()["layers"]}
        except (JSONDecodeError, KeyError):
            raise ZeroMatched(f"The requested layer does not exists in:\n{self.base_url}")

        return layers

    def get_validfields(self) -> Dict[str, str]:
        resp = self.session.get(self.base_url, {"f": "json"}).json()
        return {f["name"]: f["type"].replace("esriFieldType", "").lower() for f in resp["fields"]}

    def test_url(self) -> None:
        """Test the generated url and get the required parameters from the service."""
        try:
            resp = self.session.get(self.base_url, {"f": "json"}).json()
            try:
                self.units = resp["units"].replace("esri", "").lower()
            except KeyError:
                self.units = None
            self.maxrec_ount = resp["maxRecordCount"]
            self.query_formats = resp["supportedQueryFormats"].replace(" ", "").lower().split(",")
            self.valid_fields = list(
                set(
                    utils.traverse_json(resp, ["fields", "name"])
                    + utils.traverse_json(resp, ["fields", "alias"])
                    + ["*"]
                )
            )
        except KeyError:
            raise ServerError(self.base_url)

        self._max_nrecords = self.maxrec_ount

    @property
    def n_threads(self) -> int:
        return self._n_threads

    @n_threads.setter
    def n_threads(self, value: int) -> None:
        if not isinstance(value, int) or value < 0:
            raise InvalidInputType("n_threads", "positive int")
        self._n_threads = value

    @property
    def max_nrecords(self) -> int:
        return self._max_nrecords

    @max_nrecords.setter
    def max_nrecords(self, value: int) -> None:
        if value > self.maxrec_ount:
            raise ValueError(
                f"The server doesn't accept more than {self.maxrec_ount}" + " records per request."
            )
        if value > 0:
            self._max_nrecords = value
        else:
            raise InvalidInputType("max_nrecords", "positive int")

    @property
    def featureids(self) -> List[Tuple[str, ...]]:
        return self._featureids

    @featureids.setter
    def featureids(self, value: Union[List[int], int]) -> None:
        if not isinstance(value, (list, int)):
            raise InvalidInputType("featureids", "int or list")

        oids = [str(value)] if isinstance(value, int) else [str(v) for v in value]

        self.nfeatures = len(oids)
        if self.nfeatures == 0:
            raise ZeroMatched(self._zeromatched)

        oid_list = list(zip_longest(*[iter(oids)] * self.max_nrecords))
        oid_list[-1] = tuple(i for i in oid_list[-1] if i is not None)
        self._featureids = oid_list

    def __repr__(self) -> str:
        """Print the service configuration."""
        return (
            "Service configurations:\n"
            + f"URL: {self.base_url}\n"
            + f"Max Record Count: {self.maxrec_ount}\n"
            + f"Supported Query Formats: {self.query_formats}\n"
            + f"Units: {self.units}"
        )


@dataclass
class WMSBase:
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
        epsg:4326.
    """

    url: str
    layers: Union[str, List[str]]
    outformat: str
    version: str = "1.3.0"
    crs: str = DEF_CRS

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

    def validate_wms(self) -> None:
        """Validate input arguments with the WMS service."""
        wms = WebMapService(self.url, version=self.version)

        if not isinstance(self.layers, (str, list)):
            raise InvalidInputType("layers", "str or list")

        layers = [self.layers] if isinstance(self.layers, str) else self.layers
        valid_layers = {wms[lyr].name: wms[lyr].title for lyr in list(wms.contents)}
        if any(lyr not in valid_layers.keys() for lyr in layers):
            raise InvalidInputValue("layers", (f"{n} for {t}" for n, t in valid_layers.items()))

        valid_outformats = wms.getOperationByName("GetMap").formatOptions
        if self.outformat not in valid_outformats:
            raise InvalidInputValue("outformat", valid_outformats)

        valid_crss = {lyr: [s.lower() for s in wms[lyr].crsOptions] for lyr in layers}
        if any(self.crs not in valid_crss[lyr] for lyr in layers):
            _valid_crss = (f"{lyr}: {', '.join(cs)}\n" for lyr, cs in valid_crss.items())
            raise InvalidInputValue("CRS", _valid_crss)

    def get_validlayers(self) -> Dict[str, str]:
        """Get the layers supportted by the WMS service."""
        wms = WebMapService(self.url, version=self.version)

        return {wms[lyr].name: wms[lyr].title for lyr in list(wms.contents)}


@dataclass
class WFSBase:
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
        The WFS service version which should be either 1.1.1, 1.3.0, or 2.0.0.
        Defaults to 2.0.0.
    crs: str, optional
        The spatial reference system to be used for requesting the data, defaults to
        epsg:4326.
    """

    url: str
    layer: Optional[str] = None
    outformat: Optional[str] = None
    version: str = "2.0.0"
    crs: str = DEF_CRS

    def __repr__(self) -> str:
        """Print the services properties."""
        return (
            "Connected to the WFS service with the following properties:\n"
            + f"URL: {self.url}\n"
            + f"Version: {self.version}\n"
            + f"Layer: {self.layer}\n"
            + f"Output Format: {self.outformat}\n"
            + f"Output CRS: {self.crs}"
        )

    def validate_wfs(self):
        """Validate input arguments with the WFS service."""
        wfs = WebFeatureService(self.url, version=self.version)

        valid_layers = list(wfs.contents)
        if self.layer is None:
            raise MissingInputs(
                "The layer argument is missing."
                + " The following layers are available:\n"
                + ", ".join(valid_layers)
            )

        if self.layer not in valid_layers:
            raise InvalidInputValue("layers", valid_layers)

        valid_outformats = wfs.getOperationByName("GetFeature").parameters["outputFormat"]["values"]
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
            "outputFormat": self.outformat,
            "request": "GetFeature",
            "typeName": self.layer,
            max_features: 1,
        }

        resp = RetrySession().get(self.url, payload)
        utils.check_response(resp)

        r_json = resp.json()
        valid_fields = list(
            set(
                utils.traverse_json(r_json, ["fields", "name"])
                + utils.traverse_json(r_json, ["fields", "alias"])
                + ["*"]
            )
        )

        if None in valid_fields:
            valid_fields = list(utils.traverse_json(r_json, ["features", "properties"])[0].keys())

        return valid_fields
