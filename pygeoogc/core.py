"""Base classes and function for REST, WMS, and WMF services."""
from typing import Dict, List, Optional, Tuple, Union

import cytoolz as tlz
import pyproj
from orjson import JSONDecodeError
from owslib.wfs import WebFeatureService
from owslib.wms import WebMapService
from shapely.geometry import MultiPoint, Point, Polygon

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

        self._layer = 99
        self._outformat = outformat
        self._spatial_relation = spatial_relation
        self._outfields = outfields if isinstance(outfields, list) else [outfields]
        self._n_threads = n_threads
        self.nfeatures = 0
        self.crs = crs
        out_sr = pyproj.CRS(self.crs).to_epsg()
        if out_sr is None:
            raise InvalidInputType("crs", "a valid CRS")
        self.out_sr = out_sr
        self.valid_layers = self.get_validlayers()
        self.extent: Optional[Tuple[float, float, float, float]] = None
        self.feature_types: Optional[Dict[int, str]] = None

        self._zeromatched = "No feature ID were found within the requested region."
        self.test_url()

    @property
    def layer(self) -> int:
        """Set service layer."""
        return self._layer

    @layer.setter
    def layer(self, value: int) -> None:
        """Set service layer."""
        try:
            existing_lyr = int(self.base_url.split("/")[-1])
            self.base_url = self.base_url.replace(f"/{existing_lyr}", "")
        except ValueError:
            pass

        if f"{value}" not in self.valid_layers:
            valids = [f'"{i}" for {n}' for i, n in self.valid_layers.items()]
            raise InvalidInputValue("layer", valids)

        self._layer = value
        try:
            existing_lyr = int(self.base_url.split("/")[-1])
            self.base_url = self.base_url.replace(f"/{existing_lyr}", f"/{value}")
        except ValueError:
            self.base_url = f"{self.base_url}/{value}"

        self.test_url()

    @property
    def outformat(self) -> str:
        """Set service output format."""
        return self._outformat

    @outformat.setter
    def outformat(self, value: str) -> None:
        """Set service output format."""
        if self.base_url is not None and value.lower() not in self.query_formats:
            raise InvalidInputValue("outformat", self.query_formats)

        self._outformat = value

    @property
    def spatial_relation(self) -> str:
        """Set spatial relationship of the input geometry with the source data."""
        return self._spatial_relation

    @spatial_relation.setter
    def spatial_relation(self, value: str) -> None:
        """Set spatial relationship of the input geometry with the source data."""
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
        """Set service output field(s)."""
        return self._outfields

    @outfields.setter
    def outfields(self, value: Union[List[str], str]) -> None:
        """Set service output field(s)."""
        if not isinstance(value, (list, str)):
            raise InvalidInputType("outfields", "str or list")

        self._outfields = value if isinstance(value, list) else [value]

    @property
    def n_threads(self) -> int:
        """Set number of threads."""
        return self._n_threads

    @n_threads.setter
    def n_threads(self, value: int) -> None:
        """Set number of threads."""
        if not isinstance(value, int) or value < 0:
            raise InvalidInputType("n_threads", "positive int")
        self._n_threads = value

    @property
    def max_nrecords(self) -> int:
        """Set maximum number of features per request."""
        return self._max_nrecords

    @max_nrecords.setter
    def max_nrecords(self, value: int) -> None:
        """Set maximum number of features per request."""
        if value > self.max_nrecords:
            raise ValueError(
                f"The server doesn't accept more than {self.max_nrecords}" + " records per request."
            )
        if value < 0:
            raise InvalidInputType("max_nrecords", "positive int")

        self._max_nrecords = value

    @property
    def featureids(self) -> List[Tuple[str, ...]]:
        """Set feature ID(s)."""
        return self._featureids

    @featureids.setter
    def featureids(self, value: Union[List[int], int, None]) -> None:
        """Set feature ID(s)."""
        if value is None:
            raise ZeroMatched(self._zeromatched)

        if not isinstance(value, (list, int)):
            raise InvalidInputType("featureids", "int or list")

        oids = [str(value)] if isinstance(value, (int, str)) else [str(v) for v in value]

        self.nfeatures = len(oids)
        if self.nfeatures == 0:
            raise ZeroMatched(self._zeromatched)

        self._featureids = list(tlz.partition_all(self.max_nrecords, oids))

    def get_validlayers(self) -> Dict[str, str]:
        """Get all the valid service layer."""
        try:
            existing_lyr = int(self.base_url.split("/")[-1])
            url = self.base_url.replace(f"/{existing_lyr}", "")
        except ValueError:
            url = self.base_url

        resp = self.session.get(url, {"f": "json"})
        utils.check_response(resp)

        layers = {"No layer": ""}
        try:
            layers = {f"{lyr['id']}": lyr["name"] for lyr in resp.json()["layers"]}
        except (JSONDecodeError, KeyError):
            raise ZeroMatched(f"The requested layer does not exists in:\n{self.base_url}")

        return layers

    def get_validfields(self) -> Dict[str, str]:
        """Get all the valid service output fields."""
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
            self._max_nrecords = int(resp["maxRecordCount"])
            self.query_formats = resp["supportedQueryFormats"].replace(" ", "").lower().split(",")
            self.valid_fields = list(
                set(
                    utils.traverse_json(resp, ["fields", "name"])
                    + utils.traverse_json(resp, ["fields", "alias"])
                    + ["*"]
                )
            )
            try:
                extent = resp["extent"] if "extent" in resp else resp["fullExtent"]
                bounds = (extent["xmin"], extent["ymin"], extent["xmax"], extent["ymax"])
                crs = extent["spatialReference"]["latestWkid"]
                self.extent = utils.MatchCRS.bounds(bounds, crs, DEF_CRS)
            except KeyError:
                self.extent = None
            try:
                self.feature_types = dict(
                    zip((tlz.pluck("id", resp["types"])), tlz.pluck("name", resp["types"]))
                )
            except KeyError:
                self.feature_types = None
        except KeyError:
            raise ServerError(self.base_url)

    def _esri_query(
        self,
        geom: Union[
            Polygon,
            Point,
            MultiPoint,
            Tuple[float, float, float, float],
        ],
        geo_crs: str = DEF_CRS,
    ) -> Dict[str, Union[str, bytes]]:
        """Generate geometry queries based on ESRI template."""
        if isinstance(geom, tuple) and len(geom) == 4:
            geom = utils.MatchCRS.bounds(geom, geo_crs, self.crs)  # type: ignore
            return utils.ESRIGeomQuery(geom, self.out_sr).bbox()

        if isinstance(geom, Point):
            geom = utils.MatchCRS.geometry(geom, geo_crs, self.crs)
            return utils.ESRIGeomQuery((geom.x, geom.y), self.out_sr).point()

        if isinstance(geom, MultiPoint):
            geom = utils.MatchCRS.geometry(geom, geo_crs, self.crs)
            return utils.ESRIGeomQuery([(g.x, g.y) for g in geom], self.out_sr).multipoint()

        if isinstance(geom, Polygon):
            geom = utils.MatchCRS.geometry(geom, geo_crs, self.crs)
            return utils.ESRIGeomQuery(geom, self.out_sr).polygon()

        raise InvalidInputType("geom", "Polygon, Point, MultiPoint, tuple, or list of tuples")

    def __repr__(self) -> str:
        """Print the service configuration."""
        extent = ", ".join(f"{c:.3f}" for c in self.extent) if self.extent else None
        ftypes = ", ".join(self.feature_types.values()) if self.feature_types else None
        return "\n".join(
            [
                "Service configurations:",
                f"    URL: {self.base_url}",
                f"    Max Record Count: {self.max_nrecords}",
                f"    Supported Query Formats: {self.query_formats}",
                f"    Units: {self.units}",
                f"    Extent: ({extent})",
                f"    Feature Types: {ftypes}",
            ]
        )


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

    def __init__(
        self,
        url: str,
        layers: Union[str, List[str]],
        outformat: str,
        version: str = "1.3.0",
        crs: str = DEF_CRS,
    ) -> None:
        self.url = url
        self.layers = layers
        self.outformat = outformat
        self.version = version
        self.crs = crs

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

    def __init__(
        self,
        url: str,
        layer: Optional[str] = None,
        outformat: Optional[str] = None,
        version: str = "2.0.0",
        crs: str = DEF_CRS,
    ) -> None:
        self.url = url
        self.layer = layer
        self.outformat = outformat
        self.version = version
        self.crs = crs

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

    def validate_wfs(self) -> None:
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
