"""Base classes and function for REST, WMS, and WMF services."""
import itertools
import logging
import sys
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Dict, List, Optional, Tuple, Union

import async_retriever as ar
import cytoolz as tlz
import pyproj
import shapely.ops as ops
import yaml
from shapely.geometry import LineString, MultiPoint, MultiPolygon, Point, Polygon

from . import utils
from .core import ArcGISRESTfulBase, WFSBase, WMSBase
from .exceptions import InvalidInputType, InvalidInputValue, ZeroMatched

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
handler = logging.StreamHandler(sys.stdout)
handler.setFormatter(logging.Formatter(""))
logger.handlers = [handler]
logger.propagate = False

DEF_CRS = "epsg:4326"


class ArcGISRESTful(ArcGISRESTfulBase):
    """Access to an ArcGIS REST service.

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
        all the available fields which is the default behaviour.
    crs : str, optional
        The spatial reference of the output data, defaults to EPSG:4326
    max_workers : int, optional
        Number of simultaneous download, default to 1, i.e., no threading. Note
        that some services might face issues when several requests are sent
        simultaneously and will return the requests partially. It's recommended
        to avoid performing threading unless you are certain the web service can handle it.
    """

    def oids_bygeom(
        self,
        geom: Union[
            LineString,
            Polygon,
            Point,
            MultiPoint,
            Tuple[float, float],
            List[Tuple[float, float]],
            Tuple[float, float, float, float],
        ],
        geo_crs: str = DEF_CRS,
        spatial_relation: str = "esriSpatialRelIntersects",
        sql_clause: Optional[str] = None,
        distance: Optional[int] = None,
    ) -> None:
        """Get feature IDs within a geometry that can be combined with a SQL where clause.

        Parameters
        ----------
        geom : LineString, Polygon, Point, MultiPoint, tuple, or list of tuples
            A geometry (LineString, Polygon, Point, MultiPoint), tuple of length two
            (``(x, y)``), a list of tuples of length 2 (``[(x, y), ...]``), or bounding box
            (tuple of length 4 (``(xmin, ymin, xmax, ymax)``)).
        geo_crs : str
            The spatial reference of the input geometry, defaults to EPSG:4326.
        spatial_relation : str, optional
            The spatial relationship to be applied on the input geometry
            while performing the query. If not correct a list of available options is shown.
            It defaults to ``esriSpatialRelIntersects``. Valid predicates are:

            * ``esriSpatialRelIntersects``
            * ``esriSpatialRelContains``
            * ``esriSpatialRelCrosses``
            * ``esriSpatialRelEnvelopeIntersects``
            * ``esriSpatialRelIndexIntersects``
            * ``esriSpatialRelOverlaps``
            * ``esriSpatialRelTouches``
            * ``esriSpatialRelWithin``
            * ``esriSpatialRelRelation``

        sql_clause : str, optional
            Valid SQL 92 WHERE clause, default to None.
        distance : int, optional
            Buffer distance in meters for the input geometries, default to None.
        """
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
        if spatial_relation not in valid_spatialrels:
            raise InvalidInputValue("spatial_relation", valid_spatialrels)

        if isinstance(geom, tuple) and len(geom) == 2:
            geom = Point(geom)
        elif isinstance(geom, list) and all(len(g) == 2 for g in geom):
            geom = MultiPoint(geom)

        geom_query = self._esri_query(geom, geo_crs)

        payload = {
            **geom_query,  # type: ignore
            "spatialRel": spatial_relation,
            "returnGeometry": "false",
            "returnIdsOnly": "true",
            "f": self.outformat,
        }
        if distance:
            payload.update({"distance": f"{distance}", "units": "esriSRUnit_Meter"})

        if sql_clause:
            payload.update({"where": sql_clause})

        resp = self._get_response([payload], method="POST")[0]
        try:
            self.featureids = self.partition_oids(resp["objectIds"])
        except KeyError as ex:
            raise ZeroMatched(resp["error"]["message"]) from ex

    def oids_byfield(self, field: str, ids: Union[str, List[str]]) -> None:
        """Get Object IDs based on a list of field IDs.

        Parameters
        ----------
        field : str
            Name of the target field that IDs belong to.
        ids : str or list
            A list of target ID(s).
        """
        if field not in self.valid_fields:
            raise InvalidInputValue("field", self.valid_fields)

        ftype = self.field_types[field]
        if "string" in ftype:
            fids = ", ".join(f"'{i}'" for i in ids)
        else:
            fids = ", ".join(f"{i}" for i in ids)

        self.oids_bysql(f"{field} IN ({fids})")

    def oids_bysql(self, sql_clause: str) -> None:
        """Get feature IDs using a valid SQL 92 WHERE clause.

        Notes
        -----
        Not all web services support this type of query. For more details look
        `here <https://developers.arcgis.com/rest/services-reference/query-feature-service-.htm#ESRI_SECTION2_07DD2C5127674F6A814CE6C07D39AD46>`__.

        Parameters
        ----------
        sql_clause : str
            A valid SQL 92 WHERE clause.
        """
        if not isinstance(sql_clause, str):
            raise InvalidInputType("sql_clause", "str")

        payload = {
            "where": sql_clause,
            "returnGeometry": "false",
            "returnIdsOnly": "true",
            "f": self.outformat,
        }

        resp = self._get_response([payload])[0]
        try:
            self.featureids = self.partition_oids(resp["objectIds"])
        except KeyError as ex:
            raise ZeroMatched(resp["error"]["message"]) from ex


class WMS(WMSBase):
    """Get data from a WMS service within a geometry or bounding box.

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
    crs : str, optional
        The spatial reference system to be used for requesting the data, defaults to
        epsg:4326.
    version : str, optional
        The WMS service version which should be either 1.1.1 or 1.3.0, defaults to 1.3.0.
    validation : bool, optional
        Validate the input arguments from the WMS service, defaults to True. Set this
        to False if you are sure all the WMS settings such as layer and crs are correct
        to avoid sending extra requests.
    """

    def __init__(
        self,
        url: str,
        layers: Union[str, List[str]],
        outformat: str,
        version: str = "1.3.0",
        crs: str = DEF_CRS,
        validation: bool = True,
    ) -> None:
        super().__init__(url, layers, outformat, version, crs)

        self.layers = [self.layers] if isinstance(self.layers, str) else self.layers
        if validation:
            self.validate_wms()

    def getmap_bybox(
        self,
        bbox: Tuple[float, float, float, float],
        resolution: float,
        box_crs: str = DEF_CRS,
        always_xy: bool = False,
        max_px: int = 8000000,
        kwargs: Optional[Dict[str, Any]] = None,
    ) -> Dict[str, bytes]:
        """Get data from a WMS service within a geometry or bounding box.

        Parameters
        ----------
        bbox : tuple
            A bounding box for getting the data.
        resolution : float
            The output resolution in meters. The width and height of output are computed in pixel
            based on the geometry bounds and the given resolution.
        box_crs : str, optional
            The spatial reference system of the input bbox, defaults to
            epsg:4326.
        always_xy : bool, optional
            Whether to always use xy axis order, defaults to False. Some services change the axis
            order from xy to yx, following the latest WFS version specifications but some don't.
            If the returned value does not have any geometry, it indicates that most probably the
            axis order does not match. You can set this to True in that case.
        max_px : int, opitonal
            The maximum allowable number of pixels (width x height) for a WMS requests,
            defaults to 8 million based on some trial-and-error.
        kwargs: dict, optional
            Optional additional keywords passed as payload, defaults to None.
            For example, ``{"styles": "default"}``.

        Returns
        -------
        dict
            A dict where the keys are the layer name and values are the returned response
            from the WMS service as bytes.
        """
        utils.check_bbox(bbox)
        _bbox = utils.match_crs(bbox, box_crs, self.crs)
        bounds = utils.bbox_decompose(_bbox, resolution, self.crs, max_px)

        payload = {
            "version": self.version,
            "format": self.outformat,
            "request": "GetMap",
        }

        if not isinstance(kwargs, (dict, type(None))):
            raise InvalidInputType("kwargs", "dict or None")

        if isinstance(kwargs, dict):
            payload.update(kwargs)

        if self.version == "1.1.1":
            payload["srs"] = self.crs
        else:
            payload["crs"] = self.crs

        geographic_crs = pyproj.CRS.from_user_input(self.crs).is_geographic

        def _get_payloads(
            args: Tuple[str, Tuple[Tuple[float, float, float, float], str, int, int]]
        ) -> Tuple[str, Dict[str, str]]:
            lyr, bnds = args
            _bbox, counter, _width, _height = bnds

            if self.version != "1.1.1" and geographic_crs and not always_xy:
                _bbox = (_bbox[1], _bbox[0], _bbox[3], _bbox[2])
            _payload = payload.copy()
            _payload["bbox"] = f'{",".join(str(c) for c in _bbox)}'
            _payload["width"] = str(_width)
            _payload["height"] = str(_height)
            _payload["layers"] = lyr
            return f"{lyr}_dd_{counter}", _payload

        layers, payloads = zip(*(_get_payloads(i) for i in itertools.product(self.layers, bounds)))
        rbinary = ar.retrieve(
            [self.url] * len(payloads),
            "binary",
            [{"params": p} for p in payloads],
            max_workers=4,
        )
        return dict(zip(layers, rbinary))


class WFS(WFSBase):
    """Data from any WFS service within a geometry or by featureid.

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
        The WFS service version which should be either 1.0.0, 1.1.0, or 2.0.0.
        Defaults to 2.0.0.
    crs: str, optional
        The spatial reference system to be used for requesting the data, defaults to
        epsg:4326.
    read_method : str, optional
        Method for reading the retrieved data, defaults to ``json``. Valid options are
        ``json``, ``binary``, and ``text``.
    max_nrecords : int, optional
        The maximum number of records in a single request to be retrieved from the service,
        defaults to 1000. If the number of records requested is greater than this value,
        it will be split into multiple requests.
    validation : bool, optional
        Validate the input arguments from the WFS service, defaults to True. Set this
        to False if you are sure all the WFS settings such as layer and crs are correct
        to avoid sending extra requests.
    """

    def __init__(
        self,
        url: str,
        layer: Optional[str] = None,
        outformat: Optional[str] = None,
        version: str = "2.0.0",
        crs: str = DEF_CRS,
        read_method: str = "json",
        max_nrecords: int = 1000,
        validation: bool = True,
    ) -> None:
        super().__init__(url, layer, outformat, version, crs, read_method, max_nrecords)

        if validation:
            self.validate_wfs()

    def getfeature_bybox(
        self,
        bbox: Tuple[float, float, float, float],
        box_crs: str = DEF_CRS,
        always_xy: bool = False,
    ) -> Union[str, bytes, Dict[str, Any]]:
        """Get data from a WFS service within a bounding box.

        Parameters
        ----------
        bbox : tuple
            A bounding box for getting the data: [west, south, east, north]
        box_crs : str, optional
            The spatial reference system of the input bbox, defaults to
            epsg:4326.
        always_xy : bool, optional
            Whether to always use xy axis order, defaults to False. Some services change the axis
            order from xy to yx, following the latest WFS version specifications but some don't.
            If the returned value does not have any geometry, it indicates that most probably the
            axis order does not match. You can set this to True in that case.

        Returns
        -------
        str or bytes or dict
            WFS query response within a bounding box.
        """
        utils.check_bbox(bbox)

        if (
            self.version != "1.0.0"
            and pyproj.CRS.from_user_input(box_crs).is_geographic
            and not always_xy
        ):
            bbox = (bbox[1], bbox[0], bbox[3], bbox[2])

        payload = {
            "service": "wfs",
            "version": self.version,
            "outputFormat": self.outformat,
            "request": "GetFeature",
            "typeName": self.layer,
            "bbox": f'{",".join(str(c) for c in bbox)},{box_crs}',
            "srsName": self.crs,
        }

        return ar.retrieve([self.url], self.read_method, [{"params": payload}])[0]

    def getfeature_bygeom(
        self,
        geometry: Union[Polygon, MultiPolygon],
        geo_crs: str = DEF_CRS,
        always_xy: bool = False,
        predicate: str = "INTERSECTS",
    ) -> Union[str, bytes, Dict[str, Any]]:
        """Get features based on a geometry.

        Parameters
        ----------
        geometry : shapely.geometry
            The input geometry
        geo_crs : str, optional
            The CRS of the input geometry, default to epsg:4326.
        always_xy : bool, optional
            Whether to always use xy axis order, defaults to False. Some services change the axis
            order from xy to yx, following the latest WFS version specifications but some don't.
            If the returned value does not have any geometry, it indicates that most probably the
            axis order does not match. You can set this to True in that case.
        predicate : str, optional
            The geometric predicate to use for requesting the data, defaults to ``INTERSECTS``.
            Valid predicates are:

            * ``EQUALS``
            * ``DISJOINT``
            * ``INTERSECTS``
            * ``TOUCHES``
            * ``CROSSES``
            * ``WITHIN``
            * ``CONTAINS``
            * ``OVERLAPS``
            * ``RELATE``
            * ``BEYOND``

        Returns
        -------
        str or bytes or dict
            WFS query response based on the given geometry.
        """
        geom = utils.match_crs(geometry, geo_crs, self.crs)

        if (
            self.version != "1.0.0"
            and pyproj.CRS.from_user_input(geo_crs).is_geographic
            and not always_xy
        ):
            g_wkt = ops.transform(lambda x, y: (y, x), geom).wkt
        else:
            g_wkt = geom.wkt

        valid_predicates = [
            "EQUALS",
            "DISJOINT",
            "INTERSECTS",
            "TOUCHES",
            "CROSSES",
            "WITHIN",
            "CONTAINS",
            "OVERLAPS",
            "RELATE",
            "BEYOND",
        ]
        if predicate not in valid_predicates:
            raise InvalidInputValue("predicate", valid_predicates)

        return self.getfeature_byfilter(f"{predicate.upper()}(the_geom, {g_wkt})", method="POST")

    def getfeature_byid(
        self,
        featurename: str,
        featureids: Union[List[str], str],
    ) -> List[Union[str, bytes, Dict[str, Any]]]:
        """Get features based on feature IDs.

        Parameters
        ----------
        featurename : str
            The name of the column for searching for feature IDs.
        featureids : str or list
            The feature ID(s).

        Returns
        -------
        str or bytes or dict
            WMS query response.
        """
        valid_features = self.get_validnames()
        if featurename not in valid_features:
            raise InvalidInputValue("featurename", valid_features)

        if not isinstance(featureids, (str, int, list)):
            raise InvalidInputType("featureids", "str or list of str")

        featureids = [featureids] if isinstance(featureids, (str, int)) else featureids

        if len(featureids) == 0:
            raise InvalidInputType("featureids", "int or str or list")

        fid_list = (
            ", ".join(f"'{fid}'" for fid in fids)
            for fids in tlz.partition_all(self.max_nrecords, featureids)
        )

        return [
            self.getfeature_byfilter(f"{featurename} IN ({fids})", method="POST")
            for fids in fid_list
        ]

    def getfeature_byfilter(
        self, cql_filter: str, method: str = "GET"
    ) -> Union[str, bytes, Dict[str, Any]]:
        """Get features based on a valid CQL filter.

        Notes
        -----
        The validity of the input CQL expression is user's responsibility since
        the function does not perform any checks and just sends a request using
        the input filter.

        Parameters
        ----------
        cql_filter : str
            A valid CQL filter expression.
        method : str
            The request method, could be GET or POST (for long filters).

        Returns
        -------
        str or bytes or dict
            WFS query response
        """
        if not isinstance(cql_filter, str):
            raise InvalidInputType("cql_filter", "str")

        valid_methods = ["GET", "POST"]
        if method not in valid_methods:
            raise InvalidInputValue("method", valid_methods)

        payload = {
            "service": "wfs",
            "version": self.version,
            "outputFormat": self.outformat,
            "request": "GetFeature",
            "typeName": self.layer,
            "srsName": self.crs,
            "cql_filter": cql_filter,
        }

        if method == "GET":
            return ar.retrieve([self.url], self.read_method, [{"params": payload}])[0]

        headers = {"content-type": "application/x-www-form-urlencoded"}
        return ar.retrieve(
            [self.url], self.read_method, [{"data": payload, "headers": headers}], "POST"
        )[0]


class ServiceURL:
    """Base URLs of the supported services."""

    def __init__(self) -> None:
        fpath = Path(__file__).parent.joinpath("static/urls.yml")
        with open(fpath) as fp:
            self.urls = yaml.safe_load(fp)

    def _make_nt(self, service: str) -> SimpleNamespace:
        return SimpleNamespace(**self.urls[service])

    @property
    def restful(self) -> SimpleNamespace:
        """Read RESTful URLs from the source yml file."""
        return self._make_nt("restful")

    @property
    def wms(self) -> SimpleNamespace:
        """Read WMS URLs from the source yml file."""
        return self._make_nt("wms")

    @property
    def wfs(self) -> SimpleNamespace:
        """Read WFS URLs from the source yml file."""
        return self._make_nt("wfs")

    @property
    def http(self) -> SimpleNamespace:
        """Read HTTP URLs from the source yml file."""
        return self._make_nt("http")
