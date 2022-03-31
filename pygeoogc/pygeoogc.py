"""Base classes and function for REST, WMS, and WMF services."""
import itertools
import uuid
from ssl import SSLContext
from typing import Any, Dict, Iterator, List, NamedTuple, Optional, Sequence, Tuple, Union

import async_retriever as ar
import cytoolz as tlz
import pyproj
import shapely.ops as ops
from shapely.geometry import LineString, MultiPoint, MultiPolygon, Point, Polygon

from . import utils
from .core import DEF_CRS, ArcGISRESTfulBase, WFSBase, WMSBase
from .exceptions import InvalidInputType, InvalidInputValue, ZeroMatched


class ArcGISRESTful:
    """Access to an ArcGIS REST service.

    Notes
    -----
    By default, all retrieval methods retry to get the missing feature IDs,
    if there are any. You can disable this behavior by setting ``disable_retry``
    to ``True``. If there are any missing feature IDs after the retry,
    they are saved to a text file, path of which can be accessed by
    ``self.client.failed_path``.

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
        The spatial reference of the output data, defaults to ``epsg:4326``.
    max_workers : int, optional
        Number of simultaneous download, default to 1, i.e., no threading. Note
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
        which its path can be accessed via ``self.client.failed_path``.
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
        self.client = ArcGISRESTfulBase(
            base_url=base_url,
            layer=layer,
            outformat=outformat,
            outfields=outfields,
            crs=crs,
            max_workers=max_workers,
            verbose=verbose,
            disable_retry=disable_retry,
        )

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
        geo_crs: Union[str, pyproj.CRS] = DEF_CRS,
        spatial_relation: str = "esriSpatialRelIntersects",
        sql_clause: Optional[str] = None,
        distance: Optional[int] = None,
    ) -> Iterator[Tuple[str, ...]]:
        """Get feature IDs within a geometry that can be combined with a SQL where clause.

        Parameters
        ----------
        geom : LineString, Polygon, Point, MultiPoint, tuple, or list of tuples
            A geometry (LineString, Polygon, Point, MultiPoint), tuple of length two
            (``(x, y)``), a list of tuples of length 2 (``[(x, y), ...]``), or bounding box
            (tuple of length 4 (``(xmin, ymin, xmax, ymax)``)).
        geo_crs : str or pyproj.CRS
            The spatial reference of the input geometry.
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

        Returns
        -------
        list of tuples
            A list of feature IDs partitioned by ``self.max_nrecords``.
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

        geom_query = self.client.esri_query(geom, geo_crs)

        payload = {
            **geom_query,
            "spatialRel": spatial_relation,
            "returnGeometry": "false",
            "returnIdsOnly": "true",
            "f": self.client.outformat,
        }
        if distance:
            payload.update({"distance": f"{distance}", "units": "esriSRUnit_Meter"})

        if sql_clause:
            payload.update({"where": sql_clause})

        self.client.request_id = uuid.uuid4().hex

        resp = self.client.get_response(self.client.query_url, [payload], method="POST")[0]
        try:
            return self.partition_oids(resp["objectIds"])
        except KeyError as ex:
            raise ZeroMatched(resp["error"]["message"]) from ex

    def oids_byfield(self, field: str, ids: Union[str, List[str]]) -> Iterator[Tuple[str, ...]]:
        """Get Object IDs based on a list of field IDs.

        Parameters
        ----------
        field : str
            Name of the target field that IDs belong to.
        ids : str or list
            A list of target ID(s).

        Returns
        -------
        list of tuples
            A list of feature IDs partitioned by ``self.max_nrecords``.
        """
        if field not in self.client.valid_fields:
            raise InvalidInputValue("field", self.client.valid_fields)

        if not isinstance(ids, Sequence):
            raise InvalidInputType("ids", "str or list")

        ids_ls = [ids] if isinstance(ids, str) else ids

        ftype = self.client.field_types[field]
        if "string" in ftype:
            fids = ", ".join(f"'{i}'" for i in ids_ls)
        else:
            fids = ", ".join(f"{i}" for i in ids_ls)

        return self.oids_bysql(f"{field} IN ({fids})")

    def oids_bysql(self, sql_clause: str) -> Iterator[Tuple[str, ...]]:
        """Get feature IDs using a valid SQL 92 WHERE clause.

        Notes
        -----
        Not all web services support this type of query. For more details look
        `here <https://developers.arcgis.com/rest/services-reference/query-feature-service-.htm#ESRI_SECTION2_07DD2C5127674F6A814CE6C07D39AD46>`__.

        Parameters
        ----------
        sql_clause : str
            A valid SQL 92 WHERE clause.

        Returns
        -------
        list of tuples
            A list of feature IDs partitioned by ``self.max_nrecords``.
        """
        if not isinstance(sql_clause, str):
            raise InvalidInputType("sql_clause", "str")

        payload = {
            "where": sql_clause,
            "returnGeometry": "false",
            "returnIdsOnly": "true",
            "f": self.client.outformat,
        }
        self.client.request_id = uuid.uuid4().hex

        resp = self.client.get_response(self.client.query_url, [payload])[0]
        try:
            return self.partition_oids(resp["objectIds"])
        except KeyError as ex:
            raise ZeroMatched(resp["error"]["message"]) from ex

    def partition_oids(self, oids: Union[List[int], int]) -> Iterator[Tuple[str, ...]]:
        """Partition feature IDs based on ``self.max_nrecords``.

        Parameters
        ----------
        oids : list of int or int
            A list of feature ID(s).

        Returns
        -------
        list of tuples
            A list of feature IDs partitioned by ``self.max_nrecords``.
        """
        return self.client.partition_oids(oids)

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
        return self.client.get_features(featureids, return_m, return_geom)

    def __repr__(self) -> str:
        """Print the service configuration."""
        return self.client.__repr__()


class WMS:
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
        ``epsg:4326``.
    version : str, optional
        The WMS service version which should be either 1.1.1 or 1.3.0, defaults to 1.3.0.
    validation : bool, optional
        Validate the input arguments from the WMS service, defaults to True. Set this
        to False if you are sure all the WMS settings such as layer and crs are correct
        to avoid sending extra requests.
    ssl : bool or SSLContext, optional
        SSLContext to use for the connection, defaults to None. Set to False to disable
        SSL certification verification.
    """

    def __init__(
        self,
        url: str,
        layers: Union[str, List[str]],
        outformat: str,
        version: str = "1.3.0",
        crs: str = DEF_CRS,
        validation: bool = True,
        ssl: Union[SSLContext, bool, None] = None,
    ) -> None:
        self.client = WMSBase(
            url=url,
            layers=layers,
            outformat=outformat,
            version=version,
            crs=crs,
        )
        self.url = self.client.url
        self.outformat = self.client.outformat
        self.version = self.client.version
        self.crs = self.client.crs
        self.ssl = ssl
        self.layers = (
            [self.client.layers] if isinstance(self.client.layers, str) else self.client.layers
        )
        if validation:
            self.client.validate_wms()

    def get_validlayers(self) -> Dict[str, str]:
        """Get the layers supported by the WMS service."""
        return self.client.get_validlayers()

    def getmap_bybox(
        self,
        bbox: Tuple[float, float, float, float],
        resolution: float,
        box_crs: Union[str, pyproj.CRS] = DEF_CRS,
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
        box_crs : str, or pyproj.CRS, optional
            The spatial reference system of the input bbox, defaults to
            ``epsg:4326``.
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

        def _get_payloads(
            args: Tuple[str, Tuple[Tuple[float, float, float, float], str, int, int]]
        ) -> Tuple[str, Dict[str, str]]:
            lyr, bnds = args
            _bbox, counter, _width, _height = bnds

            if self.version != "1.1.1" and pyproj.CRS(self.crs).is_geographic and not always_xy:
                _bbox = (_bbox[1], _bbox[0], _bbox[3], _bbox[2])
            _payload = payload.copy()
            _payload["bbox"] = f'{",".join(str(c) for c in _bbox)}'
            _payload["width"] = str(_width)
            _payload["height"] = str(_height)
            _payload["layers"] = lyr
            return f"{lyr}_dd_{counter}", _payload

        layers, payloads = zip(*(_get_payloads(i) for i in itertools.product(self.layers, bounds)))
        rbinary = ar.retrieve_binary(
            [self.url] * len(payloads),
            [{"params": p} for p in payloads],
            max_workers=4,
            ssl=self.ssl,
        )
        return dict(zip(layers, rbinary))

    def __repr__(self) -> str:
        """Print the services properties."""
        return self.client.__repr__()


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
        ``epsg:4326``.
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
        super().__init__(
            url=url,
            layer=layer,
            outformat=outformat,
            version=version,
            crs=crs,
            read_method=read_method,
            max_nrecords=max_nrecords,
        )

        if validation:
            self.validate_wfs()

    def getfeature_bybox(
        self,
        bbox: Tuple[float, float, float, float],
        box_crs: Union[str, pyproj.CRS] = DEF_CRS,
        always_xy: bool = False,
    ) -> Union[str, bytes, Dict[str, Any]]:
        """Get data from a WFS service within a bounding box.

        Parameters
        ----------
        bbox : tuple
            A bounding box for getting the data: [west, south, east, north]
        box_crs : str, or pyproj.CRS, optional
            The spatial reference system of the input bbox, defaults to
            ``epsg:4326``.
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
        box_crs = pyproj.CRS(box_crs)

        if self.version != "1.0.0" and box_crs.is_geographic and not always_xy:
            bbox = (bbox[1], bbox[0], bbox[3], bbox[2])

        payload = {
            "service": "wfs",
            "version": self.version,
            "outputFormat": self.outformat,
            "request": "GetFeature",
            "typeName": self.layer,
            "bbox": f'{",".join(str(c) for c in bbox)},{box_crs.to_string()}',
            "srsName": self.crs,
        }

        return ar.retrieve([self.url], self.read_method, [{"params": payload}])[0]

    def getfeature_bygeom(
        self,
        geometry: Union[Polygon, MultiPolygon],
        geo_crs: Union[str, pyproj.CRS] = DEF_CRS,
        always_xy: bool = False,
        predicate: str = "INTERSECTS",
    ) -> Union[str, bytes, Dict[str, Any]]:
        """Get features based on a geometry.

        Parameters
        ----------
        geometry : shapely.geometry
            The input geometry
        geo_crs : str, or pyproj.CRS, optional
            The CRS of the input geometry, default to ``epsg:4326``.
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

        if self.version != "1.0.0" and pyproj.CRS(geo_crs).is_geographic and not always_xy:
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
            for fids in tlz.partition_all(self.max_nrecords, set(featureids))
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


class HttpURLs(NamedTuple):
    """URLs of the supported HTTP services."""

    ssebopeta: str = "https://edcintl.cr.usgs.gov/downloads/sciweb1/shared/uswem/web/conus/eta/modis_eta/daily/downloads"


class RESTfulURLs(NamedTuple):
    """URLs of the supported RESTful services."""

    daymet: str = "https://thredds.daac.ornl.gov/thredds/ncss/ornldaac"
    daymet_point: str = "https://daymet.ornl.gov/single-pixel/api/data"
    fema: str = "https://hazards.fema.gov/gis/nfhl/rest/services/public/NFHL/MapServer"
    fws: str = "https://www.fws.gov/wetlandsmapservice/rest/services"
    nldi: str = "https://labs.waterdata.usgs.gov/api/nldi"
    nwis: str = "https://waterservices.usgs.gov/nwis"
    wbd: str = "https://hydro.nationalmap.gov/arcgis/rest/services/wbd/MapServer"
    nhd: str = "https://hydro.nationalmap.gov/arcgis/rest/services/nhd/MapServer"
    nhdplushr: str = "https://hydro.nationalmap.gov/arcgis/rest/services/NHDPlus_HR/MapServer"
    nhdhr_edits: str = "https://edits.nationalmap.gov/arcgis/rest/services/HEM/NHDHigh/MapServer"
    nhdplushr_edits: str = (
        "https://edits.nationalmap.gov/arcgis/rest/services/NHDPlus_HR/NHDPlus_HR/MapServer"
    )
    nhdplus_epa: str = "https://watersgeo.epa.gov/arcgis/rest/services/NHDPlus/NHDPlus/MapServer"
    nhd_fabric: str = (
        "https://watersgeo.epa.gov/arcgis/rest/services/Support/CatchmentFabric/MapServer"
    )
    nid: str = "https://nid.sec.usace.army.mil/api"
    airmap: str = "https://api.airmap.com/elevation/v1/ele"
    nm_pqs: str = "https://nationalmap.gov/epqs/pqs.php"
    pygeoapi: str = "https://labs.waterdata.usgs.gov/api/nldi/pygeoapi/processes"
    nm_3dep_index: str = (
        "https://index.nationalmap.gov/arcgis/rest/services/3DEPElevationIndex/MapServer"
    )


class WFSURLs(NamedTuple):
    """URLs of the supported WFS services."""

    fema: str = "https://hazards.fema.gov/gis/nfhl/services/public/NFHL/MapServer/WFSServer"
    waterdata: str = "https://labs.waterdata.usgs.gov/geoserver/wmadata/ows"


class WMSURLs(NamedTuple):
    """URLs of the supported WMS services."""

    fema: str = "https://hazards.fema.gov/gis/nfhl/rest/services/public/NFHLWMS/MapServer/WMSServer"
    fws: str = (
        "https://www.fws.gov/wetlandsmapservice/services/Wetlands_Raster/ImageServer/WMSServer"
    )
    mrlc: str = "https://www.mrlc.gov/geoserver/mrlc_download/wms"
    nm_3dep: str = (
        "https://elevation.nationalmap.gov/arcgis/services/3DEPElevation/ImageServer/WMSServer"
    )
    gebco: str = (
        "https://www.gebco.net/data_and_products/gebco_web_services/web_map_service/mapserv"
    )


class ServiceURL(NamedTuple):
    """URLs of the supported services."""

    http: HttpURLs = HttpURLs()
    restful: RESTfulURLs = RESTfulURLs()
    wfs: WFSURLs = WFSURLs()
    wms: WMSURLs = WMSURLs()
