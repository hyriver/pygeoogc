"""Base classes and function for REST, WMS, and WMF services."""

from __future__ import annotations

import itertools
import uuid
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any, Literal, cast, overload

import cytoolz.curried as tlz
import pyproj
import shapely
from shapely import ops

import async_retriever as ar
from pygeoogc import utils
from pygeoogc.core import ArcGISRESTfulBase, WFSBase, WMSBase
from pygeoogc.exceptions import InputTypeError, InputValueError, ServiceError, ZeroMatchedError

if TYPE_CHECKING:
    from collections.abc import Iterator, Sequence

    from pyproj import CRS
    from shapely import LineString, MultiPoint, MultiPolygon, Point, Polygon

    CRSType = int | str | CRS
    Response = list[str] | list[bytes] | list[dict[str, Any]] | list[list[dict[str, Any]]]


class ArcGISRESTful:
    """Access to an ArcGIS REST service.

    Notes
    -----
    By default, all retrieval methods retry to get the missing feature IDs,
    if there are any. You can disable this behavior by setting ``disable_retry``
    to ``True``. If there are any missing feature IDs after the retry,
    they are saved to a text file, ipath of which can be accessed by
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
    crs : str, int, or pyproj.CRS, optional
        The spatial reference of the output data, defaults to ``epsg:4326``.
    verbose : bool, optional
        If True, prints information about the requests and responses,
        defaults to False.
    disable_retry : bool, optional
        If ``True`` in case there are any failed queries, no retrying attempts
        is done and object IDs of the failed requests is saved to a text file
        which its ipath can be accessed via ``self.client.failed_path``.
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
        self.client = ArcGISRESTfulBase(
            base_url=base_url,
            layer=layer,
            outformat=outformat,
            outfields=outfields,
            crs=crs,
            verbose=verbose,
            disable_retry=disable_retry,
        )

    def oids_bygeom(
        self,
        geom: (
            LineString
            | Polygon
            | Point
            | MultiPoint
            | tuple[float, float]
            | list[tuple[float, float]]
            | tuple[float, float, float, float]
        ),
        geo_crs: CRSType = 4326,
        spatial_relation: str = "esriSpatialRelIntersects",
        sql_clause: str | None = None,
        distance: int | None = None,
    ) -> Iterator[tuple[str, ...]]:
        """Get feature IDs within a geometry that can be combined with a SQL where clause.

        Parameters
        ----------
        geom : LineString, Polygon, Point, MultiPoint, tuple, or list of tuples
            A geometry (LineString, Polygon, Point, MultiPoint), tuple of length two
            (``(x, y)``), a list of tuples of length 2 (``[(x, y), ...]``), or bounding box
            (tuple of length 4 (``(xmin, ymin, xmax, ymax)``)).
        geo_crs : str, int, or pyproj.CRS, optional
            The spatial reference of the input geometry, defaults to ``epsg:4326``.
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
            raise InputValueError("spatial_relation", valid_spatialrels)

        geom_query = utils.esri_query(geom, geo_crs, self.client.crs)
        payload = {
            **geom_query,
            "spatialRel": spatial_relation,
            "returnGeometry": "false",
            "returnIdsOnly": "true",
            "f": self.client.outformat,
        }
        if distance:
            payload.update({"distance": str(int(distance)), "units": "esriSRUnit_Meter"})

        if sql_clause:
            payload.update({"where": sql_clause})

        self.client.request_id = uuid.uuid4().hex

        resp = self.client.get_response(self.client.query_url, [payload], method="POST")[0]
        try:
            return self.partition_oids(resp["objectIds"])
        except (KeyError, TypeError) as ex:
            msg = resp["error"]["message"] if "error" in resp else "No matched records"
            raise ZeroMatchedError(msg) from ex

    def oids_byfield(
        self, field: str, ids: str | int | list[str] | list[int]
    ) -> Iterator[tuple[str, ...]]:
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
            raise InputValueError("field", self.client.valid_fields)

        ids_ls = {ids} if isinstance(ids, (str, int)) else set(ids)

        if "string" in self.client.field_types.get(field, ""):
            fids = ", ".join(f"'{i}'" for i in ids_ls)
        else:
            fids = ", ".join(str(i) for i in ids_ls)

        return self.oids_bysql(f"{field} IN ({fids})")

    def oids_bysql(self, sql_clause: str) -> Iterator[tuple[str, ...]]:
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
            raise InputTypeError("sql_clause", "str")

        payload = {
            "where": sql_clause,
            "returnGeometry": "false",
            "returnIdsOnly": "true",
            "f": self.client.outformat,
        }
        self.client.request_id = uuid.uuid4().hex

        resp = self.client.get_response(self.client.query_url, [payload], method="POST")[0]
        try:
            return self.partition_oids(resp["objectIds"])
        except (KeyError, TypeError) as ex:
            msg = resp["error"]["message"] if "error" in resp else "No matched records"
            raise ZeroMatchedError(msg) from ex

    def partition_oids(self, oids: list[int] | int) -> Iterator[tuple[str, ...]]:
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
    crs : str, int, or pyproj.CRS, optional
        The spatial reference system to be used for requesting the data, defaults to
        ``epsg:4326``.
    version : str, optional
        The WMS service version which should be either 1.1.1 or 1.3.0, defaults to 1.3.0.
    validation : bool, optional
        Validate the input arguments from the WMS service, defaults to True. Set this
        to False if you are sure all the WMS settings such as layer and crs are correct
        to avoid sending extra requests.
    ssl : bool, optional
        Whether to use SSL for the connection, defaults to ``True``.
    """

    def __init__(
        self,
        url: str,
        layers: str | list[str],
        outformat: str,
        version: str = "1.3.0",
        crs: CRSType = 4326,
        validation: bool = True,
        ssl: bool = True,
    ) -> None:
        self.client = WMSBase(
            url=url,
            layers=layers,
            outformat=outformat,
            version=version,
            crs=crs,
            validation=validation,
        )
        self.url = self.client.url
        self.outformat = self.client.outformat
        self.version = self.client.version
        self.crs = self.client.crs
        self.crs_str = self.client.crs_str
        self.is_geographic = pyproj.CRS(self.crs).is_geographic
        self.ssl = ssl
        self.layers = self.client.layers

    def get_validlayers(self) -> dict[str, str]:
        """Get the layers supported by the WMS service."""
        return self.client.get_validlayers()

    @overload
    def getmap_bybox(
        self,
        bbox: tuple[float, float, float, float],
        resolution: float,
        box_crs: CRSType = ...,
        always_xy: bool = ...,
        max_px: int = ...,
        kwargs: dict[str, Any] | None = ...,
        tiff_dir: Literal[None] = None,
    ) -> dict[str, bytes]: ...

    @overload
    def getmap_bybox(
        self,
        bbox: tuple[float, float, float, float],
        resolution: float,
        box_crs: CRSType = ...,
        always_xy: bool = ...,
        max_px: int = ...,
        kwargs: dict[str, Any] | None = ...,
        tiff_dir: str | Path = ...,
    ) -> list[Path]: ...

    def getmap_bybox(
        self,
        bbox: tuple[float, float, float, float],
        resolution: float,
        box_crs: CRSType = 4326,
        always_xy: bool = False,
        max_px: int = 8000000,
        kwargs: dict[str, Any] | None = None,
        tiff_dir: str | Path | None = None,
    ) -> dict[str, bytes] | list[Path]:
        """Get data from a WMS service within a geometry or bounding box.

        Parameters
        ----------
        bbox : tuple
            A bounding box for getting the data.
        resolution : float
            The output resolution in meters. The width and height of output are computed in pixel
            based on the geometry bounds and the given resolution.
        box_crs : str, int, or pyproj.CRS, optional
            The spatial reference system of the input bbox, defaults to
            ``epsg:4326``.
        always_xy : bool, optional
            Whether to always use xy axis order, defaults to False. Some services change the axis
            order from xy to yx, following the latest WFS version specifications but some don't.
            If the returned value does not have any geometry, it indicates that most probably the
            axis order does not match. You can set this to True in that case.
        max_px : int, optional
            The maximum allowable number of pixels (width x height) for a WMS requests,
            defaults to 8 million based on some trial-and-error.
        kwargs: dict, optional
            Optional additional keywords passed as payload, defaults to None.
            For example, ``{"styles": "default"}``.
        tiff_dir : str or pathlib.Path, optional
            If given, the retrieved data will be stored on disk instead of
            returning it, defaults to ``None``, i.e., saving to memory
            and returning the data.

        Returns
        -------
        dict of bytes or list of pathlib.Path
            If ``to_disk=False``, a dict where the keys are the layer name and
            values are the returned response from the WMS service as bytes.
            If ``to_disk=True``, a list of pathlib.Path objects to the saved files.
        """
        utils.check_bbox(bbox)
        _bbox = utils.match_crs(bbox, box_crs, self.crs_str)
        bounds, sub_width, sub_height = utils.bbox_decompose(
            _bbox, resolution, self.crs_str, max_px
        )

        payload = {
            "version": self.version,
            "format": self.outformat,
            "request": "GetMap",
        }

        if kwargs is not None and not isinstance(kwargs, dict):
            raise InputTypeError("kwargs", "dict or None")

        if isinstance(kwargs, dict):
            payload.update(kwargs)

        if self.version == "1.1.1":
            payload["srs"] = self.crs_str
        else:
            payload["crs"] = self.crs_str

        precision = 2 if pyproj.CRS(self.crs).is_projected else 6

        def _get_payloads(
            lyr: str, bbox: tuple[float, float, float, float], counter: int
        ) -> tuple[str, dict[str, str]]:
            if self.version != "1.1.1" and self.is_geographic and not always_xy:
                bbox = (bbox[1], bbox[0], bbox[3], bbox[2])
            _payload = payload.copy()
            _payload["bbox"] = ",".join(str(round(c, precision)) for c in bbox)
            _payload["width"] = str(sub_width)
            _payload["height"] = str(sub_height)
            _payload["layers"] = lyr
            return f"{lyr}_dd_{counter}", _payload

        _lyr_payloads = (
            _get_payloads(lyr, box, i)
            for lyr, (i, box) in itertools.product(self.layers, enumerate(bounds))
        )
        layers, payloads = zip(*_lyr_payloads)
        layers = cast("tuple[str]", layers)
        payloads = cast("tuple[dict[str, str]]", payloads)
        if tiff_dir is not None:
            tiff_dir = Path(tiff_dir)
            tiff_dir.mkdir(parents=True, exist_ok=True)
            tiff_files = utils.streaming_download(
                [self.url] * len(payloads),
                [{"params": p} for p in payloads],
                root_dir=tiff_dir,
                file_extention="tiff",
                n_jobs=4,
                ssl=self.ssl,
            )
            tiff_files = [f for f in tiff_files if f is not None]
            return tiff_files
        rbinary = ar.retrieve_binary(
            [self.url] * len(payloads),
            [{"params": p} for p in payloads],
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
    crs : str, int, or pyproj.CRS, optional
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
        layer: str | None = None,
        outformat: str | None = None,
        version: str = "2.0.0",
        crs: CRSType = 4326,
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
            validation=validation,
        )

    def getfeature_bybox(
        self,
        bbox: tuple[float, float, float, float],
        box_crs: CRSType = 4326,
        always_xy: bool = False,
        sort_attr: str | None = None,
    ) -> Response:
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
        sort_attr : str, optional
            The column name in the database to sort request by, defaults
            to the first attribute in the schema that contains ``id`` in its name.

        Returns
        -------
        list of str or bytes or dict
            WFS query response within a bounding box.
        """
        utils.check_bbox(bbox)
        box_crs = pyproj.CRS(box_crs)

        if box_crs.is_geographic and not always_xy:
            bbox = (bbox[1], bbox[0], bbox[3], bbox[2])

        precision = 2 if box_crs.is_projected else 6
        bbox_str = f"{','.join(str(round(c, precision)) for c in bbox)},{box_crs.to_string()}"
        payload = {
            "service": "wfs",
            "version": self.version,
            "outputFormat": "text/xml",
            "request": "GetFeature",
            "typeName": self.layer,
            "bbox": bbox_str,
            "srsName": self.crs_str,
            "resultType": "hits",
        }
        resp = ar.retrieve_text([self.url], [{"params": payload}])
        try:
            nfeatures = int(resp[0].split(self.nfeat_key)[-1].split(" ")[0].strip('"'))
        except (IndexError, ValueError) as ex:
            raise ServiceError(resp[0], self.url) from ex

        payloads = [
            {
                "service": "wfs",
                "version": self.version,
                "outputFormat": self.outformat,
                "request": "GetFeature",
                "typeName": self.layer,
                "bbox": bbox_str,
                "srsName": self.crs_str,
                **self.sort_params(sort_attr, nfeatures, i),
            }
            for i in range(0, nfeatures, self.max_nrecords)
        ]

        return self.retrieve([self.url] * len(payloads), [{"params": p} for p in payloads])

    def getfeature_bygeom(
        self,
        geometry: Polygon | MultiPolygon,
        geo_crs: CRSType = 4326,
        always_xy: bool = False,
        predicate: str = "INTERSECTS",
        sort_attr: str | None = None,
    ) -> Response:
        """Get features based on a geometry.

        Parameters
        ----------
        geometry : shapely.Polygon or shapely.MultiPolygon
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

        sort_attr : str, optional
            The column name in the database to sort request by, defaults
            to the first attribute in the schema that contains ``id`` in its name.

        Returns
        -------
        str or bytes or dict
            WFS query response based on the given geometry.
        """
        geom = utils.match_crs(geometry, geo_crs, self.crs_str)
        precision = 1e-2 if pyproj.CRS(self.crs_str).is_projected else 1e-6
        geom = shapely.set_precision(geom, precision)
        geom_name = ""
        if "geometry_column" in self.schema[self.layer]:
            geom_name = self.schema[self.layer]["geometry_column"]
        elif "properties" in self.schema[self.layer]:
            geom_name = next(
                (lyr for lyr in self.schema[self.layer]["properties"] if "geom" in lyr), ""
            )
        elif "required" in self.schema[self.layer]:
            geom_name = next(
                (lyr for lyr in self.schema[self.layer]["required"] if "geom" in lyr), ""
            )

        if not geom_name:
            msg = "Cannot find the geometry column name in the schema."
            raise ValueError(msg)

        if pyproj.CRS(self.crs_str).is_geographic and not always_xy:
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
        if predicate.upper() not in valid_predicates:
            raise InputValueError("predicate", valid_predicates)

        return self.getfeature_byfilter(
            f"{predicate.upper()}({geom_name}, {g_wkt})", method="POST", sort_attr=sort_attr
        )

    def getfeature_byid(
        self,
        featurename: str,
        featureids: Sequence[int | str] | int | str,
    ) -> Response:
        """Get features based on feature IDs.

        Parameters
        ----------
        featurename : str
            The name of the column for searching for feature IDs.
        featureids : int, str, or list of them
            The feature ID(s).

        Returns
        -------
        str or bytes or dict
            WMS query response.
        """
        valid_features = self.schema[self.layer]["properties"]
        if featurename not in valid_features:
            raise InputValueError("featurename", list(valid_features))

        feat_set = {featureids} if isinstance(featureids, (str, int)) else set(featureids)

        if not feat_set:
            raise InputTypeError("featureids", "int, str, or list")

        feat_itr = tlz.partition_all(self.max_nrecords, feat_set)
        if "str" in valid_features[featurename]:
            feat_vals = (", ".join(f"'{fid}'" for fid in fids) for fids in feat_itr)
        else:
            feat_vals = (", ".join(str(fid) for fid in fids) for fids in feat_itr)

        return list(
            tlz.concat(
                self.getfeature_byfilter(
                    f"{featurename} IN ({f})", method="POST", sort_attr=featurename
                )
                for f in feat_vals
            )
        )

    def getfeature_byfilter(
        self, cql_filter: str, method: str = "GET", sort_attr: str | None = None
    ) -> Response:
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
        sort_attr : str, optional
            The column name in the database to sort request by, defaults
            to the first attribute in the schema that contains ``id`` in its name.

        Returns
        -------
        str or bytes or dict
            WFS query response
        """
        if not isinstance(cql_filter, str):
            raise InputTypeError("cql_filter", "str")

        valid_methods = ["GET", "POST"]
        if method not in valid_methods:
            raise InputValueError("method", valid_methods)

        payload = {
            "service": "wfs",
            "version": self.version,
            "outputFormat": "text/xml",
            "request": "GetFeature",
            "typeName": self.layer,
            "srsName": self.crs_str,
            "cql_filter": cql_filter,
            "resultType": "hits",
        }

        if method == "GET":
            resp = ar.retrieve_text([self.url], [{"params": payload}])
        else:
            headers = {"content-type": "application/x-www-form-urlencoded"}
            resp = ar.retrieve_text([self.url], [{"data": payload, "headers": headers}], "POST")
        try:
            nfeatures = int(resp[0].split(self.nfeat_key)[-1].split(" ")[0].strip('"'))
        except (IndexError, ValueError) as ex:
            raise ServiceError(resp[0], self.url) from ex

        payloads = [
            {
                "service": "wfs",
                "version": self.version,
                "outputFormat": self.outformat,
                "request": "GetFeature",
                "typeName": self.layer,
                "srsName": self.crs_str,
                "cql_filter": cql_filter,
                **self.sort_params(sort_attr, nfeatures, i),
            }
            for i in range(0, nfeatures, self.max_nrecords)
        ]
        if method == "GET":
            return self.retrieve([self.url] * len(payloads), [{"params": p} for p in payloads])

        headers = {"content-type": "application/x-www-form-urlencoded"}
        return self.retrieve(
            [self.url] * len(payloads),
            [{"data": p, "headers": headers} for p in payloads],
            "POST",
        )


@dataclass(frozen=True)
class HttpURLs:
    """URLs of the supported HTTP services."""

    ssebopeta: str = "https://edcintl.cr.usgs.gov/downloads/sciweb1/shared/uswem/web/conus/eta/modis_eta/daily/downloads"


@dataclass(frozen=True)
class RESTfulURLs:
    """URLs of the supported RESTful services."""

    daymet: str = "https://thredds.daac.ornl.gov/thredds/ncss/ornldaac"
    daymet_point: str = "https://daymet.ornl.gov/single-pixel/api/data"
    gridmet: str = "http://thredds.northwestknowledge.net:8080/thredds/ncss"
    fema_nfhl: str = "https://hazards.fema.gov/arcgis/rest/services/public/NFHL/MapServer"
    fema_prelim_cslf: str = (
        "https://hazards.fema.gov/arcgis/rest/services/CSLF/Prelim_CSLF/MapServer"
    )
    fema_draft_cslf: str = "https://hazards.fema.gov/arcgis/rest/services/CSLF/Draft_CSLF/MapServer"
    fema_prelim_nfhl: str = (
        "https://hazards.fema.gov/arcgis/rest/services/PrelimPending/Prelim_NFHL/MapServer"
    )
    fema_pending_nfhl: str = (
        "https://hazards.fema.gov/arcgis/rest/services/PrelimPending/Pending_NFHL/MapServer"
    )
    fema_draft_nfhl: str = (
        "https://hazards.fema.gov/arcgis/rest/services/AFHI/Draft_FIRM_DB/MapServer"
    )
    fws: str = "https://www.fws.gov/wetlandsmapservice/rest/services"
    nldi: str = "https://api.water.usgs.gov/nldi"
    nwis: str = "https://waterservices.usgs.gov/nwis"
    wbd: str = "https://hydro.nationalmap.gov/arcgis/rest/services/wbd/MapServer"
    nhd: str = "https://hydro.nationalmap.gov/arcgis/rest/services/nhd/MapServer"
    hp3d: str = "https://hydro.nationalmap.gov/arcgis/rest/services/3DHP_all/MapServer"
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
    nm_pqs: str = "https://epqs.nationalmap.gov/v1/json"
    pygeoapi: str = "https://api.water.usgs.gov/nldi/pygeoapi/processes"
    nm_3dep_index: str = (
        "https://index.nationalmap.gov/arcgis/rest/services/3DEPElevationIndex/MapServer"
    )
    geoconnex: str = "https://reference.geoconnex.us"
    stnflood: str = "https://stn.wim.usgs.gov/STNServices"
    ehydro: str = "https://services7.arcgis.com/n1YM8pTrFmm7L4hs/ArcGIS/rest/services/eHydro_Survey_Data/FeatureServer/0"
    ehydro_bins: str = "https://spatial.usace.army.mil/opjarcgis/rest/services/ehydro/SurveyJobBinsH3/FeatureServer/0"


@dataclass(frozen=True)
class WFSURLs:
    """URLs of the supported WFS services."""

    fema: str = "https://hazards.fema.gov/arcgis/rest/services/public/NFHL/MapServer/WFSServer"
    waterdata: str = "https://api.water.usgs.gov/geoserver/wmadata/ows"


@dataclass(frozen=True)
class WMSURLs:
    """URLs of the supported WMS services."""

    fema: str = "https://hazards.fema.gov/arcgis/rest/services/public/NFHLWMS/MapServer/WMSServer"
    fws: str = (
        "https://www.fws.gov/wetlandsmapservice/services/Wetlands_Raster/ImageServer/WMSServer"
    )
    mrlc: str = "https://www.mrlc.gov/geoserver/mrlc_download/wms"
    nm_3dep: str = (
        "https://elevation.nationalmap.gov/arcgis/services/3DEPElevation/ImageServer/WMSServer"
    )
    gebco: str = "https://wms.gebco.net/mapserv"


@dataclass(frozen=True)
class ServiceURL:
    """URLs of the supported services."""

    http: HttpURLs = field(default_factory=HttpURLs)
    restful: RESTfulURLs = field(default_factory=RESTfulURLs)
    wfs: WFSURLs = field(default_factory=WFSURLs)
    wms: WMSURLs = field(default_factory=WMSURLs)
