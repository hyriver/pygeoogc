"""Tests for PyGeoOGC package."""

from __future__ import annotations

import io
import sys
from dataclasses import fields

import pandas as pd
import pytest
from shapely import LineString, Polygon

import pygeoogc as ogc
from pygeoogc import WFS, WMS, ArcGISRESTful, ServiceURL, utils

DEF_CRS = 4326
ALT_CRS = 4269
GEO_NAT = Polygon(
    [[-69.77, 45.07], [-69.31, 45.07], [-69.31, 45.45], [-69.77, 45.45], [-69.77, 45.07]]
)
GEO_URB = Polygon(
    [
        [-118.72, 34.118],
        [-118.31, 34.118],
        [-118.31, 34.518],
        [-118.72, 34.518],
        [-118.72, 34.118],
    ]
)
COORDS = [(-118.72, 34.118), (-118.31, 34.518)]


class TestREST:
    wbd_url: str = ServiceURL().restful.wbd
    fab_url: str = f"{ServiceURL().restful.nhd_fabric}/1"
    epa_url: str = ServiceURL().restful.nhdplus_epa
    nhd_url: str = ServiceURL().restful.nhdplushr

    def test_byid(self):
        """RESTFul by ID."""
        wbd2 = ArcGISRESTful(self.wbd_url, 1, outfields=["huc2", "name", "areaacres"])
        assert self.wbd_url in wbd2.__repr__()
        huc2 = wbd2.get_features(wbd2.partition_oids(list(range(1, 6))))

        assert len(huc2[0]["features"]) == 5

    def test_geom_point_line(self):
        wbd2 = ArcGISRESTful(self.wbd_url, 1)
        oids = wbd2.oids_bygeom((-70.02580, 44.43280), spatial_relation="esriSpatialRelWithin")
        huc2 = wbd2.get_features(oids)
        point = len(huc2[0]["features"])

        oids = wbd2.oids_bygeom(
            LineString([(-70.02580, 44.43280), (-71.02580, 44.43280)]),
            spatial_relation="esriSpatialRelWithin",
        )
        huc2 = wbd2.get_features(oids)
        line = len(huc2[0]["features"])
        assert point == line

    def test_bygeom(self):
        """RESTFul by geometry."""
        geofab = ArcGISRESTful(self.fab_url)
        _ = geofab.oids_bygeom(GEO_NAT.bounds)
        oids = geofab.oids_bygeom(GEO_NAT)
        wb_all = geofab.get_features(oids)
        oids = geofab.oids_bygeom(GEO_NAT, sql_clause="areasqkm > 20")
        wb_large = geofab.get_features(oids)
        assert len(wb_all[0]["features"]) - len(wb_large[0]["features"]) == 915

    def test_bymultipoint(self):
        """RESTFul by geometry."""
        geom = [
            (-97.06138, 32.837),
            (-97.06133, 32.836),
            (-97.06124, 32.834),
            (-97.06127, 32.832),
        ]

        service = ArcGISRESTful(self.epa_url, 2, outformat="json")
        oids = service.oids_bygeom(
            geom, geo_crs=ALT_CRS, sql_clause="FTYPE NOT IN (420,428,566)", distance=1500
        )
        resp = service.get_features(oids, return_m=True)
        assert len(resp[0]["features"]) == 3

    @pytest.mark.skip(reason="Service is not available anymore.")
    def test_retry(self):
        rest_url = "/".join(
            [
                "http://cohgims.houstontx.gov/arcgis/rest/services",
                "PUD_Utility/StormwaterUtilities/MapServer/3",
            ]
        )
        rest = ArcGISRESTful(rest_url, verbose=True)
        assert rest_url in rest.__repr__()
        oids = [
            1006322,
            1006323,
            1006324,
            1245206,
            323351,
            324961,
            474077,
            503848,
            584849,
            585138,
        ]
        resp = rest.get_features(rest.partition_oids(oids))

        with rest.client.failed_path.open() as f:
            f_oids = [int(i) for i in f.read().splitlines()]
        assert len(resp) == 3
        assert len(f_oids) == (len(oids) - len(resp))

    def test_bysql(self):
        """RESTFul by SQL filter."""
        hr = ArcGISRESTful(self.nhd_url, 2, outformat="json")
        oids = hr.oids_bysql("nhdplusid IN (50000100239005, 50000100235192)")
        resp = hr.get_features(oids, return_m=True)

        assert len(resp[0]["features"]) == 2

    def test_byfield(self):
        """RESTFul by SQL filter."""
        hr = ArcGISRESTful(self.nhd_url, 2, outformat="json")
        oids = hr.oids_byfield("nhdplusid", ["50000100239005", "50000100235192"])
        resp = hr.get_features(oids)

        assert len(resp[0]["features"]) == 2


def test_retrysession_head():
    url = "https://httpbin.org/image/jpeg"
    with utils.RetrySession() as session:
        resp = session.head(url)
        assert resp.headers["Content-length"] == "35588"


def test_retrysession_get_length():
    url = "https://httpbin.org/image/jpeg"
    with utils.RetrySession(disable=True) as session:
        resp = session.get(url, stream=True)
        assert resp.headers["Content-length"] == "35588"


def test_retrysession_params():
    url = "https://postman-echo.com"
    params = {"params_get": "foo"}
    data = {"data_post": "foo"}
    with utils.RetrySession(disable=True) as session:
        get1 = session.get(f"{url}/get", params=params)
        get2 = session.get(f"{url}/get", payload=params)
        assert get1.json()["args"] == get2.json()["args"]

        post1 = session.post(f"{url}/post", data=data)
        post2 = session.post(f"{url}/post", payload=data)
        assert post1.json()["args"] == post2.json()["args"]


def test_stream():
    url = "http://ipv4.download.thinkbroadband.com/5MB.zip"
    fname = ogc.streaming_download(url)
    assert fname.stat().st_size == 5242880

    urls = (
        url,
        "http://ipv4.download.thinkbroadband.com:81/5MB.zip",
    )
    fname = ogc.streaming_download(urls)
    assert fname[0].stat().st_size == 5242880
    assert fname[1].stat().st_size == 5242880


@pytest.mark.filterwarnings("ignore:.*Content metadata*.")
class TestWMS:
    wms_url: str = ServiceURL().wms.gebco
    layer: str = "GEBCO_LATEST"

    def test_v111(self):
        """WMS version 1.1.1."""
        wms = WMS(
            self.wms_url,
            layers=self.layer,
            outformat="image/tiff",
            crs=DEF_CRS,
            version="1.1.1",
            validation=False,
        )
        r_dict = wms.getmap_bybox(GEO_NAT.bounds, 20, DEF_CRS)
        assert wms.get_validlayers()[self.layer] == self.layer
        assert sys.getsizeof(r_dict[f"{self.layer}_dd_0"]) == 11501067

    def test_bybox(self):
        """WMS by bounding box."""
        wms = WMS(self.wms_url, layers=self.layer, outformat="image/tiff", crs=DEF_CRS)
        assert self.wms_url in wms.__repr__()
        r_dict = wms.getmap_bybox(GEO_NAT.bounds, 20, DEF_CRS, max_px=int(3e6))
        assert sum(sys.getsizeof(r) for r in r_dict.values()) == 11498778
        flist = wms.getmap_bybox(GEO_NAT.bounds, 20, DEF_CRS, max_px=int(3e6), tiff_dir="cache")
        assert sum(f.stat().st_size for f in flist) == 11498712

    def test_valid_crs(self):
        """Get WMS valid CRSs."""
        crs = utils.valid_wms_crs(self.wms_url)
        assert sorted(crs) == ["epsg:3395", "epsg:3857", "epsg:4326"]


class TestWFS:
    wfs: WFS = WFS(
        ServiceURL().wfs.waterdata,
        layer="wmadata:gagesii",
        outformat="csv",
        read_method="text",
        version="2.0.0",
        crs=ALT_CRS,
        max_nrecords=5,
    )

    def to_df(self, resp):
        return pd.concat(pd.read_csv(io.StringIO(r)) for r in resp)

    def test_byid(self):
        """WFS by ID."""
        assert ServiceURL().wfs.waterdata in self.wfs.__repr__()
        stations = [
            "01011000",
            "01013500",
            "01015800",
            "01016500",
            "01017000",
            "01017060",
        ]
        resp = self.wfs.getfeature_byid("staid", stations)
        df = self.to_df(resp)
        assert set(df.staid.unique()) == {int(s) for s in stations}

    def test_bygeom(self):
        """WFS by geometry."""
        resp = self.wfs.getfeature_bygeom(GEO_URB, geo_crs=DEF_CRS, always_xy=False)
        assert self.to_df(resp).shape[0] == 7

    def test_bybox(self):
        """WFS by bounding box."""
        bbox = GEO_URB.bounds
        resp = self.wfs.getfeature_bybox(bbox, box_crs=DEF_CRS, always_xy=True)
        assert self.to_df(resp).shape[0] == 7

    def test_byfilter(self):
        """WFS by CQL filter."""
        resp = self.wfs.getfeature_byfilter("staid LIKE '010315%'")
        assert self.to_df(resp).shape[0] == 2

    def test_wfs110(self):
        """WFS 1.1.0 by geom."""
        resp = WFS(
            ServiceURL().wfs.waterdata,
            layer="wmadata:gagesii",
            outformat="json",
            version="1.1.0",
            crs=ALT_CRS,
        ).getfeature_bygeom(GEO_URB, geo_crs=DEF_CRS, always_xy=False)
        assert len(resp[0]["features"]) == 7


def test_decompose():
    """Bounding box decomposition."""
    bboxs = utils.bbox_decompose(GEO_URB.bounds, 10)
    assert bboxs[0][-1] == (-118.515, 34.318, -118.31, 34.518)


@pytest.mark.parametrize(
    ("geo", "gtype", "expected"),
    [
        (COORDS, "coords", -287.707),
        (GEO_URB, "geometry", -362.099),
        (GEO_URB.bounds, "bounds", -365.403),
    ],
)
def test_matchcrs(geo, gtype, expected):
    """Match CRS."""
    matched = utils.match_crs(geo, DEF_CRS, 2149)
    if gtype == "coords":
        val = matched[-1][0]
    elif gtype == "geometry":
        val = matched.centroid.x
    elif gtype == "bounds":
        val = matched[0]
    assert (val * 1e-4 - expected) < 1e-3


def test_esriquery():
    """ESRI geometry query builder."""
    point = utils.ESRIGeomQuery(COORDS[0], wkid=DEF_CRS).point()
    line = utils.ESRIGeomQuery(LineString(COORDS), wkid=DEF_CRS).polyline()
    assert point["geometryType"] == "esriGeometryPoint"
    assert line["geometryType"] == "esriGeometryPolyline"


def test_urls():
    """Service URLs."""
    total_urls = 0
    service_url_class = ServiceURL()
    for field in fields(service_url_class):
        nested_class = getattr(service_url_class, field.name)
        total_urls += sum(
            isinstance(getattr(nested_class, f.name), str) for f in fields(nested_class)
        )
    assert total_urls == 36


def test_show_versions():
    """Show versions."""
    f = io.StringIO()
    ogc.show_versions(file=f)
    assert "SYS INFO" in f.getvalue()
