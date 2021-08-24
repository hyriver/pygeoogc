"""Tests for PyGeoOGC package."""
import io
import sys
import zipfile
from unittest.mock import patch

import pytest
from shapely.geometry import LineString, Polygon

import pygeoogc as ogc
from pygeoogc import WFS, WMS, ArcGISRESTful, RetrySession, ServiceURL, utils

DEF_CRS = "epsg:4326"
ALT_CRS = "epsg:4269"
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
        """RESTFul by ID"""
        wbd2 = ArcGISRESTful(self.wbd_url, 1, outfields=["huc2", "name", "areaacres"])
        print(wbd2)
        wbd2.featureids = wbd2.partition_oids(list(range(1, 6)))
        huc2 = wbd2.get_features()

        assert len(huc2[0]["features"]) == 5

    def test_geom_point_line(self):
        wbd2 = ArcGISRESTful(self.wbd_url, 1)
        wbd2.oids_bygeom((-70.02580, 44.43280), spatial_relation="esriSpatialRelWithin")
        huc2 = wbd2.get_features()
        point = len(huc2[0]["features"])

        wbd2.oids_bygeom(
            LineString([(-70.02580, 44.43280), (-71.02580, 44.43280)]),
            spatial_relation="esriSpatialRelWithin",
        )
        huc2 = wbd2.get_features()
        line = len(huc2[0]["features"])
        assert point == line

    def test_bygeom(self):
        """RESTFul by geometry"""
        geofab = ArcGISRESTful(self.fab_url, max_workers=4)
        geofab.oids_bygeom(GEO_NAT.bounds)
        geofab.oids_bygeom(GEO_NAT)
        wb_all = geofab.get_features()
        geofab.oids_bygeom(GEO_NAT, sql_clause="areasqkm > 20")
        wb_large = geofab.get_features()
        assert len(wb_all[0]["features"]) - len(wb_large[0]["features"]) == 915

    def test_bymultipoint(self):
        """RESTFul by geometry"""
        geom = [
            (-97.06138, 32.837),
            (-97.06133, 32.836),
            (-97.06124, 32.834),
            (-97.06127, 32.832),
        ]

        service = ArcGISRESTful(self.epa_url, 2, outformat="json")
        service.oids_bygeom(
            geom, geo_crs=ALT_CRS, sql_clause="FTYPE NOT IN (420,428,566)", distance=1500
        )
        resp = service.get_features(return_m=True)
        assert len(resp[0]["features"]) == 3

    def test_retry(self):
        rest_url = "/".join(
            [
                "http://cohgims.houstontx.gov/arcgis/rest/services",
                "PUD_Utility/StormwaterUtilities/MapServer/3",
            ]
        )
        rest = ArcGISRESTful(rest_url)
        print(rest)
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
        rest.featureids = rest.partition_oids(oids)
        resp = rest.get_features()
        assert len(resp) == 3

    @pytest.mark.slow
    def test_bysql(self):
        """RESTFul by SQL filter"""
        hr = ArcGISRESTful(self.nhd_url, 2, outformat="json")
        hr.oids_bysql("PERMANENT_IDENTIFIER IN ('103455178', '103454362', '103453218')")
        resp = hr.get_features(return_m=True)

        assert len(resp[0]["features"]) == 3

    @pytest.mark.slow
    def test_byfield(self):
        """RESTFul by SQL filter"""
        hr = ArcGISRESTful(self.nhd_url, 2, outformat="json")
        hr.oids_byfield("PERMANENT_IDENTIFIER", ["103455178", "103454362", "103453218"])
        resp = hr.get_features()

        assert len(resp[0]["features"]) == 3


class TestWMS:
    wms_url: str = ServiceURL().wms.gebco
    layer: str = "GEBCO_LATEST"

    def test_v111(self):
        """WMS version 1.1.1"""
        wms = WMS(
            self.wms_url,
            layers=self.layer,
            outformat="image/tiff",
            crs=DEF_CRS,
            version="1.1.1",
            validation=False,
        )
        r_dict = wms.getmap_bybox(GEO_NAT.bounds, 20, DEF_CRS)
        assert (
            wms.get_validlayers()[self.layer] == self.layer
            and sys.getsizeof(r_dict[f"{self.layer}_dd_0_0"]) == 11485083
        )

    def test_bybox(self):
        """WMS by bounding box"""
        wms = WMS(self.wms_url, layers=self.layer, outformat="image/tiff", crs=DEF_CRS)
        print(wms)
        r_dict = wms.getmap_bybox(GEO_NAT.bounds, 20, DEF_CRS, max_px=3e6)
        assert sum(sys.getsizeof(r) for r in r_dict.values()) == 11495526


class TestWFS:
    wfs: WFS = WFS(
        ServiceURL().wfs.waterdata,
        layer="wmadata:gagesii",
        outformat="application/json",
        version="2.0.0",
        crs=ALT_CRS,
        max_nrecords=5,
    )

    def test_byid(self):
        """WFS by ID"""
        print(self.wfs)
        st = self.wfs.getfeature_byid("staid", "01031500")
        assert st[0]["numberMatched"] == 1

    def test_bygeom(self):
        """WFS by geometry"""
        r = self.wfs.getfeature_bygeom(GEO_URB, geo_crs=DEF_CRS, always_xy=False)
        assert len(r["features"]) == 7

    def test_bybox(self):
        """WFS by bounding box"""
        bbox = GEO_URB.bounds
        r = self.wfs.getfeature_bybox(bbox, box_crs=DEF_CRS, always_xy=True)
        assert len(r["features"]) == 7

    def test_byfilter(self):
        """WFS by CQL filter"""
        wb = self.wfs.getfeature_byfilter("staid LIKE '010315%'")
        assert len(utils.traverse_json(wb, ["features", "geometry", "coordinates"])) == 2

    def test_wfs111(self):
        """WFS 1.0.0 by geom"""
        wfs = WFS(
            ServiceURL().wfs.waterdata,
            layer="wmadata:gagesii",
            outformat="json",
            version="1.0.0",
            crs=ALT_CRS,
        )
        r = wfs.getfeature_bygeom(GEO_URB, geo_crs=DEF_CRS, always_xy=False)
        assert len(r["features"]) == 7


def test_decompose():
    """Bounding box decomposition"""
    bboxs = utils.bbox_decompose(GEO_URB.bounds, 10)
    assert bboxs[0][-1] == 2828


@pytest.mark.parametrize(
    "geo,gtype,expected",
    [
        (COORDS, "coords", -287.707),
        (GEO_URB, "geometry", -362.099),
        (GEO_URB.bounds, "bounds", -365.403),
    ],
)
def test_matchcrs(geo, gtype, expected):
    """Match CRS"""
    matched = utils.match_crs(geo, DEF_CRS, "epsg:2149")
    if gtype == "coords":
        val = matched[-1][0]
    elif gtype == "geometry":
        val = matched.centroid.x
    elif gtype == "bounds":
        val = matched[0]
    assert (val * 1e-4 - expected) < 1e-3


def test_esriquery():
    """ESRI geometry query builder"""
    point = utils.ESRIGeomQuery(COORDS[0], wkid=DEF_CRS).point()
    line = utils.ESRIGeomQuery(LineString(COORDS), wkid=DEF_CRS).polyline()
    assert (
        point["geometryType"] == "esriGeometryPoint"
        and line["geometryType"] == "esriGeometryPolyline"
    )


def test_ipv4():
    """Only IPv4"""
    url = f"{ServiceURL().http.ssebopeta}/det2005150.modisSSEBopETactual.zip"
    session = RetrySession()
    with patch("socket.has_ipv6", False):
        r = session.get(url)
        z = zipfile.ZipFile(io.BytesIO(r.content))
        fname = z.read(z.filelist[0].filename)

    assert sys.getsizeof(fname) == 20283088


def test_urls():
    """Service URLs"""
    assert len(ServiceURL().__dict__["urls"]) == 4


def test_show_versions():
    """Show versions"""
    f = io.StringIO()
    ogc.show_versions(file=f)
    assert "INSTALLED VERSIONS" in f.getvalue()
