"""Tests for PyGeoOGC package."""
import io
import sys
import zipfile
from datetime import datetime

import pytest
from shapely.geometry import Polygon

import pygeoogc
from pygeoogc import WFS, WMS, ArcGISRESTful, MatchCRS, RetrySession, ServiceURL, utils

DEF_CRS = "epsg:4326"
ALT_CRS = "epsg:2149"


@pytest.fixture
def geometry_nat():
    return Polygon(
        [[-69.77, 45.07], [-69.31, 45.07], [-69.31, 45.45], [-69.77, 45.45], [-69.77, 45.07]]
    )


@pytest.fixture
def geometry_urb():
    return Polygon(
        [
            [-118.72, 34.118],
            [-118.31, 34.118],
            [-118.31, 34.518],
            [-118.72, 34.518],
            [-118.72, 34.118],
        ]
    )


@pytest.fixture
def wfs():
    return WFS(
        ServiceURL().wfs.waterdata,
        layer="wmadata:gagesii",
        outformat="application/json",
        version="2.0.0",
        crs="epsg:4269",
    )


@pytest.mark.flaky(max_runs=3)
def test_restful_byid(geometry_nat):
    wbd2 = ArcGISRESTful(ServiceURL().restful.wbd)
    wbd2.layer = 1
    print(wbd2)
    wbd2.max_nrecords = 1
    wbd2.spatial_relation = "esriSpatialRelIntersects"
    wbd2.outformat = "geojson"
    wbd2.featureids = list(range(1, 6))
    fields = [f.lower() for f in wbd2.get_validfields()]
    wbd2.outfields = ["huc2", "name", "areaacres"]
    huc2 = wbd2.get_features()

    geom_type = utils.traverse_json(huc2, ["features", "geometry", "type"])

    assert (
        sum(len(h) for h in huc2) == 15 and ["MultiPolygon"] in geom_type and "areaacres" in fields
    )


@pytest.mark.flaky(max_runs=3)
def test_restful_bygeom(geometry_nat):
    geofab = ArcGISRESTful(f"{ServiceURL().restful.nhd_fabric}/1")
    geofab.n_threads = 4
    geofab.oids_bygeom(geometry_nat.bounds)
    geofab.oids_bygeom(geometry_nat)
    wb_all = geofab.get_features()
    geofab.oids_bygeom(geometry_nat, sql_clause="areasqkm > 20")
    wb_large = geofab.get_features()
    assert len(wb_all[0]["features"]) - len(wb_large[0]["features"]) == 915


@pytest.mark.flaky(max_runs=3)
def test_restful_bymultipoint():
    url = "https://watersgeo.epa.gov/arcgis/rest/services/NHDPlus/NHDPlus/MapServer"
    sql_clause = "FTYPE NOT IN (420,428,566)"
    geom = [
        (-97.06138, 32.837),
        (-97.06133, 32.836),
        (-97.06124, 32.834),
        (-97.06127, 32.832),
    ]
    geo_crs = "epsg:4269"
    distance = 1500

    service = ArcGISRESTful(url, outformat="json")
    service.layer = 2  # network flowline
    service.oids_bygeom(geom, geo_crs=geo_crs, sql_clause=sql_clause, distance=distance)
    resp = service.get_features(return_m=True)
    assert len(resp[0]["features"]) == 3


@pytest.mark.flaky(max_runs=3)
def test_restful_bysql():
    hr = ArcGISRESTful(ServiceURL().restful.nhdplushr_edits, outformat="json")
    hr.layer = 1
    hr.layer = 2
    hr.oids_bysql("NHDFlowline.PERMANENT_IDENTIFIER IN ('103455178', '103454362', '103453218')")
    resp_sql = hr.get_features(return_m=True)
    hr.oids_byfield("NHDFlowline.PERMANENT_IDENTIFIER", ["103455178", "103454362", "103453218"])
    resp_ids = hr.get_features()

    assert len(resp_sql[0]["features"]) == len(resp_ids[0]["features"])


@pytest.mark.flaky(max_runs=3)
def test_wms(geometry_nat):
    url_wms = ServiceURL().wms.fws

    wms_111 = WMS(
        url_wms, layers="0", outformat="image/tiff", crs=DEF_CRS, version="1.1.1", validation=False
    )
    r_dict_111 = wms_111.getmap_bybox(geometry_nat.bounds, 20, DEF_CRS)
    wms = WMS(url_wms, layers="0", outformat="image/tiff", crs=DEF_CRS)
    print(wms)
    r_dict = wms.getmap_bybox(geometry_nat.bounds, 20, DEF_CRS)
    assert (
        wms_111.get_validlayers()["0"] == "Wetlands_Raster"
        and sys.getsizeof(r_dict_111["0_dd_0_0"]) == 12536763
        and sys.getsizeof(r_dict["0_dd_0_0"]) == 12536763
    )


@pytest.mark.flaky(max_runs=3)
def test_wfsbyid(wfs):
    print(wfs)
    st = wfs.getfeature_byid("staid", "01031500")
    assert st.json()["numberMatched"] == 1


@pytest.mark.flaky(max_runs=3)
def test_wfsbygeom(wfs, geometry_urb):
    r = wfs.getfeature_bygeom(geometry_urb, geo_crs=DEF_CRS, always_xy=False)
    assert len(r.json()["features"]) == 7


@pytest.mark.flaky(max_runs=3)
def test_wfsbybox(wfs, geometry_urb):
    bbox = geometry_urb.bounds
    r = wfs.getfeature_bybox(bbox, box_crs=DEF_CRS, always_xy=True)
    assert len(r.json()["features"]) == 7


@pytest.mark.flaky(max_runs=3)
def test_wfsbyfilter(wfs):
    wb = wfs.getfeature_byfilter("staid LIKE '010315%'")
    assert len(utils.traverse_json(wb.json(), ["features", "geometry", "coordinates"])) == 2


def test_decompose(geometry_urb):
    bboxs = utils.bbox_decompose(geometry_urb.bounds, 10)
    assert bboxs[0][-1] == 2828


def test_matchcrs(geometry_urb):
    bounds = geometry_urb.bounds
    points = ((bounds[0], bounds[2]), (bounds[1], bounds[3]))
    coords = MatchCRS.coords(points, DEF_CRS, ALT_CRS)
    bbox = MatchCRS.bounds(geometry_urb.bounds, DEF_CRS, ALT_CRS)
    geom = MatchCRS.geometry(geometry_urb, DEF_CRS, ALT_CRS)
    assert (
        abs(geom.centroid.x * 1e-4 - (-362.099)) < 1e-3
        and abs(bbox[0] * 1e-4 - (-365.403)) < 1e-3
        and abs(coords[0][-1] * 1e-4 == (-287.707)) < 1e-3
    )


def test_esriquery():
    point = utils.ESRIGeomQuery((-118.72, 34.118), wkid=DEF_CRS).point()
    assert list(point.keys()) == ["geometryType", "geometry", "inSR"]


def test_ipv4():
    url = (
        "https://edcintl.cr.usgs.gov/downloads/sciweb1/shared/uswem/web/conus"
        + "/eta/modis_eta/daily/downloads/det2004003.modisSSEBopETactual.zip"
    )

    session = RetrySession()
    with session.onlyipv4():
        r = session.get(url)
        z = zipfile.ZipFile(io.BytesIO(r.content))
        fname = z.read(z.filelist[0].filename)

    assert sys.getsizeof(fname) == 4361682


@pytest.mark.flaky(max_runs=3)
def test_async(geometry_nat):
    west, south, east, north = geometry_nat.bounds
    base_url = "https://thredds.daac.ornl.gov/thredds/ncss/ornldaac/1299"
    url_binary = []
    dates_itr = [(datetime(y, 1, 1), datetime(y, 1, 31)) for y in range(2000, 2005)]
    url_binary = (
        (
            f"{base_url}/MCD13.A{s.year}.unaccum.nc4",
            {
                "var": "NDVI",
                "north": f"{north}",
                "west": f"{west}",
                "east": f"{east}",
                "south": f"{south}",
                "disableProjSubset": "on",
                "horizStride": "1",
                "time_start": s.strftime("%Y-%m-%dT%H:%M:%SZ"),
                "time_end": e.strftime("%Y-%m-%dT%H:%M:%SZ"),
                "timeStride": "1",
                "addLatLon": "true",
                "accept": "netcdf",
            },
        )
        for s, e in dates_itr
    )

    url_json = [
        (
            "https://labs.waterdata.usgs.gov/api/nldi/linked-data/comid/position",
            {
                "f": "json",
                "coords": "POINT(-68.325 45.0369)",
            },
        )
    ]

    url_text = [
        (
            "https://waterservices.usgs.gov/nwis/site/",
            {"format": "rdb", "sites": "01646500", "siteStatus": "all"},
        )
    ]

    r_b = pygeoogc.async_requests(url_binary, "binary")
    r_j = pygeoogc.async_requests(url_json, "json")
    r_t = pygeoogc.async_requests(url_text, "text")

    assert (
        sys.getsizeof(r_b[0]) == 986161
        and r_j[0]["features"][0]["properties"]["identifier"] == "2675320"
        and r_t[0].split("\n")[-2].split("\t")[1] == "01646500"
    )


def test_urls():
    assert len(ServiceURL().__dict__["urls"]) == 4


def test_show_versions():
    f = io.StringIO()
    pygeoogc.show_versions(file=f)
    assert "INSTALLED VERSIONS" in f.getvalue()
