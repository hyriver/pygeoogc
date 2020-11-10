"""Tests for PyGeoOGC package."""
import io
import sys
import zipfile

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
    wbd8 = ArcGISRESTful(f"{ServiceURL().restful.wbd}/4")
    wbd8.n_threads = 4
    wbd8.oids_bygeom(geometry_nat.bounds)
    wbd8.oids_bygeom(geometry_nat)
    huc8_all = wbd8.get_features()
    wbd8.oids_bygeom(geometry_nat, sql_clause="areasqkm > 5000")
    huc8_large = wbd8.get_features()
    assert len(huc8_all[0]["features"]) - len(huc8_large[0]["features"]) == 2


@pytest.mark.flaky(max_runs=3)
def test_restful_bysql():
    hr = ArcGISRESTful(ServiceURL().restful.nhdplushr, outformat="json")
    hr.layer = 1
    hr.layer = 2
    hr.oids_bysql("NHDPLUSID IN (5000500013223, 5000400039708, 5000500004825)")
    resp_sql = hr.get_features(return_m=True)
    hr.oids_byfield("NHDPLUSID", [5000500013223, 5000400039708, 5000500004825])
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


def test_urls():
    urls = ServiceURL()
    assert (
        urls.restful.nwis == "https://waterservices.usgs.gov/nwis"
        and urls.restful.nldi == "https://labs.waterdata.usgs.gov/api/nldi"
        and urls.restful.daymet_point == "https://daymet.ornl.gov/single-pixel/api/data"
        and urls.restful.daymet_grid == "https://thredds.daac.ornl.gov/thredds/ncss/ornldaac/1328"
        and urls.restful.wbd == "https://hydro.nationalmap.gov/arcgis/rest/services/wbd/MapServer"
        and urls.restful.fws == "https://www.fws.gov/wetlands/arcgis/rest/services"
        and urls.restful.fema
        == "https://hazards.fema.gov/gis/nfhl/rest/services/public/NFHL/MapServer"
        and urls.wms.mrlc == "https://www.mrlc.gov/geoserver/mrlc_download/wms"
        and urls.wms.fema
        == "https://hazards.fema.gov/gis/nfhl/rest/services/public/NFHLWMS/MapServer/WMSServer"
        and urls.wms.nm_3dep
        == "https://elevation.nationalmap.gov/arcgis/services/3DEPElevation/ImageServer/WMSServer"
        and urls.wms.fws
        == "https://www.fws.gov/wetlands/arcgis/services/Wetlands_Raster/ImageServer/WMSServer"
        and urls.wfs.waterdata == "https://labs.waterdata.usgs.gov/geoserver/wmadata/ows"
        and urls.wfs.fema
        == "https://hazards.fema.gov/gis/nfhl/services/public/NFHL/MapServer/WFSServer"
        and urls.http.ssebopeta
        == "https://edcintl.cr.usgs.gov/downloads/sciweb1/shared/uswem/web/conus/eta/modis_eta/daily/downloads"
    )


def test_show_versions():
    f = io.StringIO()
    pygeoogc.show_versions(file=f)
    assert "INSTALLED VERSIONS" in f.getvalue()
