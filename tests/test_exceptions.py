"""Tests for exceptions and requests"""
from sqlite3 import OperationalError

import pytest
from pydantic import ValidationError

from pygeoogc import (
    ArcGISRESTful,
    InvalidInputType,
    InvalidInputValue,
    MissingInputs,
    RetrySession,
    ServerError,
    ServiceURL,
    ThreadingException,
    ZeroMatched,
)


class TestRESTException:
    wbd_url: str = ServiceURL().restful.wbd
    rest_wbd: ArcGISRESTful = ArcGISRESTful(ServiceURL().restful.wbd, 1)

    def test_rest_invalid_crs(self):
        with pytest.raises(InvalidInputType) as ex:
            _ = ArcGISRESTful(f"{self.wbd_url}/1/", crs="x")
        assert "The crs argument" in str(ex.value)

    def test_rest_none_layer(self):
        with pytest.raises(ValidationError) as ex:
            _ = ArcGISRESTful(self.wbd_url)
        assert "Either layer must be passed" in str(ex.value)

    def test_rest_invalid_layer(self):
        with pytest.raises(InvalidInputValue) as ex:
            _ = ArcGISRESTful(self.wbd_url, 9999)
        assert "Given layer is invalid" in str(ex.value)

    def test_rest_invalid_max_workers(self):
        with pytest.raises(InvalidInputType) as ex:
            _ = ArcGISRESTful(f"{self.wbd_url}/1/", max_workers=-1)
        assert "positive integer" in str(ex.value)

    def test_rest_invalid_outformat(self):
        with pytest.raises(InvalidInputValue) as ex:
            _ = ArcGISRESTful(self.wbd_url, 1, outformat="png")
        assert "geojson" in str(ex.value)

    def test_rest_invalid_outfields(self):
        with pytest.raises(InvalidInputValue) as ex:
            _ = ArcGISRESTful(self.wbd_url, 1, outfields="dem")
        assert "areaacres" in str(ex.value)

    def test_rest_invalid_service_url(self):
        with pytest.raises(ServerError) as ex:
            _ = ArcGISRESTful(f"{self.wbd_url}_extra_bit", 1)
        assert "_extra_bit" in str(ex.value)

    def test_rest_oid_none(self):
        with pytest.raises(ZeroMatched) as ex:
            _ = self.rest_wbd.partition_oids(None)
        assert "Query returned no" in str(ex.value)

    def test_rest_unsupported_geometry(self):
        with pytest.raises(InvalidInputType) as ex:
            self.rest_wbd.oids_bygeom({1, 2})
        assert "The geom argument should be" in str(ex.value)

    def test_rest_unsupported_spatial_rel(self):
        with pytest.raises(InvalidInputValue) as ex:
            self.rest_wbd.oids_bygeom((-1, 1), spatial_relation="intersects")
        assert "esriSpatialRelIntersects" in str(ex.value)

    def test_rest_wrong_sql(self):
        with pytest.raises(ZeroMatched) as ex:
            self.rest_wbd.oids_bysql("NHDFlowline.PERMANENT_IDENTIFIER")
        assert "Unable to complete operation" in str(ex.value)


def server_error():
    url = "https://hydro.nationalmap.gov/arcgis/rest/services/wbd/MapServer/1000"
    raise ServerError(url)


def test_server_error():
    with pytest.raises(ServerError):
        server_error()


def threading_exception():
    raise ThreadingException(4, "Failed")


def test_threading_exception():
    with pytest.raises(ThreadingException):
        threading_exception()


def zero_matched():
    raise ZeroMatched("Query returned no matched objects.")


def test_zero_matched():
    with pytest.raises(ZeroMatched):
        zero_matched()


def invalid_value():
    raise InvalidInputValue("outFormat", ["json", "geojson"])


def test_invalid_value():
    with pytest.raises(InvalidInputValue):
        invalid_value()


def invalid_type():
    raise InvalidInputType("coords", "tuple", "(lon, lat)")


def test_invalid_type():
    with pytest.raises(InvalidInputType):
        invalid_type()


def missing_input():
    raise MissingInputs("Either coords or station_id should be provided.")


def test_missing_input():
    with pytest.raises(MissingInputs):
        missing_input()


def get_connection_error():
    url = "https://somefailedurl.com"
    s = RetrySession(retries=2)
    s.get(url)


def test_get_connection_error():
    with pytest.raises((ConnectionError, OperationalError)):
        get_connection_error()


def post_connection_error():
    url = "https://somefailedurl.com"
    s = RetrySession(retries=2)
    s.post(url)


def test_post_connection_error():
    with pytest.raises((ConnectionError, OperationalError)):
        post_connection_error()
