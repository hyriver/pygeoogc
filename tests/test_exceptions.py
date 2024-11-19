"""Tests for exceptions and requests."""

from __future__ import annotations

import pytest

import pygeoogc as ogc
from pygeoogc import ArcGISRESTful, RetrySession, ServiceURL
from pygeoogc.exceptions import (
    InputTypeError,
    InputValueError,
    MissingInputError,
    ServiceError,
    ZeroMatchedError,
)


class TestRESTException:
    wbd_url: str = ServiceURL().restful.wbd
    rest_wbd: ArcGISRESTful = ArcGISRESTful(ServiceURL().restful.wbd, 1)

    def test_rest_invalid_crs(self):
        with pytest.raises(InputTypeError) as ex:
            _ = ArcGISRESTful(f"{self.wbd_url}/1/", crs="x")
        assert "The crs argument" in str(ex.value)

    def test_rest_none_layer(self):
        with pytest.raises(MissingInputError) as ex:
            _ = ArcGISRESTful(self.wbd_url)
        assert "Either layer must be passed" in str(ex.value)

    def test_rest_invalid_layer(self):
        with pytest.raises(InputValueError) as ex:
            _ = ArcGISRESTful(self.wbd_url, 9999)
        assert "Given layer is invalid" in str(ex.value)

    def test_rest_invalid_outformat(self):
        with pytest.raises(InputValueError) as ex:
            _ = ArcGISRESTful(self.wbd_url, 1, outformat="png")
        assert "geojson" in str(ex.value)

    def test_rest_invalid_outfields(self):
        with pytest.raises(InputValueError) as ex:
            _ = ArcGISRESTful(self.wbd_url, 1, outfields="dem")
        assert "areaacres" in str(ex.value)

    def test_rest_invalid_service_url(self):
        with pytest.raises(ServiceError) as ex:
            _ = ArcGISRESTful(f"{self.wbd_url}_extra_bit", 1)
        assert "_extra_bit" in str(ex.value)

    def test_rest_zero_oid(self):
        with pytest.raises(ZeroMatchedError) as ex:
            _ = self.rest_wbd.partition_oids([])
        assert "Service returned no features" in str(ex.value)

    def test_rest_unsupported_geometry(self):
        with pytest.raises(InputTypeError) as ex:
            self.rest_wbd.oids_bygeom((1, 2, 3, 4, 5))
        assert "The geom argument should be" in str(ex.value)

    def test_rest_unsupported_spatial_rel(self):
        with pytest.raises(InputValueError) as ex:
            self.rest_wbd.oids_bygeom((-1, 1), spatial_relation="intersects")
        assert "esriSpatialRelIntersects" in str(ex.value)

    def test_rest_wrong_sql(self):
        with pytest.raises(ZeroMatchedError) as ex:
            self.rest_wbd.oids_bysql("NHDFlowline.PERMANENT_IDENTIFIER")
        assert "no features" in str(ex.value)


class TestRetrySession:
    def test_parse_err_message(self):
        url = f"{ServiceURL().restful.nwis}/dv"
        payload = {
            "format": "json",
            "sites": "01031500xx",
            "startDT": "2005-01-01",
            "endDT": "2005-01-31",
            "parameterCd": "00060",
            "statCd": "00003",
            "siteStatus": "all",
        }
        with pytest.raises(ServiceError) as ex:
            _ = RetrySession().post(url, payload)
        assert "character is not a digit" in str(ex.value)


def test_download_wrong_keyword():
    url = f"{ServiceURL().restful.nwis}/dv"
    payload = {
        "param": {
            "format": "json",
            "sites": "01031500xx",
            "startDT": "2005-01-01",
            "endDT": "2005-01-31",
            "parameterCd": "00060",
            "statCd": "00003",
            "siteStatus": "all",
        }
    }
    with pytest.raises(InputValueError) as ex:
        _ = ogc.streaming_download(url, payload)
    assert "params" in str(ex.value)


def test_download_mismatch_length():
    url = [f"{ServiceURL().restful.nwis}/dv", f"{ServiceURL().restful.nwis}/dv"]
    payload = {
        "params": {
            "format": "json",
            "sites": "01031500xx",
            "startDT": "2005-01-01",
            "endDT": "2005-01-31",
            "parameterCd": "00060",
            "statCd": "00003",
            "siteStatus": "all",
        }
    }
    with pytest.raises(InputTypeError) as ex:
        _ = ogc.streaming_download(url, payload)
    assert "list of same length" in str(ex.value)
