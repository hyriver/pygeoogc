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