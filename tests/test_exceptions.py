"""Tests for exceptions and requests"""
import pytest

from pygeoogc import (
    InvalidInputType,
    InvalidInputValue,
    MissingInputs,
    RetrySession,
    ServerError,
    ZeroMatched,
)


def server_error():
    url = "https://hydro.nationalmap.gov/arcgis/rest/services/wbd/MapServer/1000"
    raise ServerError(url)


def test_server_error():
    with pytest.raises(ServerError):
        server_error()


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
    url = "http://somefailedurl.com"
    s = RetrySession(retries=2)
    s.get(url)


def test_get_connection_error():
    with pytest.raises(ConnectionError):
        get_connection_error()


def post_connection_error():
    url = "http://somefailedurl.com"
    s = RetrySession(retries=2)
    s.post(url)


def test_post_connection_error():
    with pytest.raises(ConnectionError):
        post_connection_error()
