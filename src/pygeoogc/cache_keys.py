"""Functions for creating unique keys based on web request parameters.

This module is based on the ``aiohttp-client-cache`` package, which is
licensed under the MIT license. See the ``LICENSE`` file for more details.
"""

from __future__ import annotations

import hashlib
from collections.abc import Mapping, Sequence
from typing import Any, Literal, Union

from multidict import MultiDict
from url_normalize import url_normalize
from yarl import URL

RequestParams = Union[Mapping[Any, Any], Sequence[Any], str, None]
StrOrURL = Union[str, URL]

__all__ = ["create_request_key"]


def _normalize_url_params(url: StrOrURL, params: RequestParams = None) -> URL:
    """Normalize any combination of request parameter formats that aiohttp accepts."""
    if isinstance(url, str):
        url = URL(url)

    # Handle `params` argument, and combine with URL query string if it exists
    if params:
        norm_params = MultiDict(url.query)
        norm_params.extend(url.with_query(params).query)
        url = url.with_query(norm_params)

    # Apply additional normalization and convert back to URL object
    return URL(str(url_normalize(str(url))))


def _encode_dict(data: Any) -> bytes:
    if not data:
        return b""
    if isinstance(data, bytes):
        return data
    elif not isinstance(data, Mapping):
        return str(data).encode()
    item_pairs = [f"{k}={v}" for k, v in sorted((data or {}).items())]
    return "&".join(item_pairs).encode()


def create_request_key(
    method: Literal["GET", "POST", "get", "post"],
    url: StrOrURL,
    params: RequestParams = None,
    data: Mapping[str, Any] | None = None,
    json: Mapping[str, Any] | None = None,
) -> str:
    """Create a unique cache key based on request details.

    Parameters
    ----------
    method : str
        The HTTP method used in the request. Must be either "GET" or "POST".
    url : str or yarl.URL
        The URL of the request.
    params : dict or list or str or None, optional
        The query parameters of the request. Default is None.
    data : dict or None, optional
        The data of the request. Default is None.
    json : dict or None, optional
        The JSON data of the request. Default is None.

    Returns
    -------
    str
        The unique cache key based on the request details.
    """
    # Normalize and filter all relevant pieces of request data
    norm_url = _normalize_url_params(url, params)

    # Create a hash based on the normalized and filtered request
    return hashlib.sha256(
        method.upper().encode() + str(norm_url).encode() + _encode_dict(data) + _encode_dict(json)
    ).hexdigest()
