"""Some utilities for PyGeoOGC."""
import asyncio
from typing import Any, Callable, Dict, List, MutableMapping, Optional, Tuple, Union

import aiohttp
import cytoolz as tlz
import nest_asyncio
import orjson as json

from .exceptions import InvalidInputValue

nest_asyncio.apply()


async def _request_binary(
    url: str,
    session_req: aiohttp.ClientSession,
    **kwargs: Dict[str, Optional[MutableMapping[str, Any]]],
) -> bytes:
    """Create an async request and return the response as binary.

    Parameters
    ----------
    url : str
        URL to be retrieved
    session_req : ClientSession
        A ClientSession for sending the request
    **kwargs: dict
        Arguments to be passed to requests

    Returns
    -------
    bytes
        The retrieved response as binary
    """
    async with session_req(url, **kwargs) as response:
        return await response.read()


async def _request_json(
    url: str,
    session_req: aiohttp.ClientSession,
    **kwargs: Dict[str, Optional[MutableMapping[str, Any]]],
) -> MutableMapping[str, Any]:
    """Create an async request and return the response as json.

    Parameters
    ----------
    url : str
        URL to be retrieved
    session_req : ClientSession
        A ClientSession for sending the request
    **kwargs: dict
        Arguments to be passed to requests

    Returns
    -------
    dict
        The retrieved response as json
    """
    async with session_req(url, **kwargs) as response:
        return await response.json()


async def _request_text(
    url: str,
    session_req: aiohttp.ClientSession,
    **kwargs: Dict[str, Optional[MutableMapping[str, Any]]],
) -> str:
    """Create an async request and return the response as string.

    Parameters
    ----------
    url : str
        URL to be retrieved
    session_req : ClientSession
        A ClientSession for sending the request
    **kwargs: dict
        Arguments to be passed to requests

    Returns
    -------
    dict
        The retrieved response as string
    """
    async with session_req(url, **kwargs) as response:
        return await response.text()


async def _async_session(
    url_kwargs: Tuple[Tuple[str, MutableMapping[str, Any]], ...],
    read: str,
    request: str,
    cache_name: Optional[str],
) -> Callable:
    """Create an async session for sending requests.

    Parameters
    ----------
    url_kwargs : list of tuples of urls and payloads
        A list of URLs or URLs with their payloads to be retrieved.
    read : str
        The method for returning the request; binary, json, and text.
    request : str
        The request type; GET or POST.
    cache_name : str
        Path to a folder for caching the session.
        It is recommended to use caching when you're going to make multiple
        requests with a session. It can significantly speed up the function.

    Returns
    -------
    asyncio.gather
        An async gather function
    """
    if cache_name is None:
        client_session = aiohttp.ClientSession
        kwargs = {"json_serialize": json.dumps}
    else:
        try:
            from aiohttp_client_cache import CachedSession, SQLiteBackend
        except ImportError:
            msg = "For using cache you need to install aiohttp_client_cache and aiosqlite."
            raise ImportError(msg)

        cache = SQLiteBackend(
            cache_name=cache_name,
            expire_after=24 * 60 * 60,
            allowed_methods=("GET", "POST"),
            timeout=2.5,
        )
        kwargs = {"json_serialize": json.dumps, "cache": cache}
        client_session = CachedSession

    async with client_session(**kwargs) as session:
        read_method = {"binary": _request_binary, "json": _request_json, "text": _request_text}
        if read not in read_method:
            raise InvalidInputValue("read", list(read_method.keys()))

        request_method = {"GET": session.get, "POST": session.post}
        if request not in request_method:
            raise InvalidInputValue("method", list(request_method.keys()))

        tasks = (read_method[read](u, request_method[request], **args) for u, args in url_kwargs)
        return await asyncio.gather(*tasks, return_exceptions=True)  # type: ignore


def async_requests(
    urls: Tuple[str, ...],
    read: str,
    request_args: Optional[Tuple[MutableMapping[str, Any], ...]] = None,
    request: str = "GET",
    max_workers: int = 8,
    cache_name: Optional[str] = None,
) -> List[Union[str, MutableMapping[str, Any], bytes]]:
    """Send async requests.

    This function is based on
    `this <https://github.com/HydrologicEngineeringCenter/data-retrieval-scripts/blob/master/qpe_async_download.py>`__
    script.

    Parameters
    ----------
    urls : list of str
        A list of URLs.
    read : str
        The method for returning the request; binary, json, and text.
    request_args : list of dict, optional
        A list of requests kwargs corresponding to input URLs (1 on 1 mapping), defaults to None.
        For example, [{"params": {...}, "headers": {...}}, ...].
    request : str, optional
        The request type; GET or POST, defaults to GET.
    max_workers : int, optional
        The maximum number of async processes, defaults to 8.
    cache_name : str, optional
        Path to a folder for caching the session, default to None (no caching).
        It is recommended to use caching when you're going to make multiple
        requests with a session. It can significantly speed up the function.

    Returns
    -------
    list
        A list of responses which are not in the order of input requests.
    """
    if request_args is None:
        url_kwargs = zip(urls, len(urls) * [{"headers": None}])
    else:
        url_kwargs = zip(urls, request_args)
    chunked_urls = tlz.partition_all(max_workers, url_kwargs)

    results = (
        asyncio.get_event_loop().run_until_complete(_async_session(c, read, request, cache_name))
        for c in chunked_urls
    )
    return list(tlz.concat(results))
