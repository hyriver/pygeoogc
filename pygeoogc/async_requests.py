"""Some utilities for PyGeoOGC."""
import asyncio
from typing import (
    Any,
    Dict,
    List,
    Generator,
    MutableMapping,
    Optional,
    Tuple,
    Union,
)

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
    url_kwargs: Tuple[Tuple[str, Optional[MutableMapping[str, Any]]], ...],
    read: str,
    request: str,
    with_payload: bool,
    with_header: bool,
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
    with_payload : bool
        Specify if requests include payloads.
    with_header : bool
        Specify if requests include headers.
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
            expire_after=24*60*60,
            allowed_methods=('GET', 'POST'),
            timeout=2.5,
        )
        kwargs = {"json_serialize": json.dumps, "cache": cache}
        client_session = ClientSession
        
    async with client_session(**kwargs) as session:
        read_method = {"binary": _request_binary, "json": _request_json, "text": _request_text}
        if read not in read_method:
            raise InvalidInputValue("read", list(read_method.keys()))

        request_method = {"GET": session.get, "POST": session.post}
        if request not in request_method:
            raise InvalidInputValue("method", list(request_method.keys()))

        paylod = {"GET": "params", "POST": "data"}
        if with_payload and with_header:
            tasks = (
                read_method[read](u, request_method[request], {paylod[request]: p, "headers": h})
                for u, p, h in url_kwargs
            )
        elif with_payload:
            tasks = (
                read_method[read](u, request_method[request], {paylod[request]: p})
                for u, p in url_kwargs
            )
        elif with_header:
            tasks = (
                read_method[read](u, request_method[request], {"headers": h})
                for u, h in url_kwargs
            )
        else:
            tasks = (
                read_method[read](u, request_method[request])
                for u in url_kwargs
            )
        return await asyncio.gather(*tasks, return_exceptions=True)  # type: ignore


def async_requests(
    url_kwargs: Generator,
    read: str,
    request: str = "GET",
    max_workers: int = 8,
    with_payload: bool = False,
    with_header: bool = False,
    cache_name: Optional[str] = None,
) -> List[Union[str, MutableMapping[str, Any], bytes]]:
    """Send async requests.

    This function is based on
    `this <https://github.com/HydrologicEngineeringCenter/data-retrieval-scripts/blob/master/qpe_async_download.py>`__
    script.

    Parameters
    ----------
    url_kwargs : list of tuples
        A list of URLs. Each element of the list can be either
    read : str
        The method for returning the request; binary, json, and text.
    request : str, optional
        The request type; GET or POST, defaults to GET.
    max_workers : int, optional
        The maximum number of async processes, defaults to 8.
    with_payload : bool, optional
        Specify if requests include payloads, defaults to False.
    with_header : bool, optional
        Specify if requests include headers, defaults to False.
    cache_name : str, optional
        Path to a folder for caching the session, default to None (no caching).
        It is recommended to use caching when you're going to make multiple
        requests with a session. It can significantly speed up the function.

    Returns
    -------
    list
        A list of responses which are not in the order of input requests.
    """
    chunked_urls = tlz.partition_all(max_workers, url_kwargs)

    results = (
        asyncio.get_event_loop().run_until_complete(_async_session(c, read, request, cache_name))
        for c in chunked_urls
    )
    return list(tlz.concat(results))
