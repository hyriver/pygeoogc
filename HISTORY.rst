=======
History
=======

0.13.0 (2022-04-03)
-------------------

Breaking Changes
~~~~~~~~~~~~~~~~
- Remove caching-related arguments from all functions since now they
  can be set globally via three environmental variables:

  * ``HYRIVER_CACHE_NAME``: Path to the caching SQLite database.
  * ``HYRIVER_CACHE_EXPIRE``: Expiration time for cached requests in seconds.
  * ``HYRIVER_CACHE_DISABLE``: Disable reading/writing from/to the cache file.

  You can do this like so:

.. code-block:: python

    import os

    os.environ["HYRIVER_CACHE_NAME"] = "path/to/file.sqlite"
    os.environ["HYRIVER_CACHE_EXPIRE"] = "3600"
    os.environ["HYRIVER_CACHE_DISABLE"] = "true"

Bug Fixes
~~~~~~~~~
- In ``ArcGISRESTful.oids_byfield`` convert the input ``ids`` to a
  ``list`` if a user passes a single ``id``.

Internal Changes
~~~~~~~~~~~~~~~~
- Refactor ``ServicURL`` to hard code the supported links instead of reading
  them from a file. Also, the class now is based on ``NamedTuple`` that has a
  nicer ``__repr__``.

0.12.2 (2022-01-15)
-------------------

New Features
~~~~~~~~~~~~
- Make ``validate_crs`` public that can be accessed from the ``utils`` module.
  This is useful for checking validity of user input CRS values and getting
  its string representation.
- Add ``pygeoogc.utils.valid_wms_crs`` function for getting a list of valid
  CRS values from a WMS service.
- Add 3DEP's index WFS service for querying availability of 3DEP data within a
  bounding box.

Internal Changes
~~~~~~~~~~~~~~~~
- Add type checking with ``typeguard`` and fixed typing issues raised by
  ``typeguard``.
- Refactor ``show_versions`` to ensure getting correct versions of all
  dependencies.

0.12.1 (2021-12-31)
-------------------

Internal Changes
~~~~~~~~~~~~~~~~
- Use the three new ``ar.retrieve_*`` functions instead of the old ``ar.retrieve``
  function to improve type hinting and to make the API more consistent.

0.12.0 (2021-12-27)
-------------------

New Features
~~~~~~~~~~~~
- Add a new argument to ``ArcGISRESTful`` called ``verbose`` to turn on/off all info level logs.
- Add an option to ``ArcGISRESTful.get_features`` called ``get_geometry`` to turn on/off
  requesting the data with or without geometry.
- Now, ``ArcGISRESTful`` saves the object IDs of the features that user requested but are
  not available in the database to ``./cache/failed_request_ids.txt``.
- Add a new parameter to ``ArcGISRESTful`` called ``disable_retry`` that If ``True`` in case
  there are any failed queries, no retrying attempts is done and object IDs of the failed
  requests are saved to a text file which its path can be accessed via
  ``ArcGISRESTful.client.failed_path``.
- Set response caching expiration time to never expire, for all base classes. A new argument
  has been added to all three base classes called ``expire_after`` that can be used to set
  the expiration time.
- Add a new method to all three base classes called ``clear_cache`` that clears all cached
  responses for that specific client.

Breaking Changes
~~~~~~~~~~~~~~~~
- All ``oids_by*`` methods of ``ArcGISRESTful`` class now return a list of object IDs rather
  than setting ``self.featureids``. This makes it possible to pass the outputs of the ``oids_by*``
  functions directly to the ``get_features`` method.

Internal Changes
~~~~~~~~~~~~~~~~
- Make ``ArcGISRESTful`` less cluttered by instantiating ``ArcGISRESTfulBase`` in the
  ``init`` method of ``ArcGISRESTful`` rather than inheriting from its base class.
- Explicitly set a minimum value of 1 for the maximum number of feature IDs per request
  in ``ArcGISRESTful``, i.e., ``self.max_nrecords``.
- Add all the missing types so ``mypy --strict`` passes.

0.11.7 (2021-11-09)
-------------------

Breaking Changes
~~~~~~~~~~~~~~~~
- Remove the ``onlyipv4`` method from ``RetrySession`` since it can be easily
  be achieved using ``with unittest.mock.patch("socket.has_ipv6", False):``.

Internal Changes
~~~~~~~~~~~~~~~~
- Use the ``geoms`` method for iterating over geometries to address the
  deprecation warning of ``shapely``.
- Use ``importlib-metadata`` for getting the version instead of ``pkg_resources``
  to decrease import time as discussed in this
  `issue <https://github.com/pydata/xarray/issues/5676>`__.
- Remove unnecessary dependency on ``simplejson`` and use ``ujson`` instead.


0.11.5 (2021-09-09)
-------------------

Bug Fixes
~~~~~~~~~
- Update the code to use the latest ``requsts-cache`` API.

0.11.4 (2021-08-26)
-------------------

New Features
~~~~~~~~~~~~
- Add URL for `PyGeoAPI <https://labs.waterdata.usgs.gov/api/nldi/pygeoapi>`__ service.


0.11.3 (2021-08-21)
-------------------

Internal Changes
~~~~~~~~~~~~~~~~
- Fix a bug in ``WFS.getfeature_byid`` when the number of IDs exceeds the service's
  limit by splitting large requests into multiple smaller requests.
- Add two new arguments, ``max_nrecords`` and ``read_method``, to ``WFS`` to control
  the maximum number of records per request (defaults to 1000) and specify the response
  read method (defaults to ``json``), respectively.

0.11.2 (2021-08-19)
-------------------

Internal Changes
~~~~~~~~~~~~~~~~
- Simplify the retry logic ``ArcGISRESTFul`` by making it run four times and
  making sure that the last retry is one object ID per request.

0.11.1 (2021-07-31)
-------------------

The highlight of this release is migrating to use ``AsyncRetriever`` that can improve
the network response time significantly. Another highlight is a major refactoring of
``ArcGISRESTFul`` that improves performance and reduce code complexity.

New Features
~~~~~~~~~~~~
- Add a new method to ``ArcGISRESTFul`` class for automatically retrying the failed requests.
  This private method plucks out individual features that were in a failed request with
  several features. This happens when there are some object IDs that are not available on the
  server, and they are included in the request. In these situations the request will fail, although
  there are valid object IDs in the request. This method will pluck out the valid object IDs.
- Add support for passing additional parameters to ``WMS`` requests such as ``styles``.
- Add support for WFS version 1.0.0.

Internal Changes
~~~~~~~~~~~~~~~~
- Migrate to ``AsyncRetriever`` from ``requests-cache`` for all the web services.
- Rename ``ServiceError`` to ``ServiceUnavailable`` and ``ServerError`` to ``ServiceError``
  Since it's more representative of the intended exception.
- Raise for response status in ``RetrySession`` before the try-except block so
  ``RequestsException`` can raise and its error messaged be parsed.
- Deprecate ``utils.threading`` since all threading operations are now handled by
  ``AsyncRetriever``.
- Increase test coverage.

0.11.0 (2021-06-18)
-------------------

New Features
~~~~~~~~~~~~
- Add support for requesting ``LineString`` polygon for ``ArcGISRESTful``.
- Add a new argument called ``distance`` to ``ArcGISRESTful.oids_bygeom`` for specifying the buffer
  distance from the input geometry for getting features.

Breaking Changes
~~~~~~~~~~~~~~~~
- Drop support for Python 3.6 since many of the dependencies such as ``xarray`` and ``pandas``
  have done so.
- Remove ``async_requests`` function, since it has been packaged as a new Python library called
  `AsyncRetriever <https://github.com/cheginit/async_retriever>`__.
- Refactor ``MatchCRS``. Now, it should be instantiated by providing the in and out CRSs like so:
  ``MatchCRS(in_crs, out_crs)``. Then its methods, namely, ``geometry``, ``bounds`` and ``coords``,
  can be called. These methods now have only one input, geometry.
- Change input and output types of ``MatchCRS.coords`` from tuple of lists of coordinates
  to list of ``(x, y)`` coordinates.
- ``ArcGISRESTful`` now has a new argument, ``layer``, for specifying the layer number (int). Now,
  the target layer should either be a part of ``base_url`` or be passed with ``layer`` argument.
- Move the ``spatial_relation`` argument from ``ArcGISRESTful`` class to ``oids_bygeom`` method,
  since that's where it's applicable.

Internal Changes
~~~~~~~~~~~~~~~~
- Refactor ``ArcGISRESTfulBase`` class to reduce its code complexity and make the service
  initialization logic much simpler. The class is faster since it makes fewer requests during
  the initialization process.
- Add ``pydantic`` as a new dependency that takes care of ``ArcGISRESTfulBase`` validation.
- Use persistent caching for all send/receive requests that can significantly improve the
  network response time.
- Explicitly include all the hard dependencies in ``setup.cfg``.
- Set a default value of 1000 for ``max_nrecords`` in ``ArcGISRESTfulBase``.
- Use ``dataclass`` for ``WMSBase`` and ``WFSBase`` since support for Python 3.6 is dropped.

0.10.1 (2021-03-27)
-------------------

- Add announcement regarding the new name for the software stack, HyRiver.
- Improve ``pip`` installation and release workflow.

0.10.0 (2021-03-06)
-------------------

- The first release after renaming ``hydrodata`` to ``PyGeoHydro``.
- Fix ``extent`` property of ``ArcGISRESTful`` being set to ``None`` incorrectly.
- Add ``feature types`` property to ``ArcGISRESTFul`` for getting names and IDs of types
  of features in the database.
- Replace ``cElementTree`` with ``ElementTree`` since it's been deprecated by ``defusedxml``.
- Remove dependency on ``dataclasses`` since its benefits and usage in the code was minimal.
- Speed up CI testing by using ``mamba`` and caching.
- ``ArcGISRESTFull`` now prints number of found features before attempting to retrieve them.
- User ``logging`` module for printing information.


0.9.0 (2021-02-14)
------------------

- Bump version to the same version as PyGeoHydro.
- Add support for query by point and multi-points to ``ArcGISRESTful.bygeom``.
- Add support for buffer distance to ``ArcGISRESTful.bygeom``.
- Add support for generating ESRI-based queries for points and multi-points
  to ``ESRIGeomQuery``.
- Add all the missing type annotations.
- Update the Daymet URL to version 4. You can check the release information
  `here <https://daac.ornl.gov/DAYMET/guides/Daymet_Daily_V4.html>`_
- Use ``cytoolz`` library for improving performance of some operations.
- Add ``extent`` property to ``ArcGISRESTful`` class that get the spatial extent
  of the service.
- Add URL to ``airmap`` service for getting elevation data at 30 m resolution.

0.2.3 (2020-12-19)
-------------------

- Fix ``urlib3`` deprecation warning about using ``method_whitelist``.

0.2.2 (2020-12-05)
-------------------

- Remove unused variables in ``async_requests`` and use ``max_workers``.
- Fix the ``async_requests`` issue on Windows systems.


0.2.0 (2020-12-06)
-------------------

- Added/Renamed three class methods in ``ArcGISRESTful``: ``oids_bygeom``, ``oids_byfield``,
  and ``oids_bysql``. So you can query feature within a geometry, using specific field ID(s),
  or more generally using any valid SQL 92 WHERE clause.
- Added support for query with SQL WHERE clause to ``ArcGISRESTful``.
- Changed the NLDI's URL for migrating to its new API v3.
- Added support for CQL filter to ``WFS``, credits to `Emilio <https://github.com/emiliom>`__.
- Moved all the web services URLs to a YAML file that ``ServiceURL`` class reads. It makes
  managing the new URLs easier. The file is located at ``pygeoogc/static/urls.yml``.
- Turned off threading by default for all the services since not all web services supports it.
- Added support for setting the request method, ``GET`` or ``POST``, for ``WFS.byfilter``,
  which could be useful when the filter string is long.
- Added support for asynchronous download via the function ``async_requests``.


0.1.10 (2020-08-18)
-------------------

- Improved ``bbox_decompose`` to fix the ``WMS`` issue with high resolution requests.
- Replaces ``simplejson`` with ``orjson`` to speed up JSON operations.

0.1.8 (2020-08-12)
------------------

- Removed threading for ``WMS`` due to inconsistent behavior.
- Addressed an issue with domain decomposition for ``WMS`` where width/height becomes 0.

0.1.7 (2020-08-11)
------------------

- Renamed ``vsplit_bbox`` to ``bbox_decompose``. The function now decomposes the domain
  in both directions and return squares and rectangular.

0.1.5 (2020-07-23)
------------------

- Re-wrote ``wms_bybox`` function as a class called ``WMS`` with a similar
  interface to the ``WFS`` class.
- Added support for WMS 1.3.0 and WFS 2.0.0.
- Added a custom ``Exception`` for the threading function called ``ThreadingException``.
- Add ``always_xy`` flag to ``WMS`` and ``WFS`` which is False by default. It is useful
  for cases where a web service doesn't change the axis order from the transitional
  ``xy`` to ``yx`` for versions higher than 1.3.0.

0.1.3 (2020-07-21)
------------------

- Remove unnecessary transformation of the input bbox in WFS.
- Use ``setuptools_scm`` for versioning.

0.1.2 (2020-07-16)
------------------

- Add the missing ``max_pixel`` argument to the ``wms_bybox`` function.
- Change the ``onlyIPv4`` method of ``RetrySession`` class to ``onlyipv4``
  to conform to the ``snake_case`` convention.
- Improve docstrings.

0.1.1 (2020-07-15)
------------------

- Initial release.
