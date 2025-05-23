=======
History
=======

0.19.4 (2025-05-23)
-------------------

New Features
~~~~~~~~~~~~
- Update the links for eHydro web services to the latest addresses.

0.19.3 (2025-03-07)
-------------------

Internal Changes
~~~~~~~~~~~~~~~~
- Use ``orjson`` instead of ``ujson`` due to the package not being
  maintained anymore. The developer of ``ujson`` raised conrcerns
  about security vulnerabilities and recommended using ``orjson``
  instead.

0.19.0 (2025-01-17)
-------------------

Internal Changes
~~~~~~~~~~~~~~~~
- Update all dependencies on HyRiver libraries to the latest versions
  and modify the code to be compatible with the latest versions of
  the libraries.

0.18.0 (2024-10-05)
-------------------

New Features
~~~~~~~~~~~~
- Update the links for NLDI and PyGeoAPI web services to the latest addresses.

Bug Fixes
~~~~~~~~~
- Fix a bug in ``WFS.getfeature_bygeom`` where if the input geometry is in
  a geographic CRS, the function fails to transform it correctly to the
  web service's CRS.

Breaking Changes
~~~~~~~~~~~~~~~~
- Drop support for Python 3.8 since its end-of-life date is October 2024.
- Remove all exceptions from the main module and raise them from the
  ``exceptions`` module. This is to declutter the public API and make
  it easier to maintain.

0.17.1 (2024-09-14)
-------------------

New Features
~~~~~~~~~~~~
- Update the links for FEMA web services to the latest addresses.
- When the CRS of a WMS cannot be parsed throw a more informative error
  regarding the service being down (:issue_hydro:`122`).

Internal Changes
~~~~~~~~~~~~~~~~
- Drop support for Python 3.8 since its end-of-life date is October 2024.

0.17.0 (2024-05-16)
-------------------

Internal Changes
~~~~~~~~~~~~~~~~
- Make ``streaming_download`` more robust when encoutring request issues by
  returning ``None`` for those links that have failed to be processed,
  instead of throwing exceptions. So, for example, if only one link fails,
  the function will return a list of paths with ``None`` for the failed link.
- In all requests by geometries, set the geometry precision to 6 decimal points
  to avoid issues with large decimal points.
- Add the ``exceptions`` module to the high-level API to declutter
  the main module. In the future, all exceptions will be raised from
  this module and not from the main module. For now, the exceptions
  are raised from both modules for backward compatibility.

0.16.2 (2024-04-24)
-------------------

Internal Changes
~~~~~~~~~~~~~~~~
- Remove the deprecated AirMap URL.

0.16.1 (2024-01-15)
-------------------

Bug Fixes
~~~~~~~~~
- ``pyproj`` uses its own env variables for SSL certification. This release
  fixes the issue with ``pyproj`` not being able to download the grid database
  when using DOI SSL certification file. This release uses
  ``pyproj.network.set_ca_bundle_path`` for setting the SSL certification file
  given by the user via ``HYRIVER_SSL_CERT`` env variable.
- Fix an issue in ``WFS.getfeature_byid`` where the ``max_nrecords`` argument
  was not being used correctly, thus resulting in large requests to fail.

Internal Changes
~~~~~~~~~~~~~~~~
- For ``ServiceURL`` class, use ``dataclass`` instead for better performance
  and consistency.

0.16.0 (2024-01-03)
-------------------

New Features
~~~~~~~~~~~~
- Add a new arg to ``WMS.getmap_bybox`` called ``tiff_dir`` for storing
  the responses from a WMS request as a GeoTIFF file on disk instead of
  keeping all responses in memory. When this arg is given the function
  return a list of paths to these files. This is useful for large requests
  where the response is too large to be kept in memory. You can create
  a VRT file from these files using ``pygeoutils.gtiff2vrt`` function.

0.15.2 (2023-09-22)
-------------------

New Features
~~~~~~~~~~~~
- Added RESTfulURLs for FEMA's National Flood Hazard Layer (NFHL) service.
  Contributed by `Fernando Aristizabal <https://github.com/fernando-aristizabal>`__.
  (:pull_ogc:`62`)
- Now, ``RetrySession`` can be used as a context manager. This is useful for
  closing the session after using it. For example:

.. code-block:: python

    from pygeoogc import RetrySession

    with RetrySession() as session:
        r = session.get("https://httpbin.org/get").json()

Internal Changes
~~~~~~~~~~~~~~~~
- Improve the example in the docstring of ``traverse_json`` function.
- Improve exception handling in the ``ArcGISRESTful`` class and return
  a more informative error message.

0.15.1 (2023-08-02)
-------------------
From release 0.15 onward, all minor versions of HyRiver packages
will be pinned. This ensures that previous minor versions of HyRiver
packages cannot be installed with later minor releases. For example,
if you have ``pygeoogc==0.14.x`` installed, you cannot install
``pygeoogc==0.15.x`` series. This is to ensure that the API is
consistent across all minor versions.

New Features
~~~~~~~~~~~~
- Add the STN Flood Event Data URL to the list of RESTfuls.
  Contributed by `Fernando Aristizabal <https://github.com/fernando-aristizabal>`_.
  (:pull_ogc:`59`)
- Add the link for the eHydro's web service.

0.15.0 (2023-05-07)
-------------------
From release 0.15 onward, all minor versions of HyRiver packages
will be pinned. This ensures that previous minor versions of HyRiver
packages cannot be installed with later minor releases. For example,
if you have ``pygeoogc==0.14.x`` installed, you cannot install
``pygeoogc==0.15.x`` series. This is to ensure that the API is
consistent across all minor versions.

New Features
~~~~~~~~~~~~
- For now, retain compatibility with ``shapely<2`` while supporting
  ``shapley>=2``.

Bug Fixes
~~~~~~~~~
- Fix an issue in ``WFS`` where the ``getfeature_bygeom`` method
  fails if the requested web service does not have ``geometry_column``
  attribute in its schema. This release addresses this issue by
  trying to find the name from other attributes in the schema.
  If it fails to find, it raises a ``ValueError``.
- Catch an edge case in ``match_crs`` function where the input is
  a list of coordinates of length 4.
- Give precedence to non-default arguments for caching related arguments
  instead of directly getting them from env variables. This is to avoid
  the case where the user sets the env variables but then passes different
  arguments to the function. In this case, the function should use the
  passed arguments instead of the env variables.

Internal Changes
~~~~~~~~~~~~~~~~
- Remove ``pyyaml`` as a dependency since it is not used anymore.

0.14.0 (2023-03-05)
-------------------

Breaking Changes
~~~~~~~~~~~~~~~~
- Bump the minimum required version of ``shapely`` to 2.0,
  and use its new API.

Internal Changes
~~~~~~~~~~~~~~~~
- Sync all minor versions of HyRiver packages to 0.14.0.

0.13.12 (2023-02-10)
--------------------

New Features
~~~~~~~~~~~~
- Make ``match_crs`` less strict in terms of the input geometry type
  being ``tuple`` or ``list`` by relying on ``shapely`` and
  ``contextlib.suppress``. So, now users can pass any combination of
  ``list`` or ``tuple`` as coordinates or bounding box.
- More robust handling of inputs and outputs in ``streaming_download``.
  Now, only if input is ``str`` the function returns a single ``Path`` object.
  Previously if there was only one URL, whether ``list`` of length one or
  ``str``, the output was a single ``Path``, which could have had unintended
  consequences.

Bug Fixes
~~~~~~~~~
- In ``WFS`` when some layers have missing schema info, the class failed
  to initialize. This release fixes this issue by ignoring layers with
  missing schema info and asks the user to pass a sort parameter instead
  of trying to automatically find a sort parameter. This fix also improves
  the performance of this function by making fewer web requests.

Internal Changes
~~~~~~~~~~~~~~~~
- Fully migrate ``setup.cfg`` and ``setup.py`` to ``pyproject.toml``.
- Convert relative imports to absolute with ``absolufy-imports``.
- Sync all patch versions of HyRiver packages to x.x.12.

0.13.10 (2023-01-08)
--------------------

Bug Fixes
~~~~~~~~~
- Remove all Python 3.9 type-annotation-style in the codebase except for
  function signatures to ensure compatibility with Python 3.8.
  (:issue_ogc:`57`, :pull_ogc:`58`). Thanks to
  `Tim Cera <https://github.com/timcera>`__ for reporting and fixing the
  issue.

Internal Changes
~~~~~~~~~~~~~~~~
- Use ``pyright`` for type checking instead of ``mypy`` since it is faster
  and more accurate. Also, fix all the type errors reported by ``pyright``.
- Improve code quality by addressing issues raised by
  `DeepSource <https://deepsource.io/gh/hyriver/pygeoogc>`__.

0.13.9 (2022-12-15)
-------------------

Bug Fixes
~~~~~~~~~
- Add the missing annotation import to the ``cache_keys`` to ensure
  Python 3.8 and 3.9 work with Python 3.10 style type hinting.

0.13.8 (2022-12-09)
-------------------

New Features
~~~~~~~~~~~~
- Add a new property to ``WFS`` class called ``schema`` that contains
  information about column names and their types for all layers. It also
  the geometry type and its name for each layer.
- Automatically determine the geometry keyword that should be passed to
  ``WFS.getfeature_bygeom`` using the new ``schema`` property of ``WFS``.
- Add support for disabling SSL verification to ``RetrySession`` via ``ssl``
  parameter.
- Add support for streaming responses to ``RetrySession`` via ``stream``
  parameter to ``get`` and ``post`` methods.
- Add support for closing the session to ``RetrySession`` via ``close``
  method.
- Add support for passing ``params``, ``data``, and ``json`` to ``RetrySession``
  via ``get`` and ``post`` methods. Previously, keyword ``payload`` was used for
  ``params`` in ``get`` and ``data`` in ``post``. Now, ``params`` and ``data``
  can also be passed as keyword arguments to these methods.
- Add a new function called ``streaming_download`` for downloading large
  files in parallel and in chunks.

Bug Fixes
~~~~~~~~~
- Fix an issue in ``WFS`` class where number of requested features
  exceeds the maximum number of features allowed by the server, but
  only a portion of the features are returned. This release addresses
  this issue by first getting only the number of features and then
  requesting the features in chunks of features IDs based on the
  maximum number of features allowed by the server.

Internal Changes
~~~~~~~~~~~~~~~~
- Drop support for WFS version 1.0.0 since it does not support paging.
- Modify the codebase based on `Refurb <https://github.com/dosisod/refurb>`__
  suggestions.


Bug Fixes
~~~~~~~~~
- Fix the warning message in ``ArcGISRESTFul`` where wrong number of missing
  feature IDs were being reported.

0.13.7 (2022-11-04)
-------------------

New Features
~~~~~~~~~~~~
- Add a new method to ``RetrySession`` for getting the request head called
  ``RetrySession.head``. This is useful for getting the headers of a request
  without having to make a full request which is useful for getting the
  ``Content-Length`` header for example, i.e., download size.

Bug Fixes
~~~~~~~~~
- Fix an issue in the decompose function, ``utils.bbox_decompose``, where the generated
  bounding boxes might overlap in some cases. A new approach has been implemented based
  on finding the number of required bounding boxes from max allowable no. of pixels and
  total requested pixels without changing the input bounding box projection. This ensures
  that the decomposed bounding boxes are not overlapping so ``xarray.open_mfdataset``
  can be used without any issues.

Internal Changes
~~~~~~~~~~~~~~~~
- In the ``utils.match_crs`` function, don't perform any projection if the source
  target CRS are the same.
- Improve type hints for CRS-related arguments of all functions by including string,
  integer, and ``pyproj.CRS`` types.
- Add a new class method to ``WMSBase`` and ``WFSBase`` classes called
  ``get_service_options`` for retrieving the available layers, output formats, and
  CRSs for a given service. Here's an example:
- Use ``pyupgrade`` package to update the type hinting annotations
  to Python 3.10 style.

.. code-block:: python

    from pygeoogc.core import WMSBase

    url = "https://elevation.nationalmap.gov/arcgis/services/3DEPElevation/ImageServer/WMSServer"
    wms = WMSBase(url, validation=False)
    wms.get_service_options()
    print(wms.available_layer)

0.13.6 (2022-08-30)
-------------------

Internal Changes
~~~~~~~~~~~~~~~~
- Add the missing PyPi classifiers for the supported Python versions.

0.13.5 (2022-08-29)
-------------------

Breaking Changes
~~~~~~~~~~~~~~~~
- Append "Error" to all exception classes for conforming to PEP-8 naming conventions.

Internal Changes
~~~~~~~~~~~~~~~~
- Bump minimum version of ``owslib`` to 0.27.2 since the ``pyproj`` incompatibility issue
  has been addressed in this issue.
- Bump minimum version of ``requests-cache`` to 0.9.6 since the ``attrs`` version issue
  has been addressed.

0.13.3 (2022-07-31)
-------------------

New Features
~~~~~~~~~~~~
- Add support for disabling persistent caching in ``RetrySession``
  via an argument and also ``HYRIVER_CACHE_DISABLE`` environmental variable.

0.13.2 (2022-06-14)
-------------------

Breaking Changes
~~~~~~~~~~~~~~~~
- Set the minimum supported version of Python to 3.8 since many of the
  dependencies such as ``xarray``, ``pandas``, ``rioxarray`` have dropped support
  for Python 3.7.
- Pin ``owslib`` to version <0.26 since version 0.26 has pinned ``pyproj`` to
  version <3.3 which is not compatible with ``rasterio`` on macOS.

Internal Changes
~~~~~~~~~~~~~~~~
- Use `micromamba <https://github.com/marketplace/actions/provision-with-micromamba>`__
  for running tests
  and use `nox <https://github.com/marketplace/actions/setup-nox>`__
  for linting in CI.

0.13.1 (2022-06-11)
-------------------

New Features
~~~~~~~~~~~~
- More robust handling of errors in ``ArcGISRESTful`` by catching ``None``
  responses. Also, use the ``POST`` method for ``ArcGISRESTful.bysql`` since
  the SQL Clause could be a long string.

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
  ``RequestsException`` can raise, and its error messaged be parsed.
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
- Use ``logging`` module for printing information.


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
