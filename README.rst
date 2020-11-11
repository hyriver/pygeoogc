.. image:: https://raw.githubusercontent.com/cheginit/hydrodata/master/docs/_static/pygeoogc_logo.png
    :target: https://github.com/cheginit/pygeoogc
    :align: center

|

.. |hydrodata| image:: https://github.com/cheginit/hydrodata/workflows/build/badge.svg
    :target: https://github.com/cheginit/hydrodata/actions?query=workflow%3Abuild
    :alt: Github Actions

.. |pygeoogc| image:: https://github.com/cheginit/pygeoogc/workflows/build/badge.svg
    :target: https://github.com/cheginit/pygeoogc/actions?query=workflow%3Abuild
    :alt: Github Actions

.. |pygeoutils| image:: https://github.com/cheginit/pygeoutils/workflows/build/badge.svg
    :target: https://github.com/cheginit/pygeoutils/actions?query=workflow%3Abuild
    :alt: Github Actions

.. |pynhd| image:: https://github.com/cheginit/pynhd/workflows/build/badge.svg
    :target: https://github.com/cheginit/pynhd/actions?query=workflow%3Abuild
    :alt: Github Actions

.. |py3dep| image:: https://github.com/cheginit/py3dep/workflows/build/badge.svg
    :target: https://github.com/cheginit/py3dep/actions?query=workflow%3Abuild
    :alt: Github Actions

.. |pydaymet| image:: https://github.com/cheginit/pydaymet/workflows/build/badge.svg
    :target: https://github.com/cheginit/pydaymet/actions?query=workflow%3Abuild
    :alt: Github Actions

=========== ==================================================================== ============
Package     Description                                                          Status
=========== ==================================================================== ============
Hydrodata_  Access NWIS, HCDN 2009, NLCD, and SSEBop databases                   |hydrodata|
PyGeoOGC_   Send queries to any ArcGIS RESTful-, WMS-, and WFS-based services    |pygeoogc|
PyGeoUtils_ Convert responses from PyGeoOGC's supported web services to datasets |pygeoutils|
PyNHD_      Navigate and subset NHDPlus (MR and HR) using web services           |pynhd|
Py3DEP_     Access topographic data through National Map's 3DEP web service      |py3dep|
PyDaymet_   Access Daymet for daily climate data both single pixel and gridded   |pydaymet|
=========== ==================================================================== ============

.. _Hydrodata: https://github.com/cheginit/hydrodata
.. _PyGeoOGC: https://github.com/cheginit/pygeoogc
.. _PyGeoUtils: https://github.com/cheginit/pygeoutils
.. _PyNHD: https://github.com/cheginit/pynhd
.. _Py3DEP: https://github.com/cheginit/py3dep
.. _PyDaymet: https://github.com/cheginit/pydaymet

PyGeoOGC: Query ArcGIS RESTful, WMS, and WFS
--------------------------------------------

.. image:: https://img.shields.io/pypi/v/pygeoogc.svg
    :target: https://pypi.python.org/pypi/pygeoogc
    :alt: PyPi

.. image:: https://img.shields.io/conda/vn/conda-forge/pygeoogc.svg
    :target: https://anaconda.org/conda-forge/pygeoogc
    :alt: Conda Version

.. image:: https://codecov.io/gh/cheginit/pygeoogc/branch/master/graph/badge.svg
    :target: https://codecov.io/gh/cheginit/pygeoogc
    :alt: CodeCov

.. image:: https://mybinder.org/badge_logo.svg
    :target: https://mybinder.org/v2/gh/cheginit/hydrodata/master?filepath=docs%2Fexamples
    :alt: Binder

|

.. image:: https://img.shields.io/badge/security-bandit-green.svg
    :target: https://github.com/PyCQA/bandit
    :alt: Security Status

.. image:: https://www.codefactor.io/repository/github/cheginit/pygeoogc/badge
   :target: https://www.codefactor.io/repository/github/cheginit/pygeoogc
   :alt: CodeFactor

.. image:: https://img.shields.io/badge/code%20style-black-000000.svg
    :target: https://github.com/psf/black
    :alt: black

.. image:: https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white
    :target: https://github.com/pre-commit/pre-commit
    :alt: pre-commit

|

🚨 **This package is under heavy development and breaking changes are likely to happen.** 🚨

Features
--------

PyGeoOGC is a part of Hydrodata software stack and provides interfaces to web services
that are based on
`ArcGIS RESTful <https://en.wikipedia.org/wiki/Representational_state_transfer>`__,
`WMS <https://en.wikipedia.org/wiki/Web_Map_Service>`__, and
`WFS <https://en.wikipedia.org/wiki/Web_Feature_Service>`__. It is noted that although
all these web service have limits on the number of objects per requests (e.g., 1000
objectIDs for RESTful and 8 million pixels for WMS), PyGeoOGC divides the requests into
smaller chunks under-the-hood and then merges the returned responses.

There is also an inventory of URLs for some of these web services in form of a class called
``ServiceURL``. These URLs are in three categories: ``ServiceURL().restful``,
``ServiceURL().wms``, and ``ServiceURL().wfs``. These URLs provide you with some examples
of the services that PyGeoOGC supports. All the URLs are read from a YAML file located
`here <pygeoogc/static/urls.yml>`_. If you had success using PyGeoOGC with a web service
please consider adding its URL to this YAML file which is located at ``pygeoogc/static/urls.yml``.

There are three main classes:

* ``ArcGISRESTful``: This class can be instantiated by providing the target layer URL.
  For example, for getting Watershed Boundary Data we can use ``ServiceURL().restful.wbd``.
  By looking at the web service website
  (https://hydro.nationalmap.gov/arcgis/rest/services/wbd/MapServer) we see that there are 9
  layers; 1 for 2-digit HU (Region), 6 for 12-digit HU (Subregion), and so on. We can either
  pass the base URL or concatenate the target layer number like so
  ``f"{ServiceURL().restful.wbd}/6"``.

  If you want to change the layer you can simply set the ``layer`` property of the class.
  Afterward, we can request for the data in two steps. First, get the object IDs using
  ``oids_bygeom`` (within a geometry), ``oids_byfield`` (specific field IDs), or ``oids_bysql``
  (any valid SQL 92 WHERE clause) class methods. Second, get the actual data using ``get_features``
  class method. The returned response can be converted into a GeoDataFrame using ``json2geodf``
  function from `PyGeoOGC <https://github.com/cheginit/pygeoutils>`__ package.

* ``WMS``: Instantiation of this class requires at least 3 arguments: service URL, layer(s)
  name(s), and output format. Additionally, target CRS and the web service version can be provided.
  Upon instantiation, we could use ``getmap_bybox`` method class to get the raster data within a
  bounding box. The box can be any valid CRS and if it is different from the default EPSG:4326, it
  should be passed to the function using ``box_crs`` argumnet. The service response can be
  converted into a ``xarray.Dataset`` using ``gtiff2xarray`` function from PyGeoOGC package.

* ``WFS``: Instantiation of this class is similar to ``WMS`` and the only difference is that
  only one layer name can be passed. Upon instantiation there are three ways to get the data:

  - ``getfeature_bybox``: Get all the features within a bounding box in any valid CRS.
  - ``getfeature_byid``: Get all the features based on the IDs. Note that two arguments should be
    provided: ``featurename``, and ``featureids``. You can get a list of valid feature names using
    ``get_validnames`` class method.
  - ``getfeature_byfilter``: Get the data based on a valid
    `CQL <https://docs.geoserver.org/latest/en/user/tutorials/cql/cql_tutorial.html>`__ filter.

  You can convert the returned response to a GeoDataFrame using ``json2geodf`` function
  from PyGeoOGC package.

You can try using PyGeoOGC without installing it on you system by clicking on the binder badge
below the PyGeoOGC banner. A Jupyter notebook instance with the Hydrodata software stack
pre-installed will be launched in your web browser and you can start coding!

Moreover, requests for additional functionalities can be submitted via
`issue tracker <https://github.com/cheginit/pygeoogc/issues>`__.

Installation
------------

You can install PyGeoOGC using ``pip``:

.. code-block:: console

    $ pip install pygeoogc

Alternatively, PyGeoOGC can be installed from the ``conda-forge`` repository
using `Conda <https://docs.conda.io/en/latest/>`__:

.. code-block:: console

    $ conda install -c conda-forge pygeoogc

Quick start
-----------

We can access
`NHDPlus HR <https://edits.nationalmap.gov/arcgis/rest/services/NHDPlus_HR/NHDPlus_HR/MapServer>`__
via RESTful service,
`National Wetlands Inventory <https://www.fws.gov/wetlands/>`__ from WMS, and
`FEMA National Flood Hazard <https://www.fema.gov/national-flood-hazard-layer-nfhl>`__
via WFS. The output for these functions are of type ``requests.Response`` that
can be converted to ``GeoDataFrame`` or ``xarray.Dataset`` using
`PyGeoOGC <https://github.com/cheginit/pygeoogc>`__.

Let's start the National Map's NHDPlus HR web service. We can query the flowlines that are
within a geometry as follows:

.. code-block:: python

    from pygeoogc import ArcGISRESTful, WFS, WMS, ServiceURL
    import pygeoutils as geoutils
    from pynhd import NLDI

    basin_geom = NLDI().getfeature_byid(
        "nwissite",
        "USGS-11092450",
        basin=True
    ).geometry[0]

    hr = ArcGISRESTful(ServiceURL().restful.nhdplushr, outformat="json")
    hr.layer = 2

    hr.oids_bygeom(basin_geom, "epsg:4326")
    resp = hr.get_features()
    flowlines = geoutils.json2geodf(resp)

Note ``oids_bygeom`` has an additional argument for passing any valid SQL WHERE clause
to further filter the data on the server side.

We can also submit a query based on IDs of any valid field in the database. If the measure
property is desired you can pass ``return_m`` as ``True`` to the ``get_features`` class method:

.. code-block:: python

    hr.oids_byfield("NHDPLUSID", [5000500013223, 5000400039708, 5000500004825])
    resp = hr.get_features(return_m=True)
    flowlines = geoutils.json2geodf(resp)

Additionally, any valid SQL 92 WHERE clause can be used. For more details look
`here <https://developers.arcgis.com/rest/services-reference/query-feature-service-.htm#ESRI_SECTION2_07DD2C5127674F6A814CE6C07D39AD46>`__.

.. code-block:: python

    hr.oids_bysql("NHDPLUSID IN (5000500013223, 5000400039708, 5000500004825)")
    resp = hr.get_features()
    flowlines = geoutils.json2geodf(resp)

A WMS-based example is shown below:

.. code-block:: python

    wms = WMS(
        ServiceURL().wms.fws,
        layers="0",
        outformat="image/tiff",
        crs="epsg:3857",
    )
    r_dict = wms.getmap_bybox(
        basin_geom.bounds,
        1e3,
        box_crs="epsg:4326",
    )
    wetlands = geoutils.gtiff2xarray(r_dict, basin_geom, "epsg:4326")

Query from a WFS-based web service can be done either within a bounding box or using
any valid `CQL filter <https://docs.geoserver.org/stable/en/user/tutorials/cql/cql_tutorial.html>`__.

.. code-block:: python

    wfs = WFS(
        ServiceURL().wfs.fema,
        layer="public_NFHL:Base_Flood_Elevations",
        outformat="esrigeojson",
        crs="epsg:4269",
    )
    r = wfs.getfeature_bybox(basin_geom.bounds, box_crs="epsg:4326")
    flood = geoutils.json2geodf(r.json(), "epsg:4269", "epsg:4326")

    layer = "wmadata:huc08"
    wfs = WFS(
        ServiceURL().wfs.waterdata,
        layer=layer,
        outformat="application/json",
        version="2.0.0",
        crs="epsg:4269",
    )
    r = wfs.getfeature_byfilter(f"huc8 LIKE '13030%'")
    huc8 = geoutils.json2geodf(r.json(), "epsg:4269", "epsg:4326")


Contributing
------------

Contributions are appreciated and very welcomed. Please read
`CONTRIBUTING.rst <https://github.com/cheginit/pygeoogc/blob/master/CONTRIBUTING.rst>`__
for instructions.
