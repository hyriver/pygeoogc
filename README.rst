.. image:: https://raw.githubusercontent.com/cheginit/hydrodata/develop/docs/_static/pygeoogc_logo.png
    :target: https://github.com/cheginit/pygeoogc
    :align: center

|

.. image:: https://img.shields.io/pypi/v/pygeoogc.svg
    :target: https://pypi.python.org/pypi/pygeoogc
    :alt: PyPi

.. image:: https://img.shields.io/conda/vn/conda-forge/pygeoogc.svg
    :target: https://anaconda.org/conda-forge/pygeoogc
    :alt: Conda Version

.. image:: https://codecov.io/gh/cheginit/pygeoogc/branch/master/graph/badge.svg
    :target: https://codecov.io/gh/cheginit/pygeoogc
    :alt: CodeCov

.. image:: https://github.com/cheginit/pygeoogc/workflows/build/badge.svg
    :target: https://github.com/cheginit/pygeoogc/workflows/build/badge.svg
    :alt: Github Actions

.. image:: https://mybinder.org/badge_logo.svg
    :target: https://mybinder.org/v2/gh/cheginit/hydrodata/develop
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

PyGeoOGC is a part of `Hydrodata <https://github.com/cheginit/hydrodata>`__ software stack
and provides interfaces to web services that are based on
`ArcGIS RESTful <https://en.wikipedia.org/wiki/Representational_state_transfer>`__,
`WMS <https://en.wikipedia.org/wiki/Web_Map_Service>`__, and
`WFS <https://en.wikipedia.org/wiki/Web_Feature_Service>`__. There is also an invetory
of URLs for some of these web services in form of a class called ``ServiceURL``. These URLs
are in three categories: ``ServiceURL().resful``, ``ServiceURL().wms``, and ``ServiceURL().wfs``.

You can try using PyGeoOGC without installing it on you system by clicking on the binder badge
below the PyGeoOGC banner. A Jupyter notebook instance with the Hydrodata software stack
pre-installed will be launched in your web browser and you can start coding!

Moreover, requests for addiitonal functionalities can be submitted via
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

Quickstart
----------

We can access
`Watershed Boundary Dataset <https://hydro.nationalmap.gov/arcgis/rest/services/wbd/MapServer>`__
via RESTful service,
`National Wetlands Inventory <https://www.fws.gov/wetlands/>`__ from WMS, and
`FEMA National Flood Hazard <https://www.fema.gov/national-flood-hazard-layer-nfhl>`__
via WFS. The output for these functions are of type ``requests.Response`` that
can be converted to ``GeoDataFrame`` or ``xarray.Dataset`` using Hydrodata.

.. code-block:: python

    from pygeoogc import ArcGISRESTful, WFS, WMS, ServiceURL
    import pygeoutils as geoutils
    from hydrodata import NLDI

    basin_geom = NLDI().getfeature_byid(
        "nwissite",
        "USGS-11092450",
        basin=True
    ).geometry[0]

    wbd12 = ArcGISRESTful(f"{ServiceURL().restful.wbd}/6")
    wbd12.get_featureids(basin_geom)
    resp = wbd12.get_features()
    huc12 = geoutils.json2geodf(resp)

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

    wfs = WFS(
        ServiceURL().wfs.fema,
        layer="public_NFHL:Base_Flood_Elevations",
        outformat="esrigeojson",
        crs="epsg:4269",
    )
    r = wfs.getfeature_bybox(basin_geom.bounds, box_crs="epsg:4326")
    flood = geoutils.json2geodf(r.json(), "epsg:4269", "epsg:4326")


Contributing
------------

Contributions are very welcomed. Please read
`CONTRIBUTING.rst <https://github.com/cheginit/pygeoogc/blob/master/CONTRIBUTING.rst>`__
file for instructions.
