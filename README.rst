.. image:: https://raw.githubusercontent.com/cheginit/hydrodata/develop/docs/_static/pygeoogc_logo.png
    :target: https://github.com/cheginit/pygeoogc
    :align: center

|

.. image:: https://img.shields.io/pypi/v/pygeoogc.svg
    :target: https://pypi.python.org/pypi/pygeoogc
    :alt: PyPi

.. image:: https://codecov.io/gh/cheginit/pygeoogc/branch/develop/graph/badge.svg
    :target: https://codecov.io/gh/cheginit/pygeoogc
    :alt: CodeCov

.. image:: https://github.com/cheginit/pygeoogc/workflows/build/badge.svg
    :target: https://github.com/cheginit/pygeoogc/workflows/build/badge.svg
    :alt: Github Actions

.. image:: https://mybinder.org/badge_logo.svg
    :target: https://mybinder.org/v2/gh/cheginit/hydrodata/develop
    :alt: Binder

|

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

Features
--------

PyGeoOGC is a part of the `Hydrodata <https://github.com/cheginit/hydrodata>`__ suite that offer
interfaces to web services that are based on
`ArcGIS RESTful <https://en.wikipedia.org/wiki/Representational_state_transfer>`__,
`WMS <https://en.wikipedia.org/wiki/Web_Map_Service>`__, and
`WFS <https://en.wikipedia.org/wiki/Web_Feature_Service>`__.

You can try using PyGeoOGC without installing it by clicking on the binder badge below
the PyGeoOGC banner. A Jupyter notebook instance with Hydrodata and PyGeoOGC installed, will be
launched in yout web browser and you can start coding!

Moreover, requests functionalities can be submitted via
`issue tracker <https://github.com/cheginit/pygeoogc/issues>`__.


Installation
------------

You can install PyGeoOGC using ``pip``:

.. code-block:: console

    $ pip install pygeoogc

Quickstart
----------

We can access
`Watershed Boundary Dataset <https://hydro.nationalmap.gov/arcgis/rest/services/wbd/MapServer>`__
via RESTful service,
`3D Eleveation Program <https://www.usgs.gov/core-science-systems/ngp/3dep>`__ from WMS, and
`FEMA National Flood Hazard Layer <https://www.fema.gov/national-flood-hazard-layer-nfhl>`__
via WFS. The output fo these functions are of type ``requests.Response`` that can be converted
to ``GeoDataFrame`` or ``xarray.Dataset`` using Hydrodata.

.. code-block:: python

    from pygeoogc import ArcGISREST, WFS, wms_bybox, MatchCRS
    from hydrodata import Station, utils

    la_wshed = Station('11092450')

    wbd8 = ArcGISRESTful(base_url="https://hydro.nationalmap.gov/arcgis/rest/services/wbd/MapServer/4")
    wbd8.get_featureids(la_wshed.geometry)
    resp = wbd8.get_features()
    _huc8 = utils.json_togeodf(resp[0])
    huc8 = _huc8.append([utils.json_togeodf(r) for r in resp[1:]])

    url_wms = "https://www.fws.gov/wetlands/arcgis/services/Wetlands_Raster/ImageServer/WMSServer"
    layer = "0"
    r_dict = wms_bybox(
        url_wms,
        layer,
        la_wshed.geometry.bounds,
        1e3,
        "image/tiff",
        box_crs="epsg:4326",
        crs="epsg:3857",
    )
    geom = MatchCRS.geometry(la_wshed.geometry, "epsg:4326", "epsg:3857")
    wetlands = utils.wms_toxarray(r_dict, geom, "epsg:3857")

    url_wfs = "https://hazards.fema.gov/gis/nfhl/services/public/NFHL/MapServer/WFSServer"
    wfs = WFS(
        url_wfs,
        layer="public_NFHL:Base_Flood_Elevations",
        outformat="esrigeojson",
        crs="epsg:4269",
    )
    r = wfs.getfeature_bybox(la_wshed.geometry.bounds, box_crs="epsg:4326")
    flood = utils.json_togeodf(r.json(), "epsg:4269", "epsg:4326")

Contributing
------------

Contirbutions are very welcomed. Here's how to set up PyGeoOGC for local development.

1. Fork the PyGeoOGC repo through the GitHub website.
2. Clone your fork locally and add the main PyGeoOGC as the upstream remote:

.. code-block:: console

    $ git clone git@github.com:your_name_here/pygeoogc.git
    $ git remote add upstream git@github.com:cheginit/pygeoogc.git

3. Install your local copy into a virtualenv. Assuming you have Conda installed, this is how you
   can set up your fork for local development:

.. code-block:: console

    $ cd pygeoogc/
    $ conda env create -f ci/requirements/environment.yml
    $ conda activate pygeoogc-dev
    $ python -m pip install . --no-deps

4. Check out the ``develop`` branch and create a branch for local development:

.. code-block:: console

    $ git checkout develop
    $ git checkout -b bugfix-or-feature/name-of-your-bugfix-or-feature
    $ git push

5. Before you first commit, pre-commit hooks needs to be setup:

.. code-block:: console

    $ pre-commit install
    $ pre-commit run --all-files

6. Now you can make your changes locally, make sure to add a description of the changes to
   ``HISTORY.rst`` file and add extra tests, if applicable, to ``tests`` folder. Please
   make sure to add your name at the end of the item(s) you added to the history file like this
   ``By `Taher Chegini <https://github.com/cheginit>`_``.
   Afterwards, fetch the latest updates from the remote and resolve any merge conflicts:

.. code-block:: console

    $ git fetch upstream
    $ git merge upstream/develop

7. Then lint and test the code:

.. code-block:: console

    $ make clean
    $ make lint
    $ make install
    $ make coverage

8. If you are making breaking changes make sure to reflect them in ``docs/usage.ipynb`` and
   ``docs/quickguide.ipynb`` notebooks if necessary.

9. Commit your changes and push your branch to GitHub:

.. code-block:: console

    $ git add .
    $ git commit -m "Your detailed description of your changes."
    $ git push origin name-of-your-bugfix-or-feature

10. Submit a pull request through the GitHub website. The pull request should work for
    Python 3.6, 3.7 and 3.8. Check https://github.com/cheginit/pygeoogc/actions
    and make sure that the tests pass.
