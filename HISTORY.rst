=======
History
=======

0.1.5 (2020-07-23)
------------------

- Re-wrote ``wms_bybox`` function as a class called ``WMS`` with a similar
  interface to the ``WFS`` class.
- Added support for WMS 1.3.0 and WFS 2.0.0.
- Added a custom ``Exception`` for the threading function called ``ThreadingException``.
- Add ``always_xy`` flag to ``WMS`` and ``WFS`` which is False by default. It is useful
  for cases where a web service doesn't change the axis order from tranditional
  ``xy`` to ``yx`` for versions higher than 1.3.0.

0.1.3 (2020-07-21)
------------------

- Remove unneccassary transformation of the input bbox in WFS.
- Use ``setuptools_scm`` for versioning.

0.1.2 (2020-07-16)
------------------

- Add the missing ``max_pixel`` argument to the ``wms_bybox`` function.
- Change the ``onlyIPv4`` method of ``RetrySession`` class to ``onlyipv4``
  to conform to the ``snake_case`` convention.
- Improve docsctrings.

0.1.1 (2020-07-15)
------------------

- Initial release.
