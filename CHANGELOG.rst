2020.10
=======

* Added basic ``flex`` arrays for fixed width integer types (`#533 <https://github.com/cctbx/cctbx_project/pull/533>`_)

  * Signed types (``int8``, ``int16``, ``int32``, ``int64``)
  * Unsigned types (``uint8``, ``uint16``, ``uint32``, ``uint64``)
  * Additional functions may be wrapped in the future to support these types

* Improved building of downstream software with ``cctbx`` conda package

  * In some cases, the location of ``annlib`` is not found properly

2020.8
======

* First release on GitHub and conda-forge
