cctbx.miller package
====================

The core classes in cctbx.miller are the :ref:`set <the-miller-set>` and the
:ref:`array <the-miller-array>`.  The set (not to be confused with the built-in
Python type) contains the crystal symmetry, an array of Miller indices (h,k,l),
and a boolean flag indicating anomalous pairs.  It does not contain actual
data, although many of its methods will return an array.  The array subclasses
the Miller set and adds a flex array and many extensions for supporting a
variety
of data types.  The underlying "data" is often X-ray amplitudes or intensities,
but many other array types are also supported.

Submodules
----------

.. toctree::

   cctbx.miller.display
   cctbx.miller.reindexing

The Miller set
--------------

.. autoclass:: cctbx.miller.set
    :members:
    :undoc-members:
    :show-inheritance:

.. autofunction:: cctbx.miller.build_set

The Miller array
----------------

.. autoclass:: cctbx.miller.array
    :members:
    :undoc-members:
    :show-inheritance:

Utility classes
---------------

.. autoclass:: cctbx.miller.merge_equivalents
.. autoclass:: cctbx.miller.fft_map
.. autoclass:: cctbx.miller.array_info
.. autoclass:: cctbx.miller.normalised_amplitudes
.. autoclass:: cctbx.miller.crystal_symmetry_is_compatible_with_symmetry_from_file
.. autoclass:: cctbx.miller.binner
