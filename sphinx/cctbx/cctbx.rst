cctbx - core crystallographic objects and functions
===================================================

The :py:mod:`cctbx` module, not to be confused with the overall CCTBX project,
contains all data types and algorithms for generic crystallographic computing,
without reference to file formats, molecule types, or use cases.  The most
important subpackages, whose methods are used throughout the higher-level
modules of CCTBX, can be summarized as follows:

- :py:mod:`cctbx.sgtbx`: space group symmetry library, including handling of
  both real-space and reciprocal-space asymmetric units
- :py:mod:`cctbx.uctbx`: unit cell handling, including fractional-to-Cartesian
  coordinate conversions
- :py:mod:`cctbx.crystal`: operations using combined space group and unit cell
  information; many other classes are derived from
  :py:class:`cctbx.crystal.symmetry`
- :py:mod:`cctbx.array_family`: extensions to :py:mod:`scitbx.array_family`,
  defining several additional array types
- :py:mod:`cctbx.miller`: classes for manipulating reflections and associated
  data (stored as flex arrays)
- :py:mod:`cctbx.xray`: high-level operations on X-ray scatterers, including
  Fourier transforms between real and reciprocal space (using direct summation
  or FFT), which are at the core of any refinement program
- :py:mod:`cctbx.maptbx`: operations on real-space maps
- :py:mod:`cctbx.adptbx`: manipulation of atomic displacement parameters (ADPs,
  or equivalent B-factors), including conversion between a variety of
  conventions

Additional modules cover methods for dealing with restraints on molecular
geometry and ADPs and a complete set of scattering tables and related chemical
information.  However, these modules and data types are less likely to be
directly by high-level code.

Note that many of the data types involved, especially the
:py:class:`cctbx.miller.array` and :py:class:`cctbx.xray.structure` objects,
can either be built programatically or populated from files.  The latter
subject is covered in the documentation for :py:mod:`iotbx`.

Several dozen examples of practical use of :py:mod:`cctbx` functionality
can be found in the ``examples`` subdirectory of the module.

Submodules
----------

.. toctree::

    cctbx.adptbx
    cctbx.adp_restraints
    cctbx.array_family
    cctbx.covariance
    cctbx.crystal
    cctbx.crystal_orientation
    cctbx.dmtbx
    cctbx.eltbx
    cctbx.euclidean_model_matching
    cctbx.french_wilson
    cctbx.geometry
    cctbx.geometry_restraints
    cctbx.macro_mol
    cctbx.maptbx
    cctbx.masks
    cctbx.math_module
    cctbx.merging
    cctbx.miller
    cctbx.multipolar
    cctbx.neutron
    cctbx.r_free_utils
    cctbx.sgtbx
    cctbx.statistics
    cctbx.symmetry_search
    cctbx.translation_search
    cctbx.uctbx
    cctbx.xray

Module contents
---------------

.. automodule:: cctbx
    :members:
    :undoc-members:
    :show-inheritance:
