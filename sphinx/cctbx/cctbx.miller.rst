
.. py:module:: cctbx.miller

====================
cctbx.miller package
====================

:py:mod:`cctbx.miller` is one of the most important modules in CCTBX; it
encompasses nearly every operation performed directly on experimental data.
The core classes in cctbx.miller are the :ref:`set <the-miller-set>` and the
:ref:`array <the-miller-array>`.  The set (not to be confused with the built-in
Python type) contains the crystal symmetry, an array (type
:py:class:`cctbx.array_family.flex.miller_index`) of Miller indices (h,k,l),
and a boolean flag indicating anomalous pairs.  It does not contain actual
data, although many of its methods will return an array.  The array subclasses
the Miller set and adds a flex array containing data (and, optionally, a
flex array of experimental sigmas) and many extensions for supporting a variety
of data types.  The underlying "data" is often X-ray amplitudes or intensities,
but many other array types are also supported.

One important distinction needs to be made for developers used to working with
specific file formats or more archaic programming languages: **the Miller
arrays that you will work with do not necessarily correspond to a single
column of data in a reflection file**.  There are several major differences:

- **Friedel mates** will be stored separate if present.  This means that a
  pair of columns ``F(+)`` and ``F(-)`` from an MTZ file will become a single
  array with both ``(h,k,l)`` and ``(-h,-k,-l)`` present as distinct items.
  The same also applies to any other data type.  (Note that one consequence
  of this behavior is that the number of reflections will appear to
  double-count acentric reflections for which both Friedel mates are present.)
- For **experimental data** (amplitudes or intensities), the array will also
  store the corresponding experimental sigmas; ``array.data()`` returns the
  experimental data, while ``array.sigmas()`` returns sigmas.  In combination
  with the treatment of anomalous data, this means that a single Miller array
  can represent the combination of columns ``I(+),SIGI(+),I(-),SIGI(-)`` from
  a file.
- **Weighted map coefficients** such as ``FWT,DELFWT`` or ``2FOFCWT,PH2FOFCWT``
  will be treated as a single array with data type
  :py:class:`scitbx.array_family.flex.complex_double`.
- **Hendrickson-Lattman** phase probability coefficients are also grouped
  together, and have their own data type
  :py:class:`cctbx.array_family.flex.hendrickson_lattman`.

These conventions greatly simplify keeping track of and manipulating related
data items.
Output to various file formats will still follow the appropriate conventions.

---------------
Getting started
---------------

Miller sets (and arrays) can be created in three ways: programatically, by
reading from a file, or from a :py:class:`cctbx.xray.structure` object.  (In
practice, the latter two options almost always return an array object rather
than a set.)  Programmatic creation can be done directly, or through the
convenience method :py:func:`cctbx.miller.build_set`::

  >>> from cctbx import miller
  >>> from cctbx import crystal
  >>> ms = miller.build_set(
  ...   crystal_symmetry=crystal.symmetry(
  ...     space_group_symbol="P212121",
  ...     unit_cell=(6,6,6,90,90,90)),
  ...   anomalous_flag=False,
  ...   d_min=3.0)
  >>> print ms.size()
  7
  >>> print list(ms.indices())
  [(0, 0, 2), (0, 1, 1), (0, 2, 0), (1, 0, 1), (1, 1, 0), (1, 1, 1), (2, 0, 0)]
  >>> print ms.d_max_min()
  (4.242640687119285, 3.0)

The same set, instantiated directly::

  >>> from cctbx import miller
  >>> from cctbx import crystal
  >>> from cctbx.array_family import flex
  >>> ms = miller.set(
  ...   crystal_symmetry=crystal.symmetry(
  ...     space_group_symbol="P212121",
  ...     unit_cell=(6,6,6,90,90,90)),
  ...   anomalous_flag=False,
  ...   indices=flex.miller_index(
  ...     [(0,0,2),(0,1,1),(0,2,0),(1,0,1),(1,1,0),(1,1,1),(2,0,0)]))

From here we can retrieve a variety of information, even before we have
experimental data.  For instance, exploring systematic absences (starting
from the above example)::

  >>> from cctbx import sgtbx
  >>> point_group = ms.space_group().build_derived_point_group()
  >>> ms_base = ms.customized_copy(space_group=point_group)
  >>> ms_all = ms_base.complete_set()
  >>> print ms_all.size()
  10
  >>> sys_abs = ms_all.lone_set(other=ms_base)
  >>> print type(sys_abs)
  <class 'cctbx.miller.set'>
  >>> print list(sys_abs.indices())
  [(0, 0, 1), (0, 1, 0), (1, 0, 0)]
  >>> ms_all_p212121 = ms_all.customized_copy(
  ...   space_group_info=ms.space_group_info())
  >>> sys_abs_flags = ms_all_p212121.sys_absent_flags()
  >>> print type(sys_abs_flags)
  <class 'cctbx.miller.array'>
  >>> print list(sys_abs_flags.indices())
  [(0, 0, 1), (0, 0, 2), (0, 1, 0), (0, 1, 1), (0, 2, 0), (1, 0, 0), (1, 0, 1), (1, 1, 0), (1, 1, 1), (2, 0, 0)]
  >>> print list(sys_abs_flags.data())
  [True, False, True, False, False, True, False, False, False, False]
  >>> sys_abs = ms_all_p212121.select(sys_abs_flags.data())
  >>> print list(sys_abs.indices())
  [(0, 0, 1), (0, 1, 0), (1, 0, 0)]
  >>> not_sys_abs = ms_all_p212121.select(~sys_abs_flags.data())
  >>> print list(not_sys_abs.indices())
  [(0, 0, 2), (0, 1, 1), (0, 2, 0), (1, 0, 1), (1, 1, 0), (1, 1, 1), (2, 0, 0)]

This block of code performed the following actions:

- change the symmetry to the point group (``P 2 2 2``) corresponding to the
  original space group (``P 21 21 21``)
- generate the complete list of reflections for the new set in ``P 2 2 2``
- obtain the "lone set" of reflections missing from the original set relative
  to the new complete set; these correspond to reflections that are
  systematically absent in ``P 21 21 21`` (but not ``P 2 2 2``)
- change the symmetry for the complete set in ``P 2 2 2`` back to the original
  space group ``P 21 21 21``
- call the method ``sys_absent_flags()`` to obtain a Miller array whose data
  are a ``flex.bool`` array indicating those reflections that are
  systematically absent
- call the method ``select()`` using the resulting boolean array and its
  inverse, first to extract the set of systematic absences for
  ``P 21 21 21``, and then extract the non-absent set we started with

There are two more important details that are not immediately obvious from
the code example:

1) ``customized_copy()`` will create a new ``set`` object, but it will not
copy any underlying ``flex`` arrays (the same applies to the ``array``
class).  This means that modifications to these arrays via the new copy will
be propagated back to the original object.  If you want to avoid this
behavior, use the ``deep_copy()`` method::

  >>> ms_base = ms.customized_copy(space_group=point_group).deep_copy()

2) The comparison of sets in ``lone_set()``, and in general most other methods
that involve an ``other`` argument, will fail if the crystal symmetry is
not identical.  For instance, in the above example, if we instead tried to
call ``lone_set()`` using the original ``P 21 21 21`` set as ``other``::

  >>> sys_abs = ms_all.lone_set(other=ms)
  Traceback (most recent call last):
    File "<stdin>", line 1, in <module>
    File "/home/nat/src/cctbx_project/cctbx/miller/__init__.py", line 1047, in lone_set
      assert_is_similar_symmetry=assert_is_similar_symmetry).singles(1))
    File "/Users/nat/src/cctbx_project/cctbx/miller/__init__.py", line 1007, in match_indices
      assert self.is_similar_symmetry(other)
  AssertionError

We can prevent this if we want::

  >>> sys_abs = ms_all.lone_set(other=ms, assert_is_similar_symmetry=False)

However, you should use caution when disabling the symmetry check, as this
will also mean that comparisons between radically different crystal symmetries
(e.g. ``P 63 2 2`` versus ``I 41``) will be performed silently.

--------
File I/O
--------

Of course, if you are interested in working with actual experimental data,
additional APIs are required.  Methods for reading input files are covered in
more detail in the documentation for :py:mod:`iotbx.reflection_file_reader`
and :py:mod:`iotbx.file_reader`, but in the simplest case we can obtain
experimental data in just a couple of lines of code::

  >>> from iotbx.reflection_file_reader import any_reflection_file
  >>> hkl_in = any_reflection_file(file_name="data.sca")
  >>> miller_arrays = hkl_in.as_miller_arrays()
  >>> i_obs = miller_arrays[0]

This of course assumes that the file format includes crystal symmetry, which
is not the case for some popular formats; in these cases you will need to
obtain the symmetry information separately and pass it to
``as_miller_arrays()``.

Some of the file metadata will be preserved in the embedded ``array_info``
object; other attributes are properties of the array itself::

  >>> print i_obs.info()
  data.sca:I(+),SIGI(+),I(-),SIGI(-)
  >>> print i_obs.observation_type()
  xray.intensity
  >>> i_obs.show_summary()
  Miller array info: data.sca:I(+),SIGI(+),I(-),SIGI(-)
  Observation type: xray.intensity
  Type of data: double, size=7
  Type of sigmas: double, size=7
  Number of Miller indices: 7
  Anomalous flag: False
  Unit cell: (6.000, 6.000, 6.000, 90, 90, 90)
  Space group: P 21 21 21 (No. 19)
  <cctbx.miller.array object at 0x1071a6690>

(The final line is simply printing the Python representation of the array
itself - this is because the ``show_summary()`` method returns a reference to
``self``, which allows **chaining** of successive methods.)  Note that the
``array_info`` object returned by ``array.info()`` contains the file name
and original column labels; elsewhere, these attributes are used to select
specific arrays from multi-purpose formats such as MTZ.

From here we can quickly convert to amplitudes::

  >>> f_obs = i_obs.f_sq_as_f()
  >>> print f_obs.observation_type()
  xray.amplitude

(Note that for macromolecular data, the more sophisticated French-Wilson
treatment is recommended for dealing sensibly with weak or negative
intensities; this can be performed by calling ``array.french_wilson()``.  For
many purposes, however, the simpler and faster ``f_sq_as_f()`` will be
sufficient.)

The Miller array can also be easily output to CIF, MTZ, Scalepack (unmerged
format only), SHELX, or CNS formats, although some restrictions apply.  Some
of these methods (where the format is limited to certain data types) can
directly write to a file::

  >>> i_obs.export_as_scalepack_unmerged(file_name="data2.sca")
  >>> f = open("data.hkl", "w")
  >>> i_obs.export_as_shelx_hklf(file_object=f)
  >>> f.close()

Others require multiple steps, but this has the advantage of allowing multiple
arrays to be combined (provided that they have identical crystal symmetry)::

  >>> mtz_dataset = i_obs.as_mtz_dataset(column_root_label="I")
  >>> mtz_dataset.add_miller_array(f_obs, column_root_label="F")
  >>> mtz_dataset.add_miller_array(r_free_flags,
  ...   column_root_label="FreeR_flag")
  >>> mtz_dataset.mtz_object().write("data.mtz")

In addition to conventional formats, since all of the internal types can be
serialized as Python pickles, the same applies to set and array objects.

Processing input data - practical aspects
-----------------------------------------

In practice, preparing input arrays for the various other algorithms in CCTBX
is often significantly more complicated than implied in the previous section.
Suppose for example we have an MTZ
file with these columns (output excerpted from ``phenix.mtz.dump data.mtz``)::

    label      #valid  %valid    min     max type
    H           75612 100.00% -40.00   38.00 H: index h,k,l
    K           75612 100.00%   0.00   25.00 H: index h,k,l
    L           75612 100.00%   0.00   70.00 H: index h,k,l
    I(+)        74981  99.17% -11.10 2777.70 K: I(+) or I(-)
    SIGI(+)     74981  99.17%   0.10   89.50 M: standard deviation
    I(-)        69529  91.95% -14.00 2808.50 K: I(+) or I(-)
    SIGI(-)     69529  91.95%   0.10   64.90 M: standard deviation
    FreeR_flag  75773 100.00%   0.00    1.00 I: integer

We would eventually like to be able to use these data for calculation of
R-factors (versus some hypothetical set of structure factors derived from a
model or map).  This requires several steps to ensure that any subsequent
actions will behave sensibly: we must make sure that the input data are of
the correct type, symmetry-unique, using conventional indices, and consistent
with the R-free flags; we also want to use the R-free flags to separate out
"working" and "test" arrays.  To begin with, we read in the data::

  >>> from iotbx.reflection_file_reader import any_reflection_file
  >>> hkl_in = any_reflection_file(file_name="data.mtz")
  >>> miller_arrays = hkl_in.as_miller_arrays()
  >>> i_obs = miller_arrays[0]
  >>> flags = miller_arrays[1]
  >>> f_obs = i_obs.f_sq_as_f().map_to_asu()
  >>> flags = flags.customized_copy(data=flags.data()==1).map_to_asu()
  >>> f_obs = f_obs.merge_equivalents().array()
  >>> flags = flags.merge_equivalents().array()
  >>> flags_plus_minus = flags.generate_bijvoet_mates()
  >>> f_obs, flags_plus_minus = f_obs.common_sets(other=flags_plus_minus)
  >>> f_obs_work = f_obs.select(~flags_plus_minus.data())
  >>> f_obs_free = f_obs.select(flags_plus_minus.data())

By the end of this block of code, we have ensured that the experimental
amplitudes and R-free flags have identically sized and ordered arrays, and
are suitable for comparison with any other (similarly prepared) set of data.
A number of consistency checks built in to the various set and array methods
may raise exceptions if these steps are skipped - for instance, the separate
storage of ``(h,k,l)`` and ``(-h,-k,-l)`` requires us to expand the R-free
flags to be "anomalous", and if we skip that step::

  >>> f_obs, flags = f_obs.common_sets(other=flags)
  Traceback (most recent call last):
    File "<stdin>", line 1, in <module>
    File "/home/nat/src/cctbx_project/cctbx/miller/__init__.py", line 1032, in common_sets
      assert_is_similar_symmetry=assert_is_similar_symmetry)
    File "/home/nat/src/cctbx_project/cctbx/miller/__init__.py", line 1008, in match_indices
      assert self.anomalous_flag() == other.anomalous_flag()
  AssertionError

Or if we omit the call to ``common_sets()``::

  >>> f_obs_work = f_obs.select(~flags_plus_minus.data())
  Traceback (most recent call last):
    File "<stdin>", line 1, in <module>
    File "/home/nat/src/cctbx_project/cctbx/miller/__init__.py", line 3232, in select
      i = self.indices().select(selection)
  RuntimeError: scitbx Internal Error: /home/nat/src/cctbx_project/scitbx/array_family/selections.h(44): SCITBX_ASSERT(flags.size() == self.size()) failure.

(In practice, our use of ``common_sets()`` here is less than ideal; programs
in Phenix instead check that every reflection for which we have data also has
a corresponding R-free flag - after expansion to anomalous if necessary - and
exit with an error if this is not the case.)

Note that in the above example we are making many assumptions about the
contents and order of the input file, whereas in practice MTZ and CIF formats
may be arbitrarily complex and contain multiple arrays of each type.
(Additionally, the conventions for specifying R-free flags differ between
various software suites, and ``1`` will not necessarily be the appropriate
test flag value; fortunately, it is usually possible to guess the convention
being used.)   For
actual applications (as opposed to quick scripts and development code),
the utilities available in :py:mod:`iotbx.reflection_file_utils` enable
standardized retrieval of different array types based on a combination of
automatic behavior and user input (i.e. label strings), and the general-purpose
input wrapper in :py:mod:`mmtbx.command_line` encapsulates nearly all of these 
steps.

----------------
Comparing arrays
----------------

Here is a slightly more complex example of comparing data output by a
refinement program.  The input arrays are assumed to already be merged and
in the same ASU, but normally this would be taken care of by previous routines.

::

  def compute_r_factors (fobs, fmodel, flags) :
    fmodel, fobs = fmodel.common_sets(other=fobs)
    fmodel, flags = fmodel.common_sets(other=flags)
    fc_work = fmodel.select(~(flags.data()))
    fo_work = fobs.select(~(flags.data()))
    fc_test = fmodel.select(flags.data())
    fo_test = fobs.select(flags.data())
    r_work = fo_work.r1_factor(fc_work)
    r_free = fo_test.r1_factor(fc_test)
    print "r_work = %.4f" % r_work
    print "r_free = %.4f" % r_free
    print ""
    flags.setup_binner(n_bins=20)
    fo_work.use_binning_of(flags)
    fc_work.use_binner_of(fo_work)
    fo_test.use_binning_of(fo_work)
    fc_test.use_binning_of(fo_work)
    for i_bin in fo_work.binner().range_all() :
      sel_work = fo_work.binner().selection(i_bin)
      sel_test = fo_test.binner().selection(i_bin)
      fo_work_bin = fo_work.select(sel_work)
      fc_work_bin = fc_work.select(sel_work)
      fo_test_bin = fo_test.select(sel_test)
      fc_test_bin = fc_test.select(sel_test)
      if fc_test_bin.size() == 0 : continue
      r_work_bin = fo_work_bin.r1_factor(other=fc_work_bin,
        assume_index_matching=True)
      r_free_bin = fo_test_bin.r1_factor(other=fc_test_bin,
        assume_index_matching=True)
      cc_work_bin = fo_work_bin.correlation(fc_work_bin).coefficient()
      cc_free_bin = fo_test_bin.correlation(fc_test_bin).coefficient()
      legend = flags.binner().bin_legend(i_bin, show_counts=False)
      print "%s  %8d %8d  %.4f %.4f  %.3f %.3f" % (legend, fo_work_bin.size(),
        fo_test_bin.size(), r_work_bin, r_free_bin, cc_work_bin, cc_free_bin)

(The full source code is available in
``iotbx/examples/recalculate_phenix_refine_r_factors.py``.)

------------------------------
Working with experimental data
------------------------------

TODO

-------------------
From arrays to maps
-------------------

TODO

::

  fft_map = map_coeffs.fft_map(resolution_factor=0.25)
  fft_map.apply_sigma_scaling()
  real_map = fft_map.real_map_unpadded()
  site_map_values = flex.double()
  for site in xray_structure.sites_frac() :
    rho = real_map.eight_point_interpolation(site)
    site_map_values.append(rho)

As an extreme example of the ease of scripting repetitive actions using
``cctbx.miller``, here is a complete six-line script to convert all sets of map
coefficients in an MTZ file to sigma-scaled CCP4-format map files (covering the
unit cell)::

  from iotbx.reflection_file_reader import any_reflection_file
  hkl_in = any_reflection_file(file_name="map_coeffs.mtz")
  for i_map, array in enumerate(hkl_in.as_miller_arrays()) :
    if array.is_complex_array() :
      fft_map = array.fft_map(resolution_factor=0.25).apply_sigma_scaling()
      fft_map.as_ccp4_map(file_name="map_%d.ccp4" % (i_map+1))

See the documentation for :py:mod:`cctbx.maptbx` for details of working with
map objects.

--------------
The Miller set
--------------

.. autoclass:: cctbx.miller.set
    :members:
    :undoc-members:
    :show-inheritance:

.. autofunction:: cctbx.miller.build_set

----------------
The Miller array
----------------

.. autoclass:: cctbx.miller.array
    :members:
    :undoc-members:
    :show-inheritance:

---------------
Utility classes
---------------

.. autoclass:: cctbx.miller.merge_equivalents
    :members:

.. autoclass:: cctbx.miller.fft_map
.. autoclass:: cctbx.miller.array_info
.. autoclass:: cctbx.miller.normalised_amplitudes
.. autoclass:: cctbx.miller.crystal_symmetry_is_compatible_with_symmetry_from_file
.. autoclass:: cctbx.miller.binner

Submodules
----------

.. toctree::

   cctbx.miller.display
   cctbx.miller.reindexing
