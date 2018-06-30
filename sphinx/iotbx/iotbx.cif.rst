=========
iotbx.cif
=========

.. contents:: Table of Contents

--------
Overview
--------

:py:mod:`iotbx.cif` is a module for the development of applications that make use of
the CIF format. Comprehensive tools are provided for input, output and
validation of CIFs, as well as for interconversion with high-level cctbx
crystallographic objects. The interface to the library is written in Python,
whilst parsing is carried out using a compiled parser, combining the
performance of a compiled language (C++) with the benefits of using an
interpreted language.

If you use :py:mod:`iotbx.cif` in your work please cite:

  R. J. Gildea, L. J. Bourhis, O. V. Dolomanov, R. W. Grosse-Kunstleve,
  H. Puschmann, P. D. Adams and J. A. K. Howard:
  *iotbx.cif: a comprehensive CIF toolbox*.
  `J. Appl. Cryst. (2011). 44, 1259-1263 <https://doi.org/10.1107/S0021889811041161>`_.

See also http://www.iucr.org/resources/cif.

Reading a CIF file
------------------

First we define a Python string containing a short extract from a CIF
representation of the quartz structure::

  quartz_as_cif = """\
  data_quartz
  _space_group_name_H-M_alt         'P 62 2 2'
  _cell_length_a                    5.01
  _cell_length_b                    5.01
  _cell_length_c                    5.47
  _cell_angle_alpha                 90
  _cell_angle_beta                  90
  _cell_angle_gamma                 120
  loop_
    _atom_site_label
    _atom_site_type_symbol
    _atom_site_fract_x
    _atom_site_fract_y
    _atom_site_fract_z
    _atom_site_U_iso_or_equiv
     Si Si 0.500 0.500 0.333 0.200
     O O 0.197 -0.197 0.833 0.200
  """

Next we import the :py:mod:`iotbx.cif` module, and extract an instance of
:py:class:`cctbx.xray.structure` from the CIF string::

  import iotbx.cif
  quartz_structure = iotbx.cif.reader(
    input_string=quartz_as_cif).build_crystal_structures()["quartz"]

The :py:class:`iotbx.cif.reader` class reads the CIF string and builds a Python
representation of the CIF. We then call the
:py:meth:`~iotbx.cif.reader.build_crystal_structures`
method which constructs an instance of :py:class:`cctbx.xray.structure` for each crystal
structure it finds in the CIF. Since a single CIF file may contain multiple
crystal structures, the :py:meth:`~iotbx.cif.reader.build_crystal_structures`
method returns a dictionary
of all the crystal structures found in the given CIF, where the keys are the
name of each data block containing a crystal structure.

We then call consecutively the :py:meth:`~cctbx.xray.structure.show_summary`
and :py:meth:`~cctbx.xray.structure.show_scatterers`
methods of the :py:class:`cctbx.xray.structure` object ::

  quartz_structure.show_summary().show_scatterers()

which produces the following output::

  Number of scatterers: 2
  At special positions: 2
  Unit cell: (5.01, 5.01, 5.47, 90, 90, 120)
  Space group: P 62 2 2 (No. 180)
  Label, Scattering, Multiplicity, Coordinates, Occupancy, Uiso, Ustar as Uiso
  Si   Si     3 ( 0.5000  0.5000  0.3333) 1.00 0.2000 [ - ]
  O    O      6 ( 0.1970 -0.1970  0.8333) 1.00 0.2000 [ - ]

Next we explore the :py:obj:`quartz_structure` object further by examining the
symmetry of each site in the structure::

  for scatterer in quartz_structure.scatterers():
    print "%s:" % scatterer.label, "%8.4f %8.4f %8.4f" % scatterer.site
    site_symmetry = quartz_structure.site_symmetry(scatterer.site)
    print "  point group type:", site_symmetry.point_group_type()
    print "  special position operator:", site_symmetry.special_op_simplified()

This will produce the output::

  Si:   0.5000   0.5000   0.3333
    point group type: 222
    special position operator: 1/2,1/2,1/3
  O:   0.1970  -0.1970   0.8333
    point group type: 2
    special position operator: x,-x,5/6

Now we would like to calculate some intensities for our :py:obj:`quartz_structure` and
output them in CIF format. The cctbx provides several pretabulated  scattering
factor tables that we can use; here we choose to use those from the
International Tables::

  quartz_structure.scattering_type_registry(table="it1992")

First we calculate the structure factors from the :py:obj:`quartz_structure`
by calling the :py:meth:`~cctbx.xray.structure.structure_factors` method of
:py:class:`cctbx.xray.structure`,
passing the argument ``d_min=2`` to indicate that we only want to calculate the
structure factors for miller indices with d-spacings not less than 2 Angstroms.
Following that we call the :py:meth:`~cctbx.miller.array.as_intensity_array`
method of :py:class:`~cctbx.miller.array`
to convert the complex structure factors to intensities::

  f_calc = quartz_structure.structure_factors(d_min=2).f_calc()
  f_calc_sq = f_calc.as_intensity_array()
  f_calc_sq.show_summary().show_array()

Finally we output the calculated intensities to a CIF file::

  f_calc_sq.as_cif_simple(
    array_type="calc", data_name="quartz", out=open("quartz.hkl", "wb"))


Python representation of a CIF file
-----------------------------------

The module :py:mod:`iotbx.cif.model` contains three important classes that
represent different levels of the CIF hierarchy. The class
:py:class:`iotbx.model.cif` is equivalent to
a full CIF file, and contains zero or more CIF data blocks, which in turn are
represented by the class :py:class:`iotbx.model.block`. These two classes
behave in a very similar way to a normal Python dictionary, and values can be
set or accessed using the familiar square bracket syntax.
A :py:class:`iotbx.model.block` object consists of a set of data names and
associated data values, which may or may not be looped items. Querying a
:py:class:`iotbx.model.block` object for a given data name will return a
simple string if the item is not looped, or a
:py:class:`scitbx.flex.std_string` array for looped values.

Starting with the :py:obj:`quartz_structure` already obtained above we will create a
cif block containing the symmetry information about the structure::

  from iotbx.cif import model

  cif = model.cif()
  cif_block = model.block()

First we add the unit cell parameters to the :py:obj:`cif_block` using the square
bracket syntax to add key-value (data name, data value) pairs to the
:py:obj:`cif_block`::

  unit_cell = quartz_structure.unit_cell()
  params = unit_cell.parameters()
  cif_block["_cell_length_a"] = params[0]
  cif_block["_cell_length_b"] = params[1]
  cif_block["_cell_length_c"] = params[2]
  cif_block["_cell_angle_alpha"] = params[3]
  cif_block["_cell_angle_beta"] = params[4]
  cif_block["_cell_angle_gamma"] = params[5]
  cif_block["_cell_volume"] = unit_cell.volume()

Now we will create a CIF loop object containing the space group symmetry
operations. First we create an instance of :py:class:`iotbx.cif.model.loop`
providing the data names (``header``) for the loop, before adding data values
one row at a time::

  space_group = quartz_structure.space_group()
  symop_loop = model.loop(header=("_space_group_symop_id",
                                  "_space_group_symop_operation_xyz"))
  for symop_id, symop in enumerate(space_group):
    symop_loop.add_row((symop_id + 1, symop.as_xyz()))

Next we add the :py:obj:`symop_loop` and other space group items to the
:py:obj:`cif_block`::

  space_group_type = quartz_structure.space_group_info().type()
  cif_block["_space_group_crystal_system"] = space_group.crystal_system().lower()
  cif_block["_space_group_IT_number"] = space_group_type.number()
  cif_block["_space_group_name_H-M_alt"] = space_group_type.lookup_symbol()
  cif_block["_space_group_name_Hall"] = space_group_type.hall_symbol()
  cif_block.add_loop(symop_loop)

Finally we add :py:obj:`cif_block` to the :py:obj:`cif` object with the
data block name "quartz" and print the :py:obj:`cif` object to the standard output::

  cif["quartz"] = cif_block
  print cif


Dictionary validation
---------------------

We shall take the :py:obj:`quartz_as_cif` from above and intentionally add some
extra CIF items that contain will be interpreted as  errors when validated
against the `core CIF dictionary`_::

  from iotbx.cif import validation

  cif_model = iotbx.cif.reader(input_string=quartz_as_cif).model()

Here we call the :py:class:`iotbx.cif.reader` class as before, and then
following that we call the :py:meth:`iotbx.cif.reader.model` method to get
direct access to the Python representation of the CIF that has been
constructed by the :py:class:`~iotbx.cif.reader` class. We can then
interact with this to add new data items::

  cif_model["quartz"]["_diffrn_radiation_probe"] = "xray"
  cif_model["quartz"]["_space_group_crystal_system"] = "Monoclinic"
  cif_model["quartz"]["_space_group_IT_number"] = "one hundred and eighty"

Next we add take the :py:obj:`symop_loop` constructed above and add another column
to the loop, before adding the loop to the "quartz" CIF block::

  symop_loop.add_column("_space_group_symop_sg_id", [1]*12)
  cif_model["quartz"].add_loop(symop_loop)

We can load a cif dictionary object using the
:py:func:`iotbx.cif.validation.smart_load_dictionary` function.
This can load a dictionary from a given file path or url, or given the name of
a dictionary it will attempt to first find a local copy of the dictionary, then
alternatively look up the dictionary in the IUCr cif dictionary register and
download it from the IUCr website::

  cif_core_dic = validation.smart_load_dictionary(name="cif_core.dic")
  cif_model.validate(cif_core_dic, show_warnings=True)

This will give the following output detailing the errors found::

  Invalid loop: missing parent for loop containing '_atom_site_type_symbol': '_atom_type_symbol' required but not present in data block
  Invalid loop: missing parent for loop containing '_space_group_symop_sg_id': '_space_group_id' required but not present in data block
  Type error for _space_group_IT_number: 'one hundred and eighty' is not interpretable as type 'numb'
  Invalid enumeration value for _diffrn_radiation_probe: 'xray' not in ('x-ray', 'neutron', 'electron', 'gamma')
  Warning: case-sensitive match failure for value 'Monoclinic' for '_space_group_crystal_system'


Creating a standalone C++ CIF parser
------------------------------------

Using the :py:mod:`ucif` module it is possible to very easily integrate a CIF parser
with minimal dependencies into any C++ software. For further details see the
`ucif documentation`_.


[`Complete example script`_]

.. _`flex.std_string`: http://cctbx.sourceforge.net/current/tour.html#scitbx-array-family-flex

.. _`core CIF dictionary`: http://www.iucr.org/resources/cif/dictionaries/cif_core

.. _`Complete example script`: http://cci.lbl.gov/cctbx_sources/iotbx/examples/iotbx_cif.py

.. _`ucif documentation`: http://cctbx.sourceforge.net/ucif

-----------------
API documentation
-----------------

iotbx.cif objects
---------------------

.. autoclass:: iotbx.cif.reader
    :members:
    :undoc-members:
    :show-inheritance:

.. autoclass:: iotbx.cif.model.cif
    :members:
    :undoc-members:
    :show-inheritance:

.. autoclass:: iotbx.cif.model.save
    :members:
    :undoc-members:
    :show-inheritance:

.. autoclass:: iotbx.cif.model.block
    :members:
    :undoc-members:
    :show-inheritance:

.. autoclass:: iotbx.cif.model.loop
    :members:
    :undoc-members:
    :show-inheritance:

.. autoclass:: iotbx.cif.crystal_symmetry_as_cif_block
    :members:
    :undoc-members:
    :show-inheritance:

.. autoclass:: iotbx.cif.xray_structure_as_cif_block
    :members:
    :undoc-members:
    :show-inheritance:

.. autoclass:: iotbx.cif.miller_arrays_as_cif_block
    :members:
    :undoc-members:
    :show-inheritance:

Accessory classes and functions
-------------------------------
