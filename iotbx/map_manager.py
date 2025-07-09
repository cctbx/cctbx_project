"""
High-level manager for reading and writing 3D maps and functions to operate on maps.  This is the class to use for most map operations.
"""
from __future__ import absolute_import, division, print_function
from libtbx.utils import to_str, null_out, Sorry
from libtbx import group_args, Auto
from libtbx.test_utils import approx_equal
import sys
import io
from cctbx import miller
from iotbx.mrcfile import map_reader, write_ccp4_map
from scitbx.array_family import flex
from cctbx import maptbx
from cctbx import miller
import mmtbx.ncs.ncs
from copy import deepcopy
from scitbx.matrix import col


class map_manager(map_reader, write_ccp4_map):
  '''
   map_manager, includes map_reader and write_ccp4_map and functions to
   operate on a 3D map object.

   This class is intended to be the principal mechanism for reading
   and writing map information.  It is intended to be used by the
   iotbx.data_manager for both of these purposes.

   Use map_manager to read, write, and carry information about
   one map.  Map_manager keeps track of the origin shifts and also the
   original full unit cell and cell dimensions.  It writes out the map
   in the same place as it was read in.

   Note on wrapping:  Wrapping means that the map value outside the map
   boundaries can be obtained as the value inside the boundaries, (translated
   by some multiple of the unit cell translations.)  Normally crystallographic
   maps can be wrapped and cryo EM maps cannot.

   Wrapping should be specified on initialization if not read from a file. If
   read from a file, the value from the file labels is used if available,
   otherwise it is assumed to be wrapping = False unless specified (normal
   for a cryo-EM map. If not specified at all, it will need to be specified
   before a map_model_manager is created or the map_manager is written out.

   Map_manager also keeps track of any changes in magnification. These
   are reflected in changes in unit_cell and crystal_symmetry cell dimensions
   and angles.

   You normally create a new map_manager by initializing map_manager with a
   file name.  Then you apply the shift_origin() method and the map is
   shifted to place the origin at (0, 0, 0) and the original origin is
   recorded as self.origin_shift_grid_units.

   NOTE: do not set self.origin_shift_grid_units directly; instead use
   the method set_original_origin_and_gridding() or supply the value when
   initializing the map_manager.

   You can also create a map_manager with a map_data object (3D flex.double()
   array) along with the meta-data below.

   NOTE: MRC Maps may not represent the entire unit cell.  Normally maps that
    have an origin (corner with minimum i, j, k) that is not zero will be
    shifted at a later stage to have the origin at (0, 0, 0), along with
    any models and ncs objects (typically done with iotbx.map_and_model).
    To be able to write out a map in the same place as it was read in
    after shifting the origin and/or boxing the map, you need to keep track
    of 3 things.  These are:
    1. unit_cell_grid: grid representing one full unit cell as read in.
        Saved in map_manager as self.unit_cell_grid
    2. unit_cell_crystal_symmetry: dimensions and space group of full unit cell
        Saved in map_manager as self._unit_cell_crystal_symmetry
    3. origin_shift_grid_units: the shift in grid units to apply to the
       working map to superimpose it on the original map. When you read the
       map in this is (0, 0, 0). If you shift the map origin from (i, j, k) to
       (0, 0, 0) then the origin_shift_grid_units is (i, j, k).
         Saved in map_manager as self.origin_shift_grid_units

   Magnification (pixel size scaling) of a map: there is no general parameter
   describing magnification of an MRC map.  Changes in scaling are
   recorded in map_manager as changes in the scaling matrix/translation that
   relates grid points in a map to real-space position.

   Normal usage (NOTE: read/write should normally be done through data_manager):

     Read in a map:
       mm = map_manager('input_map.mrc')
     Summarize:
       mm.show_summary()

     Normally shift origin of map to (0, 0, 0) (you can do this here
         or you can use iotbx.map_and_model to shift models and maps together):
       mm.shift_origin()

     Get the map_data (shifted if origin was shifted above):
       map_data = mm.map_data()

     Get the crystal_symmetry of the box of data that is present:
       cs = mm.crystal_symmetry()

     Get the crystal_symmetry of the whole unit cell (even if not present):
       unit_cell_cs = mm.unit_cell_crystal_symmetry()

     Write out the map in map_data() in original location:
       mm.write_map(file_name = 'output_map.ccp4')

   CONVENTIONS
   See http://www.ccpem.ac.uk/mrc_format/mrc2014.php for MRC format
   See https://pypi.org/project/mrcfile/ for mrcfile library documentation

   Same conventions as iotbx.ccp4_map

   Default is to write maps with INTERNAL_STANDARD_ORDER of axes of [3, 2, 1]
     corresponding to columns in Z, rows in Y, sections in X to match
     flex array layout.  This can be modified by changing the values in
     output_axis_order.

   Hard-wired to convert input maps of any order to
     INTERNAL_STANDARD_ORDER = [3, 2, 1] before conversion to flex arrays
     This is not modifiable.

    INTERNAL_STANDARD_ORDER = [3, 2, 1]

  Standard limitations and associated message.
  These can be checked with: limitations = mrc.get_limitations()
    which returns a group_args object with a list of limitations and a list
    of corresponding error messages, or None if there are none
    see phenix.show_map_info for example
  These limitations can also be accessed with specific calls placed below:
   For example:
   mrc.can_be_sharpened()  returns False if "extract_unique" is present

  Map labels that are not limitations can be accessed with:
      additional_labels = mrc.get_additional_labels()

  STANDARD_LIMITATIONS_DICT = {
    "extract_unique":
     "This map is masked around unique region and not suitable for auto-sharpening.",
    "map_is_sharpened":
     "This map is auto-sharpened and not suitable for further auto-sharpening.",
    "map_is_density_modified": "This map has been density modified.",
     }


   NOTES ON ORDER OF AXES

    Phenix standard order is 3, 2, 1 (columns Z, rows Y, sections in X).
        Convert everything to this order.

    This is the order that allows direct conversion of a numpy 3D array
     with axis order (mapc, mapr, maps) to a flex array.

    For reverse = True, supply order that converts flex array to numpy 3D array
     with order (mapc, mapr, maps)

    Note that this does not mean input or output maps have to be in this order.
     It just means that before conversion of numpy to flex or vice-versa
     the array has to be in this order.

     Note that MRC standard order for input/ouput is 1, 2, 3.

     NOTE: numpy arrays indexed from 0 so this is equivalent to
      order of 2, 1, 0 in the numpy array

    NOTE:  MRC format allows data axes to be swapped using the header
      mapc mapr and maps fields. However the mrcfile library does not
      attempt to swap the axes and assigns the columns to X, rows to Y and
      sections to Z. The data array is indexed C-style, so data values can
      be accessed using mrc.data[z][y][x].

    NOTE: normal expectation is that phenix will read/write with the
      order 3, 2, 1. This means X-sections (index = 3), Y rows (index = 2),
      Z columns (index = 1). This correxponds to
       mapc (columns) =   3 or Z
       mapr (rows)    =   2 or Y
       maps (sections) =  1 or X

    In the numpy array (2, 1, 0 instead of 3, 2, 1):

    To transpose, specify i0, i1, i2 where:
        i0 = 2 means input axis 0 becomes output axis 2
        NOTE:  axes are 0, 1, 2 etc, not 1, 2, 3
      Examples:
        np.transpose(a, (0, 1, 2))  does nothing
        np.transpose(a, (1, 2, 0)) output axis 0 is input axis 1



    We want output axes to always be 2, 1, 0 and input axes for numpy array are
      (mapc-1, mapr-1, maps-1):

    For example, in typical phenix usage, the transposition is:
      i_mapc = 3    i_mapc_np = 2
      i_mapr = 2    i_mapr_np = 1
      i_maps = 1    i_maps_np = 0

   END CONVENTIONS

  '''


  def __init__(self,
     file_name = None,  # USUAL: Initialize from file and specify wrapping
     map_data = None,   # OR map_data, unit_cell_grid, unit_cell_crystal_symmetry
     unit_cell_grid = None,
     unit_cell_crystal_symmetry = None,
     origin_shift_grid_units = None, # OPTIONAL first point in map in full cell
     ncs_object = None, # OPTIONAL ncs_object with map symmetry
     wrapping = Auto,   # OPTIONAL but recommended if not read from file
     experiment_type = Auto,   # OPTIONAL can set later also
     scattering_table = Auto,   # OPTIONAL can set later also
     resolution = Auto,   # OPTIONAL can set later also
     log = None,
     ):

    '''
      Allows reading a map file or initialization with map_data

      Normally call with file_name to read map file in CCP4/MRC format.

      Alternative is initialize with map_data and metadata
       Required: specify map_data, unit_cell_grid, unit_cell_crystal_symmetry
       Optional: specify origin_shift_grid_units

      Optional in either case: supply
        ncs_object with map symmetry of full map
        wrapping (True if map repeats infinitely with repeat unit of unit cell)
        experiment_type (xray cryo_em neutron)
        scattering_table (electron n_gaussian wk1995 it1992 neutron)
        resolution (nominal resolution of map)

     Note on external origin: If input map had ORIGIN specified,
       so that the value of self.external_origin is not (0,0,0) and not None,
       and apply_external_origin_if_present is set, and shift_origin() is
       called, then:
       determine if self.external_origin is on a grid point and if so, convert
        and use negative of it as origin .
       NOTE: apply_external_origin_if_present will be ignored if
       origin_shift_grid_units is also set and is not (0,0,0).

      NOTE on "crystal_symmetry" objects
      There are two objects that are "crystal_symmetry" objects:
      A.  unit_cell_crystal_symmetry():  This is the symmetry of the
        entire unit cell. It can be any space group. The dimensions
        correspond to the dimensions of unit_cell_grid.

      B.  crystal_symmetry():  This is the symmetry of the part of the map
        that is present.  If the entire map is present, this can be any
        space group. Otherwise it is set to P 1 (no symmetry other than unity).
        The dimensions correspond to the dimensions of the map_data.all().

      NOTE: As of 2020-05-22 both map_reader and map_manager ALWAYS convert
      map_data to flex.double.

      Map_manager does not save any extra information about
      the map except the details specified in this __init__.

      After reading you can access map data with self.map_data()
        and other attributes (see class utils in ccp4_map/__init__py)
    '''

    assert (file_name is not None) or [map_data,unit_cell_grid,
        unit_cell_crystal_symmetry].count(None)==0

    assert (ncs_object is None) or isinstance(ncs_object, mmtbx.ncs.ncs.ncs)
    assert (wrapping is Auto) or isinstance(wrapping, bool)

    if origin_shift_grid_units is not None:
      origin_shift_grid_units = tuple(origin_shift_grid_units)
      assert len(origin_shift_grid_units) ==3
    else:
      origin_shift_grid_units = (0, 0, 0)

    # Initialize log filestream
    self.set_log(log)


    # NOTE: If you add anything here to be initialized, add it to the
    #  customized_copy method

    # Initialize that we don't have crystal_symmetry:
    self._crystal_symmetry = None

    # Initialize mask to be not present
    self._created_mask = None

    # Initialize that this is not a mask
    self._is_mask = False

    # Initialize that this is not a dummy map_manager
    self._is_dummy_map_manager = False

    # Initialize that there is no output_external_origin
    self.output_external_origin = None

    # Initialize program_name, limitations, labels
    self.file_name = file_name # input file (source of this manager)
    self.program_name = None  # Name of program using this manager
    self.limitations = None  # Limitations from STANDARD_LIMITATIONS_DICT
    self.labels = None  # List of labels (usually from input file) to be written

    # Initialize wrapping
    self._wrapping = None
    self._cannot_figure_out_wrapping = None

    # Initialize ncs_object
    self._ncs_object = ncs_object


    # Initialize origin shift representing position of original origin in
    #  grid units.  If map is shifted, this is updated to reflect where
    #  to place current origin to superimpose map on original map.

    # Usual initialization with a file

    if self.file_name is not None:
      self._read_map()
      # Sets self.unit_cell_grid, self._unit_cell_crystal_symmetry, self.data,
      #  self._crystal_symmetry.  Sets also self.external_origin

      # read_map does not set self.origin_shift_grid_units. Set them here:

      # Set starting values:
      self.origin_shift_grid_units = (0, 0, 0)

      # Assume this map is not wrapped unless wrapping is set or is obvious
      if isinstance(wrapping, bool):  # Take it...
        self._wrapping = wrapping
      elif self.wrapping_from_input_file() is not None:
        self._wrapping = self.wrapping_from_input_file()
      elif self.crystal_symmetry().space_group_number() > 1 and \
         self.is_full_size():  # crystal structure and full size
        self._wrapping = True
      else:
        self._wrapping = False

    else:
      '''
         Initialization with map_data object and metadata
      '''

      assert map_data and unit_cell_grid and unit_cell_crystal_symmetry
      # wrapping must be specified

      assert wrapping in [True, False]

      # Required initialization information:
      self.data = map_data
      self.unit_cell_grid = unit_cell_grid
      self._unit_cell_crystal_symmetry = unit_cell_crystal_symmetry
      self._wrapping = wrapping
      self.external_origin = (0, 0, 0)

      # Calculate values for self._crystal_symmetry
      # Must always run this method after changing
      #    self._unit_cell_crystal_symmetry  or self.unit_cell_grid
      self.set_crystal_symmetry_of_partial_map()

      # Optional initialization information
      self.origin_shift_grid_units = origin_shift_grid_units

    # Initialization steps always done:

    # Make sure map is full size if wrapping is set
    if self._wrapping:
      assert self.is_full_size()

    # make sure labels are strings
    if self.labels is not None:
      self.labels = [to_str(label, codec = 'utf8') for label in self.labels]

    # Initialize experiment type and scattering_table and set defaults
    self._experiment_type = experiment_type
    self._scattering_table = scattering_table
    self._resolution = resolution
    self._minimum_resolution = None
    self._set_up_experiment_type_and_scattering_table_and_resolution()


  # prevent pickling error in Python 3 with self.log = sys.stdout
  # unpickling is limited to restoring sys.stdout
  def __getstate__(self):
    pickle_dict = self.__dict__.copy()
    if isinstance(self.log, io.TextIOWrapper):
      pickle_dict['log'] = None
    return pickle_dict

  def __setstate__(self, pickle_dict):
    self.__dict__ = pickle_dict
    if not hasattr(self, 'log') or self.log is None:
      self.log = sys.stdout

  def __repr__(self):
    if self.is_dummy_map_manager():
      return "Dummy map_manager"
    text = "Map manager (from %s)" %(self.file_name)+\
        "\n%s, \nUnit-cell grid: %s, (present: %s), origin shift %s " %(
      str(self.unit_cell_crystal_symmetry()).replace("\n"," "),
      str(self.unit_cell_grid),
      str(self.map_data().all()),
      str(self.origin_shift_grid_units)) + "\n"+\
      "Working coordinate shift %s" %( str(self.shift_cart()))
    if self._ncs_object:
      text += "\n%s" %str(self._ncs_object)
    return text


  def set_log(self, log = sys.stdout):
    '''
       Set output log file
    '''
    if log is None:
      self.log = null_out()
    else:
      self.log = log

  def _read_map(self):
      '''
       Read map using mrcfile/__init__.py
       Sets self.unit_cell_grid, self._unit_cell_crystal_symmetry, self.data
           self._crystal_symmetry
       Does not set self.origin_shift_grid_units
       Does set self.file_name
      '''
      self._print("Reading map from %s " %(self.file_name))

      self.read_map_file(file_name = self.file_name)  # mrcfile/__init__.py

  def _print(self, m):
    if (self.log is not None) and hasattr(self.log, 'closed') and (
        not self.log.closed):
      print(m, file = self.log)

  def set_unit_cell_crystal_symmetry(self, unit_cell_crystal_symmetry):
    '''
      Specify the dimensions and space group of the full unit cell for this
      map.  This also changes the crystal_symmetry of the part that is
      present and the grid spacing. Also resets crystal_symmetry of ncs object

      Purpose is to redefine the dimensions of the map without changing values
      of the map.

      Can be used to correct the dimensions of a map where something was
      defined incorrectly.

      Can also be used to redefine the dimensions of a map after re-estimating
      the magnification of the map.

      **********************************************************************
      NOTE: Be careful using this function if the map is boxed.  If a map is
      boxed then it has a unit_cell_crystal_symmetry which describes the
      symmetry and cell dimensions of the original (full) map, and it also
      has a (different) crystal_symmetry that describes the symmetry and
      cell dimensions of the current boxed version of the map.

      This function operates on the current map to change its
      unit_cell_crystal_symmetry to the value that is supplied. This
      automatically also changes the crystal_symmetry.

      If the map is boxed, you need to supply the new value of the
      unit_cell_crystal_symmetry, NOT the new value of crystal_symmetry.
      **********************************************************************


      Does not change self.unit_cell_grid or self.origin_shift_grid_units

      Does change the result of self.shift_cart(), which is based on
        self.origin_shift_grid_units and self.unit_cell_grid
        and self.crystal_symmetry

       Fundamental parameters set:
        self._unit_cell_crystal_symmetry: dimensions of full unit cell
        self._crystal_symmetry: dimensions of part of unit cell that is present
    '''

    from cctbx import crystal
    assert isinstance(unit_cell_crystal_symmetry, crystal.symmetry)
    self._unit_cell_crystal_symmetry = unit_cell_crystal_symmetry

    # Always follow a set of _unit_cell_crystal_symmetry with this:
    self.set_crystal_symmetry_of_partial_map()

    # Now apply crystal symmetry and new shift cart to ncs object if any
    if self._ncs_object:
      self._ncs_object = self.shift_ncs_object_to_match_map_and_return_new_ncs_object(self._ncs_object)

  def set_output_external_origin(self, value):
    '''Set the value of the output external origin'''
    assert isinstance(value, tuple) or isinstance(value,list)
    self.output_external_origin = tuple(value)


  def set_original_origin_and_gridding(self,
      original_origin = None,
      gridding = None):
    '''
       Specify the location in the full unit cell grid where the origin of
       the map that is present is to be placed to match its original position.
       This is referred to here as the "original" origin, as opposed to the
       current origin of this map.

       Note that this method does not actually shift the origin of the working
       map.  It just defines where that origin is going to be placed when
       restoring the map to its original position.

       Also optionally redefine the definition of the unit cell, keeping the
       grid spacing the same.

       This allows redefining the location of the map that is present
       within the full unit cell.  It also allows redefining the
       unit cell itself.  Only use this to create a new partial map
       in a defined location.

       Previous definition of the location of the map that is present
       is discarded.

       Fundamental parameters set:
        self.origin_shift_grid_units: shift to place origin in original location
        self._unit_cell_crystal_symmetry: dimensions of full unit cell
        self.unit_cell_grid: grid units of full unit cell

       At end, recheck wrapping
    '''
    if original_origin:
      if (self.origin_shift_grid_units !=  (0, 0, 0)) or (
          not self.origin_is_zero()):
        self.shift_origin()
        self._print("Previous origin shift of %s will be discarded" %(
          str(self.origin_shift_grid_units)))

      # Set the origin
      self.origin_shift_grid_units = original_origin
      self._print("New origin shift will be %s " %(
        str(self.origin_shift_grid_units)))

    if gridding: # reset definition of full unit cell.  Keep grid spacing

       # If gridding does not match original, set space group always to P1

       current_unit_cell_parameters = self.unit_cell_crystal_symmetry(
            ).unit_cell().parameters()
       current_unit_cell_grid = self.unit_cell_grid
       new_unit_cell_parameters = []
       for a, g, new_g in zip(
          current_unit_cell_parameters[:3], current_unit_cell_grid, gridding):
         new_a = a*new_g/g
         new_unit_cell_parameters.append(new_a)

       unit_cell_parameters = \
          new_unit_cell_parameters+list(current_unit_cell_parameters[3:])

       if current_unit_cell_grid !=  gridding:
         space_group_number_use = 1
       else:
         space_group_number_use = \
            self._unit_cell_crystal_symmetry.space_group_number()
       from cctbx import crystal
       self._unit_cell_crystal_symmetry = crystal.symmetry(
          unit_cell_parameters, space_group_number_use)

       self.unit_cell_grid = gridding
       if current_unit_cell_grid !=  gridding:
         self._print ("Resetting gridding of full unit cell from %s to %s" %(
           str(current_unit_cell_grid), str(gridding)))
         self._print ("Resetting dimensions of full unit cell from %s to %s" %(
           str(current_unit_cell_parameters),
            str(new_unit_cell_parameters)))

       # Always run after setting unit_cell_grid or _unit_cell_crystal_symmetry
       # This time it should not change anything
       original_crystal_symmetry = self.crystal_symmetry()
       self.set_crystal_symmetry_of_partial_map()
       new_crystal_symmetry = self.crystal_symmetry()
       assert original_crystal_symmetry.is_similar_symmetry(
         new_crystal_symmetry)

       if not self.is_full_size():
         self.set_wrapping(False)

  def is_dummy_map_manager(self):
    ''' Is this a dummy map manager'''
    return self._is_dummy_map_manager

  def is_mask(self):
    ''' Is this a mask '''
    return self._is_mask

  def set_is_mask(self, value=True):
    ''' define if this is a mask'''
    assert isinstance(value, bool)
    self._is_mask = value

  def origin_is_zero(self):
    '''Return whether this map currently has an origin of (0,0,0)'''
    if self.map_data().origin() == (0, 0, 0):
      return True
    else:
      return False

  def shift_origin(self, desired_origin = (0, 0, 0),
     apply_external_origin_if_present = True,):
    '''
    Shift the origin of the map to desired_origin
      (normally desired_origin is (0, 0, 0), so just update
      origin_shift_grid_units

    Origin is the value of self.map_data().origin()
    origin_shift_grid_units is the shift to apply to self.map_data() to
      superimpose it on the original map.

    If you shift the origin by (i, j, k) then add -(i, j, k) to
      the current origin_shift_grid_units.

    If current origin is at (a, b, c) and
       desired origin = (d, e, f) and
       existing origin_shift_grid_units is (g, h, i):

    the shift to make is  (d, e, f) - (a, b, c)

    the new value of origin_shift_grid_units will be:
       (g, h, i)+(a, b, c)-(d, e, f)
       or new origin_shift_grid_units is: (g, h, i)- shift

    the new origin of map_data will be (d, e, f)

     Note on external origin: If input map had ORIGIN specified,
       so that the value of self.external_origin is not (0,0,0) and not None,
       and apply_external_origin_if_present is set, then:
       determine if self.external_origin is on a grid point and if so, convert
        and use negative of it as origin. Then set self.external_origin to zero.
       Does not apply if origin is already not (0,0,0).

    '''
    if(self.map_data() is None): return

    # Don't do anything and stop here if required origin is the same as in self
    alleq = all(a == b for a, b in zip(self.map_data().origin(), desired_origin))
    if alleq:
      if (not apply_external_origin_if_present):
        return
      elif tuple(self.external_origin) == (0, 0, 0):
        return

    # Figure out what the shifts are (in grid units)
    shift_info = self._get_shift_info(desired_origin = desired_origin,
      apply_external_origin_if_present = apply_external_origin_if_present)

    # Update the value of origin_shift_grid_units
    #  This is position of the origin of the new map in the full unit cell grid
    self.origin_shift_grid_units = shift_info.new_origin_shift_grid_units

    # Set external_origin to zero as it has now been used
    self.external_origin = (0, 0, 0)

    # Shift map_data if necessary
    if shift_info.shift_to_apply !=  (0, 0, 0):
      # map will start at desired_origin and have current size:
      acc = flex.grid(shift_info.desired_origin, shift_info.new_end)
      self.map_data().reshape(acc)

    # Checks
    new_current_origin = self.map_data().origin()
    assert new_current_origin == shift_info.desired_origin

    assert add_tuples_int(shift_info.current_origin, shift_info.shift_to_apply
        ) == shift_info.desired_origin

    # Original location of first element of map should agree with previous

    assert shift_info.map_corner_original_location  ==  add_tuples_int(
       new_current_origin, self.origin_shift_grid_units)

    # If there is an associated ncs_object, shift it too
    if self._ncs_object:
      self._ncs_object=self._ncs_object.coordinate_offset(
        shift_info.shift_to_apply_cart)

  def _get_shift_info(self, desired_origin = None,
    apply_external_origin_if_present = True):
    '''
      Utility to calculate the shift necessary (grid units)
      map to place the origin of the current map
      at the position defined by desired_origin.
      See definitions in shift_origin method.

    '''
    if(desired_origin is None):
      desired_origin = (0, 0, 0)
    desired_origin = tuple(desired_origin)

    if(self.origin_shift_grid_units is None):
      self.origin_shift_grid_units = (0, 0, 0)

    # Current origin and shift to apply
    current_origin = self.map_data().origin()

    self._warning_message = ""
    if apply_external_origin_if_present and \
         tuple(self.external_origin) != (0,0,0): # check for external origin
      if self.external_origin_is_compatible_with_gridding():
         external_origin_as_grid_units = self.external_origin_as_grid_units()
      else:
        external_origin_as_grid_units = (0,0,0)
        self._warning_message="External origin is not on a grid point" +\
         "...ignoring external origin" +\
         "\n***Please contact the Phenix "+\
          "developers if you need Phenix to use this external_origin***\n"
    else:
      external_origin_as_grid_units = (0,0,0)

    if self.external_origin_as_grid_units and \
        (external_origin_as_grid_units != (0,0,0)):
      if current_origin and \
          (current_origin != (0,0,0)):
        self._warning_message="Map has external origin as well as existing " +\
         "origin shift...ignoring external origin" +\
         "\n***Please contact the Phenix "+\
          "developers if you need Phenix to use this external_origin***\n"
      else:  # take it
        self._warning_message="Map has external origin " +\
         "...using external origin as origin shift after "+\
         "conversion to grid units"
        current_origin = external_origin_as_grid_units

    # Original location of first element of map
    map_corner_original_location = add_tuples_int(current_origin,
         self.origin_shift_grid_units)

    shift_to_apply = subtract_tuples_int(desired_origin, current_origin)

    assert add_tuples_int(current_origin, shift_to_apply) == desired_origin

    new_origin_shift_grid_units = subtract_tuples_int(
        self.origin_shift_grid_units, shift_to_apply)

    current_end = add_tuples_int(current_origin, self.map_data().all())
    new_end = add_tuples_int(desired_origin, self.map_data().all())
    shift_to_apply_cart = self.grid_units_to_cart(shift_to_apply)

    shift_info = group_args(
      map_corner_original_location = map_corner_original_location,
      current_origin = current_origin,
      current_end = current_end,
      current_origin_shift_grid_units = self.origin_shift_grid_units,
      shift_to_apply = shift_to_apply,
      desired_origin = desired_origin,
      new_end = new_end,
      new_origin_shift_grid_units = new_origin_shift_grid_units,
      shift_to_apply_cart = shift_to_apply_cart,
       )
    return shift_info

  def external_origin_is_compatible_with_gridding(self):
    '''Determine if external origin falls on a grid point.'''
    value = self.external_origin_as_grid_units()
    if value is not None:
      return True
    else:
      return False

  def external_origin_as_grid_units(self, as_inverse = False):
    ''' Convert external_origin to value in grid units.
        See notes on external origin. '''
    unit_cell = self.unit_cell_crystal_symmetry().unit_cell()
    unit_cell_grid = self.unit_cell_grid
    spacings = [(a/n) for a,n in zip(unit_cell.parameters()[:3],
        unit_cell_grid)]
    if self.external_origin:
      external_origin = flex.double(self.external_origin)
    else:
      external_origin = flex.double((0.,0.,0.))
    if flex.sum(flex.abs(external_origin)) > 0:
      import math
      # external_origin is the location of the external origin in xyz coords
      # origin_shift is the location in gridding space of the external_origin
      # we are hoping that it will fall on an integer grid point in this space
      # The gridding coordinate system is the unit_cell, with one unit
      #  in each direction corresponding to the spacings in that
      #   direction (e.g., a/grid_points_along_a along the x axis)

      # use unit_cell.fractionalize(origin_shift) to get fractional coords
      fractional_external_origin = unit_cell.fractionalize(col(external_origin))

      # origin_shift in grid units is fractional_external_origin multiplied
      #   by the number of grid units along the 3 axes, rounded
      origin_shift = tuple(
         [round(f * s) for f,s in zip(
            fractional_external_origin, unit_cell_grid)])

      # origin_check is position of external_origin calculated from
      #   origin_shift and the unit_cell_grid
      origin_check = unit_cell.orthogonalize(col( tuple(
        [os / a for os, a in zip(origin_shift, unit_cell_grid)])))

      origin_distance_to_grid = math.sqrt(flex.sum(
            flex.pow2(flex.double(external_origin)-flex.double(origin_check))))
      if origin_distance_to_grid > 0.001:
        return None # not compatible with grid
      else:
        if as_inverse:
          return tuple([-a for a in origin_shift])
        else: # as-is:
          return origin_shift
    else: # no shift
      return (0,0,0)

  def shift_origin_to_match_original(self):
    '''
     Shift origin by self.origin_shift_grid_units to place origin in its
     original location
    '''
    original_origin = add_tuples_int(self.map_data().origin(),
                               self.origin_shift_grid_units)

    self.shift_origin(desired_origin = original_origin)

  def set_ncs_object(self, ncs_object):
    '''
      set the ncs object for this map_manager.  Incoming ncs_object must
     be compatible (shift_cart values must match or be defined).
      Incoming ncs_object is deep_copied.
    '''
    if not ncs_object:
      return # Nothing to do
    assert isinstance(ncs_object, mmtbx.ncs.ncs.ncs)
    if (not self.is_compatible_ncs_object(ncs_object)):
      ncs_object = self.shift_ncs_object_to_match_map_and_return_new_ncs_object(ncs_object)
    self._ncs_object = deepcopy(ncs_object)

  def set_program_name(self, program_name = None):
    '''
      Set name of program doing work on this map_manager for output
      (string)
    '''
    self.program_name = program_name
    self._print("Program name of %s added" %(program_name))

  def add_limitation(self, limitation = None):
    '''
      Add a limitation from STANDARD_LIMITATIONS_DICT
      Supply the key (such as "map_is_sharpened")
    '''
    from iotbx.mrcfile import STANDARD_LIMITATIONS_DICT
    assert limitation in list(STANDARD_LIMITATIONS_DICT.keys())

    if not self.limitations:
      self.limitations = []
    self.limitations.append(limitation)
    self._print("Limitation of %s ('%s') added to map_manager" %(
      limitation, STANDARD_LIMITATIONS_DICT[limitation]))

  def remove_labels(self):
    '''
     Remove all labels
    '''
    self.labels = []

  def add_label(self, label = None, verbose = False):
    '''
     Add a label (up to 80-character string) to write to output map.
     Default is to specify the program name and date
    '''
    if not self.labels:
      self.labels = []
    if len(label)>80:  label = label[:80]
    self.labels.reverse()  # put at beginning
    self.labels.append(to_str(label, codec = 'utf8')) # make sure it is a string
    self.labels.reverse()
    if verbose:
      self._print("Label added: %s " %(label))

  def write_map(self, file_name):

    '''
      Simple version of write

      file_name is output file name
      map_data is map_data object with 3D values for map. If not supplied,
        use self.map_data()

      Normally call with file_name (file to be written)
      Output labels are generated from existing self.labels,
      self.program_name, and self.limitations

      If self.output_external_origin is specified, write that value to file


    '''

    # Make sure we have map_data
    assert self.map_data()

    map_data = self.map_data()

    assert isinstance(self.wrapping(), bool)  # need wrapping set to write file
    # remove any labels about wrapping
    for key in ["wrapping_outside_cell","no_wrapping_outside_cell"]:
      self.remove_limitation(key)
    # Add limitation on wrapping
    new_labels=[]
    if self.wrapping():
      self.add_limitation("wrapping_outside_cell")
    else:
      self.add_limitation("no_wrapping_outside_cell")


    from iotbx.mrcfile import create_output_labels
    labels = create_output_labels(
      program_name = self.program_name,
      input_file_name = self.file_name,
      input_labels = self.labels,
      limitations = self.limitations)

    crystal_symmetry = self.unit_cell_crystal_symmetry()
    unit_cell_grid = self.unit_cell_grid
    origin_shift_grid_units = self.origin_shift_grid_units

    if map_data.origin()  ==  (0, 0, 0):  # Usual
      self._print("Writing map with origin at %s and size of %s to %s" %(
        str(origin_shift_grid_units), str(map_data.all()), file_name))
      from six.moves import StringIO
      f=StringIO()
      write_ccp4_map(
        file_name   = file_name,
        crystal_symmetry = crystal_symmetry, # unit cell and space group
        map_data    = map_data,
        unit_cell_grid = unit_cell_grid,  # optional gridding of full unit cell
        origin_shift_grid_units = origin_shift_grid_units, # origin shift
        labels      = labels,
        external_origin = self.output_external_origin,


        out = f)
      self._print(f.getvalue())
    else: # map_data has not been shifted to (0, 0, 0).  Shift it and then write
          # and then shift back
      self._print("Writing map after shifting origin")
      if self.origin_shift_grid_units and origin_shift_grid_units!= (0, 0, 0):
        self._print (
          "WARNING: map_data has origin at %s " %(str(map_data.origin()))+
         " and this map_manager will apply additional origin shift of %s " %(
          str(self.origin_shift_grid_units)))

      # Save where we are
      current_origin = map_data.origin()

      # Set origin at (0, 0, 0)
      self.shift_origin(desired_origin = (0, 0, 0))
      self.write_map(file_name = file_name)
      self.shift_origin(desired_origin = current_origin)

  def create_mask_with_map_data(self, map_data):
    '''
      Set mask to be map_data

      Does not apply the mask (use apply_mask_to_map etc for that)

      Uses cctbx.maptbx.mask.create_mask_with_mask_data to do it

      Requires origin to be zero of both self and new mask
    '''

    assert isinstance(map_data, flex.double)
    assert self.map_data().all() == map_data.all()
    assert map_data.origin() == (0,0,0)
    assert self.origin_is_zero()

    from cctbx.maptbx.mask import create_mask_with_map_data as cm
    self._created_mask = cm(map_data = map_data,
      map_manager = self)


  def create_mask_around_density(self,
      resolution = None,
      molecular_mass = None,
      sequence = None,
      solvent_content = None):
    '''
      Use cctbx.maptbx.mask.create_mask_around_density to create a
       mask automatically

      Does not apply the mask (use apply_mask_to_map etc for that)

      Parameters are:
       resolution : resolution of map, taken from self.resolution() if not
          specified
       molecular_mass: optional mass (Da) of object in density
       sequence: optional sequence of object in density
       solvent_content : optional solvent_content of map


    '''

    if not resolution:
      resolution = self.resolution()

    from cctbx.maptbx.mask import create_mask_around_density as cm
    self._created_mask = cm(map_manager = self,
        resolution = resolution,
        molecular_mass = molecular_mass,
        sequence = sequence,
        solvent_content = solvent_content, )

  def create_mask_around_edges(self, boundary_radius = None):
    '''
      Use cctbx.maptbx.mask.create_mask_around_edges to create a mask around
      edges of map.  Does not make a soft mask.  For a soft mask,
      follow with soft_mask(boundary_radius =boundary_radius)
      The radius is to define the boundary around the map.

      Does not apply the mask (use apply_mask_to_map etc for that)
    '''

    if boundary_radius is None:
      boundary_radius = self.resolution()

    from cctbx.maptbx.mask import create_mask_around_edges as cm
    self._created_mask = cm(map_manager = self,
      boundary_radius = boundary_radius)

  def create_mask_around_atoms(self, model, mask_atoms_atom_radius = None,
      invert_mask = None):
    '''
      Use cctbx.maptbx.mask.create_mask_around_atoms to create a mask around
      atoms in model

      Does not apply the mask (use apply_mask_to_map etc for that)

      mask_atoms_atom_radius default is max(3, resolution)

      invert_mask makes outside 1 and inside 0
    '''

    assert model is not None
    if mask_atoms_atom_radius is None:
      mask_atoms_atom_radius = max(3, self.resolution())

    from cctbx.maptbx.mask import create_mask_around_atoms as cm
    self._created_mask = cm(map_manager = self,
      model = model,
      mask_atoms_atom_radius = mask_atoms_atom_radius,
      invert_mask = invert_mask)

  def soft_mask(self, soft_mask_radius = None):
    '''
      Make mask a soft mask. Just uses method in cctbx.maptbx.mask
      Use resolution for soft_mask radius if not specified
    '''
    assert self._created_mask is not None
    if soft_mask_radius is None:
      soft_mask_radius = self.resolution()
    self._created_mask.soft_mask(soft_mask_radius = soft_mask_radius)

  def apply_mask(self, set_outside_to_mean_inside = False):
    '''
      Replace map_data with masked version based on current mask
      Just uses method in create_mask_around_atoms
    '''

    assert self._created_mask is not None
    new_mm = self._created_mask.apply_mask_to_other_map_manager(
      other_map_manager = self,
      set_outside_to_mean_inside = set_outside_to_mean_inside)
    self.set_map_data(map_data = new_mm.map_data())  # replace map data

  def delete_mask(self):
    '''Remove working mask'''
    self._created_mask = None

  def get_mask_as_map_manager(self):
    '''Return a map_manager containing the working mask'''
    assert self._created_mask is not None
    return self._created_mask.map_manager()

  def initialize_map_data(self, map_value = 0):
    '''
      Set all values of map_data to map_value
    '''
    s = (self.map_data() != map_value )
    self.map_data().set_selected(s, map_value)

  def set_map_data(self, map_data = None):
    '''
      Replace self.data with map_data. The two maps must have same gridding

      NOTE: This uses selections to copy all the data in map_data into
      self.data.  The map_data object is not associated with self.data, the
      data is simply copied.  Also as self.data is modified in place, any
      objects that currently are just pointers to self.data are affected.
    '''
    assert self.map_data().origin() == map_data.origin()
    assert self.map_data().all() == map_data.all()
    sel = flex.bool(map_data.size(), True)
    self.data.as_1d().set_selected(sel, map_data.as_1d())

  def as_map_model_manager(self):
    '''  Return a map_model_manager'''
    from iotbx.map_model_manager import map_model_manager
    return map_model_manager(map_manager = self)

  def invert_hand(self):
    ''' Invert the hand of this map by swapping the order of all sections
        in Z'''

    map_data = self.map_data()
    # Swap all sections in Z, one pair at a time, in place
    nx,ny,nz = [ne - ns  for ne,ns in zip(map_data.all(), map_data.origin())]
    assert nx*ny*nz == map_data.size()
    nxstart,nystart,nzstart= map_data.origin()
    for k in range (int(nz//2)):
      i = k + nzstart
      ii = nzstart + nz -i -1  # sections to swap are i, ii
      data_i = maptbx.copy(map_data,
         (nxstart,nystart,i),
         (nxstart+nx-1, nystart+ny-1,i)).deep_copy()
      assert data_i.size() == nx * ny
      data_ii = maptbx.copy(map_data,
         (nxstart,nystart,ii),
         (nxstart+nx-1, nystart+ny-1,ii)).deep_copy()
      assert data_i.size() == nx * ny
      assert data_ii.size() == data_i.size()
      acc = flex.grid(data_i.all())
      data_i.reshape(acc)
      data_ii.reshape(acc)

      maptbx.set_box(
        map_data_from = data_ii,
        map_data_to   = map_data,
        start         = (nxstart,nystart,i),
        end           = (nxstart+nx, nystart+ny, i+1))

      maptbx.set_box(
        map_data_from = data_i,
        map_data_to   = map_data,
        start         = (nxstart,nystart,ii),
        end           = (nxstart+nx, nystart+ny, ii+1))

  def as_full_size_map(self):
    '''
      Create a full-size map with the current map inside it, padded by zero

      A little tricky because the starting map is going to have its origin at
      (0, 0, 0) but the map we are creating will have that point at
      self.origin_shift_grid_units.

      First use box.with_bounds to create map from -self.origin_shift_grid_units
       to -self.origin_shift_grid_units+self.unit_cell_grid-(1, 1, 1).  Then
      shift that map to place origin at (0, 0, 0)

      If the map is full size already, return the map as is
      If the map is bigger than full size stop as this is not suitable

    '''

    # Check to see if this is full size or bigger
    full_size_minus_working=subtract_tuples_int(self.unit_cell_grid,
      self.map_data().all())


    if full_size_minus_working in [(-1,-1,-1),(0, 0, 0)]: # Exactly full size already. Done
      assert self.origin_shift_grid_units == (0, 0, 0)
      assert self.map_data().origin() == (0, 0, 0)
      return self
    # Must not be bigger than full size already
    assert flex.double(full_size_minus_working).min_max_mean().min >= -1
    working_lower_bounds = self.origin_shift_grid_units
    working_upper_bounds = tuple([i+j-1 for i, j in zip(working_lower_bounds,
      self.map_data().all())])
    lower_bounds = tuple([-i for i in self.origin_shift_grid_units])
    upper_bounds = tuple([i+j-1 for i, j in zip(lower_bounds, self.unit_cell_grid)])
    new_lower_bounds = tuple([i+j for i, j in zip(
      lower_bounds, self.origin_shift_grid_units)])
    new_upper_bounds = tuple([i+j for i, j in zip(
      upper_bounds, self.origin_shift_grid_units)])
    print("Creating full-size map padding outside of current map with zero",
      file = self.log)
    print("Bounds of current map: %s to %s" %(
     str(working_lower_bounds), str(working_upper_bounds)), file = self.log)
    print("Bounds of new map: %s to %s" %(
     str(new_lower_bounds), str(new_upper_bounds)), file = self.log)

    from cctbx.maptbx.box import with_bounds
    box = with_bounds(self,
       lower_bounds = lower_bounds,
       upper_bounds = upper_bounds,
       log = self.log)
    box.map_manager().set_original_origin_and_gridding(
       original_origin = (0, 0, 0))

    box.map_manager().add_label(
       "Restored full size from box %s - %s, pad with zero" %(
     str(working_lower_bounds), str(working_upper_bounds)))
    assert box.map_manager().origin_shift_grid_units == (0, 0, 0)
    assert box.map_manager().map_data().origin() == (0, 0, 0)
    assert box.map_manager().map_data().all() == box.map_manager().unit_cell_grid
    if box.map_manager().unit_cell_crystal_symmetry().space_group_number() == 1:
      assert box.map_manager().unit_cell_crystal_symmetry().is_similar_symmetry(
        box.map_manager().crystal_symmetry())
    else:
      assert box.map_manager().crystal_symmetry().space_group_number() == 1
      from cctbx import crystal
      assert box.map_manager().crystal_symmetry().is_similar_symmetry(
        crystal.symmetry(
           box.map_manager().unit_cell_crystal_symmetry().unit_cell(),
           1))
    return box.map_manager()


  def cc_to_other_map_manager(self, other_map_manager):
    '''Get map correlation to other map manager'''
    assert self.is_similar(other_map_manager)

    return flex.linear_correlation(self.map_data().as_1d(),
     other_map_manager.map_data().as_1d()).coefficient()

  def density_at_sites_cart(self, sites_cart):
    '''
    Return flex.double list of density values corresponding to sites (Cartesian
    coordinates in A).
    Note that coordinates are relative to the current
    origin of the map (normally set to (0,0,0) before working with the map,
    see sites_cart_to_sites_cart_absolute and
    sites_cart_absolute_to_sites_cart.)
    '''
    assert isinstance(sites_cart, flex.vec3_double)

    from cctbx.maptbx import real_space_target_simple_per_site
    density_values = real_space_target_simple_per_site(
      unit_cell = self.crystal_symmetry().unit_cell(),
      density_map = self.map_data(),
      sites_cart = sites_cart)

    # Cross off anything outside the box if wrapping is false
    if not self.wrapping():
      sites_frac = self.crystal_symmetry().unit_cell().fractionalize(sites_cart)
      x,y,z = sites_frac.parts()
      s = (
         (x < 0) |
         (y < 0) |
         (z < 0) |
         (x > 1) |
         (y > 1) |
         (z > 1)
         )
      density_values.set_selected(s,0)
    return density_values



  def get_density_along_line(self,
      start_site = None,
      end_site = None,
      n_along_line = 10,
      include_ends = True):

    '''
      Return group_args object with density values and coordinates
      along a line segment from start_site to end_site
      (Cartesian coordinates in A) with n_along_line sampling points.
      Optionally include/exclude ends.
      Note that coordinates are relative to the current
      origin of the map (normally set to (0,0,0) before working with the map,
      see sites_cart_to_sites_cart_absolute and
      sites_cart_absolute_to_sites_cart.)
    '''
    along_sites = flex.vec3_double()
    if include_ends:
      start = 0
      end = n_along_line+1
    else:
      start = 1
      end = n_along_line

    for i in range(start, end):
      weight = (i/n_along_line)
      along_line_site = col(start_site)*weight+col(end_site)*(1-weight)
      along_sites.append(along_line_site)
    along_density_values = self.density_at_sites_cart(sites_cart = along_sites)
    return group_args(
     along_density_values = along_density_values,
       along_sites = along_sites)

  def apply_spectral_scaling(self, d_min = None, d_max = None,
    n_bins = 100):
    '''Apply spectral scaling to a map to approximate intensity vs
       resolution expected for a protein structure'''

    print("Applying spectral scaling", file = self.log)
    map_coeffs = self.map_as_fourier_coefficients(d_min = d_min,
      d_max = d_max)
    from iotbx.map_model_manager import get_map_coeffs_as_fp_phi
    f_array_info = get_map_coeffs_as_fp_phi(
        map_coeffs, d_min= map_coeffs.d_min(), n_bins = n_bins)

    from cctbx.development.approx_amplitude_vs_resolution import \
       approx_amplitude_vs_resolution
    aavr = approx_amplitude_vs_resolution(generate_mock_rms_fc_list=False)
    target_scale_factors = aavr.get_target_scale_factors(f_array_info.f_array)

    # Now interpolate these scale factors:

    scale_array=f_array_info.f_array.binner().interpolate(
        target_scale_factors, 1) # d_star_power=1
    scaled_f_array=f_array_info.f_array.customized_copy(
          data=f_array_info.f_array.data()*scale_array)

    new_map_coeffs = scaled_f_array.phase_transfer(
       phase_source=f_array_info.phases, deg=True)
    new_mm = self.fourier_coefficients_as_map_manager(new_map_coeffs)

    self.set_map_data(map_data = new_mm.map_data())  # replace map data


  def scale_map(self, scale = None):
    '''
      Multiply values in map by scale.
    '''
    self.set_map_data(map_data = scale * self.map_data())  # replace map data

  def resolution_filter(self, d_min = None, d_max = None):
    '''
      High- or low-pass filter the map in map_manager.
      Changes and overwrites contents of this map_manager.
      Remove all components with resolution < d_min or > d_max
      Either d_min or d_max or both can be None.
      To make a low_pass filter with cutoff at 3 A, set d_min=3
      To make a high_pass filter with cutoff at 2 A, set d_max=2

    '''
    map_coeffs = self.map_as_fourier_coefficients(d_min = d_min,
      d_max = d_max)
    mm = self.fourier_coefficients_as_map_manager( map_coeffs=map_coeffs)
    self.set_map_data(map_data = mm.map_data())  # replace map data


  def gaussian_filter(self, smoothing_radius):
    '''
      Gaussian blur the map in map_manager with given smoothing radius.
      Changes and overwrites contents of this map_manager.
    '''
    assert smoothing_radius is not None

    map_data = self.map_data()
    from cctbx.maptbx import smooth_map
    smoothed_map_data = smooth_map(
        map              = map_data,
        crystal_symmetry = self.crystal_symmetry(),
        rad_smooth       = smoothing_radius)
    self.set_map_data(map_data = smoothed_map_data)  # replace map data

  def binary_filter(self, threshold = 0.5):
    '''
      Apply a binary filter to the map (value at pixel i,j,k=1 if average
      of all 27 pixels within 1 of this one is > threshold, otherwise 0)
      Changes and overwrites contents of this map_manager.
    '''

    assert self.origin_is_zero()

    map_data=self.map_data()

    from cctbx.maptbx import binary_filter
    bf=binary_filter(map_data,threshold).result()
    self.set_map_data(map_data = bf)  # replace map data

  def randomize(self,
      d_min = None,
      low_resolution_fourier_noise_fraction=0.01,
      high_resolution_fourier_noise_fraction=2,
      low_resolution_real_space_noise_fraction=0,
      high_resolution_real_space_noise_fraction=0,
      low_resolution_noise_cutoff=None,
      random_seed = None,
         ):
    '''
      Randomize a map.

      Unique aspect of this noise generation is that it can be specified
      whether the noise is local in real space (every point in a map
      gets a random value before Fourier filtering), or local in Fourier
      space (every Fourier coefficient gets a complex random offset).
      Also the relative contribution of each type of noise vs resolution
      can be controlled.

      Parameters:

      d_min:  high-resolution limit in Fourier transformations

      low_resolution_fourier_noise_fraction (float, 0): Low-res Fourier noise
      high_resolution_fourier_noise_fraction (float, 0): High-res Fourier noise
      low_resolution_real_space_noise_fraction(float, 0): Low-res
          real-space noise
      high_resolution_real_space_noise_fraction (float, 0): High-res
          real-space noise
      low_resolution_noise_cutoff (float, None):  Low resolution where noise
          starts to be added

    '''

    assert self.origin_is_zero()

    if d_min is None:
      d_min = self.resolution()

    map_data=self.map_data()
    if random_seed is None:
      import random
      random_seed = random.randint(1,100000)
    from cctbx.development.create_models_or_maps import generate_map
    new_map_manager =generate_map(
      map_manager = self,   # gridding etc
      map_coeffs = self.map_as_fourier_coefficients(),
      d_min = d_min,
      low_resolution_fourier_noise_fraction=
         low_resolution_fourier_noise_fraction,
      high_resolution_fourier_noise_fraction=
         high_resolution_fourier_noise_fraction,
      low_resolution_real_space_noise_fraction=
         low_resolution_real_space_noise_fraction,
      high_resolution_real_space_noise_fraction=
         high_resolution_real_space_noise_fraction,
      low_resolution_noise_cutoff=
         low_resolution_noise_cutoff,
      random_seed = random_seed)

    self.set_map_data(map_data = new_map_manager.map_data())  # replace map data

  def deep_copy(self):
    '''
      Return a deep copy of this map_manager object
      Uses customized_copy to deepcopy everything including map_data

      Origin does not have to be at (0, 0, 0) to apply
    '''
    return self.customized_copy(map_data = self.map_data())

  def customized_copy(self, map_data = None, origin_shift_grid_units = None,
      use_deep_copy_for_map_data = True,
      crystal_symmetry_space_group_number = None,
      wrapping = None,):
    '''
      Return a customized deep_copy of this map_manager, replacing map_data with
      supplied map_data.

      The map_data and any _created_mask will be deep_copied before using
      them unless use_deep_copy_for_map_data = False

      Normally this customized_copy is applied with a map_manager
      that has already shifted the origin to (0, 0, 0) with shift_origin.

      Normally the new map_data will have the same dimensions of the current
      map_data. Normally new map_data will also have origin at (0, 0, 0).

      NOTE: It is permissible for map_data to have different bounds or origin
      than the current self.map_data.  In this case you must specify a new
      value of origin_shift_grid_units corresponding to this new map_data.
      This new origin_shift_grid_units specifies the original position in the
      full unit cell grid of the most-negative corner grid point of the
      new map_data. The new map_manager will still have the same unit
      cell dimensions and grid as the original.

      NOTE: It is permissible to get a customized copy before shifting the
      origin.  Applying with non-zero origin requires that:
         self.origin_shift_grid_units == (0, 0, 0)
         origin_shift_grid_units = (0, 0, 0)
         map_data.all() (size in each direction)  of current and new maps
            are the same.
         origins of current and new maps are the same

       NOTE: wrapping is normally copied from original map, but if new map is
       not full size then wrapping is always set to False.

      If crystal_symmetry_space_group_number is specified, use it
    '''

    # Make a deep_copy of map_data and _created_mask unless
    #    use_deep_copy_for_map_data = False

    if use_deep_copy_for_map_data:
      map_data = map_data.deep_copy()
      created_mask = deepcopy(self._created_mask)
    else:
      created_mask = self._created_mask

    assert map_data is not None # Require map data for the copy

    if map_data.origin() !=  (0, 0, 0):

      # Make sure all the assumptions are satisfied so we can just copy
      assert self.origin_shift_grid_units == (0, 0, 0)
      assert origin_shift_grid_units in [None, (0, 0, 0)]
      assert self.map_data().all() == map_data.all()
      assert self.map_data().origin() == map_data.origin()

      # Now just go ahead and copy using origin_shift_grid_units = (0, 0, 0)
      origin_shift_grid_units = (0, 0, 0)

    elif origin_shift_grid_units is None:  # use existing origin shift
      assert map_data.all()  ==  self.map_data().all() # bounds must be same
      origin_shift_grid_units = deepcopy(self.origin_shift_grid_units)

    # Keep track of change in shift_cart
    original_shift_cart=self.shift_cart()

    # Deepcopy this object and then set map_data and origin_shift_grid_units


    mm = map_manager(
     file_name = None,
     map_data = map_data,
     unit_cell_grid = self.unit_cell_grid,
     unit_cell_crystal_symmetry = self._unit_cell_crystal_symmetry,
     origin_shift_grid_units = origin_shift_grid_units,
     ncs_object = self._ncs_object.deep_copy() if self._ncs_object else None,
     wrapping = False,
     experiment_type = self._experiment_type,
     scattering_table = self._scattering_table,
     resolution = self._resolution,
     log = self.log,)

    mm._is_mask = self._is_mask
    mm._is_dummy_map_manager = self._is_dummy_map_manager
    # Set things that are not necessarily the same as in self:
    mm._created_mask = created_mask  # using self._created_mask or a
    mm.file_name = self.file_name
    mm.program_name = self.program_name
    mm.limitations = self.limitations
    mm.labels = self.labels


    if wrapping is not None:
      desired_wrapping = wrapping
    else:
      desired_wrapping = self._wrapping

    if mm.is_full_size():
      mm.set_wrapping(desired_wrapping)
    else: #
      mm.set_wrapping(False)

    # Set up _crystal_symmetry for the new object
    mm.set_crystal_symmetry_of_partial_map(
      space_group_number = crystal_symmetry_space_group_number)
      # Required and must be last


    # Keep track of change in shift_cart
    delta_origin_shift_grid_units = tuple([new - orig for new, orig in zip (
        origin_shift_grid_units, self.origin_shift_grid_units)])
    delta_shift_cart = tuple([-x for x in self.grid_units_to_cart(
       delta_origin_shift_grid_units)])
    new_shift_cart= tuple([
        o+d for o,d in zip(original_shift_cart,delta_shift_cart)])

    if self._ncs_object:
      mm._ncs_object = self._ncs_object.deep_copy(
        coordinate_offset=delta_shift_cart)
      assert approx_equal(mm.shift_cart(),mm._ncs_object.shift_cart())
    else:
      mm._ncs_object = None

    return mm

  def set_experiment_type(self, experiment_type):
    ''' Set the experiment type
       xray,neutron, or cryo_em
       If scattering_table is not defined, it is guessed from experiment_type
    '''
    self._experiment_type = experiment_type
    self._set_up_experiment_type_and_scattering_table_and_resolution()

  def set_scattering_table(self, scattering_table):
    ''' Set the scattering table (type of scattering)
       electron:  cryo_em
       n_gaussian x-ray (standard)
       wk1995:    x-ray (alternative)
       it1992:    x-ray (alternative)
       neutron:   neutron scattering
    '''
    self._scattering_table = scattering_table
    self._set_up_experiment_type_and_scattering_table_and_resolution()

  def set_resolution(self, resolution):
    ''' Set the nominal resolution of map
    '''
    self._resolution = resolution

  def experiment_type(self):
    '''Return the experiment type (xray or cryo_em)'''
    return self._experiment_type

  def minimum_resolution(self, set_minimum_resolution = True):
    '''
      Get minimum resolution.  If set previously, use that value
    '''
    if self._minimum_resolution:
      return self._minimum_resolution

    from cctbx.maptbx import d_min_from_map
    minimum_resolution = d_min_from_map(
           map_data=self.map_data(),
           unit_cell=self.crystal_symmetry().unit_cell())

    if set_minimum_resolution:
      self._minimum_resolution = minimum_resolution

    return minimum_resolution

  def resolution(self, force = False, method = 'd99', set_resolution = True):
    ''' Get nominal resolution
        Return existing if present unless force is True
        choices:
                  d9: resolution correlated at 0.9 with original
                  d99: resolution correlated at 0.99 with original
                  d999: resolution correlated at 0.999 with original
                  d_min: d_min_from_map
    '''
    if self._resolution is not None and (not force):
      return self._resolution


    assert method in ['d99','d9','d999','d_min']


    working_resolution = -1 # now get it

    if method in ['d99','d9','d999'] and \
        self.map_data().count(0) != self.map_data().size():
      from cctbx.maptbx import d99
      if self.origin_is_zero():
        map_data = self.map_data()
      else:
        map_data = self.map_data().deep_copy()
      d99_object = d99(
         map = map_data, crystal_symmetry = self.crystal_symmetry())

      working_resolution = getattr(d99_object.result,method,-1)

    from cctbx.maptbx import d_min_from_map  # get this to check
    d_min_estimated_from_map = self.minimum_resolution()

    if working_resolution < d_min_estimated_from_map:  # we didn't get it or want to use d_min
      working_resolution = d_min_estimated_from_map

    if set_resolution:
      self._resolution = working_resolution
    return working_resolution

  def scattering_table(self):
    '''Return the scattering table to use:
       electron:  cryo_em
       n_gaussian x-ray (standard)
       wk1995:    x-ray (alternative)
       it1992:    x-ray (alternative)
       neutron:   neutron scattering
    '''
    return self._scattering_table

  def ncs_object(self):
    ''' Return the NCS object '''
    return self._ncs_object

  def _set_up_experiment_type_and_scattering_table_and_resolution(self):
    '''Set up the experiment type, scattering table, and resolution
    '''
    default_scattering_table_dict = {
     'xray':'n_gaussian',
     'neutron':'neutron',
     'cryo_em':'electron',
     }

    if self._experiment_type not in [None, Auto]:
      assert self._experiment_type in ['xray','neutron','cryo_em']
      if self.wrapping() and self._experiment_type=='cryo_em':
        raise Sorry("Cannot use wrapping if experiment_type is 'cryo_em'")

    else:  # Try to guess experiment_type
      if self.crystal_symmetry().space_group_number() > 1:
        # Has space-group symmmetry, not cryo_em
        self._experiment_type = 'xray'  # could be neutron of course
      elif self.is_full_size() and self.wrapping() is False:
        # No space-group symmetry, full size map, no wrapping: cryo_em
        self._experiment_type = 'cryo_em'  # full size map and no wrapping
      elif self.is_full_size() and self.wrapping() is True:
        # P1 symmetry, full size map, wrapping True: xray
        self._experiment_type = 'xray'  # full size map and wrapping
      else:
        # P1 symmetry, not a full-size map...cannot tell
        self._experiment_type = None

    if self._experiment_type is not None:
      if self._scattering_table is None:
        self._scattering_table = default_scattering_table_dict[
          self._experiment_type]

    if self._scattering_table not in [None, Auto]:
      assert self._scattering_table in ['electron','n_gaussian',
       'wk1995','it1992','neutron']

    if self._scattering_table is Auto:
      self._scattering_table = None

    if self._resolution is Auto:
      self._resolution = None

    if self._experiment_type is Auto:
      self._experiment_type = None

    assert not (self._wrapping is Auto)


  def set_wrapping(self, wrapping_value):
    '''
       Set wrapping to be wrapping_value
    '''
    assert isinstance(wrapping_value, bool)
    self._wrapping = wrapping_value
    if self._wrapping:
      if not self.is_full_size():
        raise Sorry("You cannot set wrapping=True for a map that is not full size")

  def wrapping(self):
    '''
      Report if map can be wrapped

    '''
    return self._wrapping

  def is_full_size(self):
    '''
      Report if map is full unit cell
    '''
    if self.map_data().all()  ==  self.unit_cell_grid:
      return True
    else:
      return False

  def is_consistent_with_wrapping(self, relative_sd_tol = 0.01):
    '''
      Report if this map looks like it is a crystallographic map and can be
      wrapped

      If it is not full size...no wrapping
      If origin is not at zero...no wrapping
      If it is not periodic, no wrapping
      If very small or resolution_factor for map is close to 0.5...cannot tell
      If has all zeroes (or some other constant on edges) ... no wrapping

      relative_sd_tol defines how close to constant values at edges must be
      to qualify as "constant"

      Returns True, False, or None (unsure)

    '''
    if not self.is_full_size():
      return False

    if self.map_data().origin() != (0, 0, 0):
      return False
    from cctbx.maptbx import is_periodic, is_bounded_by_constant
    if is_bounded_by_constant(self.map_data(),
       relative_sd_tol = relative_sd_tol):  # Looks like a cryo-EM map
      return False

    # Go with whether it looks periodic (cell translations give similar values
    #  or transform of high-res data is mostly at edges of cell)
    return is_periodic(self.map_data())  # Can be None if unsure


  def is_similar(self, other = None,
     absolute_angle_tolerance = 0.01,
     absolute_length_tolerance = 0.01,
     ):
    '''Determine whether this map_manager is similar (symmetry, gridding,
        size) to another map_manager.
    '''
    # Check to make sure origin, gridding and symmetry are similar
    self._warning_message=""

    if tuple(self.origin_shift_grid_units) !=  tuple(
        other.origin_shift_grid_units):
      self._warning_message="Origin shift grid units "+  \
        "(%s) does not match other (%s)" %(
        str(self.origin_shift_grid_units),str(other.origin_shift_grid_units))
      return False
    if not self.unit_cell_crystal_symmetry().is_similar_symmetry(
      other.unit_cell_crystal_symmetry(),
      absolute_angle_tolerance = absolute_angle_tolerance,
      absolute_length_tolerance = absolute_length_tolerance,):
      self._warning_message="Unit cell crystal symmetry:"+ \
        "\n%s\n does not match other:\n%s\n" %(
        str(self.unit_cell_crystal_symmetry()),
         str(other.unit_cell_crystal_symmetry()))
      return False
    if not self.crystal_symmetry().is_similar_symmetry(
      other.crystal_symmetry(),
      absolute_angle_tolerance = absolute_angle_tolerance,
      absolute_length_tolerance = absolute_length_tolerance):
      self._warning_message="Crystal symmetry:"+ \
        "\n%s\ndoes not match other: \n%s\n" %(
        str(self.crystal_symmetry()),
         str(other.crystal_symmetry()))
      return False
    if self.map_data().all()!=  other.map_data().all():
      self._warning_message="Existing map gridding "+ \
        "(%s) does not match other (%s)" %(
         str(self.map_data().all()),str(other.map_data().all()))
      return False
    if self.unit_cell_grid !=  other.unit_cell_grid:
      self._warning_message="Full map gridding "+ \
        "(%s) does not match other (%s)" %(
         str(self.map_data().all()),str(other.map_data().all()))
      return False

    # Make sure wrapping is same for all
    if ( self.wrapping() !=  other.wrapping()):
      self._warning_message="Wrapping "+ "(%s) does not match other (%s)" %(
         str(self.wrapping()),
         str(other.wrapping()))
      return False

    # Make sure ncs objects are similar if both have one
    if (self.ncs_object() is not None) and (
        other.ncs_object() is not None):
      if not other.ncs_object().is_similar_ncs_object(self.ncs_object()):
        text1=self.ncs_object().as_ncs_spec_string()
        text2=other.ncs_object().as_ncs_spec_string()
        self._warning_message="NCS objects do not match"+ \
           ":\n%s\n does not match other:\n%s" %( text1,text2)
        return False

    return True

  def cart_to_grid_units(self, xyz):
    '''Convert xyz (Cartesian) to grid units'''

    return tuple([int(0.5 + x * n) for x,n in
       zip(self.crystal_symmetry().unit_cell().fractionalize(xyz),
        self.map_data().all())])

  def grid_units_to_cart(self, grid_units):
    ''' Convert grid units to Cartesian coordinates '''
    x = grid_units[0]/self.unit_cell_grid[0]
    y = grid_units[1]/self.unit_cell_grid[1]
    z = grid_units[2]/self.unit_cell_grid[2]
    return self.unit_cell().orthogonalize(tuple((x, y, z)))


  def shifted(self, eps=1.e-3):
    ''' Return True if this map has been shifted from its original
     location (e.g., by boxing the map).
     Checks self.shift_cart() to determine if map has been shifted.
    '''

    r = self.shift_cart()
    if(r is None): return False
    if(flex.max(flex.abs(flex.double(r)))<=eps): return False
    return True

  def shift_cart(self):
    '''
     Return the shift_cart of this map from its original location.

     (the negative of the origin shift ) in Cartesian coordinates
     '''
    return tuple(
       [-x for x in self.grid_units_to_cart(self.origin_shift_grid_units)])

  def shift_ncs_object_to_match_map_and_return_new_ncs_object(self,ncs_object):
    '''
      Move the ncs_object to match this map. Also sets ncs_object shift_cart
      Returns new copy of ncs_obect and does not affect the original

      Note difference from set_ncs_object_shift_cart_to_match_map which
        sets the shift_cart but does not move the object
    '''
    if ncs_object.shift_cart():
      offset = tuple(
        [s - n for s, n in zip(self.shift_cart(), ncs_object.shift_cart())])
      ncs_object = ncs_object.coordinate_offset(offset)

    else:
      ncs_object = ncs_object.coordinate_offset(self.shift_cart())
    return ncs_object

  def shift_model_to_match_map(self, model):
    '''
      Move the model to match this map.
      Note difference from set_model_symmetries_and_shift_cart_to_match_map
       which sets model symmetry and shift_cart but does not move the model
    '''
    if model.shift_cart():
      offset = tuple(
        [s - n for s, n in zip(self.shift_cart(), model.shift_cart())])
      model.shift_model_and_set_crystal_symmetry(shift_cart=offset)
    else:
      model.shift_model_and_set_crystal_symmetry(shift_cart=self.shift_cart())

  def set_ncs_object_shift_cart_to_match_map(self, ncs_object):
    '''
      Set the ncs_object shift_cart to match map

      Overwrites any information in ncs_object on shift_cart
      Modifies ncs_object in place. Does not return anything

      Do not use this to try to shift the ncs object. That is done in
      the ncs object itself with ncs_object.coordinate_offset(shift_cart),
      which returns a new ncs object

      You can use ncs_object =
       self.shift_ncs_object_to_match_map_and_return_new_ncs_object(ncs_object)
         to shift the ncs object and set its shift_cart and get a new copy.
    '''

    # Set shift_cart (shift since readin) to match shift_cart for
    #   map (shift of origin is opposite of shift applied)
    ncs_object.set_shift_cart(self.shift_cart())

  def set_crystal_symmetry_to_p1(self,
     space_group_number = 1):
    '''
      Change the working crystal symmetry to P1
      This changes map in place
      Do a deep_copy first if you do not want it changed
      (Actually can set space group number to anything so you can set it back)
    '''
    print("\nSetting working crystal symmetry to P1 so "+
       "that edges can be masked", file = self.log)

    self.set_crystal_symmetry_of_partial_map(
      space_group_number = space_group_number)


  def set_model_symmetries_and_shift_cart_to_match_map(self,model):
    '''
      Set the model original and working crystal_symmetry to match map.

      Overwrites any information in model on symmetry and shift_cart
      Modifies model in place

      NOTE: This does not shift the coordinates in model.  It is used
      to fix crystal symmetry and set shift_cart, not to actually shift
      a model.
      For shifting a model, use:
         model.shift_model_and_set_crystal_symmetry(shift_cart=shift_cart)
    '''
    # Check if we really need to do anything
    if not model:
      return # nothing to do

    if self.is_compatible_model(model,
       require_match_unit_cell_crystal_symmetry = True):
      return # already fine

    if model.shift_cart() is not None and tuple(model.shift_cart()) != (0,0,0):
      # remove shift_cart
      model.set_shift_cart((0,0,0))

    # Set crystal_symmetry to match map. This changes the xray_structure.
    model.set_crystal_symmetry(self.crystal_symmetry())

    # Set original crystal symmetry to match map unit_cell_crystal_symmetry
    # This just changes a specification in the map, nothing else changes
    model.set_unit_cell_crystal_symmetry(self.unit_cell_crystal_symmetry())

    # Set shift_cart (shift since readin) to match shift_cart for
    #   map (shift of origin is opposite of shift applied)
    model.set_shift_cart(self.shift_cart())

  def check_consistency(self, stop_on_errors = True, print_errors = True,
        absolute_angle_tolerance = None,
        absolute_length_tolerance = None,
        shift_tol = None):
    """
    Carry out overall consistency checks. Used in map_model_manager
    Note: the stop_on_errors, print_errors, and 3 tolerance kw are used in
      map_model_manager when checking consistency there
    """

    if absolute_angle_tolerance is None:
      absolute_angle_tolerance = 0.01
    if absolute_length_tolerance is None:
      absolute_length_tolerance = 0.01
    if shift_tol is None:
      shift_tol = 0.001

    # Check crystal_symmetry, unit_cell_crystal_symmetry, shift_cart
    # For now, only shift_cart and only ncs_object are relevant

    ok = True
    if self.ncs_object():
      if (not self.is_compatible_ncs_object(self.ncs_object(),
         tol = shift_tol)):
        ok = False
        text = "NCS object does not have same shift_cart as map_manager" +\
          " %s vs %s" %(self.ncs_object().shift_cart(),
           self.shift_cart())

    if (not ok):
      if print_errors:
         print("** Mismatch in model object\n%s" %(text))
      if stop_on_errors:
        raise AssertionError(text)

    return ok

  def is_compatible_ncs_object(self, ncs_object, tol = 0.001):
    '''
      ncs_object is compatible with this map_manager if shift_cart is
      the same as map
    '''

    ok=True
    text=""

    map_shift=flex.double(self.shift_cart())
    ncs_object_shift=flex.double(ncs_object.shift_cart())
    delta=map_shift-ncs_object_shift
    mmm=delta.min_max_mean()
    if mmm.min < -tol or mmm.max > tol: # shifts do not match
      text="Shift of ncs object (%s) does not match shift of map (%s)" %(
         str(ncs_object_shift),str(map_shift))
      ok=False

    self._warning_message=text
    return ok

  def is_compatible_model(self, model,
       require_match_unit_cell_crystal_symmetry = False,
        absolute_angle_tolerance = 0.01,
        absolute_length_tolerance = 0.01,
        shift_tol = 0.001):
    '''
      Model is compatible with this map_manager if it is not specified as being
      different.

      They are different if:
        1. original and current symmetries are present and different from each
          other and do not match
        2. model current symmetry does not match map original or current
        3. model has a shift_cart (shift applied) different than map shift_cart

      If require_match_unit_cell_crystal_symmetry is True, then they are
      incompatible if anything is different.

      If require_match_unit_cell_crystal_symmetry is False, original
       symmetries do not have to match.  Model crystal_symmetry can match
       the unit_cell_crystal_symmetry or crystal_symmetry of the map.
    '''

    ok=None
    text=""

    if not model:
      return None

    model_uc=model.unit_cell_crystal_symmetry()
    model_sym=model.crystal_symmetry()
    map_uc=self.unit_cell_crystal_symmetry()
    map_sym=self.crystal_symmetry()

    model_uc = model_uc if model_uc and model_uc.unit_cell() is not None else None
    model_sym = model_sym if model_sym and model_sym.unit_cell() is not None else None
    map_uc = map_uc if map_uc and map_uc.unit_cell() is not None else None
    map_sym = map_sym if map_sym and map_sym.unit_cell() is not None else None


    if (not require_match_unit_cell_crystal_symmetry) and \
        (model_uc and model_sym and model_uc.is_similar_symmetry(model_sym)):
      # Ignore the model_uc because it may or may not have come from
      # model_sym
      model_uc = None

    if (not require_match_unit_cell_crystal_symmetry) and \
        (map_uc and map_sym and map_uc.is_similar_symmetry(map_sym)):
      # Ignore the map_uc because it may or may not have come from
      # map_sym
      map_uc = None

    text_model_uc=str(model_uc).replace("\n"," ")
    text_model=str(model_sym).replace("\n"," ")
    text_map_uc=str(map_uc).replace("\n"," ")
    text_map=str(map_sym).replace("\n"," ")

    assert map_sym # map_sym should always should be there. model_sym could be missing

    if require_match_unit_cell_crystal_symmetry and map_uc and (
      not model_uc) and (
       not map_sym.is_similar_symmetry(map_uc,
        absolute_angle_tolerance = absolute_angle_tolerance,
        absolute_length_tolerance = absolute_length_tolerance,
         )):
      ok=False
      text="Model and map are different because "+\
          "require_match_unit_cell_crystal_symmetry is set and "+\
          "model does not have original_crystal_symmetry, and " +\
        "model symmetry: \n%s\n does not match map original symmetry:" %(
          model_sym) +\
        "\n%s\n. Current map symmetry is: \n%s\n " %(
         text_map_uc,text_map)

    elif model_sym and model_uc and map_uc and (
        (not map_uc.is_similar_symmetry(map_sym,
        absolute_angle_tolerance = absolute_angle_tolerance,
        absolute_length_tolerance = absolute_length_tolerance,))
         or (not model_uc.is_similar_symmetry(model_sym,
        absolute_angle_tolerance = absolute_angle_tolerance,
        absolute_length_tolerance = absolute_length_tolerance,))
         ) and (
         (not model_uc.is_similar_symmetry(map_uc,
        absolute_angle_tolerance = absolute_angle_tolerance,
        absolute_length_tolerance = absolute_length_tolerance,
        )) or
         (not model_sym.is_similar_symmetry(map_sym,
        absolute_angle_tolerance = absolute_angle_tolerance,
        absolute_length_tolerance = absolute_length_tolerance,
         ) ) ):
       ok=False # model and map_manager uc present and some symmetries do not match
       text="Model original symmetry: \n%s\n and current symmetry :\n%s\n" %(
          text_model_uc,text_model)+\
          "do not match map unit_cell symmetry:"+\
         " \n%s\n and map current symmetry: \n%s\n symmetry" %(
           text_map_uc,text_map)
    elif map_uc and model_sym and (not model_sym.is_similar_symmetry(map_uc,
        absolute_angle_tolerance = absolute_angle_tolerance,
        absolute_length_tolerance = absolute_length_tolerance,
        )) and (not
              model_sym.is_similar_symmetry(map_sym,
        absolute_angle_tolerance = absolute_angle_tolerance,
        absolute_length_tolerance = absolute_length_tolerance,
        )):
       ok=False # model does not match either map symmetry
       text="Model current symmetry: \n%s\n" %(
          text_model)+\
          " does not match map unit_cell symmetry:"+\
           " \n%s\n or map current symmetry: \n%s\n" %(
           text_map_uc,text_map)
    elif model_sym and (not model_uc) and (not map_uc) and (
           not model_sym.is_similar_symmetry(map_sym,
        absolute_angle_tolerance = absolute_angle_tolerance,
        absolute_length_tolerance = absolute_length_tolerance,
        )):
       ok=False # model has no uc and model symmetry does not match map symmetry
       text="Model current symmetry: \n%s\n" %(
          text_model)+\
          " does not match map symmetry: \n%s\n" %( text_map)

    elif require_match_unit_cell_crystal_symmetry and (
        not model_sym) and (not model_uc):
       ok=False # model does not have any symmetry so it does not match
       text="Model has no symmetry and cannot match any map"

    elif (not model_sym) and (not model_uc):
       ok=True # model does not have any symmetry so anything is ok
       text="Model has no symmetry and can match any map symmetry"

    else:  # match

       ok=True
       text="Model original symmetry: \n%s\n and current symmetry: \n%s\n" %(
          text_model_uc,text_model)+\
          "are compatible with "+\
          "map unit_cell symmetry:\n%s\n and current map symmetry:\n%s\n" %(
           text_map_uc,text_map)

    assert isinstance(ok, bool)  # must have chosen

    map_shift_cart=self.shift_cart()
    if ok and (map_shift_cart != (0,0,0)):
      if model.shift_cart() is None: # map is shifted but not model
        ok=False
        text+=" However map is shifted (shift_cart=%s) but model is not" %(
           str(map_shift_cart))
      else:
        map_shift=flex.double(map_shift_cart)
        model_shift=flex.double(model.shift_cart())
        delta=map_shift-model_shift
        mmm=delta.min_max_mean()
        if mmm.min<-shift_tol or mmm.max > shift_tol: # shifts do not match
          ok=False
          text+=" However map shift "+\
              "(shift_cart=%s) does not match model shift (%s)" %(
           str(map_shift),str(model_shift))
    self._warning_message=text
    return ok

  def warning_message(self):
    '''Return the warning message, if any'''
    if hasattr(self,'_warning_message'):
       return self._warning_message

  def set_mean_zero_sd_one(self):
    '''
     Function to normalize the map
     If the map has a constant value, do nothing
    '''
    map_data = self.map_data()
    map_data = map_data - flex.mean(map_data)
    sd = map_data.sample_standard_deviation()
    if sd is not None and sd != 0:
      map_data = map_data/sd
      self.set_map_data(map_data)

  def ncs_cc(self):
    '''Return value of NCS correlation if available'''
    if hasattr(self,'_ncs_cc'):
       return self._ncs_cc

  def absolute_center_cart(self,
       use_assumed_end = False,
       place_on_grid_point = False,
       use_unit_cell_grid = False):
    '''
     Return center of map (absolute position) in Cartesian coordinates
     A little tricky because for example the map goes from 0 to nx-1, not nx
       If use_assumed_end, go to nx
     Also map could start at non-zero origin
     If place_on_grid_point then guess the end by whether the center ends
       on a grid point
     If use_unit_cell_grid just find center of full unit cell

    '''

    if use_unit_cell_grid:  # Find center of unit cell
      return self.unit_cell_crystal_symmetry().unit_cell().orthogonalize(
        (0.5, 0.5, 0.5))

    elif place_on_grid_point:
      return tuple(col(self.crystal_symmetry().unit_cell().orthogonalize(
        tuple(col([int (0.5*n)/n + o/n for n,o in zip(
          self.map_data().all(),
          self.map_data().origin())])))) - col(self.shift_cart()))

    else:
      if use_assumed_end:
        n_end = 0
      else:
        n_end = 1
      return tuple(col(self.crystal_symmetry().unit_cell().orthogonalize(
        tuple(col([0.5*(n-n_end)/n + o/n for n,o in zip(
          self.map_data().all(),
          self.map_data().origin())])))) - col(self.shift_cart()))


  def map_map_cc(self, other_map_manager):
   ''' Return simple map correlation to other map_manager'''
   import iotbx.map_manager
   assert isinstance(other_map_manager, iotbx.map_manager.map_manager)
   return flex.linear_correlation(
      self.map_data().as_1d(), other_map_manager.map_data().as_1d()
       ).coefficient()


  def find_map_symmetry(self,
      include_helical_symmetry = False,
      symmetry_center = None,
      min_ncs_cc = None,
      symmetry = None,
      ncs_object = None,
      check_crystal_symmetry = True,
      only_proceed_if_crystal_symmetry = False,):

    '''
       Use run_get_symmetry_from_map tool in segment_and_split_map to find
       map symmetry and save it as an mmtbx.ncs.ncs.ncs object

       Here map symmetry is the reconstruction symmetry used to generate the
       map. Normally it is essentially perfect symmetry and normally the
       principal axes are aligned with x,y,z and normally the center is at
       the original center of the map.

       Sets self._warning_message if failure, sets self._ncs_object and
           self._ncs_cc if success

       This procedure may fail if the above assumptions do not hold.
       Optional center of map can be supplied, and minimum NCS correlation
       can also be supplied

       Requires that map_manager is already shifted to place origin at (0, 0, 0)

       Assumes that center of symmetry is at (1/2, 1/2, 1/2) in the full map

       It is optional to include search for helical symmetry. Reason is that
       this is much slower than other symmetries.

       symmetry (symbol such as c1, O, D7) can be supplied and search will be
       limited to that symmetry

       ncs_object can be supplied in which case it is just checked

       If check_crystal_symmetry, try to narrow down possibilities by looking
       for space-group symmetry first

       If only_proceed_if_crystal_symmetry, skip looking if nothing comes up
        with check_crystal_symmetry


    '''

    assert self.origin_is_zero()

    self._warning_message = ""
    self._ncs_cc = None

    from cctbx.maptbx.segment_and_split_map import \
       run_get_ncs_from_map, get_params

    if symmetry is None:
      symmetry = 'ALL'


    if symmetry_center is None:
      # Most likely map center is (1/2,1/2,1/2) in full grid
      symmetry_center = self.absolute_center_cart(use_assumed_end=True)
      # Our map is already shifted, so subtract off shift_cart
      symmetry_center = tuple(
        flex.double(symmetry_center) + flex.double(self.shift_cart()))

    params = get_params(args=[],
      symmetry = symmetry,
      include_helical_symmetry = include_helical_symmetry,
      symmetry_center = symmetry_center,
      min_ncs_cc = min_ncs_cc,
      return_params_only = True,
      )

    space_group_number = None
    if check_crystal_symmetry and symmetry == 'ALL' and (not ncs_object):
      # See if we can narrow it down looking at intensities at low-res
      d_min = 0.05*self.crystal_symmetry().unit_cell().volume()**0.333
      map_coeffs = self.map_as_fourier_coefficients(d_min=d_min)
      from iotbx.map_model_manager import get_map_coeffs_as_fp_phi
      f_array_info = get_map_coeffs_as_fp_phi(map_coeffs, d_min = d_min,
        n_bins = 15)
      ampl = f_array_info.f_array
      data = ampl.customized_copy(
        data = ampl.data(),sigmas = flex.double(ampl.size(),1.))
      from mmtbx.scaling.twin_analyses import symmetry_issues
      si = symmetry_issues(data)
      cs_possibility = si.xs_with_pg_choice_in_standard_setting
      space_group_number = cs_possibility.space_group_number()
      if space_group_number < 2:
        space_group_number = None
      if space_group_number is None and only_proceed_if_crystal_symmetry:
        return # skip looking further

    params.reconstruction_symmetry.\
          must_be_consistent_with_space_group_number = space_group_number
    new_ncs_obj, ncs_cc, ncs_score = run_get_ncs_from_map(params = params,
        map_data = self.map_data(),
        crystal_symmetry = self.crystal_symmetry(),
        out = self.log,
        ncs_obj = ncs_object)
    if (space_group_number) and (not new_ncs_obj):
      # try again without limits
      params.reconstruction_symmetry.\
          must_be_consistent_with_space_group_number = None
      new_ncs_obj, ncs_cc, ncs_score = run_get_ncs_from_map(params = params,
        map_data = self.map_data(),
        crystal_symmetry = self.crystal_symmetry(),
        out = self.log,
        ncs_obj = ncs_object)

    if new_ncs_obj:
      self._ncs_object = new_ncs_obj
      self._ncs_cc = ncs_cc
      self._ncs_object.set_shift_cart(self.shift_cart())
    else:
      self._warning_message = "No map symmetry found; ncs_cc cutoff of %s" %(
        min_ncs_cc)

  def _resample_on_different_grid_and_rebox(self, n_real = None,
       target_grid_spacing = None):
    '''
      Resample the boxed map on a grid of n_real and return new map_manager
      If an ncs_object is present, set its shift_cart too

      Returns new boxed map of similar size
      and location
    '''

    # Get starting lower and upper bounds (current map)
    lower_bounds_cart = self.grid_units_to_cart(self.origin_shift_grid_units)
    upper_bounds_cart = self.grid_units_to_cart(
     [o + a for o,a in zip(self.origin_shift_grid_units,
       self.map_data().all())]
     )

    # Determine if box is cubic
    box_dims = self.map_data().all()
    is_cubic_box = (box_dims[0] == box_dims[1]) and (box_dims[0] == box_dims[2])

    # Save NCS object if any and current shift_cart
    if self.ncs_object():
      working_ncs_object = self.ncs_object().deep_copy()
    else:
      working_ncs_object = None
    working_shift_cart = self.shift_cart()

    # Create full size map so that we can work easily
    mm_boxed_fs = self.as_full_size_map()
    mm_boxed_fs.set_ncs_object(working_ncs_object)

    assert tuple(mm_boxed_fs.origin_shift_grid_units) == (0,0,0)

    # Resample the full size map on new grid
    mm_boxed_fs_resample = mm_boxed_fs.resample_on_different_grid(
       n_real = n_real,
       target_grid_spacing = target_grid_spacing)

    # Now rebox the newly-gridded map

    # New bounds
    lower_bounds = mm_boxed_fs_resample.cart_to_grid_units(lower_bounds_cart)
    upper_bounds = mm_boxed_fs_resample.cart_to_grid_units(upper_bounds_cart)
    if is_cubic_box:
      box_dims = [1 + a - b for a,b in zip (upper_bounds,lower_bounds)]
      min_dim = min(box_dims)
      max_dim = max(box_dims)
      if min_dim != max_dim: # make them the same
       new_upper_bounds = []
       for lb, ub in zip(lower_bounds, upper_bounds):
         new_upper_bounds.append(lb + min_dim - 1)
       upper_bounds = new_upper_bounds
      box_dims = [1 + a - b for a,b in zip (upper_bounds,lower_bounds)]
      assert min(box_dims) == max(box_dims)
    mmm_boxed_fs_resample = mm_boxed_fs_resample.as_map_model_manager()

    mmm_boxed_fs_resample_boxed = \
      mmm_boxed_fs_resample.extract_all_maps_with_bounds(
        lower_bounds=lower_bounds, upper_bounds=upper_bounds)
    return mmm_boxed_fs_resample_boxed.map_manager()

  def resample_on_different_grid(self, n_real = None,
       target_grid_spacing = None):
    '''
      Resample the map on a grid of n_real and return new map_manager
      If an ncs_object is present, set its shift_cart too

      Allows map to be boxed; if so returns new boxed map of similar size
      and location.  If boxed map has cubic gridding, keep that after
      resampling.
    '''

    assert n_real or target_grid_spacing

    if tuple(self.origin_shift_grid_units) != (0,0,0):
      return self._resample_on_different_grid_and_rebox(n_real = n_real,
       target_grid_spacing = target_grid_spacing)

    if n_real is None:
      n_real = []
      for a, nn in zip(self.crystal_symmetry().unit_cell().parameters()[:3],
         self.map_data().all()):
        new_n = int (0.5+ a/target_grid_spacing)
        n_real.append( new_n)

    original_n_real = self.map_data().all()
    original_shift_cart = self.shift_cart()
    original_origin_shift_grid_units = self.origin_shift_grid_units

    map_coeffs = self.map_as_fourier_coefficients()
    map_data=maptbx.map_coefficients_to_map(
        map_coeffs       = map_coeffs,
        crystal_symmetry = map_coeffs.crystal_symmetry(),
        n_real           = n_real)


    new_origin_shift_grid_units = (0,0,0)

    if self.ncs_object():
      new_ncs_object = self.ncs_object().deep_copy()
      new_ncs_object.set_shift_cart(self.shift_cart())
    else:
      new_ncs_object = None
    mm = map_manager(
      map_data = map_data,
      unit_cell_grid = n_real,
      unit_cell_crystal_symmetry = map_coeffs.crystal_symmetry(),
      origin_shift_grid_units = new_origin_shift_grid_units,
      ncs_object = new_ncs_object,
      wrapping = self.wrapping(),
      experiment_type = self.experiment_type(),
      scattering_table = self.scattering_table(),
      resolution = self.resolution(),
     )
    return mm

  def get_boxes_to_tile_map(self,
     target_for_boxes = 24,
     box_cushion = 3,
     get_unique_set_for_boxes = None,
     dist_min = None,
     do_not_go_over_target = None,
     target_xyz_center_list = None,
       ):
    '''
     Return a group_args object with a list of lower_bounds and upper_bounds
     corresponding to a set of boxes that tiles the part of the map that is
     present.  The boxes may not be the same size but will tile to exactly
     cover the existing part of the map.
     Approximately target_for_boxes will be returned (may be fewer or greater)
     Also return boxes with cushion of box_cushion
     If get_unique_set_for_boxes is set, try to use map symmetry to identify
       duplicates and set ncs_object
     If target_xyz_center_list is set, use these points as centers but try
      to use standard box size.
    '''
    assert self.origin_is_zero()
    cushion_nx_ny_nz = tuple([int(0.5 + x * n) for x,n in
       zip(self.crystal_symmetry().unit_cell().fractionalize(
        (box_cushion,box_cushion,box_cushion)),
        self.map_data().all())])
    from cctbx.maptbx.box import get_boxes_to_tile_map
    box_info = get_boxes_to_tile_map(
       target_for_boxes = target_for_boxes,
       n_real = self.map_data().all(),
       crystal_symmetry = self.crystal_symmetry(),
       cushion_nx_ny_nz = cushion_nx_ny_nz,
       wrapping = self.wrapping(),
       do_not_go_over_target = do_not_go_over_target,
       target_xyz_center_list = target_xyz_center_list,
     )
    box_info.ncs_object = None
    if get_unique_set_for_boxes:
      if dist_min:
         max_distance = dist_min
      else:
         max_distance = self.resolution()
      n_before = len(box_info.lower_bounds_list)
      box_info = self._get_unique_box_info(
         box_info = box_info,
         max_distance = max_distance)

    return box_info

  def get_n_real_for_grid_spacing(self, grid_spacing = None):
    '''Identify values of gridding to match target grid spacing'''
    n_real = []
    for n,a in zip(self.map_data().all(),
       self.crystal_symmetry().unit_cell().parameters()):
      spacing = a/n
      target_n = (spacing/grid_spacing) * n
      n_real.append(int(target_n + 0.999))
    return n_real

  def sites_cart_to_sites_cart_absolute(self, sites_cart):
    """ Shift sites_cart that are relative to the boxed map to
        make them relative to the point (0,0,0) in absolute coordinates
   NOTE: sites_cart is a flex.vec3_double array
   NOTE: This is the opposite of sites_cart_absolute_to_sites_cart

   NOTE on shift_cart:

    Position of origin of boxed map:
     origin_position = self.grid_units_to_cart(self.origin_shift_grid_units)
     shift_cart = self.shift_cart() == - origin_position

    If you have Cartesian coordinates xyz for an atom relative to the boxed map,
    the absolute coordinates are:
      coords_abs =  col(xyz) + col(origin_position)
      coords_abs =  col(xyz) - col(self.shift_cart())

    """
    return  sites_cart - col(self.shift_cart())

  def sites_cart_absolute_to_sites_cart(self, sites_cart_absolute):
    """ Shift sites_cart that are in absolute coordinates to make them
        relative to the boxed map.
       NOTE: This is the opposite of sites_cart_to_sites_cart_absolute
       NOTE: sites_cart is a flex.vec3_double array
    """
    return  sites_cart_absolute + col(self.shift_cart())

  def peak_search(self,
      peak_search_level = 3,
      max_peaks = None,
      peak_cutoff            = None,
      interpolate            = True,
      min_distance_sym_equiv = 0,
      general_positions_only = False,
      min_cross_distance     = None,
      min_cubicle_edge       = 5):

    """ Run peak search on this map.
     returns group_args with:
        sites (fractional)
        sites_cart (orthogonal)
        sites_cart_absolute (orthogonal, with shift_cart applied)
        heights
        full_result (original peak_search_result object)

     Note: normally supply at least max_peaks or peak_cutoff

    position of origin of boxed map:
     origin_position = self.grid_units_to_cart(self.origin_shift_grid_units)
     shift_cart = self.shift_cart() == - origin_position

    """
    if peak_cutoff is None and max_peaks is None: # give them 1000
      max_peaks = 1000
    if min_cross_distance is None: # use half resolution
      min_cross_distance = 0.5 * self.resolution()
    if peak_search_level is None:  # this is how finely to search 1 to 3
      peak_searchlevel = 3


    map_data = self.map_data()
    cg = maptbx.crystal_gridding(
      space_group_info = self.crystal_symmetry().space_group_info(),
      symmetry_flags   = maptbx.use_space_group_symmetry,
      unit_cell        = self.crystal_symmetry().unit_cell(),
      pre_determined_n_real = map_data.all())

    # Set parameters for peak peaking and find peaks
    cgt = maptbx.crystal_gridding_tags(gridding = cg)
    peak_search_parameters = maptbx.peak_search_parameters(
      peak_search_level = peak_search_level,
      max_peaks = max_peaks,
      peak_cutoff = peak_cutoff,
      interpolate            = interpolate,
      min_distance_sym_equiv = min_distance_sym_equiv,
      general_positions_only = general_positions_only,
      min_cross_distance     = min_cross_distance,
      min_cubicle_edge       = min_cubicle_edge)
    psr = cgt.peak_search(
      parameters = peak_search_parameters,
      map        = map_data).all(max_clusters = 99999999)
    sites_cart = self.crystal_symmetry().unit_cell().orthogonalize(psr.sites())
    sites_cart_absolute = self.sites_cart_to_sites_cart_absolute(sites_cart)
    result = group_args(group_args_type = 'peak search result',
      full_result = psr,
      heights = psr.heights(),
      sites = psr.sites(),
      sites_cart = sites_cart,
      sites_cart_absolute = sites_cart_absolute,
     )

    return result

  def find_n_highest_grid_points_as_sites_cart(self, n = 0,
    n_tolerance = 0, max_tries = 100):
    '''
      Return the n highest grid points in the map as sites_cart
    '''

    # Find threshold to get exactly n points
    low_bounds = 0.
    high_bounds = 20
    mm = self.deep_copy()
    mm.set_mean_zero_sd_one() # avoid altering the working map
    tries = 0

    # Check ends
    count_high = (mm.map_data() >= high_bounds).count(True)
    count_low = (mm.map_data() >=  low_bounds).count(True)
    if count_low < n or count_high > n:
      return flex.vec3_double()

    last_threshold = None
    while tries < max_tries:
      tries += 1
      threshold = 0.5 * (low_bounds + high_bounds)
      count = (mm.map_data() >= threshold ).count(True)
      if count == n or low_bounds == high_bounds or threshold == last_threshold:
        break
      elif count > n:
        low_bounds = max(low_bounds, threshold)
      else:
        high_bounds = min(high_bounds, threshold)
      last_threshold = threshold
    if abs (count - n ) > n_tolerance:
      return flex.vec3_double()
    # Now convert to xyz and we are done
    sel = (mm.map_data() >= threshold )
    from scitbx.array_family.flex import grid
    g = grid(mm.map_data().all())
    mask_data = flex.int(mm.map_data().size(),0)
    mask_data.reshape(g)
    mask_data.set_selected(sel,1)
    mask_data.set_selected(~sel,0)

    volume_list = flex.int((sel.count(False),sel.count(True)))
    sampling_rates = flex.int((1,1))
    from cctbx.maptbx import sample_all_mask_regions
    sample_regs_obj = maptbx.sample_all_mask_regions(
      mask = mask_data,
      volumes = volume_list,
      sampling_rates = sampling_rates,
      unit_cell = mm.crystal_symmetry().unit_cell())

    return sample_regs_obj.get_array(1)

  def trace_atoms_in_map(self,
       dist_min,
       n_atoms):
     '''
       Utility to find positions where n_atoms atoms separated by
       dist_min can be placed in density in this map
     '''
     assert self.origin_is_zero()
     assert dist_min > 0.01
     assert n_atoms > 0
     n_real = self.get_n_real_for_grid_spacing(grid_spacing = dist_min)
     # temporarily remove origin shift information so we can resample
     origin_shift_grid_units_sav = tuple(self.origin_shift_grid_units)
     if self.ncs_object():
       assert tuple(self.ncs_object().shift_cart()) == tuple(
         self.shift_cart())
       self.ncs_object().set_shift_cart((0,0,0))
     self.origin_shift_grid_units = (0,0,0)
     working_map_manager = self.resample_on_different_grid(n_real = n_real)
     self.origin_shift_grid_units = origin_shift_grid_units_sav
     if self.ncs_object():
       self.ncs_object().set_shift_cart(self.shift_cart())
     return working_map_manager.find_n_highest_grid_points_as_sites_cart(
          n = n_atoms)

  def map_as_fourier_coefficients(self, d_min = None, d_max = None, box=True,
     resolution_factor=1./3, include_000 = True):
    '''
       Convert a map to Fourier coefficients to a resolution of d_min,
       if d_min is provided, otherwise box full of map coefficients
       will be created.

       Filter results with low resolution of d_max if provided

       NOTE: Fourier coefficients are relative the working origin (not
       original origin).  A map calculated from the Fourier coefficients will
       superimpose on the working (current map) without origin shifts.

       This method and fourier_coefficients_as_map_manager interconvert
       map_data and
       map_coefficients without changing origin.  Both are intended for use
       with map_data that has an origin at (0, 0, 0).

       The map coefficients are always in space group P1.
    '''
    assert self.map_data()
    assert self.map_data().origin() == (0, 0, 0)
    # Choose d_min and make sure it is bigger than smallest allowed
    if(d_min is None and not box):
      d_min = maptbx.d_min_from_map(
        map_data  = self.map_data(),
        unit_cell = self.crystal_symmetry().unit_cell(),
        resolution_factor = resolution_factor)
      print("\nResolution of map coefficients using "+\
        "resolution_factor of %.2f: %.1f A\n" %(resolution_factor, d_min),
        file=self.log)
    elif (not box):  # make sure d_min is big enough
      d_min_allowed = maptbx.d_min_from_map(
        map_data  = self.map_data(),
        unit_cell = self.crystal_symmetry().unit_cell(),
        resolution_factor = 0.5)
      if d_min < d_min_allowed:
        print("\nResolution of map coefficients allowed by gridding is %.3f " %(
          d_min_allowed),file=self.log)
        d_min=d_min_allowed
    from cctbx import crystal
    crystal_symmetry = crystal.symmetry(
      self.crystal_symmetry().unit_cell().parameters(), 1)
    ma = miller.structure_factor_box_from_map(
      crystal_symmetry = crystal_symmetry,
      include_000      = include_000,
      map              = self.map_data(),
      d_min            = d_min)
    if(d_max is not None):
      ma = ma.resolution_filter(d_max = d_max)
    return ma

  def fourier_coefficients_as_map_manager(self, map_coeffs):
    '''
       Convert Fourier coefficients into to a real-space map with gridding
        matching this existing map_manager.  Returns a map_manager object.

       Requires that this map_manager has origin at (0, 0, 0) (i.e.,
       shift_origin() has been applied if necessary)

       NOTE: Assumes that the map_coeffs are in the same frame of reference
       as this map_manager (i.e., similar to those that would be written out
       using map_as_fourier_coefficients).
    '''

    assert isinstance(map_coeffs, miller.array)
    assert isinstance(map_coeffs.data(), flex.complex_double)
    assert self.map_data() and self.map_data().origin() == (0, 0, 0)

    return self.customized_copy(
      map_data=maptbx.map_coefficients_to_map(
        map_coeffs       = map_coeffs,
        crystal_symmetry = map_coeffs.crystal_symmetry(),
        n_real           = self.map_data().all())
      )

  def shift_aware_rt(self,
     from_obj = None,
     to_obj = None,
     working_rt_info = None,
     absolute_rt_info = None):
   '''
   Returns shift_aware_rt object

   Uses rt_info objects (group_args with members of r, t).

   Simplifies keeping track of rotation/translation between two
    objects that each may have an offset from absolute coordinates.

   absolute rt is rotation/translation when everything is in original,
      absolute Cartesian coordinates.

   working_rt is rotation/translation of anything in "from_obj" object
      to anything in "to_obj" object using working coordinates in each.

   Usage:
   shift_aware_rt = self.shift_aware_rt(absolute_rt_info = rt_info)
   shift_aware_rt = self.shift_aware_rt(working_rt_info = rt_info,
      from_obj=from_obj, to_obj = to_obj)

   apply RT using working coordinates in objects
   sites_cart_to_obj = shift_aware_rt.apply_rt(sites_cart_from_obj,
      from_obj=from_obj, to_obj=to_obj)

   apply RT absolute coordinates
   sites_cart_to = shift_aware_rt.apply_rt(sites_cart_from)

   '''
   return shift_aware_rt(
     from_obj = from_obj,
     to_obj = to_obj,
     working_rt_info = working_rt_info,
     absolute_rt_info = absolute_rt_info)

  def _get_unique_box_info(self, box_info, max_distance = 1):
    if self.ncs_object() is None:
      # try to get map symmetry but do not try too hard..
      try:
        self.find_map_symmetry()
      except Exception as e:
        pass
    if not self.ncs_object() or self.ncs_object().max_operators()<2:
      return box_info # nothing to do

    box_info.ncs_object = self.ncs_object() # save it

    # Get just the unique parts of this box (apply symmetry later)
    new_lower_bounds_list = []
    new_upper_bounds_list = []
    new_lower_bounds_with_cushion_list = []
    new_upper_bounds_with_cushion_list = []
    existing_xyz_list = flex.vec3_double()
    existing_unique_xyz_list = flex.vec3_double()
    from scitbx.matrix import col
    for lower_bounds, upper_bounds,lower_bounds_with_cushion, \
      upper_bounds_with_cushion in zip (
        box_info.lower_bounds_list,
        box_info.upper_bounds_list,
        box_info.lower_bounds_with_cushion_list,
        box_info.upper_bounds_with_cushion_list,
      ):
      # NOTE: lower_bounds, upper_bounds are relative to the working
      #    map_data with origin at (0,0,0).  Our ncs_object is also
      #    relative to this same origin

      xyz = tuple([ a * 0.5*(lb+ub-1) / n for a, lb, ub, n in zip(
         self.crystal_symmetry().unit_cell().parameters()[:3],
         lower_bounds,
         upper_bounds,
         self.map_data().all())])
      target_site = flex.vec3_double((xyz,))
      ncs_object = self.ncs_object()
      if existing_xyz_list.size() > 0 :
       dist_n, id1_n, id2_n = target_site.min_distance_between_any_pair_with_id(
              existing_xyz_list)
      else:
        dist_n = 1.e+30
      if dist_n <= max_distance:  # duplicate
        pass
      else:
        ncs_sites = ncs_object.apply_ncs_to_sites( sites_cart=target_site)
        existing_xyz_list.extend(ncs_sites)
        existing_unique_xyz_list.extend(
          flex.vec3_double((xyz,)*ncs_sites.size()))
        new_lower_bounds_list.append(lower_bounds)
        new_upper_bounds_list.append(upper_bounds)
        new_lower_bounds_with_cushion_list.append(lower_bounds_with_cushion)
        new_upper_bounds_with_cushion_list.append(upper_bounds_with_cushion)

    box_info.lower_bounds_list = new_lower_bounds_list
    box_info.upper_bounds_list = new_upper_bounds_list
    box_info.lower_bounds_with_cushion_list = new_lower_bounds_with_cushion_list
    box_info.upper_bounds_with_cushion_list = new_upper_bounds_with_cushion_list

    return box_info

#   Methods for map_manager

class shift_aware_rt:
  '''
  Class to simplify keeping track of rotation/translation between two
  objects that each may have an offset from absolute coordinates.

  Basic idea:  absolute rt is rotation/translation when everything is in
  original, absolute Cartesian coordinates.

  working_rt is rotation/translation of anything in "from_obj" object to anything
   in "to_obj" object using working coordinates in each.

  The from_obj and to objects must have a shift_cart method
  '''

  def __init__(self,
     from_obj = None,
     to_obj = None,
     working_rt_info = None,
     absolute_rt_info = None):

     assert (
      (absolute_rt_info and (not from_obj) and (not to_obj) and (not working_rt_info))
      or
      (from_obj and to_obj and working_rt_info))

     if from_obj:
       assert hasattr(from_obj, 'shift_cart')
     if to_obj:
       assert hasattr(to_obj, 'shift_cart')

     if not absolute_rt_info:
       absolute_rt_info = self.get_absolute_rt_info(
         working_rt_info = working_rt_info,
         from_obj = from_obj, to_obj = to_obj)

     self._absolute_rt_info = group_args(
        r =  absolute_rt_info.r,
        t =  absolute_rt_info.t,)


  def is_similar(self, other_shift_aware_rt_info, tol = 0.001):
    '''Check whether this shift_aware_rt is similar to another one'''
    r = self._absolute_rt_info.r
    t = self._absolute_rt_info.t
    other_r = other_shift_aware_rt_info._absolute_rt_info.r
    other_t = other_shift_aware_rt_info._absolute_rt_info.t
    for x,y in zip(r,other_r):
      if (abs(x-y)) > tol:
        print(x,y,abs(x-y))
        return False
    for x,y in zip(t,other_t):
      if (abs(x-y)) > tol:
        print(x,y,abs(x-y))
        return False
    return True

  def apply_rt(self, site_cart = None, sites_cart = None,
    from_obj = None, to_obj = None):
    '''
    Apply absolute rt if from and to not specified.
    Apply relative if specified
    '''
    # get absolute if from and to not specified, otherwise working
    rt_info = self.working_rt_info(from_obj=from_obj, to_obj=to_obj)
    if site_cart:
      return rt_info.r * col(site_cart) + rt_info.t

    else:
      return rt_info.r.elems * sites_cart + rt_info.t.elems

  def get_absolute_rt_info(self, working_rt_info = None,
      from_obj = None, to_obj = None):

    '''
    working_rt_info describes how to map from_xyz -> to_xyz in local coordinates
    from_xyz is shifted from absolute by from.shift_cart()
    to_xyz is shifted from absolute by to.shift_cart()

    We have:
      r from_xyz + t = to_xyz    in working frame of reference

    We want to describe how to map:
       (from_xyz - from.shift_cart()) -> (to_xyz - to.shift_cart())
    where r is going to be the same and T will be different than t
       r ((from_xyz - from.shift_cart()) + T = (to_xyz - to.shift_cart())
       T = (to_xyz - to.shift_cart() - r from_xyz + r from.shift_cart()
         but: to_xyz -  r from_xyz = t
       T =  t - to.shift_cart() + r from.shift_cart()

    Note reverse:
       t = T + to.shift_cart() - r from.shift_cart()
    '''

    r = working_rt_info.r
    t = working_rt_info.t
    new_t =  t -  col(to_obj.shift_cart())  + r * col(from_obj.shift_cart())

    return group_args(
      r = r,
      t = new_t
    )

  def working_rt_info(self, from_obj=None, to_obj=None):
    ''' Get rt in working frame of reference
    '''
    if (not from_obj) and (not to_obj):  # as is
      return self._absolute_rt_info

    assert hasattr(from_obj, 'shift_cart')
    assert hasattr(to_obj, 'shift_cart')
    r = self._absolute_rt_info.r
    t = self._absolute_rt_info.t
    working_t =  t + col(to_obj.shift_cart()) - r * col(from_obj.shift_cart())
    return group_args(
      r = r,
      t = working_t)


  def absolute_rt_info(self):
    '''Return the absolute RT info for this shift_aware_rt object'''
    return self._absolute_rt_info


  def inverse(self):
    '''Return the inverse for this shift_aware_rt object'''
    r = self._absolute_rt_info.r
    t = self._absolute_rt_info.t

    r_inv = r.inverse()
    t_inv = - r_inv * t
    inverse_absolute_rt_info = group_args(
      r = r_inv,
      t = t_inv,)

    return shift_aware_rt(absolute_rt_info = inverse_absolute_rt_info)


def dummy_map_manager(crystal_symmetry, n_grid = 12):
  '''
   Make a map manager with crystal symmetry and unit sized map
  '''

  map_data = flex.double(n_grid*n_grid*n_grid,1)
  acc = flex.grid((n_grid, n_grid, n_grid))
  map_data.reshape(acc)
  mm = map_manager(
    map_data = map_data,
    unit_cell_grid = (n_grid, n_grid, n_grid),
    unit_cell_crystal_symmetry = crystal_symmetry,
    wrapping = False)
  mm.set_resolution(min(crystal_symmetry.unit_cell().parameters()[:3])/n_grid)
  mm._is_dummy_map_manager = True
  return mm


def get_indices_from_index(index = None, all = None):
        '''Get indices (in a 3D map) for a grid point with given 1D index'''
        #index = k+j*all[2]+i*(all[1]*all[2])
        i = index//(all[1]*all[2])
        j =  (index-i*(all[1]*all[2]))//all[2]
        k =  index-i*(all[1]*all[2])-j*all[2]
        assert k+j*all[2]+i*(all[1]*all[2]) == index
        return (i, j, k)

def get_sites_cart_from_index(
      indices_list = None,
      points = None, map_data = None, crystal_symmetry = None, all = None):
    '''  Get sites_cart from linear (1d) map indices.
       Supply either map_data or all to provide n_real
       crystal_symmetry is required
       Supply either 3D indices (i,j,k) or 1-D indices (points)
    '''

    if all is None:
      all = map_data.all()
    sites_frac = flex.vec3_double()
    if not indices_list:
      if not points: return sites_frac # nothing there
      indices_list = []
      for point in points:
        if point is None: continue
        indices_list.append(get_indices_from_index(index = point, all = all))
    for indices in indices_list:
      i, j, k = indices
      site_frac = tuple((i/all[0], j/all[1], k/all[2]))
      sites_frac.append(site_frac)
    sites_cart = crystal_symmetry.unit_cell().orthogonalize(sites_frac)
    return sites_cart

def _round_tuple_int(t):
  new_t = []
  for x in t:
    new_t.append(int(round(x)))
  return new_t

def add_tuples_int(t1, t2):
  ''' Add two tuples (can be integers or floats)'''
  try:
    return tuple(flex.int(t1)+flex.int(t2))
  except Exception as e: # not integers
    return tuple(
       flex.int(_round_tuple_int(t1)) + flex.int(_round_tuple_int(t2)))

def subtract_tuples_int(t1, t2):
  try:
    return tuple(flex.int(t1)-flex.int(t2))
  except Exception as e: # not integers
    return tuple(
       flex.int(_round_tuple_int(t1)) - flex.int(_round_tuple_int(t2)))

def remove_site_with_most_neighbors(sites_cart):
  '''Remove the site with the most neighbors'''
  useful_norms_list = []
  closest_distance = 1.e+30
  for i in range(sites_cart.size()):
    compare_xyz = flex.vec3_double(sites_cart.size(), sites_cart[i])
    delta_xyz = sites_cart - compare_xyz
    norms = delta_xyz.norms()
    useful_norms = norms[:i]
    useful_norms.extend(norms[i+1:])
    assert useful_norms.size() == sites_cart.size() -1
    useful_norms_list.append(useful_norms)
    closest_distance=min(closest_distance,useful_norms.min_max_mean().min)

  distance_list=[]
  for i in range(sites_cart.size()):
    useful_norms = useful_norms_list[i]
    s = (useful_norms <= closest_distance * 1.25)
    count = s.count(True)
    distance_list.append([count,i])
  distance_list.sort()
  distance_list.reverse()
  i = distance_list[0][1]
  new_sites_cart = sites_cart[:i]
  new_sites_cart.extend(sites_cart[i+1:])
  return new_sites_cart

def select_n_in_biggest_cluster(sites_cart,
   dist_min = None,
   n = None,
   dist_min_ratio = 1.,
   dist_min_ratio_min = 0.5,
   minimize_density_of_points = None):
  '''
    Select n of sites_cart, taking those near biggest cluster if possible
    If minimize_density_of_points, remove those with the most neighbors
  '''

  if sites_cart.size() < 1:
    return sites_cart

  if minimize_density_of_points:
    while sites_cart.size() > n:
      sites_cart = remove_site_with_most_neighbors(sites_cart)
    return sites_cart

  # Guess size of cluster (n atoms, separated by about dist_min)
  target_radius = dist_min * float(n)**0.5
  dist_list = []
  for i in range (sites_cart.size()):
    diffs = sites_cart.deep_copy() - col(sites_cart[i])
    norms = diffs.norms()
    sel = (norms <= target_radius)
    dist_list.append([sel.count(True),i])
  dist_list.sort()
  dist_list.reverse()
  i = dist_list[0][1]

  # Now take the n points close to center_point but separated from
  #  each other and we are done
  diffs = sites_cart.deep_copy() - col(sites_cart[i])
  norms = diffs.norms()  # how close each one is to the center point
  used_sites=flex.bool(sites_cart.size(), False)

  new_sites_cart = flex.vec3_double()
  unused_sites_cart = flex.vec3_double()
  for j in range(n):  # pick closest to center_point that is at least
                      # dist_min from all in new_sites_cart
    found = False
    for k in range(n):
      if found: break # go on to next
      if used_sites[k]: continue
      ok = False
      test_sites = sites_cart[k:k+1]
      if new_sites_cart.size() == 0:
        ok = True
      else:
        dist, id1, id2= new_sites_cart.min_distance_between_any_pair_with_id(
            test_sites)
        if dist >= dist_min_ratio*dist_min: # keep it
           ok = True
      if ok:
        new_sites_cart.append(test_sites[0])
        used_sites[k] = True
        found = True # go on to next one
    if not found:  # didn't get anything ... reduce dist_min_ratio
      if dist_min_ratio >= dist_min_ratio_min:
        return select_n_in_biggest_cluster(sites_cart,
          dist_min = dist_min,
          n = n,
          dist_min_ratio = dist_min_ratio * 0.9)

  return new_sites_cart
