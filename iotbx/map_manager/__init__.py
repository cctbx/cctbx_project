from __future__ import absolute_import, division, print_function
from libtbx.utils import Sorry
import sys
from libtbx import group_args
import iotbx.mtz
from iotbx.mrcfile import map_reader, write_ccp4_map
import os
from scitbx.array_family import flex
from cctbx import maptbx
from cctbx import miller

class map_manager(map_reader,write_ccp4_map):

  '''
   map_manager, includes map_reader and write_ccp4_map

   Normally use map_manager to read, write, and carry information about
   maps.  Map_manager keeps track of the origin shifts and also the
   original full unit cell and cell dimensions.  It writes out the map
   in the same place as it was read in.

   Map_manager also has high-level tools to manipulate maps and to get
   lower-level objects.

   You can create a new map_manager by initializing map_manager with an
   existing map_manager and a new map_data object. The new map_manager
   will contain the original origin shifts and full unit cell, etc, so
   the new map_data object you supply should match the existing one (i.e.,
   both origin-shifted in the same way).

   You can also get a new_map_manager with the customized_copy() method.

   NOTE: MRC Maps may not represent the entire unit cell.  Normally maps that
    have an origin (corner with minimum i,j,k) that is not zero will be
    shifted at a later stage to have the origin at (0,0,0), along with
    any models and ncs objects (typically done with iotbx.map_and_model).
    To be able to write out a map in the same place as it was read in
    after shifting the origin and/or boxing the map, you need to keep track
    of 5 things.  These are:
    1. unit_cell_grid: grid representing one full unit cell as read in.
        Available from map_manager as self.unit_cell_grid
    2. unit_cell_parameters: dimensions of one full unit cell
    3. origin_shift_grid_units: the shift in grid units to apply to the
       working map to superimpose it on the original map. When you read the
       map in this is (0,0,0). If you shift the map origin from (i,j,k) to
       (0,0,0) then the origin_shift_grid_units is (i,j,k).
       The shift you can get from map_and_model if you use it,
          is origin_shift=map_and_model.original_origin_grid_units()

    Also (4) you will want to note the gridding of the part of the map that is
      present (nx,ny,nz), available from the map_manager
      as working_map_n_xyz=self.get_working_map_n_xyz()
    And (5) you want the space-group number in case it is not 1

   Normal usage:

     Read in a map:
       mm=map_manager('input_map.mrc')
     Summarize:
       mm.show_summary()

     Optionally shift origin of map to (0,0,0) (you can do this here
         or you can use iotbx.map_and_model to shift models and maps together):
       mm.shift_origin()

     Get the map_data (shifted if origin was shifted above):
       map_data=mm.map_data()

     Get the crystal_symmetry of the box of data that is present:
       cs=mm.crystal_symmetry()

     Get the crystal_symmetry of the whole unit cell (even if not present):
       unit_cell_cs=mm.unit_cell_crystal_symmetry()

     Optional:
       Get a map_manager object that just contains info on how to shift. You
         can do anything with the map_shift_tracker you can do with the
         original map_manage except you supply the map_data separately:
       map_shift_tracker=mm.map_shift_tracker()
         NOTE: if you shift the map AFTER setting up map_shift_tracker use:
          shifted_map=map_shift_tracker.shift_origin(map_data,update_shift=True)
         so that map_shift_tracker knows about the shift.  If you shift
         another map with same map_shift_tracker, use update_shift=False

     Write out the map:
       mm.write_map(file_name='output_map.ccp4')

     Or write out the map with map_shift_tracker:
       map_shift_tracker.write_map(file_name='output_map.ccp4',
         map_data=map_data)

   See http://www.ccpem.ac.uk/mrc_format/mrc2014.php for MRC format
   See https://pypi.org/project/mrcfile/ for mrcfile library documentation

   Same conventions as iotbx.ccp4_map

   Default is to write maps with INTERNAL_STANDARD_ORDER of axes of [3,2,1]
     corresponding to columns in Z, rows in Y, sections in X to match
     flex array layout.  This can be modified by changing the values in
     output_axis_order.

   Hard-wired to convert input maps of any order to
     INTERNAL_STANDARD_ORDER = [3,2,1] before conversion to flex arrays
     This is not modifiable.

    INTERNAL_STANDARD_ORDER=[3,2,1]

  Standard limitations and associated message.
  These can be checked with: limitations=mrc.get_limitations()
    which returns a group_args object with a list of limitations and a list
    of corresponding error messages, or None if there are none
    see phenix.show_map_info for example
  These limitations can also be accessed with specific calls placed below:
   For example:
   mrc.can_be_sharpened()  returns False if "extract_unique" is present

  Map labels that are not limitations can be accessed with:
      additional_labels=mrc.get_additional_labels()

  STANDARD_LIMITATIONS_DICT={
    "extract_unique":
     "This map is masked around unique region and not suitable for auto-sharpening.",
    "map_is_sharpened":
     "This map is auto-sharpened and not suitable for further auto-sharpening.",
    "map_is_density_modified": "This map has been density modified.",
     }


   NOTES ON ORDER OF AXES

    Phenix standard order is 3,2,1 (columns Z, rows Y, sections in X).
        Convert everything to this order.

    This is the order that allows direct conversion of a numpy 3D array
     with axis order (mapc,mapr,maps) to a flex array.

    For reverse=True, supply order that converts flex array to numpy 3D array
     with order (mapc,mapr,maps)

    Note that this does not mean input or output maps have to be in this order.
     It just means that before conversion of numpy to flex or vice-versa
     the array has to be in this order.

     Note that MRC standard order for input/ouput is 1,2,3.

     NOTE: numpy arrays indexed from 0 so this is equivalent to
      order of 2,1,0 in the numpy array

    NOTE:  MRC format allows data axes to be swapped using the header
      mapc mapr and maps fields. However the mrcfile library does not
      attempt to swap the axes and assigns the columns to X, rows to Y and
      sections to Z. The data array is indexed C-style, so data values can
      be accessed using mrc.data[z][y][x].

    NOTE: normal expectation is that phenix will read/write with the
      order 3,2,1. This means X-sections (index=3), Y rows (index=2),
      Z columns (index=1). This correxponds to
       mapc (columns) =   3 or Z
       mapr (rows)    =   2 or Y
       maps (sections) =  1 or X

    In the numpy array (2,1,0 instead of 3,2,1):

    To transpose, specify i0,i1,i2 where:
        i0=2 means input axis 0 becomes output axis 2
        NOTE:  axes are 0,1,2 etc, not 1,2,3
      Examples:
        np.transpose(a,(0,1,2))  does nothing
        np.transpose(a,(1,2,0)) output axis 0 is input axis 1



    We want output axes to always be 2,1,0 and input axes for numpy array are
      (mapc-1,mapr-1,maps-1):

    For example, in typical phenix usage, the transposition is:
      i_mapc=3    i_mapc_np=2
      i_mapr=2    i_mapr_np=1
      i_maps=1    i_maps_np=0
  '''


  def __init__(self,
     file_name=None,  # Normally initialize by reading from a file
     data_manager=None, # Can associate with a data_manager XXX Not implemented
     map_manager_object=None, # Also can initialize with existing map_manager
     map_data=None,    # and optional map_data
     unit_cell_grid=None,  # Optional specification of unit cell, space group
     unit_cell_parameters=None,    # and origin shift instead of map_manager
     space_group_number=None,
     origin_grid_units=None,
     working_map_n_xyz=None,
     high_resolution=None,  # optional resolution used to create map
     log=sys.stdout,
     ):

    '''
      Allows simple version of reading a map file or blank initialization

      Can call with file_name to read map file.

      Alternative is initialize with existing map_manager object and
        optional map_data, or blank initialization

      Final alternative is to specify unit_cell_grid, unit_cell_parameters,
        space_group_number, origin_grid_units and
        map_data or working_map_n_xyz

      NOTE: Different from map_reader, map_manager ALWAYS converts
      map_data to flex.double and does not save any extra information about
      the map except the details specified in this __init__.

      After reading you can access map data with self.map_data()
        and other attributes (see class utils in ccp4_map/__init__py)
    '''

    # Initialize origin shift representing position of original origin in
    #  grid units.  If map is shifted, this is updated to reflect where
    #  to place current origin to superimpose map on original map.
    self.log = log
    self.origin_shift_grid_units=(0,0,0)
    self.data_manager=None

    # Usual initialization with a file
    if file_name is not None:
      self.read_map(file_name=file_name,log=log)
      self.data_manager=data_manager

    # Initialization with map_manager and map_data objects
    elif map_manager_object:
      self.initialize_with_map_manager(
        map_manager_object=map_manager_object,
        map_data=map_data,
        data_manager=data_manager,
        high_resolution=high_resolution,log=log)
    elif unit_cell_grid and unit_cell_parameters and \
          (space_group_number is not None) and \
          (map_data or working_map_n_xyz):
      if origin_grid_units is None:
        origin_grid_units=(0,0,0)
      manager_object=group_args(
        unit_cell_grid=unit_cell_grid,
        unit_cell_parameters=unit_cell_parameters,
        space_group_number=space_group_number,
        origin_shift_grid_units=origin_grid_units,
        _map_data=None,
        working_map_n_xyz=working_map_n_xyz)

      self.initialize_with_map_manager(
        map_manager_object=manager_object,
        map_data=map_data,
        data_manager=data_manager,
        high_resolution=high_resolution,log=log)

    else:
      # Initialize with nothing except data_manager
     print("Initializing empty map_manager ",file=log)
     self.unit_cell_grid=None
     self.unit_cell_parameters=None
     self.space_group_number=None
     self.origin_shift_grid_units=None
     self.working_map_n_xyz=None
     self._map_data=None
     self.data=None
     self._high_resolution=None
     self.data_manager=data_manager

  def initialize_with_map_manager(self,
       map_manager_object=None,
       data_manager=None,
       map_data=None,
       high_resolution=None,log=sys.stdout):

     mmo=map_manager_object

     self.unit_cell_grid= mmo.unit_cell_grid
     self.unit_cell_parameters= mmo.unit_cell_parameters
     self.space_group_number= mmo.space_group_number

     self.origin_shift_grid_units= mmo.origin_shift_grid_units
     self.working_map_n_xyz= mmo.working_map_n_xyz
     self._map_data=mmo._map_data
     self.data=None
     self._high_resolution=high_resolution

     if data_manager: # replace data_manager
       self.data_manager=data_manager

     if map_data:
       self.add_map_data(map_data,overwrite_map_gridding=True,log=log)

  def read_map(self,file_name=None,log=sys.stdout):
      print("Reading map from %s " %(file_name),file=log)
      if hasattr(self,'_map_data') and self._map_data:
        print("NOTE: Removing existing map_data",file=log)

      self.read_map_file(file_name=file_name)  # in mrcfile/__init__.py
      # Save gridding, convert to double and delete original data
      self.working_map_n_xyz=self.get_working_map_n_xyz()
      self.convert_to_double()

  def add_map_data(self,map_data,overwrite_map_gridding=None,log=sys.stdout):
       # NOTE: Used to initialize, replacing existing map data
       self.data=map_data
       self._map_data=None
       self.convert_to_double()
       if self._map_data.all() != self.working_map_n_xyz:
         if overwrite_map_gridding:
           print("Note: Map gridding (%s) " %(str(self._map_data.all())),
            " is changed from input map_manager (%s)" %(
              str(self.working_map_n_xyz)),"Assuming map has been boxed",
              file=log)
           self.working_map_n_xyz=self._map_data.all()
         else:
           raise Sorry("Map gridding (%s) " %(str(self._map_data.all()) +
            " does not match gridding specified by input map_manager (%s)" %(
              str(self.working_map_n_xyz))))

  def _print(self, m):
    if(self.log is not None and not self.log.closed):
      print(m, file=self.log)

  def read_mtz(self, file_name):
    assert os.path.isfile(file_name)
    # read in map stored as mtz dataset and convert to map and add resulting
    #  map_data.
    #  NOTE: use data_manager for this
    # Requires that self.working_map_n_xyz is defined
    assert self.working_map_n_xyz is not None
    mtz_object = iotbx.mtz.object(file_name=file_name)
    miller_arrays = mtz_object.as_miller_arrays()
    assert len(miller_arrays) == 1
    map_coeffs = mtz_object.as_miller_arrays()[0]
    assert isinstance(map_coeffs.data(), flex.complex_double)
    self._print(
      "Read %s map coefficients from %s" %(map_coeffs.size(), file_name))
    # convert to map
    map_data=self.fourier_coefficients_as_map(map_coeffs=map_coeffs)
    self.add_map_data(map_data, log=self.log)
    if self.map_data():
      self._print(
        "Converted map coefficients to map_data and replaced existing map_data")
    else:
      self._print(
        "Converted map coefficients to map_data and saved new map_data")
    self._print("Returning map coefficients")
    return map_coeffs

  def set_origin_and_gridding(self,origin_grid_units,
      gridding=None,
      allow_non_zero_previous_value=False,log=sys.stdout):
    if not allow_non_zero_previous_value:
      if (not self.origin_shift_grid_units in [None,(0,0,0)]):
              # make sure we
              #   didn't already shift it
        print("Origin shift is already set...cannot change it unless you"+
          " also specify allow_non_zero_previous_value=True",file=log)
        return

    if self.origin_shift_grid_units in [None,(0,0,0)]:
      self.origin_shift_grid_units=origin_grid_units
    else:
      new_origin=[]
      for a,b in zip(self.origin_shift_grid_units,origin_grid_units):
        new_origin.append(a+b)
      self.origin_shift_grid_units=new_origin

    if gridding: # reset definition of full unit cell.  Keep grid spacing
       current_unit_cell_parameters=self.unit_cell_parameters
       current_unit_cell_grid=self.unit_cell_grid
       new_unit_cell_parameters=[]
       for a,g,new_g in zip(
          current_unit_cell_parameters[:3],current_unit_cell_grid,gridding):
         new_a=a*new_g/g
         new_unit_cell_parameters.append(new_a)

       self.unit_cell_parameters=\
          new_unit_cell_parameters+list(current_unit_cell_parameters[3:])
       self.unit_cell_grid=gridding
       print ("Resetting gridding of full unit cell from %s to %s" %(
         str(current_unit_cell_grid),str(gridding)),file=log)
       print ("Resetting dimensions of full unit cell from %s to %s" %(
         str(current_unit_cell_parameters),
            str(new_unit_cell_parameters)),file=log)

  def already_shifted(self):
    # Check if origin is already at (0,0,0)
    if self._map_data is not None and self._map_data.origin() == (0,0,0):
      return True
    else:
      return False

  def shift_origin(self):
    '''
    Shift the origin of the map to (0,0,0) and update origin_shift_grid_units
    '''
    if(self._map_data is None): return
    if(self.origin_shift_grid_units is None):
      self.origin_shift_grid_units=(0,0,0)
    new_origin_shift = self._map_data.origin()
    if new_origin_shift != (0,0,0):
      self._map_data = self._map_data.shift_origin()
      full_shift=[]
      for a,b in zip(self.origin_shift_grid_units,new_origin_shift):
        full_shift.append(a+b)
      self.origin_shift_grid_units = full_shift

  def map_shift_tracker(self):
    # Produces map_manager with no data so it is small but still knows how
    #   to write out maps
    # PVA: Should it be a shift manager (for clarity)?
    mst=map_manager(map_manager_object=self)
    mst._map_data=None
    from copy import deepcopy
    return deepcopy(mst)

  def write_map(self,
     file_name=None,
     map_data=None,  # map_data to be written out
     crystal_symmetry=None,  # optional crystal_symmetry of region present
     unit_cell_grid=None,  # optional gridding of full unit cell
     origin_shift_grid_units=None, # optional origin shift (grid units) to apply
     labels=None,  # optional list of output labels
     verbose=None,
     use_data_manager_if_available=True,
     log=sys.stdout,
     ):

    '''
      Simple version of write

      file_name is output file name
      map_data is map_data object with 3D values for map. If not supplied,
        use self.map_data()

      Normally only file_name and map_data needs to be supplied. Optional below:
      crystal_symmetry is crystal.symmetry object describing unit cell and
        space group. If not supplied, use self.unit_cell_crystal_symmetry()
      unit_cell_grid is gridding of full unit cell, normally self.unit_cell_grid
      origin_shift_grid_units is optional origin shift (grid units) to apply
      labels are optional list of strings describing the map

      If data_manager is available and use_data_manager_if_available, do the
      write through the data_manager and only file_name and map_data are used
    '''

    if self.data_manager and use_data_manager_if_available:
      if map_data is None:
        mm=self
      else:
        mm=self.customized_copy(map_data=map_data)
      self.data_manager.write_map_with_map_manager(mm,filename=file_name)
      return

    if not file_name:
      raise Sorry("Need file_name for write_map")

    if labels and (type(labels) != type([1,2,3]) or \
         type(labels[0])!=type("ABC")):
      raise Sorry("Output labels must be a list of text strings")

    if not map_data:
      if not self.map_data():
        raise Sorry("Need map_data for write_map")
      # Take map_data from self if not supplied
      map_data=self.map_data()

    # Take other values from self unless supplied
    if not crystal_symmetry: # Use crystal symmetry of entire cell
      crystal_symmetry=self.unit_cell_crystal_symmetry()
    if not unit_cell_grid: # and gridding of entire cell
      unit_cell_grid=self.unit_cell_grid
    if not origin_shift_grid_units: # and origin shift to apply
      origin_shift_grid_units=self.origin_shift_grid_units
    # NOTE crystal_symmetry and unit_cell_grid are for entire cell

    if map_data.origin() == (0,0,0):  # Usual
      print("Writing map with origin at %s and size of %s to %s" %(
        str(origin_shift_grid_units),str(map_data.all()),file_name),file=log)
      write_ccp4_map(
        file_name   = file_name,
        crystal_symmetry = crystal_symmetry, # unit cell and space group
        map_data    = map_data,
        unit_cell_grid=unit_cell_grid,  # optional gridding of full unit cell
        origin_shift_grid_units=origin_shift_grid_units, # optional origin shift
        labels      = labels,
        verbose=verbose)
    else: # map_data has not been shifted to (0,0,0).  Shift it and then write
      print("Writing map after shifting origin",file=log)
      if self.origin_shift_grid_units and origin_shift_grid_units!=(0,0,0):
        print ("WARNING: map_data has origin at %s " %(str(map_data.origin())),
         " and this map_manager will apply additional origin shift of %s " %(
          str(self.origin_shift_grid_units)), file=log)

      new_map_manager=map_manager(map_manager_object=self,
         map_data=map_data.deep_copy())
      new_map_manager.shift_origin()
      new_map_manager.write_map(file_name=file_name,log=log)

  def deep_copy(self):
    ''' Return a deep copy of this map_manager object'''
    from copy import deepcopy
    return deepcopy(self)

  def customized_copy(self,map_data=None,
     keep_data_manager=True, log=sys.stdout):
    dc=self.deep_copy()
    dc.add_map_data(map_data,overwrite_map_gridding=False,log=log)
    if keep_data_manager:
      dc.add_data_manager(self.data_manager)
    else:
      dc.remove_data_manager()
    return dc

  def add_data_manager(self,data_manager):
    self.data_manager=data_manager

  def remove_data_manager(self):
    self.data_manager=None

  def is_similar(self,other=None,eps=0.5):
    from libtbx.test_utils import approx_equal
    # Check to make sure gridding and symmetry is similar
    if self.origin_shift_grid_units != other.origin_shift_grid_units:
      return False
    if self.space_group_number != other.space_group_number:
      return False
    if self.working_map_n_xyz != other.working_map_n_xyz:
      return False
    if self.unit_cell_grid != other.unit_cell_grid:
      return False
    if not approx_equal(tuple(self.unit_cell_parameters),
         tuple(other.unit_cell_parameters),eps=eps):
      return False
    return True


  def map_as_fourier_coefficients(self, high_resolution=None):
    '''
       Convert a map to Fourier coefficients to a resolution of high_resolution,
       if high_resolution is provided, otherwise box full of map coefficients
       will be created.
       NOTE: Fourier coefficients are relative the working origin (not
       original origin).  A map calculated from the Fourier coefficients will
       superimpose on the working (current map) without origin shifts.
    '''
    assert self.map_data()
    assert self.map_data().origin()==(0,0,0)
    return miller.structure_factor_box_from_map(
      crystal_symmetry = self.crystal_symmetry(),
      map              = self.map_data(),
      d_min            = high_resolution)

  def fourier_coefficients_as_map(self, map_coeffs):
    '''
       Convert Fourier coefficients into to a real-space map with gridding
       matching this existing map_manager.
       Requires that this map_manager has origin at (0,0,0) (i.e.,
       shift_origin())
    '''
    assert map_coeffs
    assert isinstance(map_coeffs.data(), flex.complex_double)
    assert (self.map_data() and self.map_data().origin()==(0,0,0) ) or (
       self.working_map_n_xyz)
    if self.map_data():
      n_real=self.map_data().all()
    else:
      n_real=self.working_map_n_xyz
    return maptbx.map_coefficients_to_map(
      map_coeffs       = map_coeffs,
      crystal_symmetry = self.crystal_symmetry(),
      n_real           = n_real)

  def generate_map(self,
     log=sys.stdout,
     **kw):  # Accepts keywords for generate_model, generate_map_coefficients,
             #  and generate_map and all are passed to generate_model etc.
    '''
      Simple interface to iotbx.create_models_or_maps to generate a map
      Simple use:
       mm=map_manager.generate_map()  # generates a small map (and model)
       map is in mm.map_data()
       model is in mm.model()
      Advanced use...see keywords in iotbx.create_models_or_maps.generate_map
    '''
    from iotbx.map_model_manager import map_model_manager
    mmm=map_model_manager()
    mmm.generate_map(log=log,**kw)  # generate small map
    self.__init__(map_manager_object=mmm._map_manager)  # make this map_manager the new map
