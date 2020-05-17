from __future__ import absolute_import, division, print_function
from libtbx.utils import Sorry
import sys
from libtbx import group_args

from iotbx.mrcfile import map_reader, write_ccp4_map

class map_manager(map_reader,write_ccp4_map):

  '''
   map_manager, includes map_reader and write_ccp4_map

   Normally use map_manager to read, write, and carry information about
   maps.  Map_manager keeps track of the origin shifts and also the
   original full unit cell and cell dimensions.  It writes out the map
   in the same place as it was read in.

   You can create a new map_manager by initializing map_manager with an
   existing map_manager and a new map_data object. The new map_manager
   will contain the original origin shifts and full unit cell, etc, so
   the new map_data object you supply should match the existing one (i.e.,
   both origin-shifted in the same way).

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

    All of these items are available as a single map_manager object with:
     map_shift_tracker=self.map_shift_tracker()

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

     Shift a model object to match the working origin of the map:
       shifted_model=mm.shift_model_to_match_working_map(model)

     Shift a ncs object to match the working origin of the map:
       shifted_ncs=mm.shift_ncs_to_match_working_map(ncs_object)

     NOTE: If you run map_and_model=iotbx.map_and_model separately, you can
       update your map_shift_tracker with:
         origin_shift=map_and_model.original_origin_grid_units()
         map_shift_tracker.set_origin_shift_grid_units(origin_shift)


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
     map_manager_object=None, # Also can initialize with existing map_manager
     map_data=None,    # and optional map_data
     unit_cell_grid=None,  # Optional specification of unit cell, space group
     unit_cell_parameters=None,    # and origin shift instead of map_manager
     space_group_number=None,
     origin_shift_grid_units=None,
     working_map_n_xyz=None,
     ):

    '''
      Simple version of reading a map file

      Normally just call with file_name.

      Alternative is initialize with existing map_manager object and
        optional map_data

      Final alternative is to specify unit_cell_grid, unit_cell_parameters,
        space_group_number, origin_shift_grid_units and
        map_data or working_map_n_xyz

      NOTE: Different from map_reader, map_manager ALWAYS converts
      map_data to flex.double and does not save any extra information except
      the details specified in this __init__.

      After reading you can access map data with self.map_data()
        and other attributes (see class utils in ccp4_map/__init__py)
    '''

    self.origin_shift_grid_units=(0,0,0)
    # Usual initialization with a file
    if file_name is not None:
      self.read_map(file_name=file_name)

      # Save gridding, convert to double and delete original data
      self.working_map_n_xyz=self.get_working_map_n_xyz()
      self.convert_to_double()

    # Initialization with map_manager and map_data objects
    elif map_manager_object:
      self.initialize_with_map_manager(
        map_manager_object=map_manager_object,
        map_data=map_data)
    elif unit_cell_grid and unit_cell_parameters and \
          (space_group_number is not None) and \
          (origin_shift_grid_units is not None) and \
          (map_data or working_map_n_xyz):
      manager_object=group_args(
        unit_cell_grid=unit_cell_grid,
        unit_cell_parameters=unit_cell_parameters,
        space_group_number=space_group_number,
        origin_shift_grid_units=origin_shift_grid_units,
        _map_data=None,
        working_map_n_xyz=working_map_n_xyz)

      self.initialize_with_map_manager(
        map_manager_object=manager_object,
        map_data=map_data)

    else:
      raise Sorry("You need to initialize map_manager with a file, another "+
        "map_manager, or by specifying unit_cell_grid,unit_cell_parameters, "+
        "space_group_number, origin_shift_grid_units, and "+
       "a map_data object or working_map_n_xyz")

  def initialize_with_map_manager(self,
       map_manager_object=None,
       map_data=None):

     mmo=map_manager_object

     self.unit_cell_grid= mmo.unit_cell_grid
     self.unit_cell_parameters= mmo.unit_cell_parameters
     self.space_group_number= mmo.space_group_number
     self.origin_shift_grid_units= mmo.origin_shift_grid_units
     self.working_map_n_xyz= mmo.working_map_n_xyz
     self._map_data=mmo._map_data
     self.data=None

     if map_data:
       self.data=map_data
       self.convert_to_double()
       if self._map_data.all() != self.working_map_n_xyz:
         raise Sorry("Map gridding (%s) does not match gridding specified by input map_manager (%s)" %(str(self._map_data.all()),str(self.working_map_n_xyz)))

  def set_origin_shift_grid_units(self,origin_shift_grid_units,
      allow_non_zero_previous_value=False):
    if not allow_non_zero_previous_value:
      assert self.original_origin_grid_units in [None,(0,0,0)] # make sure we
              #   didn't already shift it
    self.origin_shift_grid_units=origin_shift_grid_units

  def shift_origin(self,map_data=None,update_shift=None,log=sys.stdout):
    # Shift the origin of the map to (0,0,0) and
    #  optionally update origin_shift_grid_units
    #  if map_datais supplied, update_shift default=False
    #  if map_data is not supplied, update_shift default=True

    if update_shift is None:
      if map_data:
        update_shift=False
      else:
        update_shift=True

    if self._map_data and not map_data: # usual
      self._map_data,origin_shift_grid_units=shift_origin_of_map(
         self._map_data,self.origin_shift_grid_units)
      if update_shift:
         self.origin_shift_grid_units=origin_shift_grid_units

    elif map_data: # supplied map_data. Shift, optionally save updated shift,
         # and return shifted map

      shifted_map_data,origin_shift_grid_units=shift_origin_of_map(
        map_data,self.origin_shift_grid_units)
      if update_shift:
         self.origin_shift_grid_units=origin_shift_grid_units
      return shifted_map_data

    else:
      print ("No map to shift",file=log)


  def shift_ncs_to_match_working_map(self,ncs_object=None,reverse=False):
    # Shift an ncs object to match the working map (based
    #    on self.origin_shift_grid_units)
    shift_manager=self.get_shift_manager(reverse=reverse,
      ncs_object=ncs_object)
    return shift_manager.ncs_object # shifted ncs object

  def shift_ncs_to_match_original_map(self,ncs_object=None):
    return shift_ncs_to_match_working_map(ncs_object=ncs_object,
      reverse=True)

  def shift_model_to_match_working_map(self,model=None,reverse=False,
     log=sys.stdout):
    # Shift a model object to match the working map (based
    #    on self.origin_shift_grid_units)
    shift_manager=self.get_shift_manager(model=model,reverse=reverse)
    if(shift_manager.shift_cart is not None):
      print ("\nMap origin is not at (0,0,0): shifting the map and model.",
       file=log)
      model.set_shift_manager(shift_manager = shift_manager)
      model._process_input_model()
    return model

  def get_shift_manager(self,reverse=False,ncs_object=None,model=None):
    if not reverse: # usual
      origin_shift=self.origin_shift_grid_units
    else:  # go backwards
      a=self.origin_shift_grid_units
      origin_shift=[-a[0],-a[1],-a[2]]

    import mmtbx.utils
    shift_manager = mmtbx.utils.shift_origin(
      xray_structure=model.get_xray_structure(),
      n_xyz=self.working_map_n_xyz,
      ncs_object=ncs_object,
      origin_grid_units=origin_shift,
      crystal_symmetry = self.crystal_symmetry())
    return shift_manager


  def shift_model_to_match_original_map(self,model=None):
    # Shift a model object to match the original map (based
    #    on -self.origin_shift_grid_units)
    return shift_model_to_match_working_map(model=model,reverse=True)

  def map_shift_tracker(self):
    # Produces map_manager with no data so it is small but still knows how
    #   to write out maps
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
    '''

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
      new_map_manager.write_map(file_name=file_name)


def shift_origin_of_map(map_data,previous_origin_shift=None):

    assert map_data
    if not previous_origin_shift:
      previous_origin_shift=(0,0,0)

    new_origin_shift=map_data.origin()
    if new_origin_shift != (0,0,0):
      map_data=map_data.shift_origin() # have to use return value
      full_shift=[]
      for a,b in zip(previous_origin_shift,new_origin_shift):
        full_shift.append(a+b)
    else:
      full_shift=previous_origin_shift
    return map_data,full_shift
