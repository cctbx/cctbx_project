"""
Low-level tools to read and write maps (supersedes ccp4_map routines).
Normally use map_manager instead.
"""

from __future__ import absolute_import, division, print_function
import cctbx.array_family.flex as flex# import dependency
import os,time,sys
from libtbx.utils import Sorry
from iotbx.file_reader import splitext

import mrcfile
import warnings
from scitbx.array_family.flex import grid
from cctbx import maptbx
import numpy as np
from six.moves import range
from six.moves import zip
import math

# ====== HARD-WIRED DEFAULTS SEE DOC STRINGS BELOW FOR INFO =====

INTERNAL_STANDARD_ORDER=[3,2,1]  # THIS CANNOT BE MODIFIED

STANDARD_LIMITATIONS_DICT={  # THIS CAN BE MODIFIED AND ADDED TO
    "no_wrapping_outside_cell":
     "Values outside boundaries have no meaning",
    "wrapping_outside_cell":
     "Values outside boundaries are wrapped inside",
    "extract_unique":
     "This map is masked around unique region and not suitable for auto-sharpening.",
    "map_is_sharpened":
     "This map is auto-sharpened and not suitable for further auto-sharpening.",
    "map_is_density_modified": "This map has been density modified.",
     }
# ======= END OF HARD-WIRED DEFAULTS SEE DOC STRINGS BELOW FOR INFO ====




class map_reader:

  '''
   NOTE: See doc for iotbx.map_manager for further information on maps
   NORMALLY:   Use the map_manager class to read, write and manipulate maps
              Do not use these routines directly unless necessary

   Read an mrc/ccp4 map file
     NOTE: the values of nxstart,nystart,nzstart refer to columns, rows,
        sections, not to X,Y,Z.  The must be mapped using mapc,mapr,maps
        to get the value of "origin" that phenix uses to indicate the
        lower left corner of the map
      self.origin is self.nxstart_nystart_nzstart mapped to represent the
       origin in XYZ directions

     NOTE phenix calls "origin" the position of the lower left corner
       of the map.
     mrcfile calls "origin" the value of the field "origin" which is 3
       real numbers (A units, not grid units) indicating the placement
       of the grid point (0,0,0) relative to an external reference frame
       (typically that of a big map that this map is cut out from).

     Here we will use "external_origin" to refer to the mrc origin and
     "origin" for nxstart_nystart_nzstart.


  '''


  def __init__(self, file_name=None,
     internal_standard_order=INTERNAL_STANDARD_ORDER,
     header_only=False,
     ignore_missing_machine_stamp=True,
     print_warning_messages=False,
     ignore_all_errors=False,
     verbose=None,
     out=sys.stdout):

     self.read_map_file(file_name=file_name,
       internal_standard_order=internal_standard_order,
       header_only=header_only,
       ignore_missing_machine_stamp=ignore_missing_machine_stamp,
       print_warning_messages=print_warning_messages,
       ignore_all_errors=ignore_all_errors,
       verbose=verbose,
       out=out)

  def read_map_file(self, file_name=None,
     internal_standard_order=INTERNAL_STANDARD_ORDER,
     header_only=False,
     ignore_missing_machine_stamp=True,
     print_warning_messages=False,
     ignore_all_errors=False,
     verbose=None,
     out=sys.stdout):

    '''Read a map file with the MRC file reader'''
    # Check for file

    if not file_name:
      raise Sorry("Missing file name for map reader")
    if not os.path.isfile(file_name):
      raise Sorry("Missing file %s for map reader" %(file_name))
    base_name, file_ext, compress_ext = splitext(file_name)

    # Read the data


    with warnings.catch_warnings(record=True) as w:
      if compress_ext is None:
        mrc=mrcfile.mmap(file_name, mode='r', permissive=True)
        # Read memory-mapped for speed.

      else:
        mrc = mrcfile.open(file_name, mode='r', permissive=True)

      # Permissive to allow reading files with no machine stamp
      # Here we can deal with them
      for war in w:
        text = "\n  NOTE: WARNING message for the file '%s':\n  '%s'\n " % (
          file_name, war.message)
        if print_warning_messages:
          print(text, file=out)

        if ignore_all_errors:
          pass
        elif str(war.message).find("Unrecognised machine stamp") > -1 and \
                ignore_missing_machine_stamp:
          pass
        else:
          raise Sorry(text)


    # Note: the numpy function tolist() returns the python objects we need

    self.header_min=mrc.header.dmin.tolist() # returns float
    self.header_max=mrc.header.dmax.tolist()
    self.header_mean=mrc.header.dmean.tolist()
    self.header_rms=mrc.header.rms.tolist()

    self.unit_cell_grid=tuple((mrc.header.mx.tolist(),
                               mrc.header.my.tolist(),
                               mrc.header.mz.tolist(),))
    # NOTE: the values of nxstart,nystart,nzstart refer to columns, rows,
    #    sections, not to X,Y,Z.  The must be mapped using mapc,mapr,maps
    #    to get the value of "origin" that phenix uses to indicate the
    #    lower left corner of the map
    #  self.origin is self.nxstart_nystart_nzstart mapped to represent the
    #   origin in XYZ directions

    self.nxstart_nystart_nzstart=tuple(( mrc.header.nxstart.tolist(),
                                         mrc.header.nystart.tolist(),
                                         mrc.header.nzstart.tolist(),))
    self.origin=origin_as_xyz(
       nxstart_nystart_nzstart=self.nxstart_nystart_nzstart,
       mapc=mrc.header.mapc,mapr=mrc.header.mapr,maps=mrc.header.maps)

    # Labels
    self.labels=[]
    for i in range(mrc.header.nlabl):
      text=mrc.header.label[i]
      if not text:
        continue
      text = text.decode("utf-8")
      text = text.strip()
      self.labels.append(text)

    # NOTE phenix calls "origin" the position of the lower left corner
    #   of the map.
    # mrcfile calls "origin" the value of the field "origin" which is 3
    #   real numbers indicating the placement of the grid point (0,0,0) relative
    #   to an external reference frame (typically that of a model)

    # Here we will use "external_origin" to refer to the mrc origin and
    # "origin" for nxstart_nystart_nzstart.

    self.external_origin=tuple(mrc.header.origin.tolist())

    # Unit cell parameters (size of full unit cell in A)
    unit_cell_parameters=tuple(
          mrc.header.cella.tolist()+
          mrc.header.cellb.tolist())
    # Be permissive on cell angles. If all 0, set them to 90
    if list(unit_cell_parameters)[3:6] == [0,0,0]:
      print("Cell angles in %s are all zero...setting to 90,90,90" %(
        str(file_name)), file = out)
      unit_cell_parameters=tuple( list(unit_cell_parameters)[:3] + [90,90,90])

    # Space group number (1 for cryo-EM data)
    space_group_number=mrc.header.ispg.tolist()
    if space_group_number <= 0:
      space_group_number=1

    from cctbx import crystal
    try:
      self._unit_cell_crystal_symmetry=crystal.symmetry(
        unit_cell_parameters,
        space_group_number)
    except Exception as e:
      print("\nUnable to process cell parameters in %s:\n %s\n" %(
       str(file_name),str(
        unit_cell_parameters)), file = out)
      raise Sorry("Unable to process cell parameters : %s" %(str(
        unit_cell_parameters)))

    self._crystal_symmetry=None # Set this below after reading data

    if verbose:
      mrc.print_header(print_file=out)

    if header_only:
      return

    # Get the data. Note that the map file may have axes in any order. The
    #  order is defined by mapc, mapr, maps (columns, rows, sections).
    # Convert everything to the order 3,2,1 (X-sections, Y-rows, Z-columns).
    #  self.data is going to be a  flex double array (It originally was
    # float, same as ccp4_map.reader.data but now it is converted to double)

    self.data=numpy_map_as_flex_standard_order(np_array=mrc.data,
       mapc=mrc.header.mapc,mapr=mrc.header.mapr,maps=mrc.header.maps,
       internal_standard_order=internal_standard_order)

    # Shift the origin of this map to nxstart,nystart,nzstart
    if self.nxstart_nystart_nzstart != (0,0,0):
      grid_start=self.origin
      grid_end=tuple(add_list(grid_start,self.data.all()))
      g=grid(grid_start,grid_end)
      self.data.reshape(g)

    # Convert to double from float (NEW 2020-05-22)
    self.data=self.data.as_double()

    # Ready with self.data as flex.double() array.
    # Get the data with map_data=self.data or map_data=self.map_data()

    map_all=self.data.all()

    # Check if self.data contains 'nan'
    if math.isnan(flex.sum(self.data)):
       raise Sorry("Input map file contains 'nan':  not a valid map file")

    # Set self._crystal_symmetry which is crystal symmetry of part of map
    #  that is present

    self.set_crystal_symmetry_of_partial_map()

    # Set self.unit_cell_parameters self.space_group_number in case anyone
    #  still uses them
    self.set_legacy_parameters()

  def set_legacy_parameters(self):

    '''Set legacy unit_cell_parameters. Do not use these now.'''
    if self.unit_cell_crystal_symmetry():
      self.unit_cell_parameters=self.unit_cell().parameters()
      self.space_group_number=self.unit_cell_crystal_symmetry(
         ).space_group_number()

  def set_crystal_symmetry_of_partial_map(self,
     space_group_number  = None):
    '''
      This sets the crystal_symmetry of a partial map based on the
      gridding of the part of the map that is present.
      If space_group_number is specified, use it. Otherwise,
        if exactly the entire map is present, use space group of
        current self._crystal_symmetry, otherwise use space group P1

      Note that self._crystal_symmetry space group might not match
      self._unit_cell_crystal_symmetry (space group of entire map) because
      map may already have been boxed.

    '''
    map_all=self.map_data().all()

    a,b,c, al,be,ga = self.unit_cell().parameters()
    a = a * map_all[0]/self.unit_cell_grid[0]
    b = b * map_all[1]/self.unit_cell_grid[1]
    c = c * map_all[2]/self.unit_cell_grid[2]

    if tuple(map_all) == tuple(self.unit_cell_grid[:3]):
      if space_group_number:
        space_group_number_use=space_group_number
      # Take space group we have...
      elif self.crystal_symmetry():  # use sg of existing box crystal symmetry
        space_group_number_use=self.crystal_symmetry(
          ).space_group_number()
      else:  # otherwise take the sg of crystal symmetry of full cell
        space_group_number_use=self.unit_cell_crystal_symmetry(
           ).space_group_number()
    else: # boxed, use P1 always
      space_group_number_use=1  # use Space group 1 (P 1) for any partial cell

    from cctbx import crystal
    self._crystal_symmetry=crystal.symmetry((a,b,c, al,be,ga),
         space_group_number_use)

  # Code to check for specific text in map labels limiting the use of the map

  def cannot_be_sharpened(self):
    '''Determine whether map has limitations that prevent sharpening'''
    if self.is_in_limitations("extract_unique"):
      return True
    if self.is_in_limitations("map_is_sharpened"):
      return True
    return False

  def wrapping_from_input_file(self):
    '''Guess wrapping of map from limitations in input file'''
    if self.is_in_limitations("no_wrapping_outside_cell"):
      return False # cannot wrap outside cell
    elif self.is_in_limitations("wrapping_outside_cell"):
      return True # can wrap outside cell
    else:
      return None # no information


  def remove_limitation(self,text):
    '''Remove specified limitations'''
    limitations=self.get_limitations()
    new_labels=[]
    if not self.labels:
      return new_labels
    for label in self.labels:
      text_matches_key=False
      for key,message in zip(limitations.limitations,
         limitations.messages):
        if label==message and key==text:
          text_matches_key=True
          break
      if not text_matches_key:
        new_labels.append(label)
    return new_labels

  def is_in_limitations(self,text):
    '''Check whether a limitation is in limitations of this map'''
    limitations=self.get_limitations()
    if not limitations:
      return False
    elif text in limitations.limitations:
      return True
    else:
      return False

  def get_labels(self):
    '''Get the labels for this map'''
    return self.labels

  def get_additional_labels(self):
    '''Get all labels that are not limitations'''
    limitations=self.get_limitations()
    if not limitations:
      return self.labels
    else:
      additional_labels=[]
      for label in self.labels:
        if not label.strip() in limitations.limitations:
          additional_labels.append(label.strip())
    return additional_labels

  def get_limitations(self):
    '''Get all limitations of this map'''
    limitations=[]
    limitation_messages=[]
    if self.labels:
      for label in self.labels:
        limitation_message=self.get_limitation(label.strip())
        if limitation_message:
          limitations.append(label.strip())
          limitation_messages.append(limitation_message.strip())
    from libtbx import group_args
    return group_args(
       limitations=limitations,
       messages=limitation_messages,
      )

  def get_limitation(self,label):
    '''Get limitation corresponding to this label'''
    return STANDARD_LIMITATIONS_DICT.get(label,None)

  def show_summary(self, out=sys.stdout, prefix=""):
    '''Summarize map information'''
    data=self.map_data()

    if hasattr(self,'header_min'):
      print(prefix + "header_min: ", self.header_min, file=out)
      print(prefix + "header_max: ", self.header_max, file=out)
      print(prefix + "header_mean:", self.header_mean, file=out)
      print(prefix + "header_rms: ", self.header_rms, file=out)
    print("\n"+prefix + "Information about FULL UNIT CELL:",file=out)
    print(prefix + "unit cell grid:", self.unit_cell_grid, file=out)
    print(prefix + "unit cell parameters:", self.unit_cell().parameters(), file=out)
    print(prefix + "space group number:  ", self.unit_cell_crystal_symmetry().space_group_number(), file=out)

    if not data:
      print("No map data available")
    else:
      print("\n"+prefix + "Information about the PART OF MAP THAT IS PRESENT:",
       file=out)
      print(prefix + "map cell grid:", data.all(), file=out)
      print(prefix + "map cell parameters:",
        self.crystal_symmetry().unit_cell().parameters(), file=out)
      print(prefix + "map origin:", data.origin(), file=out)
      print(prefix + "pixel size: (%.4f, %.4f, %.4f) " %(
        self.pixel_sizes()), file=out)
    if hasattr(self,'origin_shift_grid_units'):
      print(prefix + "Shift (grid units) to place origin at original position:",
          self.origin_shift_grid_units, file=out)

    if hasattr(self,'wrapping'):
      print(prefix +
       "Wrapping (using unit_cell_translations to get map values) allowed:",
          self.wrapping(), file=out) # don't try too hard


    if hasattr(self,'_ncs_object') and self._ncs_object:
      print (prefix + "Associated ncs_object with",
          self._ncs_object.max_operators(),"operators",
           file=out)

    if hasattr(self,'_model') and self._model:
      print (prefix + "Associated model with",
          self._model.get_hierarchy().overall_counts().n_residues,"residues",
           file=out)

    if self.high_resolution():
      print (prefix + "High-resolution limit of map: ",self.high_resolution(),
            file=out)

  def pixel_sizes(self):
    '''Return tuple with pixel size in each direction (normally all the same)'''
    data=self.map_data()
    if not data:
      return None
    cs=self.crystal_symmetry()
    cell_params=cs.unit_cell().parameters()[:3]
    map_all=data.all()
    pa=cell_params[0]/map_all[0]
    pb=cell_params[1]/map_all[1]
    pc=cell_params[2]/map_all[2]
    return (pa,pb,pc)

  def crystal_symmetry(self):
    '''
      This is "crystal_symmetry" of a box the size of the map that is present
    '''
    return self._crystal_symmetry


  def unit_cell(self):
    '''
     This is the cell dimensions and angles of the full unit_cell
    '''
    cs=self.unit_cell_crystal_symmetry()
    if cs:
      return cs.unit_cell()
    else:
      return None

  def unit_cell_crystal_symmetry(self):
    '''
     This is the cell dimensions and angles of the full unit_cell
    '''
    return self._unit_cell_crystal_symmetry


  def statistics(self):
    '''Calculate statistics for this map'''
    from cctbx import maptbx
    return maptbx.statistics(self.map_data())

  def get_origin(self):
    '''Return the origin of this map'''
    data=self.map_data()
    if data:
      return data.origin()
    else:
      return None

  def map_data(self):

    '''
       Input data is converted to double and stored in self.data
       self.map_data() always returns self.data
    '''

    return self.data

  def high_resolution(self):
    '''Return high resolution limit, if specified'''
    if hasattr(self,'_high_resolution'):
      return self._high_resolution
    else:
      return None

  def is_similar_map(self, other):
    '''Check whether this map is similar (crystal symmetry, size) to another'''
    f1 = self.crystal_symmetry().is_similar_symmetry(other.crystal_symmetry())
    s = self.map_data()
    o = other.map_data()
    if not s or not o:
      return None

    f2 = s.focus()  == o.focus()
    f3 = s.origin() == o.origin()
    f4 = s.all()    == o.all()
    if([f1,f2,f3,f4].count(False)>0): return False
    else: return True

  def grid_unit_cell(self):
    """
    This is a unit cell describing one pixel of the map.
    It is used in maptbx.non_crystallographic_eight_point_interpolation.
    This grid_unit_cell is the original unit cell divided by the original
    grid size.
    """
    from cctbx import uctbx
    unit_cell_parameters=self.unit_cell().parameters()
    a = unit_cell_parameters[0] / self.unit_cell_grid[0]
    b = unit_cell_parameters[1] / self.unit_cell_grid[1]
    c = unit_cell_parameters[2] / self.unit_cell_grid[2]
    alpha,beta,gamma = unit_cell_parameters[3:6]
    return uctbx.unit_cell((a,b,c,alpha,beta,gamma))

class write_ccp4_map:
  '''
     NORMALLY:   Use the map_manager class to read, write and manipulate maps
                Do not use these routines directly unless necessary

     Write an mrc/CCP4 map file

     Always writes with column,row,section (mapc,mapr,maps) of (3,2,1)
      Columns are Z, rows are Y, sections in X.  This matches the shape of
      flex arrays.  Could be changed without difficulty by transposing the
      numpy array as is done in the map_reader class.

     Class parallel to iotbx.ccp4_map.write_ccp4_map

       The parameter map_data should be a flex array (normally flex double)

       Options:

       Specify unit_cell_grid, or gridding_first and gridding_last,
          or take grid values from map_data.
       If gridding_first and last supplied, input map origin must be at (0,0,0)

       Specify crystal_symmetry or unit_cell and space_group

       Specify labels (80-character information lines) or not.  Normally
        these should be specified to retain information from previous maps.

       Notes on grid values:

       All input grid values, cell dimensions, etc are along X,Y,Z

       Definitions:

         unit_cell_grid (grid units along axes of full unit cell)

         nxyz_start (starting point of map to be written out, in grid units,
         nxyz_end (ending point of map to be written out, in grid units,
          optionally set with gridding_first and gridding_last if input map
          has origin at (0,0,0))

         mapc,mapr,maps (assignment of the axes to X, Y and Z in numpy array
           on input or output using mrcfile.  Columns are mapc (X=1, Y=2, Z=3),
           rows are mapr (X=1, Y=2, Z=3), sections are maps (X=1, Y=2, Z=3).
         Phenix writes maps with (mapc,mapr,maps)=(3,2,1).  Many other
           programs use (1,2,3). This routine can read any order.

        origin of flex array (grid point of lower left point in the map)
        all of flex array (number of grid points in each direction)

        Examples:
        For map_data with origin=(0,0,0) and all=(3,3,3) the map runs from
          0 to 2 in each direction, total elements = 27.  nxyz_start=(0,0,0)
          nxyz_end=(2,2,2).

        For map_data with origin=(5,5,5) and all=(3,3,3) the map runs from
          5 to 7 in each direction, total elements = 27.  nxyz_start=(5,5,5)
          nxyz_end=(7,7,7).  map_data.all()=(3,3,3). map_data.focus()=(8,8,8)

      Make sure map is unpadded

  '''

  def __init__(self,
      file_name=None,
      crystal_symmetry=None,
      unit_cell=None,
      space_group=None,
      map_data=None,
      labels=None,
      gridding_first=None,
      gridding_last=None,
      unit_cell_grid=None,
      origin_shift_grid_units=None,  # origin shift (grid units) to be applied
      external_origin=None, # Do not use this unless required
      output_axis_order=INTERNAL_STANDARD_ORDER,
      internal_standard_order=INTERNAL_STANDARD_ORDER,
      replace_backslash = True,
      verbose=None,
      out=sys.stdout,
      ):

    '''
     Write mrc/ccp4 format file
     if replace_backslash then replace backslashes in labels with forward slash
    '''


    assert map_data  # should never be called without map_data

    if map_data.is_padded(): # copied from cctbx/miller/__init__.py
      map_data=maptbx.copy(map_data, flex.grid(map_data.focus()))

    # Get unit_cell and space_group if crystal_symmetry is supplied:
    if unit_cell is None:
      assert crystal_symmetry is not None
      unit_cell=crystal_symmetry.unit_cell()

    if space_group is None:  # use P1 if not supplied
      if crystal_symmetry is not None:
        space_group=crystal_symmetry.space_group()
      else:
        from cctbx import sgtbx
        space_group=sgtbx.space_group_info("P1").group()

    # Get empty labels if not supplied
    if not labels:
      labels=create_output_labels()


    if gridding_first is not None and gridding_last is not None:

      # Writing box from gridding_first to gridding_last (inclusive)
      # NOTE: This is not a common use of this routine

      assert len(gridding_first)==3 and len(gridding_last)==3
      assert unit_cell_grid is None
      assert map_data.origin()==(0,0,0)  # unshifted maps only

      nxyz_start=gridding_first
      nxyz_end=gridding_last
      unit_cell_grid=map_data.all()

      new_map_data = maptbx.copy(map_data, tuple(nxyz_start), tuple(nxyz_end))
      #   NOTE: end point of map is nxyz_end, so size of map (new all()) is
      #   (nxyz_end-nxyz_start+ (1,1,1))
      assert tuple(new_map_data.origin()) == tuple(nxyz_start)
      assert new_map_data.all()==tuple(
         add_list((1,1,1),subtract_list(gridding_last,gridding_first)))

      new_map_data=new_map_data.shift_origin() # this is the map we will pass in

    else:

      # This is the recommended way to use this writer
      # Optionally uses supplied unit_cell_grid.
      # Allows writing just a part of a map to a file, but retaining
      #   information on size of whole map if unit_cell_grid is supplied.
      # Takes origin and size of map from the map_data array
      # Optionally shifts origin based on input "origin_shift_grid_units"

      assert gridding_first is None and gridding_last is None

      if unit_cell_grid is None:
        # Assumes map_data.all() is the entire unit cell
        unit_cell_grid=map_data.all()
        # Note: if map_data.origin()==(0,0,0) this grid corresponds to the
        #  box of density that is present.
      else:
        assert len(unit_cell_grid)==3
        # Assumes unit_cell_grid is the entire unit cell

      if origin_shift_grid_units:  # Supplied non-zero origin. Input map must be
        #   at (0,0,0).  Use supplied origin_shift_grid_units as origin
        assert map_data.origin()==(0,0,0)
        nxyz_start=origin_shift_grid_units
        new_map_data=map_data
      else:
        nxyz_start=map_data.origin()
        new_map_data=map_data.deep_copy().shift_origin()

    # Ready to write the map

    assert new_map_data.origin()==(0,0,0) # must not be shifted at this point
    assert space_group is not None
    assert unit_cell is not None

    # Open file for writing
    mrc=mrcfile.new(file_name,overwrite=True)

    # Convert flex array to the numpy array required for mrcfile
    if hasattr(new_map_data,'as_float'):
      numpy_data=new_map_data.as_float().as_numpy_array()
    else: # was float
      numpy_data=new_map_data.as_numpy_array()

    # This numpy_data array is always in the order (3,2,1): columns are Z,
    #  rows are Y, sections in X.  This comes from the shape of flex arrays.

    # To write with another order, call with values for output_axis_order

    if output_axis_order!=internal_standard_order:
      i_order=get_standard_order(output_axis_order[0],
        output_axis_order[1],output_axis_order[2],
        internal_standard_order=internal_standard_order,reverse=True)
      numpy_data_output_axis_order=np.transpose(numpy_data,i_order)
    else:
      numpy_data_output_axis_order=numpy_data
    mrc.header.mapc=output_axis_order[0]
    mrc.header.mapr=output_axis_order[1]
    mrc.header.maps=output_axis_order[2]

    mrc.set_data(numpy_data_output_axis_order) # numpy array

    # Labels
    # Keep all limitations labels and other labels up to total of 10 or fewer
    output_labels=select_output_labels(labels,
       replace_backslash = replace_backslash)

    mrc.header.nlabl=len(output_labels)
    for i in range(min(10,len(output_labels))):
      mrc.header.label[i]=output_labels[i]
    mrc.update_header_from_data() # don't move this later as we overwrite values

    # Unit cell parameters and space group
    abc=unit_cell.parameters()[:3]
    angles=unit_cell.parameters()[3:]
    mrc.header.cella=abc
    mrc.header.cellb=angles
    space_group_number=space_group.info().type().number()
    mrc.header.ispg=space_group_number

    # Start point of the supplied map in grid units

    # nxyz_start is the origin (grid units) in XYZ coordinate system.
    # The mrc header needs the origin along columns, rows, sections which
    #  is represented here as nxstart_nystart_nzstart

    nxstart_nystart_nzstart=origin_as_crs(origin=nxyz_start,
       mapc=mrc.header.mapc,mapr=mrc.header.mapr,maps=mrc.header.maps)

    mrc.header.nxstart=nxstart_nystart_nzstart[0]
    mrc.header.nystart=nxstart_nystart_nzstart[1]
    mrc.header.nzstart=nxstart_nystart_nzstart[2]

    # Size of entire unit cell in grid units
    # This is ALWAYS along X,Y,Z regardless of the sectioning of the map
    # Note that this can cause confusion as mrc.header.nx may be a different
    #   axis than mrc.header.mx (mx is always X, nx is whatever axis is
    #   specified by mapc)
    mrc.header.mx=unit_cell_grid[0]
    mrc.header.my=unit_cell_grid[1]
    mrc.header.mz=unit_cell_grid[2]

    # External origin
    # NOTE: External origin should rarely be used.  It is a poorly-defined
    #   element that refers to the position of a non-defined external model
    #   (PDB file).   For origin purposes use instead "origin".

    if external_origin is None:
      external_origin=(0.,0.,0.,)
    # This also is ALWAYS along X,Y,Z regardless of the sectioning of the map
    mrc.header.origin.x=external_origin[0]
    mrc.header.origin.y=external_origin[1]
    mrc.header.origin.z=external_origin[2]

    # Update header
    mrc.update_header_stats()

    if verbose:
      mrc.print_header()

    # Write the file
    mrc.close()

def apply_origin_shift(map_data=None,origin_shift=None):
  from scitbx.matrix import col
  new_last=list(col(map_data.all())+col(origin_shift))
  new_map_data=map_data.deep_copy()
  new_map_data.resize(flex.grid(origin_shift,new_last))
  return new_map_data

def create_output_labels(
    program_name=None,  # e.g., auto_sharpen
    input_file_name=None, # mymrc_file.mrc
    input_labels=None,    # input labels from mymrc_file.mrc
    output_labels=None,  # any specific output labels
    limitations=None,    # any standard limitations
     ):

  output_map_labels=[]

  # get standard label from program_name, input_file_name and date
  text=""
  if program_name:
    text+="%s" %(program_name)
  if input_file_name:
    text+=' on %s' %(input_file_name)
  text+=' %s' %(time.asctime())
  text=text[:80]
  output_map_labels.append(text)

  # any specific limitations
  if limitations:
    for limitation in limitations:
      if not limitation in STANDARD_LIMITATIONS_DICT:
        print("The limitation '%s' is not in STANDARD_LIMITATIONS_DICT: '%s'" %(
       limitation,str(list(STANDARD_LIMITATIONS_DICT))))
      assert limitation in STANDARD_LIMITATIONS_DICT
      output_map_labels.append(limitation)

  # any specific labels given
  if output_labels:
    output_map_labels+=output_labels

  # any input labels to pass on
  if input_labels:
    output_map_labels+=input_labels

  # Now write out up to 10 unique labels
  final_labels=[]
  for label in output_map_labels:
    if not label in final_labels:
      final_labels.append(label)
  final_labels=select_output_labels(final_labels)
  return final_labels

def select_output_labels(labels,max_labels=10, replace_backslash = None):
  n_limitations=0
  used_labels=[]
  for label in labels:
    if label in STANDARD_LIMITATIONS_DICT and not label in used_labels:
      n_limitations+=1
      used_labels.append(label)
  output_labels=[]
  n_available=max_labels-n_limitations
  n_general=0
  for label in labels:
    if label in output_labels: continue
    if len(output_labels)>=max_labels:
      continue
    if label in STANDARD_LIMITATIONS_DICT:
      output_labels.append(label)
    elif n_general < n_available:
      n_general+=1
      output_labels.append(label)

  if replace_backslash:
    new_labels = []
    for label in output_labels:
      new_label = label
      if isinstance(label, str):
        new_label = label.replace("\\","/")
      new_labels.append(new_label)
    output_labels = new_labels
  return output_labels

def get_standard_order(mapc,mapr,maps,internal_standard_order=None,
    reverse=False):
  '''

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

  assert internal_standard_order==INTERNAL_STANDARD_ORDER  # This is hard-wired for flex array

  standard_order_np=[
    internal_standard_order[0]-1,
    internal_standard_order[1]-1,
    internal_standard_order[2]-1]

  mapc_np=mapc-1
  mapr_np=mapr-1
  maps_np=maps-1

  # Set up ordering for transposition
  # if mapc,mapr,maps=(3,2,1) then with internal_standard_order=(3,2,1),
  #   i_order=(1,2,3)
  #  (0,1,2) and (2,1,0) for axis number 0-2)

  #  from above:
  #   i_order[0]=2 means input axis 0 becomes output axis 2

  i_order=[None,None,None]
  i_order[mapc_np]=standard_order_np[0]
  i_order[mapr_np]=standard_order_np[1]
  i_order[maps_np]=standard_order_np[2]

  i_order=tuple(i_order)

  if not reverse:
    assert i_order.count(None)==0
    for i in range(3):
      assert i_order.count(i)==1
    return i_order
  else:
    # To reverse it:
    #   if standard is: input axis 0 -> output 2 :i0=2
    #   then reversed:  input 2 -> output 0      :i2=0
    i_order_reverse=[None,None,None]
    for i in range(3):
      i_order_reverse[i_order[i]]=i
    assert i_order_reverse.count(None)==0
    for i in range(3):
      assert i_order_reverse.count(i)==1
    return i_order_reverse

def origin_as_crs(origin=None,mapc=None,mapr=None,maps=None):

  # convert  origin (origin along x,y,z) to nxstart,nystart,nzstart (origin
  #   along columns,rows,sections)

  # mapc is (1,2,or 3) indicating that columns are the (X,Y,or Z) axis
  #  same for mapr (rows) and maps (sections).

  # So if mapc=1, origin along columns is origin[0]
  #    if mapc=3, origin along columns is origin[3]

  order=(mapc-1,mapr-1,maps-1)

  assert origin is not None and mapc is not None and \
    mapr is not None and maps is not None

  # order[0]=2 means mapc-1=2, or columns are Z, or that nxstart (column start)
  #   is the origin along Z

  nxstart_nystart_nzstart=[None,None,None]
  for i in range(3):
    nxstart_nystart_nzstart[order[i]]=origin[i]
  return nxstart_nystart_nzstart

def origin_as_xyz(nxstart_nystart_nzstart=None,mapc=None,mapr=None,maps=None):

  # convert  nxstart,nystart,nzstart (origin along columns,rows,sections) to
  #   origin (origin along x,y,z)

  # mapc is (1,2,or 3) indicating that columns are the (X,Y,or Z) axis
  #  same for mapr (rows) and maps (sections).

  # So if mapc=1, origin along x is nxstart.
  #    if mapc=3, origin along x is nzstart

  order=(mapc-1,mapr-1,maps-1)

  assert nxstart_nystart_nzstart is not None and mapc is not None and \
    mapr is not None and maps is not None

  # order[0]=2 means mapc-1=2, or columns are Z, or that nxstart (column start)
  #   is the origin along Z

  origin=[None,None,None]
  for i in range(3):
    origin[order[i]]=nxstart_nystart_nzstart[i]
  return origin

def numpy_map_as_flex_standard_order(np_array=None,
  mapc=None,mapr=None,maps=None,internal_standard_order=None):


  # Convert numpy version of map (from CCP4-EM mrcfile) to flex.float array

  assert np_array is not None and mapc is not None and \
    mapr is not None and maps is not None

  # Verify that the numpy array is float.  If it isn't then this conversion
  #   will not work below when we convert to a flex array (but it will not
  #   give an error message).

  assert type(np_array[0,0,0].tolist())==type(float(1))

  # Input map can have columns, rows, sections corresponding to any of (X,Y,Z)

  # Get the order for transposing input map (columns are mapc (X=1, Y=2, Z=3),
  #   rows are mapr (X=1, Y=2, Z=3), sections are maps (X=1, Y=2, Z=3).
  i_order=get_standard_order(mapc,mapr,maps,
     internal_standard_order=internal_standard_order)

  # Transpose the input numpy array to match Phenix expected order of (3,2,1)
  np_array_standard_order=np.transpose(np_array,i_order)

  # this np array may have any actual memory layout. We have to
  #  flatten it out, read into flex array and reshape

  # Save the shape
  shape=tuple(np_array_standard_order.shape)

  # Flatten it out
  np_array_standard_order_1d=np_array_standard_order.flatten()

  # Read in to flex array
  flex_array=flex.float(np_array_standard_order_1d)

  # set up new shape (same as was in the numpy array after transposing it)
  flex_grid=grid(shape)

  # Reshape the flex array
  flex_array.reshape(flex_grid)

  # All done. Returning float array
  return flex_array


def add_list(list_a,list_b):
  assert len(list_a)==len(list_b)
  new_list=[]
  for a,b in zip(list_a,list_b):
    new_list.append(a+b)
  return new_list

def subtract_list(list_a,list_b):
  assert len(list_a)==len(list_b)
  new_list=[]
  for a,b in zip(list_a,list_b):
    new_list.append(a-b)
  return new_list
