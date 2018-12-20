from __future__ import division
import cctbx.array_family.flex as flex# import dependency
import os
from libtbx.utils import Sorry
from iotbx.ccp4_map import utils  # utilities in common with ccp4_map

#  mrcfile

#  Access to ccp4/mrc maps in a fashion parallel to iotbx.ccp4_map except
#    using mrcfile libraries instead of cmap libraries from ccp4-em

#  See http://www.ccpem.ac.uk/mrc_format/mrc2014.php for MRC format
#  See https://pypi.org/project/mrcfile/ for mrcfile library documentation


class map_reader(utils):

  # Read an mrc/ccp4 map file

  def __init__(self, file_name=None, header_only=False, verbose=False):

    # Check for file

    if not file_name:
      raise Sorry("Missing file name for map reader")
    if not os.path.isfile(file_name):
      raise Sorry("Missing file %s for map reader" %(file_name))

    # Read the data

    import mrcfile

    mrc=mrcfile.mmap(file_name, mode='r')  # Read memory-mapped for speed

    # Set the same attributes that are set in ccp4_map.reader()

    # Note: the numpy function tolist() returns the python objects we need

    self.header_min=mrc.header.dmin.tolist() # returns float
    self.header_max=mrc.header.dmax.tolist()
    self.header_mean=mrc.header.dmean.tolist()
    self.header_rms=mrc.header.rms.tolist()

    self.unit_cell_grid=tuple((mrc.header.mx.tolist(),
                               mrc.header.my.tolist(),
                               mrc.header.mz.tolist(),))

    self.nxstart_nystart_nzstart=tuple(( mrc.header.nxstart.tolist(),
                                         mrc.header.nystart.tolist(),
                                         mrc.header.nzstart.tolist(),))

    # NOTE phenix calls "origin" the value of nxstart_nystart_nzstart
    # mrcfile calls "origin" the value of the field "origin" which is 3
    #   real numbers indicating the placement of the grid point (0,0,0) relative
    #   to an external reference frame (typically that of a model)

    # Here we will use "external_origin" to refer to the mrc origin and
    # "origin" for nxstart_nystart_nzstart.

    self.external_origin=tuple(mrc.header.origin.tolist())

    # Unit cell parameters (size of full unit cell in A)
    self.unit_cell_parameters=tuple(
          mrc.header.cella.tolist()+
          mrc.header.cellb.tolist())

    # Space group number (1 for cryo-EM data)
    self.space_group_number=mrc.header.ispg.tolist()

    if verbose:
      mrc.print_header()

    if header_only:
      return

    # Get the data. Note that the map file may have axes in any order. The
    #  order is defined by mapc, mapr, maps (columns, rows, sections).
    # Convert everything to the order 3,2,1 (Z-sections, Y-rows, X-columns).
    #  self.data is a flex float array (same as ccp4_map.reader.data)

    self.data=numpy_map_as_flex_standard_order(np_array=mrc.data,
     mapc=mrc.header.mapc,mapr=mrc.header.mapr,maps=mrc.header.maps)

    # Shift the origin of this map to nxstart,nystart,nzstart
    if self.nxstart_nystart_nzstart != (0,0,0):
      from scitbx.array_family.flex import grid
      grid_start=self.nxstart_nystart_nzstart
      grid_end=tuple(add_list(grid_start,self.data.all()))
      g=grid(grid_start,grid_end)
      self.data.reshape(g)

    # Ready with self.data as float flex array. For normal use convert to
    # double with map_data=self.data.as_double()

class write_ccp4_map:

    # Write an mrc/CCP4 map file

    # Class parallel to iotbx.ccp4_map.write_ccp4_map

  def __init__(self,
      file_name=None,
      unit_cell=None,
      space_group=None,
      map_data=None,
      labels=None,
      gridding_first=None,
      gridding_last=None,
      unit_cell_grid=None,
      external_origin=None,
      standard_order=[3,2,1],
      verbose=False,
      ):

    #  The parameter map_data should be a flex array (normally flex.vec3_double)

    #  Options:  specify unit_cell_grid, or gridding_first and gridding_last,
    #     or take grid values from map_data.

    #  If gridding_first and last supplied, input map origin must be at (0,0,0)

    #  Notes on grid values:
    #    unit_cell_grid (grid units along axes of full unit cell)
    #    nxyz_start (starting point of map that is present, in grid units)
    #    nxyz_end   (ending point of map that is present, in grid units)

    #   For map_data with origin=(0,0,0) and all=(3,3,3) the map runs from
    #     0 to 2 in each direction, total elements = 27.  nxyz_start=(0,0,0)
    #     nxyz_end=(2,2,2).

    #   For map_data with origin=(5,5,5) and all=(3,3,3) the map runs from
    #     5 to 7 in each direction, total elements = 27.  nxyz_start=(5,5,5)
    #     nxyz_end=(7,7,7).  map_data.all()=(3,3,3). map_data.focus()=(8,8,8)

    if gridding_first is not None and gridding_last is not None:

      # Writing box from gridding_first to gridding_last (inclusive)
      # NOTE: This is not a common use of this routine

      assert len(gridding_first)==3 and len(gridding_last)==3
      assert unit_cell_grid is None
      assert map_data.origin()==(0,0,0)  # unshifted maps only

      nxyz_start=gridding_first
      nxyz_end=gridding_last
      unit_cell_grid=map_data.all()

      from cctbx import maptbx
      new_map_data = maptbx.copy(map_data, tuple(nxyz_start), tuple(nxyz_end))
      #   NOTE: end point of map is nxyz_end, so size of map (new all()) is
      #   (nxyz_end-nxyz_start+ (1,1,1))
      assert new_map_data.origin() == nxyz_start
      assert new_map_data.all()==tuple(
         add_list((1,1,1),subtract_list(gridding_last,gridding_first)))

      new_map_data=new_map_data.shift_origin() # this is the map we will pass in

    elif unit_cell_grid is not None:
      # This is the recommended way to use this writer
      # Uses supplied unit_cell_grid, allowing writing just a part of
      #   a map to a file, but retaining information on size of whole map
      # Takes origin and size of map to write from the map_data array

      assert len(unit_cell_grid)==3
      assert gridding_first is None and gridding_last is None

      nxyz_start=map_data.origin()
      nxyz_end=map_data.last(False)  # last points where data is present
      new_map_data=map_data.deep_copy().shift_origin()

    else: # common case, but not recommended as unit_cell_grid not supplied
      # Use information from map_data
      # Assumes entire map is the entire unit cell but does not shift it
      unit_cell_grid=map_data.all()
      nxyz_start=map_data.origin()
      nxyz_end=map_data.last(False) # last points where data is present
      new_map_data=map_data.deep_copy().shift_origin()

    # Ready to write the map

    assert new_map_data.origin()==(0,0,0) # must not be shifted at this point

    # Open file for writing
    import mrcfile
    mrc=mrcfile.new(file_name,overwrite=True)

    # Data.  Convert to numpy array required for mrcfile
    numpy_data=new_map_data.as_float().as_numpy_array()
    mrc.set_data(numpy_data) # numpy array

    # Labels
    mrc.header.nlabl=labels.size()
    for i in xrange(min(10,labels.size())):
      mrc.header.label[i]=labels[i]
    mrc.update_header_from_data() # don't move later as we overwrite values

    # Unit cell parameters
    abc=unit_cell.parameters()[:3]
    angles=unit_cell.parameters()[3:]
    mrc.header.cella=abc
    mrc.header.cellb=angles

    # Axis order. Always is (3,2,1)  (see standard_order)
    mrc.header.mapc=standard_order[0]
    mrc.header.mapr=standard_order[1]
    mrc.header.maps=standard_order[2]

    # Start point of the supplied map in grid units
    mrc.header.nxstart=nxyz_start[0]
    mrc.header.nystart=nxyz_start[1]
    mrc.header.nzstart=nxyz_start[2]

    # Size of entire unit cell in grid units
    mrc.header.mx=unit_cell_grid[0]
    mrc.header.my=unit_cell_grid[1]
    mrc.header.mz=unit_cell_grid[2]

    # External origin
    if external_origin is None:
      external_origin=(0.,0.,0.,)
    mrc.header.origin=external_origin

    # Update header
    mrc.update_header_stats()

    if verbose:
      mrc.print_header()

    # Write the file
    mrc.close()

def get_standard_order(mapc,mapr,maps,
    standard_order=[3,2,1]):

  # MRC standard order is 3,2,1.  Convert everything to this order.

  #  NOTE: numpy arrays indexed from 0 so this is equivalent to
  #   order of 2,1,0 in the numpy array

  # NOTE:  MRC format allows data axes to be swapped using the header
  #   mapc mapr and maps fields. However the mrcfile library does not
  #   attempt to swap the axes and assigns the columns to X, rows to Y and
  #   sections to Z. The data array is indexed C-style, so data values can
  #   be accessed using mrc.data[z][y][x].

  # NOTE: normal expectation is that phenix will read/write with the
  #   order 3,2,1. This means Z-sections (index=3), Y rows (index=2),
  #   X columns (index=1). This correxponds to
  #    mapc (columns) =   1 or X
  #    mapr (rows)    =   2 or Y
  #    maps (sections) =  3 or Z

  # In the numpy array (2,1,0 instead of 3,2,1):

  # To transpose, specify i0,i1,i2 where:
  #     i0=2 means input axis 0 becomes output axis 2
  #     NOTE:  axes are 0,1,2 etc, not 1,2,3
  #   Examples:
  #     np.transpose(a,(0,1,2))  does nothing
  #     np.transpose(a,(1,2,0)) output axis 0 is input axis 1
  #
  # We want output axes to always be 2,1,0 and input axes for numpy array are
  #   (mapc-1,mapr-1,maps-1):

  # For hard-wired mrcfile library, The transposition is always:
  #   i_mapc=3    i_mapc_np=2
  #   i_mapr=2    i_mapr_np=1
  #   i_maps=1    i_maps_np=0

  assert standard_order==[3,2,1]  # Take this out if you want to change it

  standard_order_np=[
    standard_order[0]-1,
    standard_order[1]-1,
    standard_order[2]-1]

  mapc_np=mapc-1
  mapr_np=mapr-1
  maps_np=maps-1

  # Set up ordering for transposition

  i_order=[None,None,None]
  i_order[mapc_np]=standard_order_np[0]
  i_order[mapr_np]=standard_order_np[1]
  i_order[maps_np]=standard_order_np[2]

  i_order=tuple(i_order)

  assert i_order.count(None)==0
  for i in xrange(3):
    assert i_order.count(i)==1

  return i_order

def numpy_map_as_flex_standard_order(np_array=None,
  mapc=None,mapr=None,maps=None):


  # Convert numpy version of map (from CCP4-EM mrcfile) to flex.float array

  assert np_array is not None and mapc is not None and \
    mapr is not None and maps is not None

  # Verify that the numpy array is float.  If it isn't then this conversion
  #   will not work below when we convert to a flex array (but it will not
  #   give an error message).

  assert type(np_array[0,0,0].tolist())==type(float(1))

  i_order=get_standard_order(mapc,mapr,maps)

  import numpy as np
  np_array_standard_order=np.transpose(np_array,i_order)

  # this np array may have any actual memory layout. We have to
  #  flatten it out, read into flex array and reshape

  # Save the shape
  shape=tuple(np_array_standard_order.shape)

  # Flatten it out
  np_array_standard_order_1d=np_array_standard_order.flatten().tolist()

  # Read in to flex array
  flex_array=flex.float(np_array_standard_order_1d)

  # set up new shape (same as was in the numpy array after transposing it)
  from scitbx.array_family.flex import grid
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

