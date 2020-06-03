from __future__ import absolute_import, division, print_function
from libtbx.math_utils import ifloor, iceil
from cctbx import crystal, maptbx, uctbx
from scitbx.array_family import flex
from libtbx.utils import null_out
from libtbx import group_args
import iotbx.map_manager
import mmtbx.model

class around_model(object):
  """
  Extract map box around atomic model. Box is in P1 and has origin at (0,0,0).

  Input map_manager must have origin (working position) at (0,0,0)
  Input model coordinates must correspond to working position of map_manager

  Wrapping must be specified on initialization
    wrapping=True means that grid points outside of the unit cell can be
    mapped inside with unit translations and box can effectively come
    from anywhere in space.
    If wrapping=True, the supplied box must be a complete unit cell so that
    map_data.all() must be the same as unit_cell_grid.

   wrapping=False means that grid points outside of the unit cell are
    undefined.  If a box is specified that uses points outside the defined
    region, those points are set to zero.

  """
  def __init__(self, map_manager=None, model=None, cushion=None,
      wrapping=None):
    self.map_manager = map_manager
    self.model = model
    self.wrapping = wrapping

    self.basis_for_boxing_string='using model, wrapping=%s'  %(wrapping)

    # safeguards
    assert wrapping is not None
    assert isinstance(map_manager, iotbx.map_manager.map_manager)
    assert isinstance(model, mmtbx.model.manager)
    assert self.map_manager.map_data().accessor().origin() == (0,0,0)
    assert map_manager.crystal_symmetry().is_similar_symmetry(
      model.crystal_symmetry())
    assert cushion >= 0

    if wrapping:
      assert map_manager.unit_cell_grid==map_manager.map_data().all()

    # get items needed to do the shift
    cs = map_manager.crystal_symmetry()
    uc = cs.unit_cell()
    sites_frac = model.get_sites_frac()
    map_data = map_manager.map_data()
    # convert cushion into fractional vector
    cushion_frac = flex.double(uc.fractionalize((cushion,)*3))
    # find fractional corners
    frac_min = sites_frac.min()
    frac_max = sites_frac.max()
    frac_max = list(flex.double(frac_max)+cushion_frac)
    frac_min = list(flex.double(frac_min)-cushion_frac)

    # find corner grid nodes
    all_orig = map_data.all()
    self.gridding_first = [ifloor(f*n) for f,n in zip(frac_min, all_orig)]
    self.gridding_last  = [ iceil(f*n) for f,n in zip(frac_max, all_orig)]

    # Ready with gridding...set up shifts and box crystal_symmetry
    self.set_shifts_and_crystal_symmetry()

  def set_shifts_and_crystal_symmetry(self):
    # set items needed to do the shift

    cs = self.map_manager.crystal_symmetry()
    uc = cs.unit_cell()

    self.box_all=[j-i+1 for i,j in zip(self.gridding_first,self.gridding_last)]
    # get shift vector as result of boxing
    all_orig = self.map_manager.map_data().all()
    self.shift_frac = [-self.gridding_first[i]/all_orig[i] for i in range(3)]
    self.shift_cart = cs.unit_cell().orthogonalize(self.shift_frac)
    # get crystal symmetry of the box
    p = uc.parameters()
    abc = [p[i] * self.box_all[i]/all_orig[i] for i in range(3)]
    box_uc = uctbx.unit_cell(parameters=(abc[0],abc[1],abc[2],p[3],p[4],p[5]))
    self.crystal_symmetry = crystal.symmetry(
      unit_cell = box_uc, space_group="P1")

  def apply_to_map(self, map_manager=None):
    if(map_manager is None):
      # Apply to data already present
      map_manager = self.map_manager
    else:
      # Apply to a map_manager that is similar to the one used to generate
      #   this around_model object
      assert map_manager.crystal_symmetry().is_similar_symmetry(
        self.map_manager.crystal_symmetry())

      ma1 = map_manager.map_data().accessor()
      ma2 = self.map_manager.map_data().accessor()
      assert ma1.all()    == ma2.all()
      assert ma1.origin() == ma2.origin()
      assert ma1.focus()  == ma2.focus()

    map_data = map_manager.map_data()

    # Check if map is all valid
    bounds_info=get_bounds_of_valid_region(map_data=map_data,
      gridding_first=self.gridding_first,
      gridding_last=self.gridding_last)

    if self.wrapping or bounds_info.inside_allowed_bounds:
      # Just copy everything
      map_box = maptbx.copy(map_data, self.gridding_first, self.gridding_last)
      # Note: map_box gridding is self.gridding_first to self.gridding_last

    else: # Need to copy and then zero outside of defined region
      map_box=copy_and_zero_map_outside_bounds(map_data=map_data,
         bounds_info=bounds_info)

    #  Now reshape map_box to put origin at (0,0,0)
    map_box.reshape(flex.grid(self.box_all))

    # Create new map_manager object:
    #   Use original values for:
    #     unit_cell_grid    (gridding of original full unit cell)
    #     unit_cell_crystal_symmetry  (symmetry of original full unit cell)
    #     input_file_name

    #   Use new (boxed) values for:
    #     map_data
    #     crystal_symmetry   (symmetry of the part of the map that is present)

    #   Update:
    #     origin_shift_grid_units  (position in the original map of the
    #                                 (0,0,0) grid point in map_box)
    #     labels  (add label specifying boxing operation)

    # New origin_shift_grid_units:
    origin_shift_grid_units=[
       self.gridding_first[i]+map_manager.origin_shift_grid_units[i]
        for i in range(3)]

    # New labels:
    new_label="Boxed %s to %s %s" %(
      str(tuple(self.gridding_first)),str(tuple(self.gridding_last)),
      self.basis_for_boxing_string)

    #  Set up new map_manager.
    #  NOTE: origin_shift_grid_units is required as bounds have changed
    new_map_manager=map_manager.customized_copy(map_data=map_box,
      origin_shift_grid_units=origin_shift_grid_units)

    # Add the label
    new_map_manager.add_label(new_label)

    return new_map_manager

  def apply_to_model(self, model=None):
    if(model is None):
      # Apply to model already present
      model = self.model
    else:
      # Apply to similar model
      assert model.crystal_symmetry().is_similar_symmetry(
        self.map_manager.crystal_symmetry())

    box_uc = self.crystal_symmetry.unit_cell()
    uc     = model.crystal_symmetry().unit_cell()
    sites_frac = model.get_sites_frac()
    sites_frac_new = box_uc.fractionalize(
      uc.orthogonalize(sites_frac+self.shift_frac))
    sites_cart_new = box_uc.orthogonalize(sites_frac_new)
    pdb_hierarchy = model.get_hierarchy().deep_copy()
    pdb_hierarchy.atoms().set_xyz(sites_cart_new)
    return mmtbx.model.manager(
      model_input      = None,
      pdb_hierarchy    = pdb_hierarchy,
      crystal_symmetry = self.crystal_symmetry,
      log              = null_out())

class with_bounds(around_model):
  """
  Extract map box using specified lower_bounds and upper_bounds

  Box is in P1 and has origin at (0,0,0).

  Input map_manager must have origin (working position) at (0,0,0)
  Input model coordinates must correspond to working position of map_manager

  Wrapping must be specified on initialization
    wrapping=True means that grid points outside of the unit cell can be
    mapped inside with unit translations and box can effectively come
    from anywhere in space.
    If wrapping=True, the supplied box must be a complete unit cell so that
    map_data.all() must be the same as unit_cell_grid.

   wrapping=False means that grid points outside of the unit cell are
    undefined.  If a box is specified that uses points outside the defined
    region, those points are set to zero.

  """
  def __init__(self, map_manager=None, lower_bounds=None, upper_bounds=None,
      wrapping=None):
    self.map_manager = map_manager
    self.wrapping = wrapping
    self.basis_for_boxing_string='supplied bounds, wrapping=%s' %(wrapping)

    # safeguards
    assert wrapping is not None
    assert isinstance(map_manager, iotbx.map_manager.map_manager)
    assert self.map_manager.map_data().accessor().origin() == (0,0,0)
    assert lower_bounds is not None
    assert upper_bounds is not None

    if wrapping:
      assert map_manager.unit_cell_grid==map_manager.map_data().all()

    self.gridding_first = lower_bounds
    self.gridding_last  = upper_bounds

    # Ready with gridding...set up shifts and box crystal_symmetry
    self.set_shifts_and_crystal_symmetry()


def get_bounds_of_valid_region(map_data=None,
    gridding_first=None,
    gridding_last=None):
  '''
    If map_data is sampled from gridding_first to gridding_last with
    maptbx.copy, (1) does the sampling go outside of map_data?
    (2) What are the lower and upper bounds of the valid region of the resulting
    map?
  '''

  lower_allowed_bounds=[]
  upper_allowed_bounds=[]
  inside_allowed_bounds=True
  for o,a,f,l in zip(map_data.origin(),map_data.all(),
     gridding_first,gridding_last):
    # Available map goes from o to o+a
    # Requested map goes from f to l

    # If f is less than o, first valid grid point is o.
    # Otherwise, first valid grid point is f
    # so first valid grid point is max(f,o)

    # If f is less than o
    #   After shifting origin to (0,0,0) first valid grid point is o-f
    # Otherwise, after shifting origin, first valid grid point is 0.

    #  If l is >= a+o, last valid grid point is a+o-1
    #  Otherwise last valid grid point is l
    # So last valid grid point is min(l,a+o-1)

    #  If l is >= a+o
    #    After shifting origin to (0,0,0), last valid grid point is a+o-1-f
    #  Otherwise, after shifting origin, last valid grid point is l-f


    lower_allowed_bounds.append(max(o,f)) # first valid grid point
    upper_allowed_bounds.append(min(l,o+a-1)) # last valid grid point
    if f < o or l >= a+o:
      inside_allowed_bounds=False

  return group_args(
     lower_allowed_bounds=lower_allowed_bounds,
     upper_allowed_bounds=upper_allowed_bounds,
     gridding_first=gridding_first,
     gridding_last=gridding_last,
     inside_allowed_bounds=inside_allowed_bounds)

def copy_and_zero_map_outside_bounds(map_data=None,bounds_info=None):
  '''
     Copy part of a map and zero outside valid region

     Goes with get_bounds_of_valid_region

     First copy requested part of map, wrapping if request goes outside
     of supplied map.

     Then zero out everything that was from outside the available region.

     Returns map with bounds from bounds_info.gridding_first to
     bounds_info.gridding_last, but everything outside of
     lower_allowed_bounds to upper_allowed_bounds is zeroed out.
  '''

  # First copy the entire requested region:

  map_copy=maptbx.copy(map_data,
      bounds_info.gridding_first, bounds_info.gridding_last)
  # Note: the origin of map_copy is at gridding_first and goes to last

  # Make sure we are working with a flex.double array
  if type(map_copy) != type(flex.double()): # must be double
    map_copy=map_copy.as_double()

  # Now zero out everything outside of the valid region


  # Make sure this map matches the bounds_info
  assert tuple(map_copy.origin())==tuple(bounds_info.gridding_first)

  # We are going to shift the origin in this copy, so shift the bounds to match

  lower_bounds_after_shift=[]
  upper_bounds_after_shift=[]
  for l,u,f in zip(bounds_info.lower_allowed_bounds,
      bounds_info.upper_allowed_bounds,
      bounds_info.gridding_first):
    lower_bounds_after_shift.append(l-f)
    upper_bounds_after_shift.append(u-f)

  acc=map_copy.accessor() # save where the origin is

  map_copy=map_copy.shift_origin()  # put origin at (0,0,0)
  map_copy_all=map_copy.all() # save size of map
  # XXX work-around for set_box does not allow offset origin
  map_copy.resize(flex.grid(map_copy_all))
  new_map=maptbx.set_box_copy_inside(0,  # copies inside, zero outside bounds
    map_data_to   = map_copy,
    start         = tuple(lower_bounds_after_shift),
    end           = tuple(upper_bounds_after_shift))
  # XXX and shift map back
  new_map=new_map.as_1d()
  new_map.reshape(acc)
  return new_map
