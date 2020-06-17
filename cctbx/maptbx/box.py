from __future__ import absolute_import, division, print_function
from libtbx.math_utils import ifloor, iceil
from libtbx.utils import Sorry
from cctbx import crystal, maptbx, uctbx
from scitbx.array_family import flex
from libtbx import group_args
import iotbx.map_manager
import mmtbx.model
import sys

class with_bounds(object):
  """
  Extract map box using specified lower_bounds and upper_bounds

  NOTE: changes supplied model and map_manager in place. Call with deep_copy()
  versions if you do not want them modified.

  Output versions of map_manager and model are in P1 and have origin at (0,0,0).

  Bounds refer to grid position in this box with origin at (0,0,0)

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
  def __init__(self, map_manager, lower_bounds, upper_bounds, wrapping,
     model=None,
     ncs_object=None,
     log=sys.stdout):
    self.lower_bounds=lower_bounds
    self.upper_bounds=upper_bounds
    self.wrapping=wrapping
    self._map_manager=map_manager
    self._model=model
    self._ncs_object=ncs_object

    # safeguards
    assert lower_bounds is not None
    assert upper_bounds is not None
    assert upper_bounds is not None
    assert isinstance(wrapping, bool)
    assert isinstance(map_manager, iotbx.map_manager.map_manager)
    assert self._map_manager.map_data().accessor().origin() == (0,0,0)
    if model:
      assert isinstance(model, mmtbx.model.manager)
      assert map_manager.crystal_symmetry().is_similar_symmetry(
        model.crystal_symmetry())
    if ncs_object:
      assert isinstance(ncs_object, mmtbx.ncs.ncs.ncs)
    if wrapping:
      assert map_manager.unit_cell_grid==map_manager.map_data().all()

    self.basis_for_boxing_string='supplied bounds, wrapping=%s' %(wrapping)

    # These are lower and upper bounds of map with origin at (0,0,0)
    #   (not the original map)

    self.gridding_first = lower_bounds
    self.gridding_last  = upper_bounds

    # Ready with gridding...set up shifts and box crystal_symmetry
    self.set_shifts_and_crystal_symmetry()

    # Apply boxing to model, ncs, and map (if available)
    self.apply_to_model_ncs_and_map()

  def as_map_and_model(self):
    '''
      Return map_and_model object with contents of this class (not a deepcopy)

      Sets wrapping=False always because this class does not return a full
      unit cell.

      Sets box=False or recursion will occur.

    '''
    import iotbx.map_and_model
    mam=iotbx.map_and_model.input(
        map_manager=self.map_manager(),
        model=self.model(),
        ncs_object=self.ncs_object(),
        wrapping=False, # Always
        box=False # Always
        )
    # Keep track of the gridding and solvent_content (if used) in this boxing.
    mam.set_gridding_first(self.gridding_first)
    mam.set_gridding_last(self.gridding_last)
    if hasattr(self,'solvent_content') and self.solvent_content is not None:
      mam.set_solvent_content(self.solvent_content)
    return mam

  def model(self):
    return getattr(self,'_model',None)

  def map_manager(self):
    return getattr(self,'_map_manager',None)

  def ncs_object(self):
    return getattr(self,'_ncs_object',None)

  def set_shifts_and_crystal_symmetry(self):
    '''
      Set items needed to do the shift

      Here self.gridding_first and self.gridding_last are the grid points
      marking the start and end, in the map with origin at (0,0,0), of the
      region to be kept.
    '''

    # Save to check later if another map_manager is used as input
    self.accessor_at_initialization=self._map_manager.map_data().accessor()

    full_cs = self._map_manager.unit_cell_crystal_symmetry()
    full_uc = full_cs.unit_cell()
    self.box_all=[j-i+1 for i,j in zip(self.gridding_first,self.gridding_last)]
    # get shift vector as result of boxing
    full_all_orig = self._map_manager.unit_cell_grid
    self.shift_frac = \
        [-self.gridding_first[i]/full_all_orig[i] for i in range(3)]
    self.shift_cart = full_cs.unit_cell().orthogonalize(self.shift_frac)
    # get crystal symmetry of the box
    p = full_uc.parameters()
    abc = [p[i] * self.box_all[i]/full_all_orig[i] for i in range(3)]
    box_uc = uctbx.unit_cell(parameters=(abc[0],abc[1],abc[2],p[3],p[4],p[5]))
    self.crystal_symmetry = crystal.symmetry(
      unit_cell = box_uc, space_group="P1")

  def apply_to_model_ncs_and_map(self):
    '''
    Apply boxing to to self._model, self._ncs_object, self._map_manager
    so all are boxed

    Apply to self._map_manager last (so self._map_manager.crystal symmetry()
      is not changed until here and is used to compare to symmetries of other
      objects
    '''

    if self._model:
      self._model=self.apply_to_model(self._model)

    if self._ncs_object:
      self._ncs_object=self.apply_to_ncs_object(self._ncs_object)

    # Must be last
    self._map_manager=self.apply_to_map(self._map_manager)


  def apply_to_map(self, map_manager):
    '''
     Apply boxing to a map_manager that is similar to the one used to generate
       this around_model object
    '''
    assert isinstance(map_manager, iotbx.map_manager.map_manager)
    assert map_manager.unit_cell_crystal_symmetry().is_similar_symmetry(
      self._map_manager.unit_cell_crystal_symmetry())

    ma1 = map_manager.map_data().accessor()
    ma2 = self.accessor_at_initialization
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
    elif not bounds_info.some_valid_points:
      # Just copy everything and zero
      map_box = maptbx.copy(map_data, self.gridding_first, self.gridding_last)
      map_box = map_box * 0.

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
    #
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

  def apply_to_model(self, model):
    '''
       Apply boxing to a model that is similar to the one used to generate
       this around_model object
    '''

    assert isinstance(model, mmtbx.model.manager)

    # Allow models where either original or current symmetry
    #  match this object's original symmetry

    assert self._map_manager.is_compatible_model(model)

    # Shift the model and add self.shift_cart on to whatever shift was there
    model.shift_model_and_set_crystal_symmetry(
       shift_cart=self.shift_cart, # shift to apply
       crystal_symmetry=self.crystal_symmetry, # new crystal_symmetry
       )

    return model

  def apply_to_ncs_object(self, ncs_object):
    '''
      Apply shifts from this boxing to an ncs_object
       ncs does keep track of shifts
    '''

    return ncs_object.coordinate_offset(coordinate_offset=self.shift_cart)


class around_model(with_bounds):
  """
  Extract map box around atomic model. Box is in P1 and has origin at (0,0,0).

  NOTE: changes supplied model and map_manager in place. Call with deep_copies
  if you do not want them modified.

  Input map_manager must have origin (working position) at (0,0,0)
  Input model coordinates must correspond to working position of map_manager

  On initialization new bounds and crystal_symmetry are identified.
  Then map_manager is replaced with boxed version and model is shifted in place.

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
  def __init__(self, map_manager, model, cushion, wrapping,
      ncs_object=None,
      log=sys.stdout):

    self.wrapping=wrapping
    self._map_manager=map_manager
    self._model=model
    self._ncs_object=ncs_object

    self.basis_for_boxing_string='using model, wrapping=%s'  %(wrapping)

    # safeguards
    assert isinstance(wrapping, bool)
    assert isinstance(map_manager, iotbx.map_manager.map_manager)
    assert isinstance(model, mmtbx.model.manager)
    assert self._map_manager.map_data().accessor().origin() == (0,0,0)
    if ncs_object:
      assert isinstance(ncs_object, mmtbx.ncs.ncs.ncs)

    # Make sure working model and map_manager crystal_symmetry match

    assert map_manager.crystal_symmetry().is_similar_symmetry(
       model.crystal_symmetry())

    assert cushion >= 0

    if wrapping:  # map must be entire unit cell
      assert map_manager.unit_cell_grid==map_manager.map_data().all()

    # NOTE: We are going to use crystal_symmetry and sites_frac based on
    #   the map_manager (the model could still have different crystal_symmetry)

    # get items needed to do the shift
    cs = map_manager.crystal_symmetry()
    uc = cs.unit_cell()
    sites_cart = model.get_sites_cart()
    sites_frac = uc.fractionalize(sites_cart)
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

    # Apply boxing to model, ncs, and map (if available)
    self.apply_to_model_ncs_and_map()

class extract_unique(with_bounds):

  '''
  Identify unique part of density in a map (using ncs object) and create a
  new map_manager containing this box of density, masked around regions
  containing density.  Note: the map may be masked between nearby
  density regions so this map could have many discontinuities.

  NOTE: changes supplied model and map_manager in place. Call with deep_copies
  if you do not want them modified.
  '''

  def __init__(self, map_manager,
    wrapping=None,
    ncs_object=None,
    model=None,
    target_ncs_au_model=None,
    regions_to_keep=None,
    solvent_content=None,
    resolution=None,
    sequence=None,
    molecular_mass=None,
    symmetry=None,
    chain_type='PROTEIN',
    keep_low_density=True,  # default from map_box
    box_buffer=5,
    soft_mask_extract_unique=True,
    mask_expand_ratio=1,
    log=sys.stdout
    ):

    self.wrapping=wrapping
    self._map_manager=map_manager
    self._model=model
    self._ncs_object=ncs_object

    assert isinstance(wrapping, bool)
    assert isinstance(map_manager, iotbx.map_manager.map_manager)
    assert self._map_manager.map_data().accessor().origin() == (0,0,0)
    assert resolution
    if ncs_object:
      assert isinstance(ncs_object, mmtbx.ncs.ncs.ncs)
    if model:
      assert isinstance(model, mmtbx.model.manager)
      assert map_manager.crystal_symmetry().is_similar_symmetry(
        model.crystal_symmetry())
    if wrapping:  # map must be entire unit cell
      assert map_manager.unit_cell_grid==map_manager.map_data().all()

    # Get crystal_symmetry
    self.crystal_symmetry=map_manager.crystal_symmetry()
    # Convert to map_data
    self.map_data=map_manager.map_data()

    from cctbx.maptbx.segment_and_split_map import run as segment_and_split_map
    assert self.map_data.origin()==(0,0,0)
    args=[]

    ncs_group_obj,remainder_ncs_group_obj,tracking_data =\
      segment_and_split_map(args,
        map_data=self._map_manager.map_data(),
        crystal_symmetry=self.crystal_symmetry,
        ncs_obj=self._ncs_object,
        target_model=target_ncs_au_model,
        write_files=False,
        auto_sharpen=False,
        add_neighbors=False,
        density_select=False,
        save_box_map_ncs_au=True,
        resolution=resolution,
        solvent_content=solvent_content,
        chain_type=chain_type,
        sequence=sequence,
        molecular_mass=molecular_mass,
        symmetry=symmetry,
        keep_low_density=keep_low_density,
        regions_to_keep=regions_to_keep,
        box_buffer=box_buffer,
        soft_mask_extract_unique=soft_mask_extract_unique,
        mask_expand_ratio=mask_expand_ratio,
        out=log)

    from scitbx.matrix import col

    if not hasattr(tracking_data,'box_map_ncs_au_map_data'):
      raise Sorry(" Extraction of unique part of map failed...")

    ncs_au_map_data=tracking_data.box_map_ncs_au_map_data
    ncs_au_crystal_symmetry=tracking_data.box_map_ncs_au_crystal_symmetry

    lower_bounds=ncs_au_map_data.origin()
    upper_bounds=tuple(
      col(ncs_au_map_data.focus())-col((1,1,1)))
    print("\nBounds for unique part of map: %s to %s " %(
     str(lower_bounds),str(upper_bounds)),file=log)

    # shift the map so it is in the same position as the box map will be in
    ncs_au_map_data.reshape(flex.grid(ncs_au_map_data.all()))
    assert col(ncs_au_map_data.all())==\
        col(upper_bounds)-col(lower_bounds)+col((1,1,1))

    self.basis_for_boxing_string='extract unique, wrapping=%s' %(self.wrapping)

    self.gridding_first = lower_bounds
    self.gridding_last  = upper_bounds

    # Ready with gridding...set up shifts and box crystal_symmetry
    self.set_shifts_and_crystal_symmetry()

    # Apply boxing to model, ncs, and map (if available)
    self.apply_to_model_ncs_and_map()

    # And replace map data with boxed asymmetric unit
    self._map_manager.set_map_data(map_data=ncs_au_map_data)
    self._map_manager.add_limitation("extract_unique")

class around_mask(with_bounds):
  """
  Extract map box around masked region of a map that represents a mask
  You need to supply the mask as a map_manager object

  Box is in P1 and has origin at (0,0,0).

  Input map_manager must have origin (working position) at (0,0,0)

  Returns boxed version of map_manager supplied.  Object will contain
  the boxed version of mask as self.mask_as_map_manager.

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
  def __init__(self, map_manager, wrapping,
     mask_as_map_manager=None,
     model=None,
     ncs_object=None,
     box_cushion=3,
     log=sys.stdout):

    self.wrapping=wrapping
    self._map_manager=map_manager
    self._model=model
    self._ncs_object=ncs_object

    # safeguards
    assert isinstance(wrapping, bool)
    assert isinstance(map_manager, iotbx.map_manager.map_manager)
    assert isinstance(mask_as_map_manager, iotbx.map_manager.map_manager)
    assert self._map_manager.map_data().accessor().origin() == (0,0,0)
    assert map_manager.is_similar(mask_as_map_manager)
    if wrapping:
      assert map_manager.unit_cell_grid==map_manager.map_data().all()

    self.basis_for_boxing_string='around_mask bounds, wrapping=%s' %(wrapping)

    # Get a connectivity object that marks all the connected regions in map

    from cctbx.maptbx.segment_and_split_map import get_co
    co,sorted_by_volume,min_b,max_b=get_co(
       map_data=mask_as_map_manager.map_data(),
       threshold=0.5,
       wrapping=False)

    if len(sorted_by_volume)<2:  # didn't work
      raise Sorry("No mask obtained...")

    # Get the biggest connected region in the map

    original_id_from_id={}
    for i in range(1,len(sorted_by_volume)):
      v,id=sorted_by_volume[i]
      original_id_from_id[i]=id
    id=1
    orig_id=original_id_from_id[id]

    # Get lower and upper bounds of this region in grid units

    self.gridding_first = min_b[orig_id]
    self.gridding_last  = max_b[orig_id]

    # Increase range of bounds by box_cushion
    cs=map_manager.crystal_symmetry()
    cushion=flex.double(cs.unit_cell().fractionalize((box_cushion,)*3))
    all_orig = map_manager.map_data().all()
    self.gridding_first = [max(0,ifloor((gf/n-c)*n)) for c,gf,n in zip(
       cushion, self.gridding_first, all_orig)]
    self.gridding_last  = [ min(n-1,iceil((gf+c)*n)) for c,gf,n in zip(
       cushion, self.gridding_last, all_orig)]


    # Ready with gridding...set up shifts and box crystal_symmetry
    self.set_shifts_and_crystal_symmetry()

    self.apply_to_model_ncs_and_map()

    # Also apply to mask_as_map_manager so that mask_as_map_manager is boxed
    mask_as_map_manager=self.apply_to_map(mask_as_map_manager)
    self.mask_as_map_manager=mask_as_map_manager # save it


class around_density(with_bounds):
  """
  Extract map box around region of a map containing density

  Box is in P1 and has origin at (0,0,0).

  Input map_manager must have origin (working position) at (0,0,0)

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
  def __init__(self, map_manager, wrapping,
     threshold=0.05,
     box_cushion=3.,
     get_half_height_width=True,
     model=None,
     ncs_object=None,
     log=sys.stdout):

    self.wrapping=wrapping
    self._map_manager=map_manager
    self._model=model
    self._ncs_object=ncs_object

    # safeguards
    assert threshold is not None
    assert box_cushion is not None
    assert isinstance(wrapping, bool)
    assert isinstance(map_manager, iotbx.map_manager.map_manager)
    assert self._map_manager.map_data().accessor().origin() == (0,0,0)
    if wrapping:
      assert map_manager.unit_cell_grid==map_manager.map_data().all()

    self.basis_for_boxing_string='around_density, wrapping=%s' %(wrapping)


    # Select box where data are positive (> threshold*max)
    map_data=map_manager.map_data()
    origin=list(map_data.origin())
    assert origin==[0,0,0]
    all=list(map_data.all())
    # Get max value vs x,y,z
    value_list=flex.double()
    for i in range(0,all[0]):
      new_map_data = maptbx.copy(map_data,
         tuple((i,0,0)),
         tuple((i,all[1],all[2]))
       )
      value_list.append(new_map_data.as_1d().as_double().min_max_mean().max)
    ii=0
    for z in value_list:
      ii+=1
    x_min,x_max=get_range(value_list,threshold=threshold,
      get_half_height_width=get_half_height_width)

    value_list=flex.double()
    for j in range(0,all[1]):
      new_map_data = maptbx.copy(map_data,
         tuple((0,j,0)),
         tuple((all[0],j,all[2]))
       )
      value_list.append(new_map_data.as_1d().as_double().min_max_mean().max)
    y_min,y_max=get_range(value_list,threshold=threshold,
      get_half_height_width=get_half_height_width)

    value_list=flex.double()
    for j in range(0,all[1]):
      new_map_data = maptbx.copy(map_data,
         tuple((0,j,0)),
         tuple((all[0],j,all[2]))
       )
      value_list.append(new_map_data.as_1d().as_double().min_max_mean().max)
    y_min,y_max=get_range(value_list,threshold=threshold,
      get_half_height_width=get_half_height_width)

    value_list=flex.double()
    for k in range(0,all[2]):
      new_map_data = maptbx.copy(map_data,
         tuple((0,0,k)),
         tuple((all[0],all[1],k))
       )
      value_list.append(new_map_data.as_1d().as_double().min_max_mean().max)
    z_min,z_max=get_range(value_list,threshold=threshold,
      get_half_height_width=get_half_height_width)

    # Get lower and upper bounds of this region in grid units
    frac_min=(x_min,y_min,z_min)
    frac_max=(x_max,y_max,z_max)
    cs=map_manager.crystal_symmetry()
    cushion=flex.double(cs.unit_cell().fractionalize((box_cushion,)*3))
    all_orig = map_data.all()
    self.gridding_first = [max(0,ifloor((f-c)*n)) for c,f,n in zip(
       cushion, frac_min, all_orig)]
    self.gridding_last  = [ min(n-1,iceil((f+c)*n)) for c,f,n in zip(
       cushion, frac_max, all_orig)]

    # Ready with gridding...set up shifts and box crystal_symmetry
    self.set_shifts_and_crystal_symmetry()

    # Apply boxing to model, ncs, and map (if available)
    self.apply_to_model_ncs_and_map()

def get_range(value_list, threshold=None, ignore_ends=True,
     keep_near_ends_frac=0.02, half_height_width=2., get_half_height_width=None,
     cutoff_ratio=4,ratio_max=0.5): # XXX May need to set cutoff_ratio and
    #  ratio_max lower.
    # ignore ends allows ignoring the first and last points which may be off
    # if get_half_height_width, find width at half max hieght, go
    #  half_height_width times this width out in either direction, use that as
    #  baseline instead of full cell. Don't do it if the height at this point
    #  is over cutoff_ratio times threshold above original baseline.
    if get_half_height_width:
      z_min,z_max=get_range(value_list,threshold=0.5,
        ignore_ends=ignore_ends,keep_near_ends_frac=keep_near_ends_frac,
        get_half_height_width=False)
      z_mid=0.5*(z_min+z_max)
      z_width=0.5*(z_max-z_min)
      z_low=z_mid-2*z_width
      z_high=z_mid+2*z_width
      if ignore_ends:
        i_max=value_list.size()-2
        i_min=1
      else:
        i_max=value_list.size()-1
        i_min=0

      i_low= max(i_min,min(i_max,int(0.5+z_low* value_list.size())))
      i_high=max(i_min,min(i_max,int(0.5+z_high*value_list.size())))
      min_value=value_list.min_max_mean().min
      max_value=value_list.min_max_mean().max
      ratio_low=(value_list[i_low]-min_value)/max(
         1.e-10,(max_value-min_value))
      ratio_high=(value_list[i_high]-min_value)/max(
         1.e-10,(max_value-min_value))
      if ratio_low <= cutoff_ratio*threshold and ratio_low >0 \
           and ratio_low<ratio_max\
           and ratio_high <= cutoff_ratio*threshold and ratio_high > 0 \
           and ratio_high < ratio_max:
        ratio=min(ratio_low,ratio_high)
        z_min,z_max=get_range(
          value_list,threshold=threshold+ratio,
          ignore_ends=ignore_ends,keep_near_ends_frac=keep_near_ends_frac,
          get_half_height_width=False)
        return z_min,z_max
      else:
        z_min,z_max=get_range(value_list,threshold=threshold,
          ignore_ends=ignore_ends,keep_near_ends_frac=keep_near_ends_frac,
          get_half_height_width=False)
        return z_min,z_max

    if threshold is None: threshold=0
    n_tot=value_list.size()
    assert n_tot>0
    min_value=value_list.min_max_mean().min
    max_value=value_list.min_max_mean().max
    cutoff=min_value+(max_value-min_value)*threshold
    if ignore_ends:
      i_off=1
    else:
      i_off=0
    i_low=None
    for i in range(i_off,n_tot-i_off):
      if value_list[i]>cutoff:
        i_low=max(i_off,i-1)
        break
    i_high=None
    for i in range(i_off,n_tot-i_off):
      ii=n_tot-1-i
      if value_list[ii]>cutoff:
        i_high=min(n_tot-1-i_off,ii+1)
        break
    if i_low is None or i_high is None:
      raise Sorry("Cannot auto-select region...")
    if i_low/n_tot<keep_near_ends_frac: i_low=0
    if (n_tot-1-i_high)/n_tot<keep_near_ends_frac: i_high=n_tot-1
    return i_low/n_tot,i_high/n_tot


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
  some_valid_points=True
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


    first_valid=max(o,f) # first valid grid point
    last_valid=min(l,o+a-1) # last valid grid point
    lower_allowed_bounds.append(first_valid)
    upper_allowed_bounds.append(last_valid)
    if f < o or l >= a+o:
      inside_allowed_bounds=False
    if last_valid <= first_valid:
      some_valid_points=False

  return group_args(
     lower_allowed_bounds=lower_allowed_bounds,
     upper_allowed_bounds=upper_allowed_bounds,
     gridding_first=gridding_first,
     gridding_last=gridding_last,
     inside_allowed_bounds=inside_allowed_bounds,
     some_valid_points=some_valid_points)

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
