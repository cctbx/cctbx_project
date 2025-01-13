from __future__ import absolute_import, division, print_function
from libtbx.math_utils import ifloor, iceil
from libtbx.utils import Sorry, null_out
from cctbx import crystal, maptbx, uctbx
from scitbx.array_family import flex
from libtbx import group_args
import iotbx.map_manager
import mmtbx.model
import sys

class with_bounds(object):
  """
  Extract map box using specified lower_bounds and upper_bounds

  Creates new map_manager and modifies model in place

  Input map_manager must have origin (working position) at (0, 0, 0)
  Input model coordinates must correspond to working position of map_manager

  On initialization new bounds and crystal_symmetry are identified.
  Then map_manager is replaced with boxed version and shifted model is created

  Output versions of map_manager and model are in P1 and have origin at (0, 0, 0).
  Bounds refer to grid position in this box with origin at (0, 0, 0)

  Wrapping:
    wrapping = True means that grid points outside of the unit cell can be
    mapped inside with unit translations and box can effectively come
    from anywhere in space.
    If wrapping = True, the supplied box must be a complete unit cell so that
    map_data.all() must be the same as unit_cell_grid.

   wrapping = False means that grid points outside of the unit cell are
    undefined.  If a box is specified that uses points outside the defined
    region, those points are set to zero.

  if use_cubic_boxing, adjust bounds to make a cubic box by making box bigger.
      if this conflicts with stay_inside_current_map, make box smaller. Normally
      this option should not be used with with_bounds because the bounds should
      be directly specified.
  Note: stay_inside_current_map applies only when using cubic boxing
  If require_match_unit_cell_crystal_symmetry:  require that unit_cell
    crystal symmetry of map and model match

  """
  def __init__(self,
     map_manager,
     lower_bounds,
     upper_bounds,
     model = None,
     wrapping = None,
     model_can_be_outside_bounds = False,
     stay_inside_current_map = True,
     use_cubic_boxing = False,
     require_match_unit_cell_crystal_symmetry = False,
     log = sys.stdout):

    self.require_match_unit_cell_crystal_symmetry = \
       require_match_unit_cell_crystal_symmetry
    self.lower_bounds = lower_bounds
    self.upper_bounds = upper_bounds
    self._map_manager = map_manager
    self._model = model
    self.model_can_be_outside_bounds = model_can_be_outside_bounds
    self.use_cubic_boxing = use_cubic_boxing
    self._info = None


    # safeguards
    assert lower_bounds is not None
    assert upper_bounds is not None
    assert len(tuple(lower_bounds))==3
    assert len(tuple(upper_bounds))==3
    for i in range(3):
      assert list(upper_bounds)[i] > list(lower_bounds)[i]

    assert isinstance(map_manager, iotbx.map_manager.map_manager)

    assert self._map_manager.map_data().accessor().origin()  ==  (0, 0, 0)
    if model is not None:
      assert isinstance(model, mmtbx.model.manager)
      assert map_manager.is_compatible_model(model,
        require_match_unit_cell_crystal_symmetry=True)

    self._force_wrapping = wrapping
    if wrapping is None:
      wrapping = self.map_manager().wrapping()
    self.basis_for_boxing_string = 'supplied bounds, wrapping = %s' %(
      wrapping)

    # These are lower and upper bounds of map with origin at (0, 0, 0)
    #   (not the original map)

    self.gridding_first = lower_bounds
    self.gridding_last  = upper_bounds

    # Adjust gridding to make a cubic box if desired
    if self.use_cubic_boxing:
      self.set_cubic_boxing(stay_inside_current_map = stay_inside_current_map)

    # Ready with gridding...set up shifts and box crystal_symmetry
    self.set_shifts_and_crystal_symmetry()

    # Apply boxing to model, ncs, and map (if available)
    self.apply_to_model_ncs_and_map()

  def as_map_model_manager(self):
    '''
      Return map_model_manager object with contents of this class (not a deepcopy)

    '''
    from iotbx.map_model_manager import map_model_manager
    mmm = map_model_manager(
        map_manager = self.map_manager(),
        model = self.model(),
        )
    # Keep track of the gridding in this boxing.
    mmm.set_gridding_first(self.gridding_first)
    mmm.set_gridding_last(self.gridding_last)
    return mmm

  def model(self):
    return self._model

  def map_manager(self):
    return self._map_manager

  def info(self):
    return self._info

  def set_info(self, info):
    ''' Holder for anything you want.
      Usually set it to a group_args object:
         self.set_info(group_args(group_args_type='my_group_args', value = value))
    '''

    self._info = info

  def ncs_object(self):
    if self.map_manager():
      return self.map_manager().ncs_object()
    else:
      return None

  def set_cubic_boxing(self, stay_inside_current_map = None,
    require_even_gridding = True,
    log = sys.stdout):
    ''' Adjust bounds to make a cubic box.
      Adjust bounds to make a cubic box by making box bigger.
      if this conflicts with stay_inside_current_map, make box smaller. Normally
      this option should not be used with with_bounds because the bounds should
      be directly specified.
      Normally creates a box with an even number of grid points '''

    if not self.use_cubic_boxing:
      return  # nothing to do

    print("\nSetting up cubic box", file = log)
    map_all = self._map_manager.map_data().all()
    lmn = [1 + b - a for a,b in zip(self.gridding_first, self.gridding_last)]
    if lmn[0] == lmn[1] and lmn[0] == lmn[2] and (
       is_even(lmn[0]) or (not require_even_gridding)):
      print("Box is already cubic with dimensions ",lmn, file = log)
      return # all set already
    # Maximum dimension of box
    max_dim = max(lmn)
    if not is_even(max_dim):
      max_dim += 1

    # How many grid points to add in each direction
    dlmn = [max_dim - a for a in lmn]
    # How many to add before
    dlmn1= [a//2 for a in dlmn]
    # How many to add after
    dlmn2 = [a - b for a,b in zip( dlmn, dlmn1)]
    # New start, end
    new_first = [b - a for a, b in zip(dlmn1,self.gridding_first)]
    new_last = [b + a for a, b in zip(dlmn2,self.gridding_last)]
    new_lmn = [1 +b - a for a,b in zip(new_first, new_last)]

    print("Original map size: ",map_all, file = log)
    print("Original box: ",lmn, "cubic:", max_dim, file = log)
    print("Original start: ",self.gridding_first, file = log)
    print("Original end: ",self.gridding_last, file = log)
    print("New box: ",new_lmn, file = log)
    print("New start: ",new_first, file = log)
    print("New end: ",new_last, file = log)

    # Now make sure we are inside map if requested
    lowest_value,highest_value = cube_relative_to_box( # 0 or neg, 0 or pos
      new_first, new_last, map_all,
       require_even_gridding = require_even_gridding)
    if stay_inside_current_map and (lowest_value != 0 or highest_value != 0):
      print("Reboxing cubic map to stay inside current map", file = log)
      new_first = [a - lowest_value for a in new_first]
      new_last = [a - highest_value for a in new_last]
      lowest_value,highest_value = cube_relative_to_box(
        new_first, new_last, map_all,
        require_even_gridding = require_even_gridding)
      assert [lowest_value,highest_value] == [0,0]
    print("Final start: ",new_first, file = log)
    print("Final end: ",new_last, file = log)
    print("Final box: ",[1 + b - a for a,b in zip(new_first, new_last)],
       file = log)
    self.gridding_first = new_first
    self.gridding_last = new_last

  def set_shifts_and_crystal_symmetry(self):
    '''
      Set items needed to do the shift

      Here self.gridding_first and self.gridding_last are the grid points
      marking the start and end, in the map with origin at (0, 0, 0), of the
      region to be kept.

      Also save the input crystal_symmetry for comparison later
    '''

    # Save to check later if another map_manager is used as input
    self.accessor_at_initialization = self._map_manager.map_data().accessor()
    self.map_crystal_symmetry_at_initialization = self._map_manager.crystal_symmetry()

    full_cs = self._map_manager.unit_cell_crystal_symmetry()
    full_uc = full_cs.unit_cell()
    self.box_all = [j-i+1 for i, j in zip(self.gridding_first, self.gridding_last)]
    assert min(self.box_all) >= 1 # box size must be greater than zero. Check stay_inside_current_map and bounds
    # get shift vector as result of boxing
    full_all_orig = self._map_manager.unit_cell_grid
    self.shift_frac = \
        [-self.gridding_first[i]/full_all_orig[i] for i in range(3)]
    self.shift_cart = full_cs.unit_cell().orthogonalize(self.shift_frac)
    # get crystal symmetry of the box
    p = full_uc.parameters()
    abc = [p[i] * self.box_all[i]/full_all_orig[i] for i in range(3)]
    box_uc = uctbx.unit_cell(parameters = (abc[0], abc[1], abc[2], p[3], p[4], p[5]))
    self.crystal_symmetry = crystal.symmetry(
      unit_cell = box_uc, space_group = "P1")

    self._warning_message = ""

  def warning_message(self):
    return self._warning_message

  def apply_to_model_ncs_and_map(self):
    '''
    Apply boxing to to self._model,  self._map_manager
    so all are boxed

    '''

    if self._model:
      self._model = self.apply_to_model(self._model)

    self._map_manager = self.apply_to_map(self._map_manager)


  def apply_to_map(self, map_manager):
    '''
     Apply boxing to a map_manager that is similar to the one used to generate
       this around_model object

     Also apply to its ncs_object, if any

    '''
    assert isinstance(map_manager, iotbx.map_manager.map_manager)

    # This one should just have similar unit_cell_crystal_symmetry
    # crystal_symmetry should match self.map_crystal_symmetry_at_initialization

    assert map_manager.unit_cell_crystal_symmetry().is_similar_symmetry(
      self._map_manager.unit_cell_crystal_symmetry())
    assert map_manager.crystal_symmetry().is_similar_symmetry(
      self.map_crystal_symmetry_at_initialization)

    ma1 = map_manager.map_data().accessor()
    ma2 = self.accessor_at_initialization
    assert ma1.all()     ==  ma2.all()
    assert ma1.origin()  ==  ma2.origin()
    assert ma1.focus()   ==  ma2.focus()
    map_data = map_manager.map_data()
    # Check if map is all valid
    bounds_info = get_bounds_of_valid_region(map_data,
      self.gridding_first,
      self.gridding_last)
    # Allow override of wrapping
    if isinstance(self._force_wrapping, bool):
      wrapping = self._force_wrapping
    else:
      # Get wrapping from map_manager. If it is not defined and
      #  bounds are outside allowed, try to get the wrapping
      wrapping = map_manager.wrapping()

    if wrapping or bounds_info.inside_allowed_bounds:
      # Just copy everything
      map_box = maptbx.copy(map_data, self.gridding_first, self.gridding_last)
      # Note: map_box gridding is self.gridding_first to self.gridding_last
    elif not bounds_info.some_valid_points:
      # No valid points, Just copy everything and zero
      map_box = maptbx.copy(map_data, self.gridding_first, self.gridding_last)
      map_box = map_box * 0.
      self._warning_message += "\nWARNING: boxed map is entirely outside map"+\
         " and wrapping=%s\n...setting all values to zero" %(wrapping)

    else: # Need to copy and then zero outside of defined region
      map_box = copy_and_zero_map_outside_bounds(map_data, bounds_info)
      self._warning_message += \
            "\nWARNING: boxed map goes outside original map"+\
         " and wrapping=%s\n...setting unknown values to zero" %(wrapping)
    #  Now reshape map_box to put origin at (0, 0, 0)
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
    #                                 (0, 0, 0) grid point in map_box)
    #     labels  (add label specifying boxing operation)
    #
    # New origin_shift_grid_units:
    origin_shift_grid_units = [
      self.gridding_first[i]+map_manager.origin_shift_grid_units[i]
        for i in range(3)]
    # New labels:
    new_label = "Boxed %s to %s %s" %(
      str(tuple(self.gridding_first)), str(tuple(self.gridding_last)),
      self.basis_for_boxing_string)
    #  Set up new map_manager. This will contain new data and not overwrite
    #   original
    #  NOTE: origin_shift_grid_units is required as bounds have changed

    # Crystal symmetry is now always P1 and wrapping is False
    new_map_manager = map_manager.customized_copy(map_data = map_box,
      origin_shift_grid_units = origin_shift_grid_units,
      crystal_symmetry_space_group_number = 1,
      wrapping = False)
    if self._force_wrapping:
      # Set the wrapping of the new map if it is possible
      if (self._force_wrapping and (new_map_manager.is_full_size())) or \
       ( (not self._force_wrapping) and (not new_map_manager.is_full_size())):
        new_map_manager.set_wrapping(self._force_wrapping)

    # Add the label
    new_map_manager.add_label(new_label)
    return new_map_manager

  def apply_to_model(self, model):
    '''
       Apply boxing to a model that is similar to the one used to generate
       this around_model object

       Changes the model in place

       Allow relaxed check if require_match_unit_cell_crystal_symmetry=False
    '''

    assert isinstance(model, mmtbx.model.manager)

    # This one should have similar unit_cell_crystal_symmetry for map and
    #  model and model original_crystal_symmetry should match
    #   self.map_crystal_symmetry_at_initialization

    # Allow relaxed check if require_match_unit_cell_crystal_symmetry=False

    s = getattr(self,'require_match_unit_cell_crystal_symmetry', None)
    require_uc_crystal_symmetry = (s in (None, True))
    if require_uc_crystal_symmetry:
      if model.shift_cart() is None:
        # model not yet initialized for shifts
        assert self.map_manager().unit_cell_crystal_symmetry(
          ).is_similar_symmetry( model.crystal_symmetry())
      else:  # model is initialized: should match unless not requiring it
        assert model.unit_cell_crystal_symmetry()
        assert self.map_manager().unit_cell_crystal_symmetry(
           ).is_similar_symmetry( model.unit_cell_crystal_symmetry())

    # Shift the model and add self.shift_cart on to whatever shift was there
    model.shift_model_and_set_crystal_symmetry(
       shift_cart = self.shift_cart, # shift to apply
       crystal_symmetry = self.crystal_symmetry, # new crystal_symmetry
       )

    # if wrapping is False, check to see if model is outside the box
    if (not self.map_manager().wrapping()) and (
        not self.model_can_be_outside_bounds):
      if not model.is_inside_working_cell():
        self._warning_message += "\nWARNING: Model is not entirely "+\
          "inside working cell and wrapping is False"
    return model

  def apply_to_ncs_object(self, ncs_object):
    '''
      Apply shifts from this boxing to an ncs_object
       ncs does keep track of shifts
    '''

    return ncs_object.coordinate_offset(coordinate_offset = self.shift_cart)


class around_model(with_bounds):
  """
  Extract map box around atomic model. Box is in P1 and has origin at (0, 0, 0).

  Creates new map_manager and modifies model in place

  Input map_manager must have origin (working position) at (0, 0, 0)
  Input model coordinates must correspond to working position of map_manager

  On initialization new bounds and crystal_symmetry are identified.
  Then map_manager is replaced with boxed version and shifted model is created

  Output versions of map_manager and model are in P1 and have origin
    at (0, 0, 0).
  Bounds refer to grid position in this box with origin at (0, 0, 0)

  Wrapping:
    wrapping = True means that grid points outside of the unit cell can be
    mapped inside with unit translations and box can effectively come
    from anywhere in space.
    If wrapping = True, the supplied box must be a complete unit cell so that
    map_data.all() must be the same as unit_cell_grid.

   wrapping = False means that grid points outside of the unit cell are
    undefined.  If a box is specified that uses points outside the defined
    region, those points are set to zero.

  Bounds:
    if model_can_be_outside_bounds, allow model to be outside the bounds
    if stay_inside_current_map, adjust bounds to not go outside current map
      in the case that bounds are entirely outside current map, use current map
    Note: stay_inside_current_map applies to all boxing in this method
      because it is a normal part of the boxing process (in other boxing it
        only applies to cubic boxing which can cause a box to go outside the
        current map)
    if use_cubic_boxing, adjust bounds to make a cubic box by making box bigger.
      if this conflicts with stay_inside_current_map, make box smaller
  If require_match_unit_cell_crystal_symmetry:  require that unit_cell
    crystal symmetry of map and model match
  """
  def __init__(self, map_manager, model, box_cushion,
      wrapping = None,
      model_can_be_outside_bounds = False,
      stay_inside_current_map = None, # Note that this default is different
      use_cubic_boxing = False,
      require_match_unit_cell_crystal_symmetry = False,
      log = sys.stdout):

    self._map_manager = map_manager
    self._model = model
    self.model_can_be_outside_bounds = model_can_be_outside_bounds
    self.use_cubic_boxing = use_cubic_boxing
    self.require_match_unit_cell_crystal_symmetry = \
       require_match_unit_cell_crystal_symmetry

    s = getattr(self,'require_match_unit_cell_crystal_symmetry', None)
    self._force_wrapping = wrapping
    if wrapping is None:
      wrapping = self.map_manager().wrapping()
    self.basis_for_boxing_string = 'using_model, wrapping = %s' %(
      wrapping)

    # safeguards
    assert isinstance(map_manager, iotbx.map_manager.map_manager)
    assert isinstance(model, mmtbx.model.manager)
    assert self._map_manager.map_data().accessor().origin()  ==  (0, 0, 0)

    # Do not work with dummy map_manager
    assert not map_manager.is_dummy_map_manager()

    # Make sure working model and map_manager crystal_symmetry match
    assert map_manager.is_compatible_model(model,
      require_match_unit_cell_crystal_symmetry =
         self.require_match_unit_cell_crystal_symmetry)

    assert box_cushion >=  0

    if self.map_manager().wrapping():  # map must be entire unit cell
      assert map_manager.unit_cell_grid == map_manager.map_data().all()

    # NOTE: We are going to use crystal_symmetry and sites_frac based on
    #   the map_manager (the model could still have different crystal_symmetry)

    info = get_bounds_around_model(
      map_manager = map_manager,
      model = model,
      box_cushion = box_cushion,
      stay_inside_current_map = stay_inside_current_map)
    from scitbx.matrix import col
    if flex.double(col(info.upper_bounds) - col(info.lower_bounds)
        ).min_max_mean().min < -0.5:  # nothing there
      raise AssertionError("Sorry, model is entirely outside box,"+
        " so boxing around model staying inside current map is not possible")
      self.gridding_first = map_manager.data().origin()
      self.gridding_last = map_manager.data().all()
    else: # usual
      self.gridding_first = info.lower_bounds
      self.gridding_last = info.upper_bounds

    # Adjust gridding to make a cubic box if desired
    if self.use_cubic_boxing:
      self.set_cubic_boxing(stay_inside_current_map = stay_inside_current_map)

    # Ready with gridding...set up shifts and box crystal_symmetry
    self.set_shifts_and_crystal_symmetry()

    # Apply boxing to model, ncs, and map (if available)
    self.apply_to_model_ncs_and_map()

class around_unique(with_bounds):

  '''
  Identify unique part of density in a map (using ncs object if present)
  and create a new map_manager containing this box of density, masked
  around regions containing density.  Note: the map may be masked between
  nearby density regions so this map could have many discontinuities.

  NOTE: This method carries out both boxing and masking. Its effect is
  similar to create_box_with_bounds, where the bounds are defined by the
  asymmetric part of the map, followed by masking around that asymmetric part
  of the map.

  NOTE: ncs_object from map_manager will be used to identify the unique
  part of the map.

  Creates new map_manager and modifies model in place

  Input map_manager must have origin (working position) at (0, 0, 0)
  Input model coordinates must correspond to working position of map_manager

  On initialization new bounds and crystal_symmetry are identified.
  Then map_manager is replaced with boxed version and shifted model is created

  Output versions of map_manager and model are in P1 and have origin at (0, 0, 0).
  Bounds refer to grid position in this box with origin at (0, 0, 0)

  Wrapping:
    wrapping = True means that grid points outside of the unit cell can be
    mapped inside with unit translations and box can effectively come
    from anywhere in space.
    If wrapping = True, the supplied box must be a complete unit cell so that
    map_data.all() must be the same as unit_cell_grid.

   wrapping = False means that grid points outside of the unit cell are
    undefined.  If a box is specified that uses points outside the defined
    region, those points are set to zero.

    if use_cubic_boxing, adjust bounds to make a cubic box by making box bigger.
      if this conflicts with stay_inside_current_map, make box smaller
    Note: stay_inside_current_map applies only when using cubic boxing

      Additional parameters:
         mask_expand_ratio:   allows increasing masking radius beyond default at
                              final stage of masking
         solvent_content:  fraction of cell not occupied by macromolecule. May
                           be None in which case it is estimated from the map
         sequence:        one-letter code of sequence of unique part of molecule
         chain_type:       PROTEIN or RNA or DNA. Used with sequence to estimate
                            molecular_mass
         molecular_mass:    Molecular mass (Da) of entire molecule used to
                            estimate solvent_content
         symmetry:         Alternative way to specify symmetry (as a character
                            string like D2, T, C3)
         use_symmetry_in_extract_unique:   Use symmetry in identification
                            of unique part of map
         target_ncs_au_model: model marking center of location to choose as
                              unique
         box_cushion:        buffer around unique region to be boxed
         soft_mask:  use soft mask
         keep_low_density:  keep low density regions
         regions_to_keep:   Allows choosing just highest-density contiguous
                            region (regions_to_keep=1) or a few
         keep_this_region_only: Allows choosing any specific region (first region is 0 not 1)
         residues_per_region: Allows setting threshold to try and get about this many
                              residues in each region. Default is 50.
        If require_match_unit_cell_crystal_symmetry:  require that unit_cell
           crystal symmetry of map and model match

  '''

  def __init__(self, map_manager,
    model = None,
    target_ncs_au_model = None,
    regions_to_keep = None,
    residues_per_region = None,
    keep_this_region_only = None,
    solvent_content = None,
    resolution = None,
    sequence = None,
    molecular_mass = None,
    symmetry = None,
    use_symmetry_in_extract_unique = True,
    chain_type = 'PROTEIN',
    keep_low_density = True,  # default from map_box
    box_cushion= 5,
    soft_mask = True,
    mask_expand_ratio = 1,
    wrapping = None,
    use_cubic_boxing = False,
    stay_inside_current_map = True,
    require_match_unit_cell_crystal_symmetry = False,
    log = None):


    self.require_match_unit_cell_crystal_symmetry = \
       require_match_unit_cell_crystal_symmetry

    self.model_can_be_outside_bounds = None  # not used but required to be set
    self.use_cubic_boxing = use_cubic_boxing
    self._map_manager = map_manager
    self._model = model

    self._mask_data = None

    self._force_wrapping = wrapping
    if wrapping is None:
      wrapping = self.map_manager().wrapping()
    self.basis_for_boxing_string = 'around_unique, wrapping = %s' %(
      wrapping)

    if log is None:
      log = null_out() # Print only if a log is supplied

    assert isinstance(map_manager, iotbx.map_manager.map_manager)
    assert self._map_manager.map_data().accessor().origin()  ==  (0, 0, 0)
    assert resolution is not None
    if model is not None:
      assert isinstance(model, mmtbx.model.manager)
      assert map_manager.is_compatible_model(model,
        require_match_unit_cell_crystal_symmetry=True)
    if self.map_manager().wrapping():  # map must be entire unit cell
      assert map_manager.unit_cell_grid == map_manager.map_data().all()

    # Get crystal_symmetry
    crystal_symmetry = map_manager.crystal_symmetry()
    # Convert to map_data

    from cctbx.maptbx.segment_and_split_map import run as segment_and_split_map
    assert self._map_manager.map_data().origin() == (0, 0, 0)

    args = []
    if residues_per_region:
      args.append("residues_per_region=%s" %(residues_per_region))

    if keep_this_region_only is not None:
      regions_to_keep = -1 * keep_this_region_only


    if solvent_content is None and sequence is None and molecular_mass is None:
      from cctbx.maptbx.segment_and_split_map \
          import get_iterated_solvent_fraction
      solvent_content = get_iterated_solvent_fraction(
          crystal_symmetry = crystal_symmetry,
          mask_resolution = resolution,
          map = self._map_manager.map_data(),
          out = log)


    ncs_group_obj, remainder_ncs_group_obj, tracking_data  = \
      segment_and_split_map(args,
        map_data = self._map_manager.map_data(),
        crystal_symmetry = crystal_symmetry,
        ncs_obj = self._map_manager.ncs_object() if \
          use_symmetry_in_extract_unique else None,
        target_model = target_ncs_au_model,
        write_files = False,
        auto_sharpen = False,
        add_neighbors = False,
        density_select = False,
        save_box_map_ncs_au = True,
        resolution = resolution,
        solvent_content = solvent_content,
        chain_type = chain_type,
        sequence = sequence,
        molecular_mass = molecular_mass,
        symmetry = symmetry if use_symmetry_in_extract_unique else None,
        keep_low_density = keep_low_density,
        regions_to_keep = regions_to_keep,
        box_buffer = box_cushion,
        soft_mask_extract_unique = soft_mask,
        mask_expand_ratio = mask_expand_ratio,
        out = log)

    from scitbx.matrix import col

    # Note number of selected regions used.
    if hasattr(tracking_data, 'available_selected_regions'):
      self.set_info(group_args(
        group_args_type = 'available selected regions from around_unique',
        available_selected_regions = tracking_data.available_selected_regions,
        ))

    if not hasattr(tracking_data, 'box_mask_ncs_au_map_data'):
      raise Sorry(" Extraction of unique part of map failed...")

    ncs_au_mask_data = tracking_data.box_mask_ncs_au_map_data

    lower_bounds = ncs_au_mask_data.origin()
    upper_bounds = tuple(
      col(ncs_au_mask_data.focus())-col((1, 1, 1)))

    print("\nBounds for unique part of map: %s to %s " %(
     str(lower_bounds), str(upper_bounds)), file = log)

    # shift the map so it is in the same position as the box map will be in
    ncs_au_mask_data.reshape(flex.grid(ncs_au_mask_data.all()))
    assert col(ncs_au_mask_data.all()) == \
        col(upper_bounds)-col(lower_bounds)+col((1, 1, 1))

    self.gridding_first = lower_bounds
    self.gridding_last  = upper_bounds

    # Adjust gridding to make a cubic box if desired
    if self.use_cubic_boxing:
      self.set_cubic_boxing(stay_inside_current_map = stay_inside_current_map)

    # Ready with gridding...set up shifts and box crystal_symmetry
    self.set_shifts_and_crystal_symmetry()

    # Apply boxing to model, ncs, and map (if available)
    self.apply_to_model_ncs_and_map()

    # Note that at this point, self._map_manager has been boxed
    assert ncs_au_mask_data.all() == self._map_manager.map_data().all()
    self._mask_data = ncs_au_mask_data

    # Now separately apply the mask to the boxed map
    self.apply_around_unique_mask(
       self._map_manager,
       resolution = resolution,
       soft_mask = soft_mask)

  def apply_around_unique_mask(self,
      map_manager,
      resolution,
      soft_mask):
    '''
      This procedure matches what is done in segment_and_split_map
      It comes at the end of around_unique and can be applied to additional
      map_manager objects if desired.
    '''
    assert self._mask_data is not None

    map_manager.create_mask_with_map_data(map_data = self._mask_data)

    if soft_mask: # Make the mask a soft mask if requested
      map_manager.soft_mask(soft_mask_radius = resolution)
      map_manager.apply_mask()
      # Now mask around edges
      map_manager.create_mask_around_edges(boundary_radius = resolution)
      map_manager.soft_mask(soft_mask_radius = resolution)
      map_manager.apply_mask()

    else:  # just apply the mask
      map_manager.apply_mask()

    # And add limitation to map
    map_manager.add_limitation("extract_unique")

class around_mask(with_bounds):
  """
  Extract map box around masked region of a map that represents a mask
  You need to supply the mask as a map_manager object

  Box is in P1 and has origin at (0, 0, 0).

  Input map_manager must have origin (working position) at (0, 0, 0)

  Returns boxed version of map_manager supplied.  Object will contain
  the boxed version of mask as self.mask_as_map_manager.

  Wrapping:
    wrapping = True means that grid points outside of the unit cell can be
    mapped inside with unit translations and box can effectively come
    from anywhere in space.
    If wrapping = True, the supplied box must be a complete unit cell so that
    map_data.all() must be the same as unit_cell_grid.

   wrapping = False means that grid points outside of the unit cell are
    undefined.  If a box is specified that uses points outside the defined
    region, those points are set to zero.

  Note: stay_inside_current_map applies only when using cubic boxing
  if use_cubic_boxing, adjust bounds to make a cubic box by making box bigger.
      if this conflicts with stay_inside_current_map, make box smaller
  If require_match_unit_cell_crystal_symmetry:  require that unit_cell
    crystal symmetry of map and model match

  """
  def __init__(self, map_manager,
     mask_as_map_manager,
     model = None,
     box_cushion = 3,
     wrapping = None,
     use_cubic_boxing = False,
     model_can_be_outside_bounds = False,
     stay_inside_current_map = True,
     require_match_unit_cell_crystal_symmetry = False,
     log = sys.stdout):

    self.require_match_unit_cell_crystal_symmetry = \
       require_match_unit_cell_crystal_symmetry
    self._map_manager = map_manager
    self._model = model
    self.model_can_be_outside_bounds = model_can_be_outside_bounds
    self.use_cubic_boxing = use_cubic_boxing
    assert map_manager.shift_cart()==mask_as_map_manager.shift_cart()

    # safeguards
    assert isinstance(map_manager, iotbx.map_manager.map_manager)
    assert isinstance(mask_as_map_manager, iotbx.map_manager.map_manager)
    assert self._map_manager.map_data().accessor().origin()  ==  (0, 0, 0)
    assert map_manager.is_similar(mask_as_map_manager)
    if self.map_manager().wrapping():
      assert map_manager.unit_cell_grid == map_manager.map_data().all()

    self._force_wrapping = wrapping
    if wrapping is None:
      wrapping = self.map_manager().wrapping()
    self.basis_for_boxing_string = 'around_mask bounds, wrapping = %s' %(
      wrapping)

    # Make sure the map goes from 0 to 1
    map_data = mask_as_map_manager.map_data()
    mmm = map_data.as_1d().min_max_mean()
    minimum = mmm.min
    range_of_values = mmm.max - mmm.min
    map_data = (map_data - minimum ) / max(1.e-10,range_of_values)


    # Get a connectivity object that marks all the connected regions in map

    from cctbx.maptbx.segment_and_split_map import get_co
    co, sorted_by_volume, min_b, max_b = get_co(
       map_data = map_data,
       threshold = 0.5,
       wrapping = False)


    if len(sorted_by_volume)<2:  # didn't work
      raise Sorry("No mask obtained...")

    # Get the biggest connected region in the map

    original_id_from_id = {}
    for i in range(1, len(sorted_by_volume)):
      v, id = sorted_by_volume[i]
      original_id_from_id[i] = id
    id = 1
    orig_id = original_id_from_id[id]

    # Get lower and upper bounds of this region in grid units

    self.gridding_first = min_b[orig_id]
    self.gridding_last  = max_b[orig_id]

    # Increase range of bounds by box_cushion
    cs = map_manager.crystal_symmetry()
    cushion = flex.double(cs.unit_cell().fractionalize((box_cushion, )*3))
    all_orig = map_manager.map_data().all()
    self.gridding_first = [max(0, ifloor(gf-c*n)) for c, gf, n in zip(
       cushion, self.gridding_first, all_orig)]
    self.gridding_last  = [min(n-1, iceil(gl+c*n)) for c, gl, n in zip(
       cushion, self.gridding_last, all_orig)]

    # Adjust gridding to make a cubic box if desired
    if self.use_cubic_boxing:
      self.set_cubic_boxing(stay_inside_current_map = stay_inside_current_map)

    # Ready with gridding...set up shifts and box crystal_symmetry
    self.set_shifts_and_crystal_symmetry()

    self.apply_to_model_ncs_and_map()

    # Also apply to mask_as_map_manager so that mask_as_map_manager is boxed
    mask_as_map_manager = self.apply_to_map(mask_as_map_manager)
    self.mask_as_map_manager = mask_as_map_manager # save it


class around_density(with_bounds):
  """
  Extract map box around region of a map containing density

  Box is in P1 and has origin at (0, 0, 0).

  Input map_manager must have origin (working position) at (0, 0, 0)

  Wrapping:
    wrapping = True means that grid points outside of the unit cell can be
    mapped inside with unit translations and box can effectively come
    from anywhere in space.
    If wrapping = True, the supplied box must be a complete unit cell so that
    map_data.all() must be the same as unit_cell_grid.

   wrapping = False means that grid points outside of the unit cell are
    undefined.  If a box is specified that uses points outside the defined
    region, those points are set to zero.

   if use_cubic_boxing, adjust bounds to make a cubic box by making box bigger.
   Note: stay_inside_current_map applies only when using cubic boxing
   if this conflicts with stay_inside_current_map, make box smaller
  If require_match_unit_cell_crystal_symmetry:  require that unit_cell
    crystal symmetry of map and model match

  """
  def __init__(self, map_manager,
     threshold = 0.05,
     box_cushion = 3.,
     get_half_height_width = True,
     model = None,
     wrapping = None,
     model_can_be_outside_bounds = False,
     use_cubic_boxing = False,
     stay_inside_current_map = True,
     require_match_unit_cell_crystal_symmetry = False,
     log = sys.stdout):

    self.require_match_unit_cell_crystal_symmetry = \
       require_match_unit_cell_crystal_symmetry
    self._map_manager = map_manager
    self._model = model
    self.model_can_be_outside_bounds = model_can_be_outside_bounds
    self.use_cubic_boxing = use_cubic_boxing

    # safeguards
    assert threshold is not None
    assert box_cushion is not None
    assert isinstance(map_manager, iotbx.map_manager.map_manager)
    assert self._map_manager.map_data().accessor().origin()  ==  (0, 0, 0)
    if self.map_manager().wrapping():
      assert map_manager.unit_cell_grid == map_manager.map_data().all()

    self._force_wrapping = wrapping
    if wrapping is None:
      wrapping = self.map_manager().wrapping()
    self.basis_for_boxing_string = 'around_density, wrapping = %s' %(
      wrapping)

    # Select box where data are positive (> threshold*max)
    map_data = map_manager.map_data()
    origin = list(map_data.origin())
    assert origin == [0, 0, 0]
    all = list(map_data.all())

    edge_values = flex.double()
    ux = all[0]-1
    uy = all[1]-1
    uz = all[2]-1
    for lb, ub in (
      [(0,0,0), (0,uy,uz)],
      [(ux,0,0), (ux,uy,uz)],
      [(0,0,0), (ux,0,uz)],
      [(0,uy,0), (ux,uy,uz)],
      [(0,0,0), (ux,uy,0)],
      [(0,0,uz), (ux,uy,uz)],
      ):
      new_map_data = maptbx.copy(map_data,lb,ub)
      edge_values.append(new_map_data.as_1d().as_double().min_max_mean().max)
    edge_value = edge_values.min_max_mean().mean

    # Get max value vs x, y, z
    value_list = flex.double()
    for i in range(0, all[0]):
      new_map_data = maptbx.copy(map_data,
         tuple((i, 0, 0)),
         tuple((i, all[1], all[2]))
       )
      value_list.append(
        new_map_data.as_1d().as_double().min_max_mean().max - edge_value)
    ii = 0
    for z in value_list:
      ii+= 1
    x_min, x_max = get_range(value_list, threshold = threshold,
      get_half_height_width = get_half_height_width)

    value_list = flex.double()
    for j in range(0, all[1]):
      new_map_data = maptbx.copy(map_data,
         tuple((0, j, 0)),
         tuple((all[0], j, all[2]))
       )
      value_list.append(
        new_map_data.as_1d().as_double().min_max_mean().max - edge_value)
    y_min, y_max = get_range(value_list, threshold = threshold,
      get_half_height_width = get_half_height_width)

    value_list = flex.double()
    for k in range(0, all[2]):
      new_map_data = maptbx.copy(map_data,
         tuple((0, 0, k)),
         tuple((all[0], all[1], k))
       )
      value_list.append(
        new_map_data.as_1d().as_double().min_max_mean().max - edge_value)
    z_min, z_max = get_range(value_list, threshold = threshold,
      get_half_height_width = get_half_height_width)

    # Get lower and upper bounds of this region in grid units
    frac_min = (x_min, y_min, z_min)
    frac_max = (x_max, y_max, z_max)
    cs = map_manager.crystal_symmetry()
    cushion = flex.double(cs.unit_cell().fractionalize((box_cushion, )*3))
    all_orig = map_data.all()
    self.gridding_first = [max(0, ifloor((f-c)*n)) for c, f, n in zip(
       cushion, frac_min, all_orig)]
    self.gridding_last  = [ min(n-1, iceil((f+c)*n)) for c, f, n in zip(
       cushion, frac_max, all_orig)]
    # Adjust gridding to make a cubic box if desired
    if self.use_cubic_boxing:
      self.set_cubic_boxing(stay_inside_current_map = stay_inside_current_map)

    # Ready with gridding...set up shifts and box crystal_symmetry
    self.set_shifts_and_crystal_symmetry()

    # Apply boxing to model, ncs, and map (if available)
    self.apply_to_model_ncs_and_map()

def is_even(n):
  if n < 0:
    n = -n
  if 2 * (n//2) == n:
    return True
  else:
    return False

def cube_relative_to_box(new_first, new_last, map_all,
       require_even_gridding = None):
    ''' returns zero or negative number for lowest_value of new_first and
       zero or positive number for highest value of new_last - map_all
      If require_even_gridding, make the lowest and highest even by making
       them further from zero'''

    lowest_value = min(0, min([a for a in new_first]))
    if require_even_gridding and (not is_even(lowest_value)):
      lowest_value -= 1
    highest_value = max(0, max([a - b for a, b in zip(new_last,map_all)]))
    if require_even_gridding and (not is_even(highest_value)):
      highest_value += 1
    print("lowest, highest out of box:",lowest_value, highest_value)
    return lowest_value,highest_value

def get_range(value_list, threshold = None, ignore_ends = True,
   keep_near_ends_frac = 0.02, half_height_width = 2.,
   get_half_height_width = None,
   cutoff_ratio = 4, ratio_max = 0.5,
   smooth_list_first = True,
   smooth_units = 20,
   max_allowed_outside_of_box = 0.33 ): # XXX May need to set cutoff_ratio and
  #  ratio_max lower.
  # ignore ends allows ignoring the first and last points which may be off
  # if get_half_height_width, find width at half max hieght, go
  #  half_height_width times this width out in either direction, use that as
  #  baseline instead of full cell. Don't do it if the height at this point
  #  is over cutoff_ratio times threshold above original baseline.
  #  If value outside range is more than max_allowed_outside_of_box of max,
  #    make it bigger.
  n_tot = value_list.size()
  assert n_tot>0
  if smooth_list_first:
   value_list = flex.double(smooth_list(value_list,
     smooth_range = max(1, value_list.size()//smooth_units)))

  if get_half_height_width:
    z_min, z_max = get_range(value_list, threshold = 0.5,
      ignore_ends = ignore_ends, keep_near_ends_frac = keep_near_ends_frac,
      get_half_height_width = False, smooth_list_first = False)
    z_mid = 0.5*(z_min+z_max)
    z_width = 0.5*(z_max-z_min)
    z_low = z_mid-2*z_width
    z_high = z_mid+2*z_width
    if ignore_ends:
      i_max = value_list.size()-2
      i_min = 1
    else:
      i_max = value_list.size()-1
      i_min = 0

    i_low =  max(i_min, min(i_max, int(0.5+z_low* value_list.size())))
    i_high = max(i_min, min(i_max, int(0.5+z_high*value_list.size())))
    min_value = value_list.min_max_mean().min
    max_value = value_list.min_max_mean().max
    ratio_low = (value_list[i_low]-min_value)/max(
       1.e-10, (max_value-min_value))
    ratio_high = (value_list[i_high]-min_value)/max(
       1.e-10, (max_value-min_value))
    if ratio_low <=  cutoff_ratio*threshold and ratio_low >0 \
         and ratio_low<ratio_max\
         and ratio_high <=  cutoff_ratio*threshold and ratio_high > 0 \
         and ratio_high < ratio_max:
      ratio = min(ratio_low, ratio_high)
      z_min, z_max = get_range(
        value_list, threshold = threshold+ratio,
        ignore_ends = ignore_ends, keep_near_ends_frac = keep_near_ends_frac,
        get_half_height_width = False, smooth_list_first = False)
    else:
      z_min, z_max = get_range(value_list, threshold = threshold,
        ignore_ends = ignore_ends, keep_near_ends_frac = keep_near_ends_frac,
        get_half_height_width = False, smooth_list_first = False)
    if not too_high_outside_range(value_list, int(0.5+ n_tot * z_min),
         int(0.5+ n_tot *z_max),
        max_allowed_outside_of_box): # ok
      return z_min, z_max

  if threshold is None: threshold = 0
  min_value = value_list.min_max_mean().min
  max_value = value_list.min_max_mean().max
  cutoff = min_value+(max_value-min_value)*threshold

  # Find lowest point to left and right of highest point
  vl = list(value_list)
  max_i = vl.index(max_value)
  if max_i == 0:
    min_i_to_left = 0
  else: # usual
    min_to_left = value_list[:max_i].min_max_mean().min
    min_i_to_left = vl.index(min_to_left,0,max_i)
  if max_i == n_tot - 1:
    min_i_to_right = n_tot - 1
  else: # usual
    min_to_right = value_list[max_i:].min_max_mean().min
    min_i_to_right = vl.index(min_to_right,min(n_tot-1,max_i+1),n_tot)
  if ignore_ends:
    i_off = 1
  else:
    i_off = 0
  i_low = None
  for i in range(max(min_i_to_left,i_off), min(min_i_to_right,n_tot-i_off)):
    if value_list[i]>cutoff:
      i_low = max(i_off, i-1)
      break
  i_high = None
  for ii in range(
       min(min_i_to_right,n_tot-i_off),
       max(min_i_to_left,i_off),
        -1):
    if value_list[ii]>cutoff:
      i_high = min(n_tot-1-i_off, ii+1)
      break
  if i_low is None or i_high is None:
    raise Sorry("Cannot auto-select region...")
  if i_low/n_tot<keep_near_ends_frac: i_low = 0
  if (n_tot-1-i_high)/n_tot<keep_near_ends_frac: i_high = n_tot-1
  if not too_high_outside_range(value_list, i_low, i_high,
        max_allowed_outside_of_box): # ok
    return i_low/n_tot, i_high/n_tot

  # Failed to include high density...try again not using low point
  if threshold is None: threshold = 0
  min_value = value_list.min_max_mean().min
  max_value = value_list.min_max_mean().max
  cutoff = min_value+(max_value-min_value)*threshold
  if ignore_ends:
    i_off = 1
  else:
    i_off = 0
  i_low = None
  for i in range(i_off, n_tot-i_off):
    if value_list[i]>cutoff:
      i_low = max(i_off, i-1)
      break
  i_high = None
  for i in range(i_off, n_tot-i_off):
    ii = n_tot-1-i
    if value_list[ii]>cutoff:
      i_high = min(n_tot-1-i_off, ii+1)
      break
  if i_low is None or i_high is None:
    raise Sorry("Cannot auto-select region...")
  if i_low/n_tot<keep_near_ends_frac: i_low = 0
  if (n_tot-1-i_high)/n_tot<keep_near_ends_frac: i_high = n_tot-1
  if not too_high_outside_range(value_list, i_low, i_high,
        max_allowed_outside_of_box): # ok
    return i_low/n_tot, i_high/n_tot
  else:  # give up and take the whole thing
    return i_off/n_tot, (n_tot - i_off - 1)/n_tot



def too_high_outside_range(value_list, i_start, i_end,
       max_allowed_outside_of_box):
  inside_values = flex.double(value_list[i_start:i_end])
  outside_values = flex.double(value_list[:i_start])
  outside_values.extend(flex.double(value_list[i_end:]))
  if outside_values.min_max_mean().max > \
      max_allowed_outside_of_box * inside_values.min_max_mean().max:
    return True
  else:
    return False

def smooth_list(working_list,smooth_range = None): # smooth this list of numbers
    assert smooth_range is not None
    new_list=[]
    delta=smooth_range//2
    for i in range(len(working_list)):
      sum=0.
      sumn=0.
      for j in range(-delta,delta+1):
        jj=j+i
        if jj >=0 and jj < len(working_list):
          sum+=working_list[jj]
          sumn+=1.
      if sumn>0.:
        sum=sum/sumn
      new_list.append(sum)
    return new_list

def get_bounds_of_valid_region(map_data,
    gridding_first,
    gridding_last):

  '''
    If map_data is sampled from gridding_first to gridding_last with
    maptbx.copy, (1) does the sampling go outside of map_data?
    (2) What are the lower and upper bounds of the valid region of the resulting
    map?
  '''

  lower_allowed_bounds = []
  upper_allowed_bounds = []
  inside_allowed_bounds = True
  some_valid_points = True
  for o, a, f, l in zip(map_data.origin(), map_data.all(),
     gridding_first, gridding_last):
    # Available map goes from o to o+a-1
    # Requested map goes from f to l

    # If f is less than o, first valid grid point is o.
    # Otherwise, first valid grid point is f
    # so first valid grid point is max(f, o)

    # If f is less than o
    #   After shifting origin to (0, 0, 0) first valid grid point is o-f
    # Otherwise, after shifting origin, first valid grid point is 0.

    #  If l is >=  a+o, last valid grid point is a+o-1
    #  Otherwise last valid grid point is l
    # So last valid grid point is min(l, a+o-1)

    #  If l is >=  a+o
    #    After shifting origin to (0, 0, 0), last valid grid point is a+o-1-f
    #  Otherwise, after shifting origin, last valid grid point is l-f


    first_valid = max(o, f) # first valid grid point
    last_valid = min(l, o+a-1) # last valid grid point
    lower_allowed_bounds.append(first_valid)
    upper_allowed_bounds.append(last_valid)
    if f < o or l >=  a+o:
      inside_allowed_bounds = False
    if last_valid <=  first_valid:
      some_valid_points = False

  return group_args(
     lower_allowed_bounds = lower_allowed_bounds,
     upper_allowed_bounds = upper_allowed_bounds,
     gridding_first = gridding_first,
     gridding_last = gridding_last,
     inside_allowed_bounds = inside_allowed_bounds,
     some_valid_points = some_valid_points)

def copy_and_zero_map_outside_bounds(map_data, bounds_info):
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

  map_copy = maptbx.copy(map_data,
      bounds_info.gridding_first, bounds_info.gridding_last)
  # Note: the origin of map_copy is at gridding_first and goes to last

  # Make sure we are working with a flex.double array
  if type(map_copy) !=  type(flex.double()): # must be double
    map_copy = map_copy.as_double()

  # Make sure this map matches the bounds_info
  assert tuple(map_copy.origin()) == tuple(bounds_info.gridding_first)

  # We are going to shift the origin in this copy, so shift the bounds to match

  lower_bounds_after_shift = []
  upper_bounds_after_shift = []
  for l, u, f in zip(bounds_info.lower_allowed_bounds,
      bounds_info.upper_allowed_bounds,
      bounds_info.gridding_first):
    lower_bounds_after_shift.append(l-f)
    upper_bounds_after_shift.append(u-f)

  acc = map_copy.accessor() # save where the origin is

  map_copy = map_copy.shift_origin()  # put origin at (0, 0, 0)
  map_copy_all = map_copy.all() # save size of map
  # XXX work-around for set_box does not allow offset origin
  map_copy.resize(flex.grid(map_copy_all))
  new_map = maptbx.set_box_copy_inside(0,  # copies inside, zero outside bounds
    map_data_to   = map_copy,
    start         = tuple(lower_bounds_after_shift),
    end           = tuple(upper_bounds_after_shift))
  # XXX and shift map back
  new_map = new_map.as_1d()
  new_map.reshape(acc)
  return new_map

def shift_and_box_model(model = None,
    box_cushion = 5, shift_model = True,
    crystal_symmetry = None):
  '''
    Shift a model near the origin and box around it
    Use crystal_symmetry if supplied
    Keeps input model unchanged.
  '''
  from mmtbx.model import manager as model_manager
  from scitbx.matrix import col
  from cctbx import crystal

  ph=model.get_hierarchy()
  sites_cart=ph.atoms().extract_xyz()
  if shift_model:
    sites_cart=sites_cart-col(sites_cart.min())+col(
      (box_cushion,box_cushion,box_cushion))

  box_start=col(sites_cart.min())-col((box_cushion,box_cushion,box_cushion))
  box_end=col(sites_cart.max())+col((box_cushion,box_cushion,box_cushion))
  if not crystal_symmetry:
    a,b,c = box_end - box_start
    crystal_symmetry=crystal.symmetry((a,b,c, 90,90,90),1)
  phc = ph.deep_copy()
  phc.atoms().set_xyz(sites_cart)
  phc.atoms().reset_serial()
  mm = model_manager(
     model_input = None,
     pdb_hierarchy = phc,
     crystal_symmetry = crystal_symmetry,
     restraint_objects = model.get_restraint_objects(),
     monomer_parameters = model.get_monomer_parameters(),
     log = null_out())
  return mm

def get_boxes_to_tile_map(target_for_boxes = 24,
      n_real = None,
      crystal_symmetry = None,
      cushion_nx_ny_nz = None,
      wrapping = False,
      do_not_go_over_target = None,
      target_xyz_center_list = None,
     ):

    '''
      Get a set of boxes that tile the map
      If cushion_nx_ny_nz is set ... create a second set of boxes that are
        expanded by cushion_nx_ny_nz in each direction
      Try to make boxes symmetrical in full map
      If target_xyz_center_list is set, try to use them as centers but keep
       size the same as would otherwise be used
    '''
    nx,ny,nz = n_real
    smallest = min(nx,ny,nz)
    largest = max(nx,ny,nz)
    target_volume_per_box = (nx*ny*nz)/target_for_boxes
    target_length = target_volume_per_box**0.33
    if target_xyz_center_list:
      lower_bounds_list = []
      upper_bounds_list = []
      uc = crystal_symmetry.unit_cell()
      for site_frac in uc.fractionalize(target_xyz_center_list):
        center_ijk = tuple([ int(0.5+x * n) for x,n in zip(site_frac, n_real)])
        lower_bounds_list.append(
          tuple( [
             int(max(1,min(n-2,(i - (1+target_length)//2)))) for i,n in
            zip(center_ijk,n_real)
             ]
          ))
        upper_bounds_list.append(
          tuple( [
             int(max(1,min(n-2,(i + (1+target_length)//2)))) for i,n in
            zip(center_ijk,n_real)
             ]
          ))
    elif target_for_boxes == 1:
      lower_bounds_list = [(0,0,0)]
      upper_bounds_list = [tuple([i - 1 for i in n_real])]
    else:
      lower_bounds_list = []
      upper_bounds_list = []
      for x_info in get_bounds_list(nx, target_length,
        do_not_go_over_target = do_not_go_over_target):
        for y_info in get_bounds_list(ny, target_length,
           do_not_go_over_target = do_not_go_over_target):
          for z_info in get_bounds_list(nz, target_length,
             do_not_go_over_target = do_not_go_over_target):
            lower_bounds_list.append(
               [x_info.lower_bound,
                y_info.lower_bound,
                z_info.lower_bound])
            upper_bounds_list.append(
               [x_info.upper_bound,
                y_info.upper_bound,
                z_info.upper_bound])

    # Now make a set of boxes with a cushion if requested
    lower_bounds_with_cushion_list = []
    upper_bounds_with_cushion_list = []
    if cushion_nx_ny_nz:

      for lb,ub in zip (lower_bounds_list,upper_bounds_list):
        if (wrapping):
          new_lb = tuple([b - c for b,c in zip(lb, cushion_nx_ny_nz)])
          new_ub = tuple([u + c for u,c in zip(ub, cushion_nx_ny_nz)])
        else:
          new_lb = tuple([max(0,b - c) for b,c in zip(lb, cushion_nx_ny_nz)])
          new_ub = tuple([min(n-1,u + c) for u,c,n in zip(ub, cushion_nx_ny_nz,
            n_real)])
        lower_bounds_with_cushion_list.append(new_lb)
        upper_bounds_with_cushion_list.append(new_ub)
    else:
      lower_bounds_with_cushion_list = lower_bounds_list
      upper_bounds_with_cushion_list = upper_bounds_list

    # Now remove any duplicates
    lb_ub_list = []
    new_lb_list = []
    new_ub_list = []
    for lb,ub in zip (
         lower_bounds_with_cushion_list,upper_bounds_with_cushion_list):
       if [lb,ub] in lb_ub_list: continue
       lb_ub_list.append([lb,ub])
       new_lb_list.append(lb)
       new_ub_list.append(ub)
    lower_bounds_with_cushion_list = new_lb_list
    upper_bounds_with_cushion_list = new_ub_list

    return group_args(
      lower_bounds_list = lower_bounds_list,
      upper_bounds_list = upper_bounds_list,
      lower_bounds_with_cushion_list = lower_bounds_with_cushion_list,
      upper_bounds_with_cushion_list = upper_bounds_with_cushion_list,
      n_real = n_real,
      crystal_symmetry = crystal_symmetry,
     )

def get_bounds_list(nx, target_length,
     do_not_go_over_target = None,):
  '''
    Return start, end that are about the length target_length and that
    collectively cover exactly nx grid units.
    Try to make bounds on the ends match target_length
    Try to make bounds symmetrical
  '''

  bounds_list = []
  if nx < (3 * target_length)/2:  # take one only
    bounds_list.append(
      group_args(
       lower_bound = 0,
       upper_bound = nx - 1,)
      )
  else:  # take as many as fit
    n_target = nx/target_length  # how many we want (float)
    if do_not_go_over_target:
      n = max(1,int(n_target)) # int ... how many can fit
    else: # usual
      n = max(1,int(0.5 + n_target)) # int ... how many can fit
    exact_target_length =  nx/n  # float length of each group

    last_end_point = -1
    length_list = []
    for i in range(n):
      target_end_point = exact_target_length * (i+1)
      actual_end_point = min (nx -1, max(0, int(0.5 + target_end_point)))
      length_list.append(actual_end_point - last_end_point)
      last_end_point = actual_end_point

    # Now try and make length_list symmetric
    length_list = make_list_symmetric(length_list)
    last_end_point = -1
    for i in range(n):
      actual_end_point = last_end_point + length_list[i]
      bounds_list.append(
          group_args(
       lower_bound = last_end_point + 1,
       upper_bound = actual_end_point,)
      )
      last_end_point = actual_end_point

  return bounds_list

def make_list_symmetric(length_list):
  '''
   adjust entries in length_list to make it symmetric but same total
  '''
  from copy import deepcopy
  length_list = deepcopy(length_list)
  unused_length = 0
  n=len(length_list)

  from scitbx.array_family import flex
  total = flex.double(tuple(length_list)).min_max_mean().mean*len(length_list)
  for i_from_end in range (n//2): # may leave out middle one if present
    i = i_from_end
    n_bigger = length_list[i] - length_list[n-i-1]
    n_bigger_abs = abs(n_bigger)
    n_shift = n_bigger_abs//2
    n_bigger_even = 2*n_shift
    if n_bigger > 0:
      # move n_shift to n-i-1 and save remainder
      length_list[n-i-1] += n_shift
      length_list[i] -= n_shift
      delta = length_list[i] - length_list[n-i-1]
      assert delta >= 0
      length_list[i] -= delta
      unused_length += delta
    elif n_bigger < 0:
      # move n_shift to i and save remainder
      length_list[i] += n_shift
      length_list[n-i-1] -= n_shift
      delta = length_list[n-i-1] - length_list[i]
      assert delta >= 0
      length_list[n-i-1] -= delta
      unused_length += delta
    if unused_length//2 > 0:
      length_list[i] += unused_length//2
      length_list[n-i-1] += unused_length//2
      unused_length -= 2* (unused_length//2)
  if unused_length:
    length_list[(n+1)//2] += unused_length
  return length_list


def get_bounds_around_model(
      map_manager = None,
      model = None,
      box_cushion = None,
      stay_inside_current_map = None,
     ):
    '''
      Calculate the lower and upper bounds to box around a model
      Allow bounds to go outside the available box unless
        stay_inside_current_map (this has to be dealt with at the boxing stage)
    '''

    # get items needed to do the shift
    cs = map_manager.crystal_symmetry()
    uc = cs.unit_cell()
    sites_cart = model.get_sites_cart()
    sites_frac = uc.fractionalize(sites_cart)
    map_data = map_manager.map_data()
    # convert box_cushion into fractional vector
    cushion_frac = flex.double(uc.fractionalize((box_cushion, )*3))
    # find fractional corners
    frac_min = sites_frac.min()
    frac_max = sites_frac.max()
    frac_max = list(flex.double(frac_max)+cushion_frac)
    frac_min = list(flex.double(frac_min)-cushion_frac)
    # find corner grid nodes
    all_orig = map_data.all()

    lower_bounds = [ifloor(f*n) for f, n in zip(frac_min, all_orig)]
    upper_bounds = [ iceil(f*n) for f, n in zip(frac_max, all_orig)]
    n = all_orig[-1]
    if stay_inside_current_map:
      lower_bounds = [ min(n-1,max (0,lb)) for lb in lower_bounds]
      upper_bounds = [ min (ub, n-1) for ub,n in zip(upper_bounds,all_orig)]
    return group_args(
      lower_bounds = lower_bounds,
      upper_bounds = upper_bounds,
    )
