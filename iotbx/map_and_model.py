from __future__ import absolute_import, division, print_function
from libtbx.utils import Sorry
from cctbx import maptbx
from libtbx import group_args
from scitbx.array_family import flex
from iotbx.map_model_manager import map_model_manager
from mmtbx.model import manager as model_manager
from libtbx.utils import null_out
from libtbx.test_utils import approx_equal

class input(object):

  '''
    Class for shifting origin of map(s) and model to (0, 0, 0) and keeping
    track of the shifts.

    Typical use:
    mam = map_and_model.input(
      model = model,
      map_manager = map_manager,
      ncs_object = ncs_object)

    mam.box_around_model(wrapping=False, box_cushion=3)

    shifted_model = mam.model()  # at (0, 0, 0), knows about shifts
    shifted_map_manager = mam.map_manager() # also at (0, 0, 0) knows shifts
    shifted_ncs_object = mam.ncs_object() # also at (0, 0, 0) and knows shifts

    NOTE: Expects symmetry of model and map_manager to match unless
      ignore_symmetry_conflicts = True

    Optional:  apply soft mask to map (requires resolution)
  '''
  def __init__(self,
               model            = None,
               map_manager      = None,  # replaces map_data
               map_manager_1    = None,  # replaces map_data_1
               map_manager_2    = None,  # replaces map_data_2
               map_manager_list = None,  # replaces map_data_list
               ncs_object       = None,
               ignore_symmetry_conflicts = False):

    self._model = model
    self._map_manager = map_manager
    self._map_manager_1 = map_manager_1
    self._map_manager_2 = map_manager_2
    self._map_manager_list = map_manager_list
    self._shift_manager = None
    self._ncs_object = ncs_object
    # CHECKS


    # Make sure that map_manager is either already shifted to (0, 0, 0) or has
    #   origin_shift_grid_unit of (0, 0, 0).
    if self._map_manager:
       assert self._map_manager.origin_is_zero() or \
          self._map_manager.origin_shift_grid_units == (0, 0, 0)

    # Normally map_manager unit_cell_crystal_symmetry should match
    #  model original_crystal_symmetry (and also usually model.crystal_symmetry)

    assert self._map_manager # always need one
    # ZZZ remove crystal_symmetry input

    # Check symmetry
    original_map_cs = self._map_manager.unit_cell_crystal_symmetry()
    working_map_cs = self._map_manager.crystal_symmetry()

    if self._model:
      ok = False
      original_model_cs = None
      if self._model.get_shift_manager() and\
          self._model.get_shift_manager().original_crystal_symmetry:
        original_model_cs = \
           self._model.get_shift_manager().original_crystal_symmetry
        if original_map_cs.is_similar_symmetry(original_model_cs): # fine
          ok = True # original model matches original map
      if (not ok) and original_map_cs.is_similar_symmetry(
           self._model.crystal_symmetry()):
        ok = True # model matches original map
      if (not ok) and working_map_cs.is_similar_symmetry(
           self._model.crystal_symmetry()):
        ok = True # model matches working map

      if (not ok) and ignore_symmetry_conflicts:
        self._model.set_crystal_symmetry(original_map_cs)
      elif (not ok): # stop
        raise Sorry("Map and model symmetry do not match"+
         "\nMap original symmetry: %s\n" %(str(original_map_cs))+
         "\nMap working symmetry: %s\n" %(str(working_map_cs))+
         "\nModel original symmetry: %s\n" %(str(original_model_cs))+
       "\nModel working symmetry: %s\n" %(str(self._model.crystal_symmetry())))

    # Make sure we have what is expected: optional model, mm,
    # self._map_manager_1 and self._map_manager_2 or neither,
    #   optional list of self._map_manager_list

    if not self._map_manager_list:
      self._map_manager_list = []

    if(not [self._map_manager_1, self._map_manager_2].count(None) in [0, 2]):
      raise Sorry("None or two half-maps are required.")
    if(not self._map_manager):
      raise Sorry("A map is required.")

    # Make sure all map_managers have same gridding and symmetry
    for m in [self._map_manager_1, self._map_manager_2]+ \
         self._map_manager_list:
      if m:
        assert self._map_manager.is_similar(m)

    # READY

    # If model, make a map_model_manager with model and mm and
    #  let it check symmetry against this one
    mmm = map_model_manager()
    mmm.add_map_manager(self._map_manager)
    if self._model:
      mmm.add_model(self._model, set_model_log_to_null = False) # keep the log
    if self._ncs_object:
      mmm.add_ncs_object(self._ncs_object)
    # ALL OK if it does not stop

    # Shift origin of model and map_manager to (0, 0, 0) with
    #    mmm which knows about both
    mmm.shift_origin(log = null_out())
    self._model = mmm.model()  # this model knows about shift so far
                             # NOTE: NO SHIFTS ALLOWED COMING IN
    self._map_manager = mmm.map_manager()  # map_manager also knows about shift
    self._ncs_object = mmm.ncs_object()  # ncs object also knows about shift
    self._crystal_symmetry = self._map_manager.crystal_symmetry()

    if self._model:
      self._shift_manager = self._model.get_shift_manager() # XXX save manager
      # Make sure model shift manager agrees with map_manager shift
      if self._shift_manager and self._map_manager:
        assert approx_equal(
          tuple([-a for a in self._shift_manager.shift_cart]),
          self._map_manager.origin_shift_cart())

    # Shift origins of all other maps
    for m in [self._map_manager_1, self._map_manager_2]+\
         self._map_manager_list:
      if m:
        m.shift_origin()

    # Make sure all really match:
    for m in [self._map_manager_1, self._map_manager_2]+\
        self._map_manager_list:
      if m:
        assert self._map_manager.is_similar(m)

    # Save origin after origin shift but before any boxing
    #    so they can be accessed easily later

    self._original_origin_grid_units = self._map_manager.origin_shift_grid_units
    self._original_origin_cart = self._map_manager.origin_shift_cart()

    #  Save gridding of this original map (after shifting, whole thing):
    self._gridding_first = (0, 0, 0)
    self._gridding_last = self._map_manager.map_data().all()

    # Holder for solvent content used in boxing and transferred to box_object
    self._solvent_content = None

  def box_around_model(self,
     wrapping = None,
     box_cushion = 5.):

    '''
       Box all maps around the model, shift origin of maps, model, ncs_object

       wrapping must be specified. Wrapping means map is infinite and repeats
       outside unit cell. Requires a full unit cell in the maps.
    '''
    assert isinstance(self._model, model_manager)
    assert isinstance(wrapping, bool) # must be decided by programmer
    assert box_cushion is not None

    from cctbx.maptbx.box import around_model
    if(self._map_manager_1 is not None):
      tmp_box = around_model(
        map_manager = self._map_manager_1,
        model = self._model.deep_copy(),
        cushion = box_cushion,
        wrapping = wrapping)
      self._map_manager_1 = tmp_box.map_manager()
      tmp_box = around_model(
        map_manager = self._map_manager_2,
        model = self._model.deep_copy(),
        cushion = box_cushion,
        wrapping = wrapping)
      self._map_manager_2 = tmp_box.map_manager()
    if self._map_manager_list:
      new_list = []
      for x in self._map_manager_list:
        tmp_box = around_model(
          map_manager = x,
          model = self._model.deep_copy(),
          cushion = box_cushion,
          wrapping = wrapping)
        new_list.append(tmp_box.map_manager())
      self._map_manager_list = new_list

    # Make box around model
    box = around_model(
      map_manager = self._map_manager,
      model = self._model,
      ncs_object = self._ncs_object,
      cushion = box_cushion,
      wrapping = wrapping)

    box_as_mam = box.as_map_and_model()

    # New map_manager and model know about cumulative shifts (original
    #   shift to move origin to (0, 0, 0) plus shift from boxing
    self._map_manager = box_as_mam.map_manager()
    self._model = box_as_mam.model()
    self._ncs_object = box_as_mam.ncs_object()

    self._shift_manager = self._model.get_shift_manager().deep_copy()

    # Update self._crystal_symmetry
    self._crystal_symmetry = self._model.crystal_symmetry()
    assert self._crystal_symmetry.is_similar_symmetry(
      self._map_manager.crystal_symmetry())


  def soft_mask_all_maps_around_edges(self,
      resolution = None,
      soft_mask_radius = None):

    # Apply a soft mask around edges of all maps. Overwrites values in maps

    for mm in self.all_map_managers():
      if not mm: continue
      mm.create_mask_around_edges(
        soft_mask_radius = soft_mask_radius)
      mm.apply_mask()

  def mask_all_maps_around_model(self,
      mask_atoms_atom_radius = None,
      set_outside_to_mean_inside = None,
      soft_mask = None,
      soft_mask_radius = None):
    assert mask_atoms_atom_radius is not None
    assert (not soft_mask) or (soft_mask_radius is not None)
    assert self.model() is not None

    # Apply a mask to all maps. Overwrites values in these maps

    for mm in self.all_map_managers():
      if not mm: continue
      mm.create_mask_around_atoms(
         model = self.model(),
         mask_atoms_atom_radius = mask_atoms_atom_radius)
      if soft_mask:
        mm.soft_mask(soft_mask_radius = soft_mask_radius)
      mm.apply_mask(
         set_outside_to_mean_inside = \
           set_outside_to_mean_inside)

  def original_origin_cart(self):
    assert self._original_origin_cart is not None
    return self._original_origin_cart

  def original_origin_grid_units(self):
    assert self._original_origin_grid_units is not None
    return self._original_origin_grid_units

  def shift_manager(self):
    return self._shift_manager

  def map_data(self):
    return self.map_manager().map_data()

  def map_data_1(self):
    if self.map_manager_1():
      return self.map_manager_1().map_data()

  def map_data_2(self):
    if self.map_manager_2():
      return self.map_manager_2().map_data()

  def all_map_managers(self):
    all_map_managers_list = []
    for x in [self.map_manager()]+[self.map_manager_1()]+\
        [self.map_manager_2()]+ self.map_manager_list():
      if x: all_map_managers_list.append(x)
    return all_map_managers_list

  def map_data_list(self):
    map_data_list = []
    for mm in self.map_manager_list():
      map_data_list.append(mm.map_data())
    return map_data_list

  def map_manager(self):
     return self._map_manager

  def map_manager_1(self):
     return self._map_manager_1

  def map_manager_2(self):
     return self._map_manager_2

  def map_manager_list(self):
     if self._map_manager_list:
       return self._map_manager_list
     else:
       return []

  def model(self): return self._model

  def ncs_object(self): return self._ncs_object

  def crystal_symmetry(self): return self._crystal_symmetry

  def xray_structure(self):
    if(self.model() is not None):
      return self.model().get_xray_structure()
    else:
      return None

  def hierarchy(self): return self._model.get_hierarchy()

  def set_gridding_first(self, gridding_first):
    self._gridding_first = tuple(gridding_first)

  def set_gridding_last(self, gridding_last):
    self._gridding_last = tuple(gridding_last)

  def set_solvent_content(self, solvent_content):
    self._solvent_content = solvent_content

  def get_counts_and_histograms(self):
    self._counts = get_map_counts(
      map_data         = self.map_data(),
      crystal_symmetry = self.crystal_symmetry())
    self._map_histograms = get_map_histograms(
        data    = self.map_data(),
        n_slots = 20,
        data_1  = self.map_data_1(),
        data_2  = self.map_data_2())

  def counts(self):
    if not hasattr(self, '_counts'):
      self.get_counts_and_histograms()
    return self._counts

  def histograms(self):
    if not hasattr(self, '_map_histograms'):
      self.get_counts_and_histograms()
    return self._map_histograms

  def as_map_and_model(self):
    '''
      Return this object (allows using .as_map_and_model() on both
      map_and_model objects and others including box.around_model() etc.
    '''
    return self

  def as_map_model_manager(self):
    '''
      Return this object as a map_model_manager
    '''
    from iotbx.map_model_manager import map_model_manager
    mmm = map_model_manager()
    if self.map_manager():
      mmm.add_map(self.map_manager())
    if self.model():
      mmm.add_model(self.model())
    if self.ncs_object():
      mmm.add_ncs_object(self.ncs_object())
    return mmm

  def as_box_object(self,
        original_map_data = None,
        solvent_content = None):
    '''
      Create a box_object for backwards compatibility with methods that used
       extract_box_around_model_and_map
    '''

    if solvent_content:
      self.set_solvent_content(solvent_content)

    if self.model():
       xray_structure_box = self.model().get_xray_structure()
       hierarchy = self.model().get_hierarchy()
    else:
       xray_structure_box = None
       hierarchy = None

    output_box = box_object(
      shift_cart = tuple([-x for x in self.map_manager().origin_shift_cart()]),
      xray_structure_box = xray_structure_box,
      hierarchy = hierarchy,
      ncs_object = self.ncs_object(),
      map_box = self.map_manager().map_data(),
      map_data = original_map_data,
      map_box_half_map_list = None,
      box_crystal_symmetry = self.map_manager().crystal_symmetry(),
      pdb_outside_box_msg = "",
      gridding_first = self._gridding_first,
      gridding_last = self._gridding_last,
      solvent_content = self._solvent_content,
      origin_shift_grid_units = [
         -x for x in self.map_manager().origin_shift_grid_units],
      )
    return output_box

def get_map_histograms(data, n_slots = 20, data_1 = None, data_2 = None):
  h0, h1, h2 = None, None, None
  data_min = None
  hmhcc = None
  if(data_1 is None):
    h0 = flex.histogram(data = data.as_1d(), n_slots = n_slots)
  else:
    data_min = min(flex.min(data_1), flex.min(data_2))
    data_max = max(flex.max(data_1), flex.max(data_2))
    h0 = flex.histogram(data = data.as_1d(), n_slots = n_slots)
    h1 = flex.histogram(data = data_1.as_1d(), data_min = data_min,
      data_max = data_max, n_slots = n_slots)
    h2 = flex.histogram(data = data_2.as_1d(), data_min = data_min,
      data_max = data_max, n_slots = n_slots)
    hmhcc = flex.linear_correlation(
      x = h1.slots().as_double(),
      y = h2.slots().as_double()).coefficient()
  return group_args(h_map = h0, h_half_map_1 = h1, h_half_map_2 = h2,
    _data_min = data_min, half_map_histogram_cc = hmhcc)

def get_map_counts(map_data, crystal_symmetry = None):
  a = map_data.accessor()
  map_counts = group_args(
    origin       = a.origin(),
    last         = a.last(),
    focus        = a.focus(),
    all          = a.all(),
    min_max_mean = map_data.as_1d().min_max_mean().as_tuple(),
    d_min_corner = maptbx.d_min_corner(map_data = map_data,
      unit_cell = crystal_symmetry.unit_cell()))
  return map_counts

def add_tuples(t1, t2):
  new_list = []
  for a, b in zip(t1, t2):
    new_list.append(a+b)
  return tuple(new_list)
