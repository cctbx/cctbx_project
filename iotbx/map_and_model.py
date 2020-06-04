from __future__ import absolute_import, division, print_function
import mmtbx.model
from libtbx.utils import Sorry
from cctbx import maptbx
from libtbx import group_args
from scitbx.array_family import flex
from iotbx.map_model_manager import map_model_manager

class input(object):
  def __init__(self,
               map_manager      = None,  # replaces map_data
               map_manager_1    = None,  # replaces map_data_1
               map_manager_2    = None,  # replaces map_data_2
               map_manager_list = None,  # replaces map_data_list
               model            = None,
               crystal_symmetry = None,  # optional, used only to check
               box              = True,
               soft_mask        = None,
               resolution       = None, # required for soft_mask
               ignore_symmetry_conflicts = False):

    self._model=model
    self._map_manager=map_manager
    self._map_manager_1=map_manager_1
    self._map_manager_2=map_manager_2
    self._map_manager_list=map_manager_list
    self._crystal_symmetry=crystal_symmetry

    # Take cs from model or self.map_manager if not specified
    if self._map_manager and not self._crystal_symmetry:
      self._crystal_symmetry=self._map_manager.crystal_symmetry()
    if self._model and not self._crystal_symmetry:
      self._crystal_symmetry=self._model.crystal_symmetry()

    # Decide what to do if conflicting crystal_symmetry
    if self._crystal_symmetry and self._model and (
       not self._crystal_symmetry.is_similar_symmetry(
          self._model.crystal_symmetry())):
      if ignore_symmetry_conflicts: # take crystal_symmetry overwrite model
        from libtbx.utils import null_out
        self._model=mmtbx.model.manager(
             model_input = self._model.get_hierarchy().as_pdb_input(),
             crystal_symmetry = self._crystal_symmetry,
            log = null_out())
      else: # stop
        assert self._crystal_symmetry.is_similar_symmetry(
         self._model.crystal_symmetry())

    # Make sure we have what is expected: optional model, mm,
    # self._map_manager_1 and self._map_manager_2 or neither,
    #   optional list of self._map_manager_list

    if not self._map_manager_list:
      self._map_manager_list=[]

    if(not [self._map_manager_1, self._map_manager_2].count(None) in [0,2]):
      raise Sorry("None or two half-maps are required.")
    if(not self._map_manager):
      raise Sorry("A map is required.")

    # Make sure all map_managers have same gridding and symmetry
    for m in [self._map_manager_1,self._map_manager_2]+ \
         self._map_manager_list:
      if m:
        assert self._map_manager.is_similar(m)


    # If model, make a map_model_manager with model and mm and
    #  let it check symmetry against this one
    mmm=map_model_manager()
    mmm.add_map_manager(self._map_manager)
    if self._model:
      mmm.add_model(self._model)
    # ALL OK if it does not stop

    # Shift origin of model and map_manager to (0,0,0) with
    #    mmm which knows about both
    mmm.shift_origin()
    self._model=mmm.model()
    self._map_manager=mmm.map_manager()

    # Shift origins of all other maps
    for m in [self._map_manager_1,self._map_manager_2]+\
         self._map_manager_list:
      if m:
        m.shift_origin()

    # Make sure all really match:
    for m in [self._map_manager_1,self._map_manager_2]+\
        self._map_manager_list:
      if m:
        assert self._map_manager.is_similar(m)

    # Save things so they can be accessed easily later

    self._original_origin_grid_units=self._map_manager.origin_shift_grid_units
    self._original_origin_cart=self._map_manager.origin_shift_cart()

    box_result = None

    if(box and model is not None):

      if(self._map_manager_1 is not None):
        info_1 = mmtbx.utils.extract_box_around_model_and_map(
          model          = self._model,
          mm             = self._map_manager_1,
          soft_mask      = soft_mask,
          resolution      = resolution,
          box_cushion    = 5.0)
        self._map_manager_1=self._map_manager_1.customized_copy(
           map_data=info_1.map_box,
           origin_shift_grid_units=info_1.origin_shift_grid_units())

        info_2 = mmtbx.utils.extract_box_around_model_and_map(
          model          = self._model,
          mm             = self._map_manager_2,
          soft_mask      = soft_mask,
          resolution      = resolution,
          box_cushion    = 5.0)
        self._map_manager_2=self._map_manager_2.customized_copy(
           map_data=info_2.map_box,
           origin_shift_grid_units=info_2.origin_shift_grid_units())

      if self._map_manager_list:
        new_list=[]
        for x in self._map_manager_list:
          info = mmtbx.utils.extract_box_around_model_and_map(
            model          = self._model,
            mm             = x,
            soft_mask      = soft_mask,
            resolution      = resolution,
            box_cushion    = 5.0)
          new_list.append(x.customized_copy(
            map_data=info.map_box,
            origin_shift_grid_units=info.origin_shift_grid_units()))
        self._map_manager_list=new_list
      box_result = mmtbx.utils.extract_box_around_model_and_map(
        model          = self._model,
        mm             = self._map_manager,
        soft_mask      = soft_mask,
        resolution     = resolution,
        box_cushion    = 5.0)
      self._map_manager=self._map_manager.customized_copy(
           map_data=box_result.map_box,
           origin_shift_grid_units=box_result.origin_shift_grid_units())

      # Update model and crystal_symmetry with new values
      self._model.set_shift_manager(shift_manager= box_result)
      self._crystal_symmetry = self._model.crystal_symmetry()
      assert self._crystal_symmetry.is_similar_symmetry(
        self._map_manager.crystal_symmetry())

  def original_origin_cart(self):
    assert self._original_origin_cart is not None
    return self._original_origin_cart

  def original_origin_grid_units(self):
    assert self._original_origin_grid_units is not None
    return self._original_origin_grid_units

  def map_data(self):
    return self.map_manager().map_data()

  def map_data_1(self):
    if self.map_manager_1():
      return self.map_manager_1().map_data()

  def map_data_2(self):
    if self.map_manager_2():
      return self.map_manager_2().map_data()

  def map_data_list(self):
    map_data_list=[]
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

  def crystal_symmetry(self): return self._crystal_symmetry

  def xray_structure(self):
    if(self.model() is not None):
      return self.model().get_xray_structure()
    else:
      return None

  def hierarchy(self): return self._model.get_hierarchy()

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
    if not hasattr(self,'_counts'):
      self.get_counts_and_histograms()
    return self._counts

  def histograms(self):
    if not hasattr(self,'_map_histograms'):
      self.get_counts_and_histograms()
    return self._map_histograms

def get_map_histograms(data, n_slots=20, data_1=None, data_2=None):
  h0, h1, h2 = None, None, None
  data_min = None
  hmhcc = None
  if(data_1 is None):
    h0 = flex.histogram(data = data.as_1d(), n_slots = n_slots)
  else:
    data_min = min(flex.min(data_1), flex.min(data_2))
    data_max = max(flex.max(data_1), flex.max(data_2))
    h0 = flex.histogram(data = data.as_1d(), n_slots = n_slots)
    h1 = flex.histogram(data = data_1.as_1d(), data_min=data_min,
      data_max=data_max, n_slots = n_slots)
    h2 = flex.histogram(data = data_2.as_1d(), data_min=data_min,
      data_max=data_max, n_slots = n_slots)
    hmhcc = flex.linear_correlation(
      x=h1.slots().as_double(),
      y=h2.slots().as_double()).coefficient()
  return group_args(h_map = h0, h_half_map_1 = h1, h_half_map_2 = h2,
    _data_min = data_min, half_map_histogram_cc = hmhcc)

def get_map_counts(map_data, crystal_symmetry=None):
  a = map_data.accessor()
  map_counts = group_args(
    origin       = a.origin(),
    last         = a.last(),
    focus        = a.focus(),
    all          = a.all(),
    min_max_mean = map_data.as_1d().min_max_mean().as_tuple(),
    d_min_corner = maptbx.d_min_corner(map_data=map_data,
      unit_cell = crystal_symmetry.unit_cell()))
  return map_counts


