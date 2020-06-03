from __future__ import absolute_import, division, print_function
import mmtbx.model
from libtbx.utils import Sorry
from cctbx import maptbx
from libtbx import group_args
from scitbx.array_family import flex
from iotbx.map_manager import map_manager
from iotbx.map_model_manager import map_model_manager

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

# I suggest to rewrite it as a plain function returning crystal_symmetry.
# That's all what is needed. Everything else is being changed in place.
class input(object):
  def __init__(self,
               # I suggest we use map_input with symmetries here instead of
               # raw map_data and call something like mtriage.py:check_and_set_crystal_symmetry
               # to set CS consistently for maps and model.
               # Especially because DataManager is supposed to provide these
               # map_input objects.
               # consider reusing/replacing crystal.select_crystal_symmetry()

               # Warning! Model and all map_data are being changed in place.
               # This warning should remain here

               map_data         = None, # whole_map_input would be a better name?
               map_data_1       = None, # half_map_input_1 would be a better name?
               map_data_2       = None, # half_map_input_2 would be a better name?
               map_data_list    = None,  # optional list of any map data files
               mm               = None,  # replaces map_data
               mm_1             = None,  # replaces map_data_1
               mm_2             = None,  # replaces map_data_2
               mm_list          = None,  # replaces map_data_list
               model            = None,
               crystal_symmetry = None,  # optional, used only to check
               box              = True,
               soft_mask        = None,
               resolution       = None, # required for soft_mask
               ignore_symmetry_conflicts = False):

    # Take cs from model or mm if not specified
    if model and not crystal_symmetry:
      crystal_symmetry=model.crystal_symmetry()
    if mm and not crystal_symmetry:
      crystal_symmetry=mm.crystal_symmetry()

    # Decide what to do if conflicting crystal_symmetry
    if crystal_symmetry and model and (
       not crystal_symmetry.is_similar_symmetry(model.crystal_symmetry())):
      if ignore_symmetry_conflicts: # take crystal_symmetry
        from libtbx.utils import null_out
        model=mmtbx.model.manager(
             model_input = model.get_hierarchy().as_pdb_input(),
             crystal_symmetry = crystal_symmetry,
            log = null_out())
      else: # stop
        assert crystal_symmetry.is_similar_symmetry(model.crystal_symmetry())

    # Switch from map_data to mm inputs XXX IN PROGRESS

    if map_data and not mm:
       # Assumes crystal_symmetry model.crystal_symmetry and
       #  unit_cell_crystal_symmetry are all the same thing
       if crystal_symmetry and model:
         assert crystal_symmetry.is_similar_symmetry(model.crystal_symmetry())
       mm=map_manager(map_data=map_data,unit_cell_grid=map_data.all(),
         unit_cell_crystal_symmetry=crystal_symmetry)
    if map_data_1 and not mm_1:
       mm_1=map_manager(map_data=map_data_1,unit_cell_grid=map_data_1.all(),
         unit_cell_crystal_symmetry=crystal_symmetry)
    if map_data_2 and not mm_2:
       mm_2=map_manager(map_data=map_data_2,unit_cell_grid=map_data_2.all(),
         unit_cell_crystal_symmetry=crystal_symmetry)
    if map_data_list and not mm_list:
      mm_list=[]
      for md in map_data_list:
        mm_list.append(map_manager(map_data=md,unit_cell_grid=md.all(),
          unit_cell_crystal_symmetry=crystal_symmetry))
    if not mm_list:
      mm_list=[]

    # Make sure we have what is expected: optional model, mm,
    # mm_1 and mm_2 or neither,  optional list of mm_list

    if(not [mm_1, mm_2].count(None) in [0,2]):
      raise Sorry("None or two half-maps are required.")
    if(not mm):
      raise Sorry("A map is required.")

    # Make sure all map_managers have same gridding and symmmetry
    for m in [mm_1,mm_2]+mm_list:
      if m:
        assert mm.is_similar(m)


    # If model, make a map_model_manager with model and mm and
    #  let it check symmetry against this one
    mmm=map_model_manager()
    mmm.add_map_manager(mm)
    if model:
      mmm.add_model(model)
    # ALL OK

    # Shift origin of model and mm with mmm which knows about both
    mmm.shift_origin()
    model=mmm.model()
    mm=mmm.map_manager()

    # Shift origins of all other maps
    for m in [mm_1,mm_2]+mm_list:
      if m:
        m.shift_origin()

    # Make sure all really match:
    for m in [mm_1,mm_2]+mm_list:
      if m:
        assert mm.is_similar(m)

    # Save things so they can be accessed easily later

    self._original_origin_grid_units=mm.origin_shift_grid_units
    self._original_origin_cart=mm.origin_shift_cart()

    self.box = None
    self._model=model
    self._map_data=mm.map_data()
    self._crystal_symmetry = crystal_symmetry
    if mm_1:
      self._half_map_data_1 = mm_1.map_data()
      self._half_map_data_2 = mm_2.map_data()
    else:
      self._half_map_data_1 = None
      self._half_map_data_2 = None
    self._map_data_list=None
    if mm_list:
      self._map_data_list=[]
      for x in mm_list:
        self._map_data_list.append(x.map_data())

    # Move elsewhere XXX
    self.get_counts_and_histograms(mm=mm,mm_1=mm_1,mm_2=mm_2,
      crystal_symmetry=crystal_symmetry) # move XXX

    if(box and model is not None):

      if(mm_1 is not None):
        info_1 = mmtbx.utils.extract_box_around_model_and_map(
          model          = model,
          mm             = mm_1,
          soft_mask      = soft_mask,
          resolution      = resolution,
          box_cushion    = 5.0)
        mm_1=mm_1.customized_copy(map_data=info_1.map_box,
           origin_shift_grid_units=info_1.origin_shift_grid_units())
        info_2 = mmtbx.utils.extract_box_around_model_and_map(
          model          = model,
          mm             = mm_2,
          soft_mask      = soft_mask,
          resolution      = resolution,
          box_cushion    = 5.0)
        mm_2=mm_2.customized_copy(map_data=info_2.map_box,
           origin_shift_grid_units=info_2.origin_shift_grid_units())

      if mm_list:
        new_list=[]
        for x in mm_list:
          info = mmtbx.utils.extract_box_around_model_and_map(
            model          = model,
            mm             = x,
            soft_mask      = soft_mask,
            resolution      = resolution,
            box_cushion    = 5.0)
          new_list.append(x.customized_copy(map_data=info.map_box,
            origin_shift_grid_units=info_2.origin_shift_grid_units()))
        mm_list=new_list
      self.box = mmtbx.utils.extract_box_around_model_and_map(
        model          = model,
        mm             = mm,
        soft_mask      = soft_mask,
        resolution     = resolution,
        box_cushion    = 5.0)

      # Update everything with new values
      self._model=model
      self._model.set_xray_structure(xray_structure = self.box.xray_structure_box)

      self._crystal_symmetry = self._model.crystal_symmetry()
      self._map_data = self.box.map_box
      if mm_1:
        self._half_map_data_1 = mm_1.map_data()
        self._half_map_data_2 = mm_2.map_data()
      self._map_data_list=None
      if mm_list:
        self._map_data_list=[]
        for x in mm_list:
          self._map_data_list.append(x.map_data())

  def get_counts_and_histograms(self,
    mm=None,
    mm_1=None,
    mm_2=None,
    crystal_symmetry=None):
    self._counts = get_map_counts(
      map_data         = mm.map_data(),
      crystal_symmetry = crystal_symmetry)
    if mm_1:
      mm_1_data=mm_1.map_data()
      mm_2_data=mm_2.map_data()
    else:
      mm_1_data=None
      mm_2_data=None
    self._map_histograms = get_map_histograms(
        data    = mm.map_data(),
        n_slots = 20,
        data_1  = mm_1_data,
        data_2  = mm_2_data)

  def original_origin_cart(self):
    assert self._original_origin_cart is not None
    return self._original_origin_cart

  def original_origin_grid_units(self):
    assert self._original_origin_grid_units is not None
    return self._original_origin_grid_units

  def counts(self): return self._counts

  def histograms(self): return self._map_histograms

  def map_data(self): return self._map_data

  def map_data_1(self): return self._half_map_data_1

  def map_data_2(self): return self._half_map_data_2

  def map_data_list(self): return self._map_data_list

  def model(self): return self._model

  def crystal_symmetry(self): return self._crystal_symmetry

  def xray_structure(self):
    if(self.model() is not None):
      return self.model().get_xray_structure()
    else:
      return None

  def hierarchy(self): return self._model.get_hierarchy()
