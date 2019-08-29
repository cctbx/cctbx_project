from __future__ import absolute_import, division, print_function
import mmtbx.model
from libtbx.utils import Sorry
from cctbx import maptbx
from libtbx import group_args
from scitbx.array_family import flex

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
               model            = None,
               # where this CS is supposed to come from? After agreing on picking
               # CS here, the only thing it could be useful - to pass CS
               # obtained from command-line args or from parameters. Consider
               # renaming parameter accordingly.
               crystal_symmetry = None,
               box              = True,
               ignore_symmetry_conflicts = False):
    #
    # We should be able to work without symmetry at all. Why not just box
    # model?
    assert [model, crystal_symmetry].count(None) != 2
    if(crystal_symmetry is None and model is not None):
      crystal_symmetry = model.crystal_symmetry()
    if([model, crystal_symmetry].count(None)==0):
      if ignore_symmetry_conflicts: # Take crystal_symmetry if necessary
        if not (model.crystal_symmetry().is_similar_symmetry(crystal_symmetry)):
          model = mmtbx.model.manager(
            model_input = model.get_hierarchy().as_pdb_input(),
            crystal_symmetry = crystal_symmetry)
      else:
        assert model.crystal_symmetry().is_similar_symmetry(crystal_symmetry)
    if(not [map_data_1, map_data_2].count(None) in [0,2]):
      raise Sorry("None or two half-maps are required.")
    #

    # Suggest to get rid of self._model, self._map_data etc to make crystal
    # clear that they are changed in place. Therefore getter functions
    # at the bottom are useless and confusing.
    self._map_data         = map_data
    self._half_map_data_1  = map_data_1
    self._half_map_data_2  = map_data_2
    self._model            = model
    self._crystal_symmetry = crystal_symmetry
    #
    # I don't see any connection between _counts, map_histograms and main
    # purpose of this class (actually, it is function written using class syntax)
    # - shifting origins, cutting boxes, figuring out crystal symmetries.
    # This can be easily done just before calling this and totally separate.
    # I suggest to remove it from here --->
    self._counts = get_map_counts(
      map_data         = self._map_data,
      crystal_symmetry = crystal_symmetry)
    self._map_histograms = get_map_histograms(
      data    = self._map_data,
      n_slots = 20,
      data_1  = self._half_map_data_1,
      data_2  = self._half_map_data_2)
    # <---- End of removing suggestion.
    # Shift origin
    sites_cart = None
    if(self._model is not None):
      sites_cart = self._model.get_sites_cart()
    self.soin = maptbx.shift_origin_if_needed(
      map_data         = self._map_data,
      sites_cart       = sites_cart,
      crystal_symmetry = crystal_symmetry)
    self._original_origin_cart=self.soin.original_origin_cart
    self._original_origin_grid_units=self.soin.original_origin_grid_units
    self.box = None
    self._map_data = self.soin.map_data
    if(self._model is not None):
      self._model.set_sites_cart(sites_cart = self.soin.sites_cart)
    if(self._half_map_data_1 is not None):
      self._half_map_data_1 = maptbx.shift_origin_if_needed(
        map_data         = self._half_map_data_1,
        sites_cart       = None,
        crystal_symmetry = None).map_data
      self._half_map_data_2 = maptbx.shift_origin_if_needed(
        map_data         = self._half_map_data_2,
        sites_cart       = None,
        crystal_symmetry = None).map_data
    # Box
    if(self._model is not None and box):
      xrs = self._model.get_xray_structure()
      if(self._half_map_data_1 is not None):
        self._half_map_data_1 = mmtbx.utils.extract_box_around_model_and_map(
          xray_structure = xrs,
          map_data       = self._half_map_data_1,
          box_cushion    = 5.0).map_box
        self._half_map_data_2 = mmtbx.utils.extract_box_around_model_and_map(
          xray_structure = xrs,
          map_data       = self._half_map_data_2,
          box_cushion    = 5.0).map_box
      self.box = mmtbx.utils.extract_box_around_model_and_map(
        xray_structure = xrs,
        map_data       = self._map_data,
        box_cushion    = 5.0)
      # This should be changed to call model.set_shift_manager(shift_manager=box)
      # For now just call _model.unset_restraints_manager() afterwards.
      self._model.set_xray_structure(xray_structure = self.box.xray_structure_box)
      self._crystal_symmetry = self._model.crystal_symmetry()
      self._map_data = self.box.map_box

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

  def model(self): return self._model

  def xray_structure(self):
    if(self.model() is not None):
      return self.model().get_xray_structure()
    else:
      return None

  def crystal_symmetry(self): return self._crystal_symmetry

  def hierarchy(self): return self._model.get_hierarchy()
