from __future__ import division
import mmtbx.model
from libtbx.utils import null_out
from libtbx.utils import Sorry
from cctbx import maptbx
from libtbx import group_args
from scitbx.array_family import flex

def check_and_set_crystal_symmetry(models=[], map_inps=[], miller_arrays=[]):
  # XXX This should go into a central place
  # XXX Check map gridding here!
  for it in [models, map_inps, miller_arrays]:
    assert isinstance(it, (list, tuple))
  crystal_symmetry = None
  css = []
  all_inputs = models+map_inps+miller_arrays
  for it in all_inputs:
    if(it is not None):
      it = it.crystal_symmetry()
      if(it is None): continue
      if(not [it.unit_cell(), it.space_group()].count(None) in [0,2]):
        raise Sorry("Inconsistent box (aka crystal symmetry) info.")
      if([it.unit_cell(), it.space_group()].count(None)==0):
        css.append(it)
  if(len(css)>1):
    cs0 = css[0]
    for cs in css[1:]:
      if(not cs0.is_similar_symmetry(cs)):
        raise Sorry("Box info (aka crystal symmetry) mismatch across inputs.")
  if(len(css)==0):
    raise Sorry("No box info (aka crystal symmetry) available.")
  crystal_symmetry = css[0]
  for model in models:
    if(model is None): continue
    cs = model.crystal_symmetry()
    if(cs is None or [cs.unit_cell(), cs.space_group()].count(None)==2):
      model.set_crystal_symmetry_if_undefined(crystal_symmetry)
  if(len(map_inps)>1):
    m0 = map_inps[0].map_data()
    for m in map_inps[1:]:
      if(m is None): continue
      maptbx.assert_same_gridding(map_1=m0, map_2=m.map_data())

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

def get_map_counts(map_data, crystal_symmetry):
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

class input(object):
  def __init__(self,
               map_inp,
               map_inp_1 = None,
               map_inp_2 = None,
               pdb_inp   = None,
               box       = True):
    #
    if(not [map_inp_1, map_inp_2].count(None) in [0,2]):
      raise Sorry("None or two half-maps are required.")
    #
    self._map_data        = None
    self._half_map_data_1 = None
    self._half_map_data_2 = None
    self._model = None
    #
    if(pdb_inp is not None):
      self._model = mmtbx.model.manager(model_input = pdb_inp, log = null_out())
    check_and_set_crystal_symmetry(
      models   = [self._model],
      map_inps = [map_inp, map_inp_1, map_inp_2])
    self._map_data = map_inp.map_data()
    self._counts = get_map_counts(
      map_data         = self._map_data,
      crystal_symmetry = map_inp.crystal_symmetry())
    if(map_inp_1 is not None): self._half_map_data_1 = map_inp_1.map_data()
    if(map_inp_2 is not None): self._half_map_data_2 = map_inp_2.map_data()
    self._map_histograms = get_map_histograms(
      data    = self._map_data,
      n_slots = 20,
      data_1  = self._half_map_data_1,
      data_2  = self._half_map_data_2)
    # Shift origin
    sites_cart = None
    if(self._model is not None):
      sites_cart = self._model.get_sites_cart()
    soin = maptbx.shift_origin_if_needed(
      map_data         = self._map_data,
      sites_cart       = sites_cart,
      crystal_symmetry = self._model.crystal_symmetry())
    self._map_data = soin.map_data
    if(self._model is not None):
      self._model.set_sites_cart(sites_cart = soin.sites_cart)
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
      box = mmtbx.utils.extract_box_around_model_and_map(
        xray_structure = xrs,
        map_data       = self._map_data,
        box_cushion    = 5.0)
      self._model.set_xray_structure(xray_structure = box.xray_structure_box)
      self._map_data       = box.map_box

  def counts(self): return self._counts

  def histograms(self): return self._map_histograms

  def map_data(self): return self._map_data

  def half_map_data_1(self): return self._half_map_data_1

  def half_map_data_2(self): return self._half_map_data_2

  def model(self): return self._model

  def xray_structure(self): return self.model().get_xray_structure()

  def crystal_symmetry(self): return self.model().crystal_symmetry()

  def hierarchy(self): return self._model.get_hierarchy()

  def update_maps(self,map_data=None,half_map_data_1=None,half_map_data_2=None):
    # XXX check gridding
    # XXX do we need this?
    if(map_data is not None): self._map_data = map_data
    if(half_map_data_1 is not None): self._half_map_data_1 = half_map_data_1
    if(half_map_data_2 is not None): self._half_map_data_2 = half_map_data_2
