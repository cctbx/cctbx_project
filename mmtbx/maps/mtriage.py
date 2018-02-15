from __future__ import division
import sys
from libtbx import group_args
from cctbx import maptbx
import iotbx.phil
from libtbx import adopt_init_args
from cctbx.maptbx import resolution_from_map_and_model
import mmtbx.utils
from libtbx import group_args
from cctbx import miller, adptbx
from mmtbx.maps import correlation
from mmtbx import masks
from scitbx.array_family import flex
import time
from libtbx.utils import Sorry
from libtbx import introspection
from libtbx.str_utils import size_as_string_with_commas

def show_process_info(out):
  print >> out, "\\/"*39
  introspection.virtual_memory_info().show_if_available(out=out, show_max=True)
  print >> out, "/\\"*39
  out.flush()

master_params_str = """
  scattering_table = wk1995  it1992  n_gaussian  neutron *electron
    .type = choice
    .help = Scattering table (X-ray, neutron or electron)
  compute {
    map_counts = True
      .type = bool
      .help = Compute map counts
    fsc_curve_model = True
      .type = bool
      .help = Compute model-map CC in reciprocal space: FSC(model map, data map)
    d_fsc_model_05 = True
      .type = bool
      .help = Compute d_fsc_model (FSC=0.5)
    d_fsc_model_0 = True
      .type = bool
      .help = Compute d_fsc_model (FSC=0)
    d_fsc_model_0143 = True
      .type = bool
      .help = Compute d_fsc_model (FSC=0.143)
    d_model = True
      .type = bool
      .help = Resolution estimate using model and map
    d_model_b0 = True
      .type = bool
      .help = Resolution estimate using model and map, assumin all atoms B=0
    d99 = True
      .type = bool
      .help = Resolution estimate d99
  }
  resolution = None
    .type = float
    .help = Map resolution (d_FSC)
  mask_maps = None
    .type = bool
    .help = Mask out region outside molecule
  radius_smooth = None
    .type = float
    .help = Mask smoothing radius
    .short_caption = Mask smoothing radius
  nproc = 1
    .type = int
    .help = Number of processors to use
  show_time = False
    .type = bool
    .help = Show individual run times for each step
"""

def get_atom_radius(xray_structure=None, resolution=None, radius=None):
  if(radius is not None): return radius
  radii = []
  if(resolution is not None):
    radii.append(resolution)
  if(xray_structure is not None and resolution is not None):
    b_iso = adptbx.u_as_b(
      flex.mean(xray_structure.extract_u_iso_or_u_equiv()))
    o = maptbx.atom_curves(scattering_type="C", scattering_table="electron")
    rad_image = o.image(d_min=resolution, b_iso=b_iso,
      radius_max=max(15.,resolution), radius_step=0.01).radius
    radii.append(rad_image)
  return max(3, min(10, max(radii)))

def master_params():
  return iotbx.phil.parse(master_params_str, process_includes=False)

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

class base(object):
  def __init__(self,
               map_data,
               crystal_symmetry,
               half_map_data_1=None,
               half_map_data_2=None,
               pdb_hierarchy=None):
    self._map_data         = map_data
    self._crystal_symmetry = crystal_symmetry
    self._half_map_data_1  = half_map_data_1
    self._half_map_data_2  = half_map_data_2
    self._pdb_hierarchy    = pdb_hierarchy
    self._xray_structure = None
    #
    self._validate()
    self._counts = get_map_counts(
      map_data = self._map_data, crystal_symmetry = self._crystal_symmetry)
    self._map_histograms = get_map_histograms(
      data    = self._map_data,
      n_slots = 20,
      data_1  = self._half_map_data_1,
      data_2  = self._half_map_data_2)
    # Shift origin if needed
    sites_cart = None
    if(pdb_hierarchy is not None):
      sites_cart = self._pdb_hierarchy.atoms().extract_xyz()
    soin = maptbx.shift_origin_if_needed(
      map_data         = self._map_data,
      sites_cart       = sites_cart,
      crystal_symmetry = self._crystal_symmetry)
    self._map_data = soin.map_data
    if(pdb_hierarchy is not None):
      self._pdb_hierarchy.atoms().set_xyz(soin.sites_cart)
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
    if(self._pdb_hierarchy is not None):
      self._xray_structure = self._pdb_hierarchy.extract_xray_structure(
        crystal_symmetry = self._crystal_symmetry)
      if(self._half_map_data_1 is not None):
        self._half_map_data_1 = mmtbx.utils.extract_box_around_model_and_map(
          xray_structure = self._xray_structure,
          map_data       = self._half_map_data_1,
          box_cushion    = 5.0).map_box
        self._half_map_data_2 = mmtbx.utils.extract_box_around_model_and_map(
          xray_structure = self._xray_structure,
          map_data       = self._half_map_data_2,
          box_cushion    = 5.0).map_box
      box = mmtbx.utils.extract_box_around_model_and_map(
        xray_structure = self._xray_structure,
        map_data       = self._map_data,
        box_cushion    = 5.0)
      self._pdb_hierarchy.adopt_xray_structure(box.xray_structure_box)
      self._map_data       = box.map_box
      self._xray_structure = box.xray_structure_box
      self._crystal_symmetry = self._xray_structure.crystal_symmetry()

  def counts(self): return self._counts

  def histograms(self): return self._map_histograms

  def map_data(self): return self._map_data

  def half_map_data_1(self): return self._half_map_data_1

  def half_map_data_2(self): return self._half_map_data_2

  def xray_structure(self): return self._xray_structure

  def crystal_symmetry(self): return self._crystal_symmetry

  def pdb_hierarchy(self): return self._pdb_hierarchy

  def update_maps(self,map_data=None,half_map_data_1=None,half_map_data_2=None):
    if(map_data is not None): self._map_data = map_data
    if(half_map_data_1 is not None): self._half_map_data_1 = half_map_data_1
    if(half_map_data_2 is not None): self._half_map_data_2 = half_map_data_2

  def _validate(self):
    if(not [self._half_map_data_1, self._half_map_data_2].count(None) in [0,2]):
      raise Sorry("None or two half-maps are required.")
    if(self._half_map_data_1 is not None):
      correlation.assert_same_gridding(
        map_1 = self._half_map_data_1,
        map_2 = self._half_map_data_2,
        Sorry_message="Half-maps have different gridding.")
      correlation.assert_same_gridding(
        map_1 = self._map_data,
        map_2 = self._half_map_data_2,
        Sorry_message="Half-maps and full map have different gridding.")
    if(self._crystal_symmetry.space_group().type().number()!=1):
      raise Sorry("Symmetry must be P1")
    return self

class caller(object):
  def __init__(self, show=False):
    self.time_cumulative = 0
    self.show=show

  def call(self, f, msg):
    t0 = time.time()
    f()
    sa=size_as_string_with_commas(
      introspection.virtual_memory_info().current_max_sizes().virtual_memory)
    if(self.show):
      delta = time.time()-t0
      self.time_cumulative += delta
      print "%6.2f %8.2f %15s:"%(delta, self.time_cumulative, sa), msg
      sys.stdout.flush()

class mtriage(object):
  def __init__(self,
               map_data,
               crystal_symmetry,
               params          = None,
               half_map_data_1 = None,
               half_map_data_2 = None,
               pdb_hierarchy   = None):
    #adopt_init_args(self, locals())
    self.crystal_symmetry = crystal_symmetry
    self.params           = params
    self.map_data         = map_data.deep_copy()
    self.half_map_data_1 = None
    if(half_map_data_1 is not None):
      self.half_map_data_1 = half_map_data_1.deep_copy()
    self.half_map_data_2 = None
    if(half_map_data_2 is not None):
      self.half_map_data_2 = half_map_data_2.deep_copy()
    self.pdb_hierarchy = None
    if(pdb_hierarchy is not None):
      self.pdb_hierarchy = pdb_hierarchy.deep_copy()
    #
    self.results_masked   = None
    self.results_unmasked = None
    self.time_cumulative  = 0
    if(self.params is None):
      self.params = master_params().extract()
    self.caller = caller(show=self.params.show_time)

  def _create_base(self):
    return base(
      map_data         = self.map_data,
      crystal_symmetry = self.crystal_symmetry,
      half_map_data_1  = self.half_map_data_1,
      half_map_data_2  = self.half_map_data_2,
      pdb_hierarchy    = self.pdb_hierarchy)

  def call(self, func, prefix):
    t0 = time.time()
    result = func()
    sa=size_as_string_with_commas(
      introspection.virtual_memory_info().current_max_sizes().virtual_memory)
    if(self.params.show_time):
      delta = time.time()-t0
      self.time_cumulative += delta
      print "%6.2f %8.2f %15s:"%(delta, self.time_cumulative, sa), prefix
      sys.stdout.flush()
    return result

  def _run(self, base, include_curves, include_mask):
    if(self.params.mask_maps is None):
      # No masking
      self.params.mask_maps = False
      self.results_unmasked = _mtriage(
        base   = base,
        caller = self.caller,
        params = self.params,
      ).run().get_results(
        include_curves = include_curves,
        include_mask   = include_mask)
      # Masking
      if(self.params.radius_smooth is None):
        self.params.radius_smooth = self.results_unmasked.d99
      self.params.mask_maps = True
      self.results_masked = _mtriage(
        base   = base,
        caller = self.caller,
        params = self.params,
      ).run().get_results(
        include_curves = include_curves,
        include_mask   = include_mask)
    else:
      result = _mtriage(
        base   = base,
        caller = self.caller,
        params = self.params,
      ).run().get_results(
        include_curves = include_curves,
        include_mask   = include_mask)
      if(self.params.mask_maps): self.results_masked = result
      else:                      self.results_unmasked = result

  def get_results(self, include_curves, include_mask):
    _base = self.call(func=self._create_base, prefix="Create base")
    self._run(base           = _base,
              include_curves = include_curves,
              include_mask   = include_mask)
    return group_args(
      crystal_symmetry = _base.crystal_symmetry(),
      counts           = _base.counts(),
      histograms       = _base.histograms(),
      masked           = self.results_masked,
      unmasked         = self.results_unmasked)

class _mtriage(object):
  def __init__(self, base, caller, params):
    adopt_init_args(self, locals())
    self.call = self.caller.call
    self.resolution = self.params.resolution
    # Results
    self.d99              = None
    self.d999             = None
    self.d9999            = None
    self.d99_1            = None
    self.d99_2            = None
    self.d_model          = None
    self.d_model_b0       = None
    self.b_iso_overall    = None
    self.d_fsc            = None
    self.d_fsc_model_05   = None
    self.d_fsc_model_0    = None
    self.d_fsc_model_0143 = None
    self.fsc_curve        = None
    self.fsc_curve_model  = None
    self.mask_object      = None
    self.radius_smooth    = self.params.radius_smooth
    self.d_corner         = None
    self.d9999           = None
    # Info (results)
    self.map_counts        = None
    self.half_map_1_counts = None
    self.half_map_2_counts = None
    self.map_histograms    = None
    # Internal work objects
    self.f_map    = None
    self.f_map_1  = None
    self.f_map_2  = None
    self.f_calc   = None

  def run(self):
    # Compute radius
    self.call(f=self._compute_radius, msg="Compute radius")
    # Compute and apply mask
    self.call(f=self._compute_and_apply_mask, msg="Masking")
    # Compute F_maps
    self.call(f=self._compute_f_maps, msg="Compute F_maps")
    # Compute d99
    self.call(f=self._compute_d99, msg="Compute d99")
#    # Compute d_model at B=0
#    self.call(f=self._compute_d_model_b0, msg="Compute d_model_b0")
    # Compute half-map FSC
    self.call(f=self._compute_half_map_fsc, msg="Compute half-map FSC")
    # Compute Fcalc
    self.call(f=self._compute_f_calc, msg="Compute Fcalc")
    # Map-model FSC curve
    self.call(f=self._compute_fsc_curve_model, msg="Compute fsc_curve_model")
    # d_fsc_model_0
    self.call(f=self._compute_f_fsc_model_0, msg="Compute d_fsc_model_0")
    # d_fsc_model_0143
    self.call(f=self._compute_f_fsc_model_0143, msg="Compute d_fsc_model_0143")
    # d_fsc_model_05
    self.call(f=self._compute_f_fsc_model_05, msg="Compute d_fsc_model_05")
    # Compute d_model
    self.call(f=self._compute_d_model, msg="Compute d_model")
    return self

  def _shift_origin(self):
    sites_cart = None
    if(self.xray_structure is not None):
      sites_cart = self.xray_structure.sites_cart()
    soin = maptbx.shift_origin_if_needed(
      map_data         = self.map_data,
      sites_cart       = sites_cart,
      crystal_symmetry = self.crystal_symmetry)
    self.map_data = soin.map_data
    if(self.xray_structure is not None):
      self.xray_structure.set_sites_cart(soin.sites_cart)
      self.pdb_hierarchy.atoms().set_xyz(soin.sites_cart)
    if(self.half_map_data_1 is not None):
      soin = maptbx.shift_origin_if_needed(
        map_data         = self.half_map_data_1,
        sites_cart       = None,
        crystal_symmetry = None)
      self.half_map_data_1 = soin.map_data
      soin = maptbx.shift_origin_if_needed(
        map_data         = self.half_map_data_2,
        sites_cart       = None,
        crystal_symmetry = None)
      self.half_map_data_2 = soin.map_data

  def _compute_radius(self):
    if(not self.params.mask_maps): return
    if(self.base.pdb_hierarchy() is None): return
    self.radius_smooth = get_atom_radius(
      xray_structure   = self.base.xray_structure(),
      radius           = self.radius_smooth,
      resolution       = self.resolution)

  def _compute_and_apply_mask(self):
    if(not self.params.mask_maps): return
    if(self.base.pdb_hierarchy() is None): return
    mask_smooth = masks.smooth_mask(
      xray_structure = self.base.xray_structure(),
      n_real         = self.base.map_data().all(),
      rad_smooth     = self.radius_smooth).mask_smooth
    self.base.update_maps(map_data = self.base.map_data()*mask_smooth)
    if(self.base.half_map_data_1() is not None):
      self.base.update_maps(
        half_map_data_1 = self.base.half_map_data_1()*mask_smooth,
        half_map_data_2 = self.base.half_map_data_2()*mask_smooth)

  def _compute_f_maps(self):
    self.f_map = miller.structure_factor_box_from_map(
      map              = self.base.map_data(),
      crystal_symmetry = self.base.crystal_symmetry())
    if(self.base.half_map_data_1() is not None):
      self.f_map_1 = miller.structure_factor_box_from_map(
        map              = self.base.half_map_data_1(),
        crystal_symmetry = self.base.crystal_symmetry())
      self.f_map_2 = miller.structure_factor_box_from_map(
        map              = self.base.half_map_data_2(),
        crystal_symmetry = self.base.crystal_symmetry())

  def _compute_d99(self):
    if(not self.params.compute.d99): return
    d99 = maptbx.d99(f_map = self.f_map)
    self.d99    = d99.result.d99
    self.d999   = d99.result.d999
    self.d9999  = d99.result.d9999
    self.f_map = self.f_map.resolution_filter(d_min = self.d9999-0.1) # TRUNCATED!
    d99_obj_1, d99_obj_2 = None,None
    if(self.base.half_map_data_1() is not None):
      d99_1 = maptbx.d99(
        map              = self.base.half_map_data_1(),
        crystal_symmetry = self.base.crystal_symmetry())
      d99_2 = maptbx.d99(
        map              = self.base.half_map_data_2(),
        crystal_symmetry = self.base.crystal_symmetry())
      self.d99_1 = d99_1.result.d99
      self.d99_2 = d99_2.result.d99
      self.f_map_1 = d99_1.f_map
      self.f_map_2 = d99_2.f_map

  def _compute_f_calc(self):
    self.f_calc = self.f_map.structure_factors_from_scatterers(
      xray_structure = self.base.xray_structure()).f_calc()

  def _compute_d_model(self):
    if(not self.params.compute.d_model): return
    if(self.base.pdb_hierarchy() is not None):
      o = resolution_from_map_and_model.run(
        f_map            = self.f_map,
        d_fsc_model      = self.d_fsc_model_0,
        xray_structure   = self.base.xray_structure())
      self.d_model       = o.d_min
      self.b_iso_overall = o.b_iso
      self.d_model_b0    = o.d_model_b0

  def _compute_half_map_fsc(self):
    if(self.base.half_map_data_1() is not None):
      self.fsc_curve = self.f_map_1.d_min_from_fsc(
        other = self.f_map_2, bin_width=100, fsc_cutoff=0.143)
      self.d_fsc = self.fsc_curve.d_min

  def _compute_fsc_curve_model(self):
    if(not self.params.compute.fsc_curve_model): return
    if(self.base.pdb_hierarchy() is not None):
      self.fsc_curve_model = self.f_calc.fsc(
        other=self.f_map, bin_width=100)

  def _compute_f_fsc_model_0(self):
    if(not self.params.compute.d_fsc_model_0): return
    assert self.fsc_curve_model is not None
    self.d_fsc_model_0 = self.f_calc.d_min_from_fsc(
      fsc_curve=self.fsc_curve_model, fsc_cutoff=0.).d_min

  def _compute_f_fsc_model_0143(self):
    if(not self.params.compute.d_fsc_model_0143): return
    assert self.fsc_curve_model is not None
    self.d_fsc_model_0143 = self.f_calc.d_min_from_fsc(
      fsc_curve=self.fsc_curve_model, fsc_cutoff=0.143).d_min

  def _compute_f_fsc_model_05(self):
    if(not self.params.compute.d_fsc_model_05): return
    assert self.fsc_curve_model is not None
    self.d_fsc_model_05 = self.f_calc.d_min_from_fsc(
      fsc_curve=self.fsc_curve_model, fsc_cutoff=0.5).d_min

  def get_results(self, include_curves, include_mask):
    mask = None
    if(self.mask_object is not None and include_mask):
      mask = self.mask_object.mask_smooth
    map_histograms  = None
    fsc_curve       = None
    fsc_curve_model = None
    if(include_curves):
      map_histograms  = self.map_histograms
      fsc_curve       = self.fsc_curve
      fsc_curve_model = self.fsc_curve_model
    return group_args(
      d99               = self.d99,
      d999              = self.d999,
      d9999             = self.d9999,
      d99_1             = self.d99_1,
      d99_2             = self.d99_2,
      d_model           = self.d_model,
      d_model_b0        = self.d_model_b0,
      b_iso_overall     = self.b_iso_overall,
      d_fsc             = self.d_fsc,
      d_fsc_model_05    = self.d_fsc_model_05,
      d_fsc_model_0     = self.d_fsc_model_0,
      d_fsc_model_0143  = self.d_fsc_model_0143,
      fsc_curve         = fsc_curve,
      fsc_curve_model   = fsc_curve_model,
      mask              = mask,
      radius_smooth     = self.radius_smooth)

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
