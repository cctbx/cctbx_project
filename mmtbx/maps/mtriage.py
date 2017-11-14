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

master_params_str = """
  scattering_table = wk1995  it1992  n_gaussian  neutron *electron
    .type = choice
    .help = Scattering table (X-ray, neutron or electron)
  compute_fsc_curve_model = True
    .type = bool
    .help = Compute model-map CC in reciprocal space: FSC(model map, data map)
  compute_d_model = True
    .type = bool
    .help = Resolution estimate using model and map
  compute_d99 = True
    .type = bool
    .help = Resolution estimate d99
  mask_maps = True
    .type = bool
    .help = Mask out region outside molecule
  radius_smooth = None
    .type = float
    .help = Mask smoothing radius (by default set to max(10,d99)
    .short_caption = Mask smoothing radius
  nproc = 1
    .type = int
    .help = Number of processors to use
"""

def get_atom_radius(xray_structure=None, d_min=None, map_data=None,
                    crystal_symmetry=None, radius=None):
  if(radius is not None): return radius
  radii = []
  if(d_min is not None):
    radii.append(d_min)
  if([xray_structure, crystal_symmetry].count(None)==0):
    assert crystal_symmetry.is_similar_symmetry(
      xray_structure.crystal_symmetry())
  if([map_data, crystal_symmetry].count(None)==0):
    d99 = maptbx.d99(
      map              = map_data,
      crystal_symmetry = crystal_symmetry).result.d99
    radii.append(d99)
  if(xray_structure is not None and d_min is not None):
    b_iso = adptbx.u_as_b(
      flex.mean(xray_structure.extract_u_iso_or_u_equiv()))
    o = maptbx.atom_curves(scattering_type="C", scattering_table="electron")
    rad_image = o.image(d_min=d_min, b_iso=b_iso,
      radius_max=max(15.,d_min), radius_step=0.01).radius
    radii.append(rad_image)
  return max(3, min(10, max(radii)))

def master_params():
  return iotbx.phil.parse(master_params_str, process_includes=False)

def get_box(map_data, pdb_hierarchy, xray_structure):
  if(pdb_hierarchy is not None):
    box = mmtbx.utils.extract_box_around_model_and_map(
      xray_structure         = xray_structure,
      map_data               = map_data,
      box_cushion            = 5.0,
      selection              = None,
      density_select         = None,
      threshold              = None)
    pdb_hierarchy.adopt_xray_structure(box.xray_structure_box)
    return group_args(
      map_data       = box.map_box,
      xray_structure = box.xray_structure_box,
      pdb_hierarchy  = pdb_hierarchy)
  else:
    return None

def get_map_histograms(data, n_slots=20, data_1=None, data_2=None):
  h0, h1, h2 = None, None, None
  data_min = None
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

def get_map_counts(map_data):
  if(map_data is None): return None
  a = map_data.accessor()
  map_counts = group_args(
    origin       = a.origin(),
    last         = a.last(),
    focus        = a.focus(),
    all          = a.all(),
    min_max_mean = map_data.as_1d().min_max_mean().as_tuple())
  return map_counts

class mtriage(object):
  def __init__(self,
               map_data,
               crystal_symmetry,
               params=master_params().extract(),
               half_map_data_1=None,
               half_map_data_2=None,
               pdb_hierarchy=None):
    adopt_init_args(self, locals())
    # Objects may be altered inside
    self.map_data = self.map_data.deep_copy()
    if(self.half_map_data_1 is not None):
      self.half_map_data_1 = self.half_map_data_1.deep_copy()
      self.half_map_data_2 = self.half_map_data_2.deep_copy()
    if(self.pdb_hierarchy is not None):
      self.pdb_hierarchy = self.pdb_hierarchy.deep_copy()
    #
    assert [half_map_data_1, half_map_data_2].count(None) in [0,2]
    # Results
    self.d9               = None
    self.d99              = None
    self.d999             = None
    self.d99_1            = None
    self.d99_2            = None
    self.d_model          = None
    self.d_model_b0       = None
    self.b_iso_overall    = None
    self.d_fsc            = None
    self.d_fsc_model      = None
    self.d_fsc_model_0    = None
    self.d_fsc_model_0143 = None
    self.fsc_curve        = None
    self.fsc_curve_model  = None
    self.mask_object      = None
    self.radius_smooth    = self.params.radius_smooth
    # Info (results)
    self.crystal_symmetry = crystal_symmetry
    self.map_counts        = get_map_counts(map_data = self.map_data)
    self.half_map_1_counts = get_map_counts(map_data = self.half_map_data_1)
    self.half_map_2_counts = get_map_counts(map_data = self.half_map_data_2)
    self.map_histograms = get_map_histograms(
      data    = self.map_data,
      n_slots = 20,
      data_1  = self.half_map_data_1,
      data_2  = self.half_map_data_2)
    # Internal work objects
    self.f   = None
    self.f1  = None
    self.f2  = None
    self.box = None
    self.xray_structure = None

  def validate(self):
    if(not [self.half_map_data_1, self.half_map_data_2].count(None) in [0,2]):
      raise Sorry("None or two half-maps are required.")
    if(self.half_map_data_1 is not None):
      correlation.assert_same_gridding(
        map_1 = self.half_map_data_1,
        map_2 = self.half_map_data_2,
        Sorry_message="Half-maps have different gridding.")
      correlation.assert_same_gridding(
        map_1 = self.map_data,
        map_2 = self.half_map_data_2,
        Sorry_message="Half-maps and full map have different gridding.")
    if(self.crystal_symmetry.space_group().type().number()!=1):
      raise Sorry("Symmetry must be P1")

  def call(self, func, prefix, show_time=False):
    t0 = time.time()
    func()
    if(show_time):
      print prefix, ":", time.time()-t0
      sys.stdout.flush()

  def run(self):
    # Extract xrs from pdb_hierarchy
    self.call(func=self._get_xray_structure, prefix="xrs from pdb_hierarchy")
    # Shift origin if needed
    self.call(func=self._shift_origin, prefix="Shift origin if needed")
    # Compute mask
    self.call(func=self._compute_mask, prefix="Compute mask")
    # Apply mask to map data
    self.call(func=self._apply_mask, prefix="Apply mask to map data")
    # Extract box around model with map
    self.call(func=self._get_box, prefix="Extract box around model with map")
    # Compute d99
    self.call(func=self._compute_d99, prefix="Compute d99")
    # Compute d_model at B=0
    self.call(func=self._compute_d_model_b0, prefix="Compute d_model_b0")
    # Compute d_model
    self.call(func=self._compute_d_model, prefix="Compute d_model")
    # Compute half-map FSC
    self.call(func=self._compute_half_map_fsc, prefix="Compute half-map FSC")
    # Map-model FSC and d_fsc_model
    self.call(func=self._compute_model_map_fsc, prefix="Map-model FSC and d_fsc_model")
    return self

  def _get_xray_structure(self):
    if(self.pdb_hierarchy is not None):
      self.pdb_hierarchy.atoms().reset_i_seq()
      self.xray_structure = self.pdb_hierarchy.extract_xray_structure(
        crystal_symmetry = self.crystal_symmetry)

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

  def _compute_mask(self):
    if(not self.params.mask_maps): return
    if(self.pdb_hierarchy is None): return
    self.radius_smooth = get_atom_radius(
      xray_structure   = self.xray_structure,
      map_data         = self.map_data,
      crystal_symmetry = self.crystal_symmetry,
      radius           = self.radius_smooth)
    self.mask_object = masks.smooth_mask(
      xray_structure = self.xray_structure,
      n_real         = self.map_data.all(),
      rad_smooth     = self.radius_smooth)

  def _apply_mask(self):
    if(not self.params.mask_maps): return
    if(self.mask_object is None): return
    self.map_data = self.map_data * self.mask_object.mask_smooth
    if(self.half_map_data_1 is not None):
      self.half_map_data_1 = self.half_map_data_1 * self.mask_object.mask_smooth
      self.half_map_data_2 = self.half_map_data_2 * self.mask_object.mask_smooth

  def _get_box(self):
    if(self.pdb_hierarchy is not None):
      self.box = get_box(
        map_data       = self.map_data,
        pdb_hierarchy  = self.pdb_hierarchy,
        xray_structure = self.xray_structure)

  def _compute_d99(self):
    if(not self.params.compute_d99): return
    d99_obj = maptbx.d99(
      map              = self.map_data,
      crystal_symmetry = self.crystal_symmetry)
    self.d9   = d99_obj.result.d9
    self.d99  = d99_obj.result.d99
    self.d999 = d99_obj.result.d999
    self.f = d99_obj.f
    d99_obj_1, d99_obj_2 = None,None
    if(self.half_map_data_1 is not None):
      d99_obj_1 = maptbx.d99(
        map              = self.half_map_data_1,
        crystal_symmetry = self.crystal_symmetry)
      d99_obj_2 = maptbx.d99(
        map              = self.half_map_data_2,
        crystal_symmetry = self.crystal_symmetry)
      self.d99_1 = d99_obj_1.result.d99
      self.d99_2 = d99_obj_2.result.d99
      self.f1 = d99_obj_1.f
      self.f2 = d99_obj_2.f

  def _compute_d_model_b0(self):
    o = resolution_from_map_and_model.run_at_b0(
      map_data         = self.box.map_data,
      xray_structure   = self.box.xray_structure,
      d_min_min        = 1.7)
    self.d_model_b0 = o.d_min

  def _compute_d_model(self):
    if(not self.params.compute_d_model): return
    if(self.pdb_hierarchy is not None):
      o = resolution_from_map_and_model.run(
        map_data         = self.box.map_data,
        xray_structure   = self.box.xray_structure,
        pdb_hierarchy    = self.box.pdb_hierarchy,
        d_min_min        = 1.7,
        nproc            = self.params.nproc)
      self.d_model       = o.d_min
      self.b_iso_overall = o.b_iso

  def _compute_half_map_fsc(self):
    if(self.half_map_data_1 is not None):
      self.fsc_curve = self.f1.d_min_from_fsc(
        other = self.f2, bin_width=100, fsc_cutoff=0.143)
      self.d_fsc = self.fsc_curve.d_min

  def _compute_model_map_fsc(self):
    if(not self.params.compute_fsc_curve_model): return
    if(self.pdb_hierarchy is not None):
      f_obs = miller.structure_factor_box_from_map(
        map              = self.box.map_data,
        crystal_symmetry = self.box.xray_structure.crystal_symmetry())
      f_calc = f_obs.structure_factors_from_scatterers(
        xray_structure = self.box.xray_structure).f_calc()
      self.fsc_curve_model = f_calc.d_min_from_fsc(
        other=f_obs, bin_width=100, fsc_cutoff=0.5)
      self.d_fsc_model = self.fsc_curve_model.d_min
      self.d_fsc_model_0 = f_calc.d_min_from_fsc(
        other=f_obs, bin_width=100, fsc_cutoff=0.).d_min
      self.d_fsc_model_0143 = f_calc.d_min_from_fsc(
        other=f_obs, bin_width=100, fsc_cutoff=0.143).d_min

  def write_fsc_curve_model_plot_data(self, file_name):
    if(self.fsc_curve_model is not None):
      of = open(file_name,"w")
      for a,b in zip(self.fsc_curve_model.fsc.d_inv,
                     self.fsc_curve_model.fsc.fsc):
        print >> of, "%15.9f %15.9f"%(a,b)
      of.close()

  def write_fsc_curve_plot_data(self, file_name):
    if(self.fsc_curve is not None):
      of = open(file_name,"w")
      for a,b in zip(self.fsc_curve.fsc.d_inv, self.fsc_curve.fsc.fsc):
        print >> of, "%15.9f %15.9f"%(a,b)
      of.close()

  def show_summary(self, log=None, fsc_file_prefix="fsc_curve"):
    if(log is None): log = sys.stdout
    r = self.get_results()
    print >> log, "d99                    : ", r.d99
    print >> log, "d99_1                  : ", r.d99_1
    print >> log, "d99_2                  : ", r.d99_2
    print >> log, "d_model                : ", r.d_model
    print >> log, "b_iso_overall          : ", r.b_iso_overall
    print >> log, "d_fsc                  : ", r.d_fsc
    print >> log, "d_fsc_model (FSC=0.5)  : ", r.d_fsc_model
    print >> log, "d_fsc_model (FSC=0.143): ", r.d_fsc_model_0143
    print >> log, "d_fsc_model (FSC=0)    : ", r.d_fsc_model_0
    print >> log, "CC(half_map1,half_map2 : ", r.map_histograms.half_map_histogram_cc
    #
    of = open("%s_model"%fsc_file_prefix,"w")
    for a,b in zip(r.fsc_curve_model.fsc.d_inv, r.fsc_curve_model.fsc.fsc):
      print >>of, "%15.9f %15.9f"%(a,b)
    of.close()
    #
    if(r.fsc_curve is not None):
      of = open("%s"%fsc_file_prefix,"w")
      for a,b in zip(r.fsc_curve.fsc.d_inv, r.fsc_curve.fsc.fsc):
        print >>of, "%15.9f %15.9f"%(a,b)
      of.close()

  def get_results(self, slim=False):
    mask = None
    if(not slim):
      if(self.mask_object is not None):
        mask = self.mask_object.mask_smooth
    return group_args(
      d9                = self.d9,
      d99               = self.d99,
      d999              = self.d999,
      d99_1             = self.d99_1,
      d99_2             = self.d99_2,
      d_model           = self.d_model,
      d_model_b0        = self.d_model_b0,
      b_iso_overall     = self.b_iso_overall,
      d_fsc             = self.d_fsc,
      d_fsc_model       = self.d_fsc_model,
      d_fsc_model_0     = self.d_fsc_model_0,
      d_fsc_model_0143  = self.d_fsc_model_0143,
      fsc_curve         = self.fsc_curve,
      fsc_curve_model   = self.fsc_curve_model,
      crystal_symmetry  = self.crystal_symmetry,
      map_counts        = self.map_counts,
      half_map_1_counts = self.half_map_1_counts,
      half_map_2_counts = self.half_map_2_counts,
      map_histograms    = self.map_histograms,
      mask              = mask,
      radius_smooth     = self.radius_smooth)

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
