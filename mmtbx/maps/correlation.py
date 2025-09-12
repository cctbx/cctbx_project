from __future__ import absolute_import, division, print_function
from libtbx import adopt_init_args
from scitbx.array_family import flex
from cctbx import maptbx
from cctbx import miller
from cctbx import adptbx
from mmtbx import masks
import sys
import boost_adaptbx.boost.python as bp
from six.moves import range
cctbx_maptbx_ext = bp.import_ext("cctbx_maptbx_ext")
from libtbx import group_args
from mmtbx.maps import map_tools

def get_selection_above_cutoff(m, n):
  return m>=m.as_1d().select(flex.sort_permutation(m.as_1d(), reverse=True))[n]

class five_cc(object):
  def __init__(self,
               map,
               xray_structure,
               d_min,
               map_calc=None,
               box=None,
               keep_map_calc=False,
               compute_cc_box=False,
               compute_cc_image=False,
               compute_cc_mask=True,
               compute_cc_volume=True,
               compute_cc_peaks=True):
    adopt_init_args(self, locals())
    if(box is None and
       xray_structure.crystal_symmetry().space_group().type().number()==1):
      box=True
    else:
      box=False
    #
    cc_box    = None
    cc_image  = None
    cc_mask   = None
    cc_volume = None
    cc_peaks  = None
    #
    self.crystal_gridding = maptbx.crystal_gridding(
      unit_cell             = xray_structure.unit_cell(),
      space_group_info      = xray_structure.space_group_info(),
      pre_determined_n_real = map.accessor().all(),
      symmetry_flags        = maptbx.use_space_group_symmetry)
    self.atom_radius = None
    bs_mask = masks.mask_from_xray_structure(
      xray_structure        = self.xray_structure,
      p1                    = True,
      for_structure_factors = False,
      n_real                = self.map.accessor().all()).mask_data
    maptbx.unpad_in_place(map=bs_mask)
    self.sel_inside = (bs_mask==0.).iselection()
    self.n_nodes_inside = self.sel_inside.size()
    del bs_mask
    #
    if map_calc is None:
      map_calc = self._get_map_calc()
    else:
      maptbx.assert_same_gridding(map_calc, self.map)
    #
    if(compute_cc_mask):
      cc_mask = from_map_map_selection(map_1=self.map, map_2=map_calc,
        selection = self.sel_inside)
    del self.sel_inside
    if(compute_cc_box):
      cc_box = from_map_map(map_1=self.map, map_2=map_calc)
    if(compute_cc_image):
      self.atom_radius = self._atom_radius()
      cc_image = self._cc_image(map_calc = map_calc)
    if(compute_cc_volume):
      cc_volume = self._cc_volume(map_calc=map_calc)
    if(compute_cc_peaks):
      cc_peaks = self._cc_peaks(map_calc=map_calc)
    #
    if(not keep_map_calc): map_calc=None
    self.result = group_args(
      cc_box      = cc_box,
      cc_image    = cc_image,
      cc_mask     = cc_mask,
      cc_volume   = cc_volume,
      cc_peaks    = cc_peaks,
      map_calc    = map_calc,
      atom_radius = self.atom_radius)

  def _get_map_calc(self):
    if(self.box is True):
      f_calc = miller.structure_factor_box_from_map(
        crystal_symmetry = self.xray_structure.crystal_symmetry(),
        n_real           = self.map.focus()).structure_factors_from_scatterers(
          xray_structure = self.xray_structure).f_calc()
      d_spacings = f_calc.d_spacings().data()
      sel = d_spacings > d_min
      f_calc = f_calc.select(sel)
    else:
      f_calc = self.xray_structure.structure_factors(d_min=self.d_min).f_calc()
    fft_map = miller.fft_map(
      crystal_gridding     = self.crystal_gridding,
      fourier_coefficients = f_calc)
    return fft_map.real_map_unpadded()

  def _atom_radius(self):
    b_iso = adptbx.u_as_b(
      flex.mean(self.xray_structure.extract_u_iso_or_u_equiv()))
    o = maptbx.atom_curves(scattering_type="C", scattering_table="electron")
    return o.image(d_min=self.d_min, b_iso=b_iso,
      radius_max=max(15.,self.d_min), radius_step=0.01).radius

  def _cc_image(self, map_calc):
    return from_map_map_atoms(
      map_1      = self.map,
      map_2      = map_calc,
      sites_cart = self.xray_structure.sites_cart(),
      unit_cell  = self.xray_structure.unit_cell(),
      radius     = self.atom_radius)

  def _cc_peaks(self, map_calc):
    s1 = get_selection_above_cutoff(m=self.map, n=self.n_nodes_inside)
    s2 = get_selection_above_cutoff(m=map_calc, n=self.n_nodes_inside)
    s = (s1 | s2).iselection()
    del s1
    del s2
    #G = flex.double(flex.grid(self.map.all()), 0)
    #G = G.set_selected(s, 1)
    #ccp4_map(cg=self.crystal_gridding, file_name="m1.ccp4", map_data=self.map)
    #ccp4_map(cg=self.crystal_gridding, file_name="m2.ccp4", map_data=map_calc)
    #ccp4_map(cg=self.crystal_gridding, file_name="m3.ccp4", map_data=G)
    return flex.linear_correlation(
      x=self.map.select(s).as_1d(),
      y=map_calc.select(s).as_1d()).coefficient()

  def _cc_volume(self, map_calc):
    s = get_selection_above_cutoff(m=map_calc, n=self.n_nodes_inside).iselection()
    #G = flex.double(flex.grid(self.map.all()), 0)
    #G = G.set_selected(s, 1)
    #ccp4_map(cg=self.crystal_gridding, file_name="m1.ccp4", map_data=self.map)
    #ccp4_map(cg=self.crystal_gridding, file_name="m2.ccp4", map_data=map_calc)
    #ccp4_map(cg=self.crystal_gridding, file_name="m3.ccp4", map_data=G)
    return flex.linear_correlation(
      x=self.map.select(s).as_1d(),
      y=map_calc.select(s).as_1d()).coefficient()

class fsc_model_vs_map(object):
  def __init__(self,
               xray_structure,
               map,
               atom_radius,
               d_min):
    assert atom_radius is not None # otherwise C++ traceback from flex.double(sites_frac.size(), atom_radius)
    self.atom_radius = atom_radius
    sgn = xray_structure.crystal_symmetry().space_group().type().number()
    assert sgn == 1 # P1 only
    def compute_mc(map_coeffs, map_data):
      return map_coeffs.structure_factors_from_map(
        map            = map_data,
        use_scale      = True,
        anomalous_flag = False,
        use_sg         = False)
    # Crystal gridding
    crystal_gridding = maptbx.crystal_gridding(
      unit_cell             = xray_structure.unit_cell(),
      space_group_info      = xray_structure.space_group_info(),
      pre_determined_n_real = map.accessor().all())
    # Compute mask
    sites_frac = xray_structure.sites_frac()
    mask = cctbx_maptbx_ext.mask(
      sites_frac                  = sites_frac,
      unit_cell                   = xray_structure.unit_cell(),
      n_real                      = map.all(),
      mask_value_inside_molecule  = 1,
      mask_value_outside_molecule = 0,
      radii                       = flex.double(sites_frac.size(), atom_radius))
    # Compute Fcalc
    f_calc = xray_structure.structure_factors(
      d_min=min(d_min, 3.0)).f_calc().resolution_filter(d_min=d_min)
    # Compute Fcalc masked inplace
    fft_map = miller.fft_map(
      crystal_gridding     = crystal_gridding,
      fourier_coefficients = f_calc)
    map_calc = fft_map.real_map_unpadded()
    del fft_map
    map_calc = map_calc * mask
    map = map*mask
    del mask
    f_calc = compute_mc(map_coeffs = f_calc, map_data = map_calc)
    del map_calc
    # Compute Fobs masked
    f_obs = compute_mc(map_coeffs = f_calc, map_data = map)
    # Binning
    n_bins=min(30, f_obs.data().size()//500)
    dsd = f_obs.d_spacings().data()
    if(dsd.size()>1500):
      f_obs.setup_binner(n_bins = n_bins)
    else:
      f_obs.setup_binner(reflections_per_bin = dsd.size())
    # Compute FSC
    self.result = []
    for i_bin in f_obs.binner().range_used():
      sel       = f_obs.binner().selection(i_bin)
      d         = dsd.select(sel)
      d_min     = flex.min(d)
      d_max     = flex.max(d)
      n         = d.size()
      fc        = f_calc.select(sel)
      fo        = f_obs.select(sel)
      cc        = fc.map_correlation(other = fo)
      bin = group_args(
        d_min = d_min, d_max = d_max, n = d.size(), cc = cc)
      self.result.append(bin)

  def show_lines(self, prefix=""):
    lines = []
    msg = prefix+"Atom radius used for masking: %s A"%(
      str("%5.2f"%self.atom_radius).strip())
    lines.append(msg)
    lines.append(prefix+"Bin  Resolution      CC(FSC)  N_coeffs")
    fmt="%2d: %7.3f-%-7.3f %7.4f %5d"
    for i, bin in enumerate(self.result):
      lines.append(prefix+fmt%(i, bin.d_max, bin.d_min, bin.cc, bin.n))
    return lines

  def show(self, log=None, prefix="", allcaps=False):
    if(log is None): log = sys.stdout
    lines = self.show_lines(prefix=prefix)
    for l in lines:
      if(allcaps): l = l.upper()
      print(l, file=log)

def assert_same_gridding(map_1, map_2):
  # XXX remove it!
  maptbx.assert_same_gridding(map_1, map_2)

def from_map_map(map_1, map_2):
  assert_same_gridding(map_1, map_2)
  return flex.linear_correlation(
    x=map_1.as_1d(),
    y=map_2.as_1d()).coefficient()

def from_map_map_atom(map_1, map_2, site_cart, unit_cell, radius):
  assert_same_gridding(map_1, map_2)
  sel = maptbx.grid_indices_around_sites(
    unit_cell  = unit_cell,
    fft_n_real = map_1.focus(),
    fft_m_real = map_1.all(),
    sites_cart = flex.vec3_double([site_cart]),
    site_radii = flex.double(1, radius))
  return flex.linear_correlation(
    x=map_1.select(sel).as_1d(),
    y=map_2.select(sel).as_1d()).coefficient()

def from_map_map_atoms_optimal_radius(map_1, map_2, sites_cart, unit_cell, d_min):
  radii_coarse = [r/100. for r in range(150,550,50)]+[d_min,]
  cc_best = -999
  r_best = None
  for r in radii_coarse:
    cc = from_map_map_atoms(map_1=map_1, map_2=map_2, sites_cart=sites_cart,
      unit_cell=unit_cell, radius=r)
    if(cc>cc_best):
      r_best = r
      cc_best = cc
  r_fine = [r/100. for r in range(max(150,int(r_best*100)-50), int(r_best*100)+50,10)]
  for r in r_fine:
    cc = from_map_map_atoms(map_1=map_1, map_2=map_2, sites_cart=sites_cart,
      unit_cell=unit_cell, radius=r)
    if(cc>cc_best):
      r_best = r
      cc_best = cc
  return cc_best, r_best

def from_map_map_atoms(map_1, map_2, sites_cart, unit_cell, radius):
  assert_same_gridding(map_1, map_2)
  sel = maptbx.grid_indices_around_sites(
    unit_cell  = unit_cell,
    fft_n_real = map_1.focus(),
    fft_m_real = map_1.all(),
    sites_cart = sites_cart,
    site_radii = flex.double(sites_cart.size(), radius))
  return flex.linear_correlation(
    x=map_1.select(sel).as_1d(),
    y=map_2.select(sel).as_1d()).coefficient()

def from_map_map_selection(map_1, map_2, selection):
  assert_same_gridding(map_1, map_2)
  return flex.linear_correlation(
    x=map_1.select(selection).as_1d(),
    y=map_2.select(selection).as_1d()).coefficient()

def from_map_map_atoms_per_atom(map_1, map_2, sites_cart, unit_cell, radius):
  assert_same_gridding(map_1, map_2)
  result = flex.double()
  for site_cart in sites_cart:
    cc = from_map_map_atom(map_1=map_1, map_2=map_2, site_cart=site_cart,
      unit_cell=unit_cell, radius=radius)
    result.append(cc)
  return result

class histogram_per_atom(object):
  def __init__(self, map_1, map_2, sites_cart, unit_cell, radius, n_slots):
    assert_same_gridding(map_1, map_2)
    self.ccs = from_map_map_atoms_per_atom(
      map_1      = map_1,
      map_2      = map_2,
      sites_cart = sites_cart,
      unit_cell  = unit_cell,
      radius     = radius)
    self.hist = flex.histogram(data = self.ccs, n_slots = n_slots)

  def format(self, prefix=""):
    h = self.hist
    lc_1 = h.data_min()
    s_1 = enumerate(h.slots())
    lines = []
    for (i_1,n_1) in s_1:
      hc_1 = h.data_min() + h.slot_width() * (i_1+1)
      line = "%s %7.4f - %-7.4f: %6.2f %s" % (prefix, lc_1, hc_1,
        n_1/self.ccs.size()*100, "%")
      lines.append(line)
      lc_1 = hc_1
    return "\n".join(lines)

class from_map_and_xray_structure_or_fmodel(object):

  def __init__(self,
        xray_structure    = None,
        fmodel            = None,
        map_data          = None,
        d_min             = None,
        resolution_factor = 0.25,
        map_type          = "2mFo-DFc"):
    """
    Utility to calculate correlation between two maps:
      CC(xray_structure, map_data), xray_structure are map_data inputs
    or
      CC(2mFo-DFc, Fc), 2mFo-DFc and Fc are from input fmodel .
    """
    assert [fmodel, map_data].count(None) == 1
    assert [xray_structure, map_data].count(None) in [0, 2]
    assert [fmodel, xray_structure].count(None) == 1
    assert [d_min, fmodel].count(None) == 1
    adopt_init_args(self, locals())
    if(fmodel is not None): self.xray_structure = fmodel.xray_structure
    # get map_data defined
    if(self.fmodel is not None):
      e_map_obj = fmodel.electron_density_map()
      isotropize = True
      if(fmodel.is_twin_fmodel_manager()): isotropize = False
      mc = e_map_obj.map_coefficients(
        map_type           = map_type,
        fill_missing       = False,
        isotropize         = isotropize)
      crystal_gridding = self.fmodel.f_obs().crystal_gridding(
        d_min              = self.fmodel.f_obs().d_min(),
        resolution_factor  = resolution_factor)
      fft_map = miller.fft_map(
        crystal_gridding     = crystal_gridding,
        fourier_coefficients = mc)
      self.map_data = fft_map.real_map_unpadded()
    # get model map
    if(self.fmodel is not None):
      if(fmodel.is_twin_fmodel_manager()):
        f_model = self.fmodel.f_model()
      else:
        f_model = self.fmodel.f_model_scaled_with_k1()
      fft_map = miller.fft_map(
        crystal_gridding     = crystal_gridding,
        fourier_coefficients = f_model)
      fft_map.apply_sigma_scaling()
      self.map_model = fft_map.real_map_unpadded()
    else:
      crystal_gridding = maptbx.crystal_gridding(
        unit_cell             = self.xray_structure.unit_cell(),
        space_group_info      = self.xray_structure.space_group_info(),
        pre_determined_n_real = self.map_data.accessor().all())
      f_model = self.xray_structure.structure_factors(d_min=self.d_min).f_calc()
      fft_map = miller.fft_map(
        crystal_gridding     = crystal_gridding,
        fourier_coefficients = f_model)
      del f_model
      fft_map.apply_sigma_scaling()
      self.map_model = fft_map.real_map_unpadded()
      del fft_map
    if(self.fmodel is not None):
      self.sites_cart = self.fmodel.xray_structure.sites_cart()
      self.sites_frac = self.fmodel.xray_structure.sites_frac()
    else:
      self.sites_cart = self.xray_structure.sites_cart()
      self.sites_frac = self.xray_structure.sites_frac()

  def cc(self, selections=None, selection=None, atom_radius=2.0, per_atom=None):
    def compute(sites_cart):
      return from_map_map_atoms(
        map_1      = self.map_data,
        map_2      = self.map_model,
        sites_cart = sites_cart,
        unit_cell  = self.xray_structure.unit_cell(),
        radius     = atom_radius)
    if(selections is not None):
      result = []
      for s in selections:
        result.append(compute(sites_cart=self.sites_cart.select(s)))
      return result
    elif(selection is not None):
      return compute(sites_cart=self.sites_cart.select(selection))
    elif(per_atom):
      return from_map_map_atoms_per_atom(
        map_1      = self.map_data,
        map_2      = self.map_model,
        sites_cart = self.sites_cart,
        unit_cell  = self.fmodel.xray_structure.unit_cell(),
        radius     = atom_radius)
    elif([selection, selections].count(None)==2):
      return from_map_map(
        map_1 = self.map_data,
        map_2 = self.map_model)
    else:
      raise RuntimeError

def map_model_cc_and_vals_per_atom_xtal(
      fmodel,
      map_obs_type  = "2mFo-DFc",
      map_calc_type = "Fmodel",
      grid_step     = 0.6, # From mosaic work
      atom_radius   = 2.0):
  """
  Calculate CC(map_obs_type, map_calc_type) and map values at atom center of
  map_obs_type. Crystallography specific since it uses fmodel.
  fmodel is expected to have all scales set and up to date.
  """
  def get_map(map_type):
    map_coeffs = map_tools.electron_density_map(
       fmodel = fmodel).map_coefficients(
         map_type     = map_type,
         isotropize   = True,
         fill_missing = False)
    fft_map = miller.fft_map(
      crystal_gridding     = crystal_gridding,
      fourier_coefficients = map_coeffs)
    fft_map.apply_sigma_scaling()
    return fft_map.real_map_unpadded()
  xrs = fmodel.xray_structure
  uc  = xrs.unit_cell()
  crystal_gridding = maptbx.crystal_gridding(
    unit_cell        = uc,
    space_group_info = xrs.space_group_info(),
    symmetry_flags   = maptbx.use_space_group_symmetry,
    step             = grid_step)
  map_obs  = get_map(map_type = map_obs_type)
  map_calc = get_map(map_type = map_calc_type)
  ccs  = flex.double()
  vals = flex.double()
  for site_cart, site_frac in zip(xrs.sites_cart(), xrs.sites_frac()):
    cc = from_map_map_atom(
      map_1     = map_obs,
      map_2     = map_calc,
      site_cart = site_cart,
      unit_cell = uc,
      radius    = atom_radius)
    mv = map_obs.tricubic_interpolation(site_frac)
    ccs .append(cc)
    vals.append(mv)
  return group_args(ccs = ccs, vals = vals)
