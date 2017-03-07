from __future__ import division
from libtbx import adopt_init_args
from scitbx.array_family import flex
from cctbx import maptbx
from cctbx import miller
from cctbx import adptbx
from mmtbx import masks
import sys
import boost.python
cctbx_maptbx_ext = boost.python.import_ext("cctbx_maptbx_ext")

class d99(object):
  def __init__(self, map, crystal_symmetry):
    # First thing first: shift origin if needed
    shift_needed = not \
      (map.focus_size_1d() > 0 and map.nd() == 3 and map.is_0_based())
    if(shift_needed):
      if(not crystal_symmetry.space_group().type().number() in [0,1]):
        raise RuntimeError("Not implemented")
    map = map.shift_origin()
    #
    n_real = map.focus()
    max_index = [(i-1)//2 for i in n_real]
    complete_set = miller.build_set(
      crystal_symmetry = crystal_symmetry,
      anomalous_flag   = False,
      max_index        = max_index)
    self.f = complete_set.structure_factors_from_map(
      map            = map,
      use_scale      = True,
      anomalous_flag = False,
      use_sg         = False)
    ds = self.f.d_spacings().data()
    d_max, d_min = flex.max(ds), flex.min(ds)
    o = maptbx.d99(f=self.f.data(), d_spacings=ds,
      hkl=self.f.indices(), d_min=d_min, d_max=d_max)
    self.d_min_cc9, self.d_min_cc99, self.d_min_cc999=\
      o.d_min_cc9(), o.d_min_cc99(), o.d_min_cc999()
    self.ccs=o.ccs()
    self.d_mins=o.d_mins()

  def write_mtz_file(self, file_name):
    mtz_dataset = self.f.as_mtz_dataset(column_root_label="F")
    mtz_object = mtz_dataset.mtz_object()
    mtz_object.write(file_name = file_name)

  def write_log(self, file_name):
    of=open(file_name,"w")
    fmt = "%12.6f %8.6f"
    for d_min, cc in zip(self.d_mins, self.ccs):
      print >> of, fmt%(d_min, cc)
    of.close()

def get_selection_above_cutoff(m, n):
  m_ = m.as_1d()
  s = flex.sort_permutation(m_, reverse=True)
  return m>=m_.select(s)[n]

class five_cc(object):
  def __init__(self, map, xray_structure, d_min):
    adopt_init_args(self, locals())
    self.crystal_gridding = maptbx.crystal_gridding(
      unit_cell             = xray_structure.unit_cell(),
      space_group_info      = xray_structure.space_group_info(),
      pre_determined_n_real = map.accessor().all(),
      symmetry_flags        = maptbx.use_space_group_symmetry)
    f_calc = xray_structure.structure_factors(d_min=d_min).f_calc()
    fft_map = miller.fft_map(
      crystal_gridding     = self.crystal_gridding,
      fourier_coefficients = f_calc)
    self.atom_radius = self._atom_radius()
    self.map_calc = fft_map.real_map_unpadded()
    self.bs_mask = masks.mask_from_xray_structure(
      xray_structure        = self.xray_structure,
      p1                    = True,
      for_structure_factors = False,
      n_real                = self.map.accessor().all()).mask_data
    maptbx.unpad_in_place(map=self.bs_mask)
    self.sel_inside = (self.bs_mask==0.).iselection()
    #
    self.cc_overall = from_map_map(map_1=self.map, map_2=self.map_calc)
    self.cc_image = self._cc_image()
    self.cc_mask = from_map_map_selection(map_1=self.map, map_2=self.map_calc,
      selection = self.sel_inside)
    self.cc_peaks = self._cc_peaks()
    self.cc_volume = self._cc_volume()
    # Free memory
    del self.bs_mask, self.sel_inside

  def _atom_radius(self):
    b_iso = adptbx.u_as_b(
      flex.mean(self.xray_structure.extract_u_iso_or_u_equiv()))
    o = maptbx.atom_curves(scattering_type="C", scattering_table="electron")
    return o.image(d_min=self.d_min, b_iso=b_iso,
      radius_max=max(10.,self.d_min), radius_step=0.01).radius

  def _cc_image(self):
    return from_map_map_atoms(
      map_1      = self.map,
      map_2      = self.map_calc,
      sites_cart = self.xray_structure.sites_cart(),
      unit_cell  = self.xray_structure.unit_cell(),
      radius     = self.atom_radius)

  def _cc_peaks(self):
    n_nodes_inside = self.sel_inside.size()
    s1 = get_selection_above_cutoff(m=self.map, n=n_nodes_inside)
    s2 = get_selection_above_cutoff(m=self.map_calc, n=n_nodes_inside)
    s = s1 | s2
    #G = flex.double(flex.grid(self.map.all()), 0)
    #G = G.set_selected(s, 1)
    #ccp4_map(cg=self.crystal_gridding, file_name="m1.ccp4", map_data=self.map)
    #ccp4_map(cg=self.crystal_gridding, file_name="m2.ccp4", map_data=self.map_calc)
    #ccp4_map(cg=self.crystal_gridding, file_name="m3.ccp4", map_data=G)
    return flex.linear_correlation(
      x=self.map.select(s.iselection()).as_1d(),
      y=self.map_calc.select(s.iselection()).as_1d()).coefficient()

  def _cc_volume(self):
    n_nodes_inside = self.sel_inside.size()
    s = get_selection_above_cutoff(m=self.map_calc, n=n_nodes_inside)
    #G = flex.double(flex.grid(self.map.all()), 0)
    #G = G.set_selected(s, 1)
    #ccp4_map(cg=self.crystal_gridding, file_name="m1.ccp4", map_data=self.map)
    #ccp4_map(cg=self.crystal_gridding, file_name="m2.ccp4", map_data=self.map_calc)
    #ccp4_map(cg=self.crystal_gridding, file_name="m3.ccp4", map_data=G)
    return flex.linear_correlation(
      x=self.map.select(s.iselection()).as_1d(),
      y=self.map_calc.select(s.iselection()).as_1d()).coefficient()

def fsc_model_map(xray_structure, map, d_min, log=sys.stdout, radius=2.,
                  prefix=""):
  # XXX make radius dynamic
  sgn = xray_structure.crystal_symmetry().space_group().type().number()
  f_calc = xray_structure.structure_factors(d_min=d_min).f_calc()
  n_bins=min(30, f_calc.data().size()//500)
  def compute_mc(f_calc, map):
    return f_calc.structure_factors_from_map(
      map            = map,
      use_scale      = True,
      anomalous_flag = False,
      use_sg         = False)
  sites_frac = xray_structure.sites_frac()
  if(sgn==1):
    mask = cctbx_maptbx_ext.mask(
      sites_frac                  = sites_frac,
      unit_cell                   = xray_structure.unit_cell(),
      n_real                      = map.all(),
      mask_value_inside_molecule  = 1,
      mask_value_outside_molecule = 0,
      radii                       = flex.double(sites_frac.size(), radius))
  mc = compute_mc(f_calc=f_calc, map=map)
  if(sgn==1):
    mc_masked = compute_mc(f_calc=f_calc, map=map*mask)
    del mask
  print >> log, prefix, "Overall (entire box):  %7.4f"%\
    f_calc.map_correlation(other = mc)
  if(sgn==1):
    cc = f_calc.map_correlation(other = mc_masked)
    if(cc is not None): print >> log, prefix, "Around atoms (masked): %7.4f"%cc
  dsd = f_calc.d_spacings().data()
  if(dsd.size()>1500):
    f_calc.setup_binner(n_bins = n_bins)
  else:
    f_calc.setup_binner(reflections_per_bin = dsd.size())
  if(sgn==1):
    print >> log, prefix, "Bin# Resolution (A)     CC   CC(masked)"
  else:
    print >> log, prefix, "Bin# Resolution (A)     CC"
  fmt1="%2d: %7.3f-%-7.3f %7.4f"
  for i_bin in f_calc.binner().range_used():
    sel       = f_calc.binner().selection(i_bin)
    d         = dsd.select(sel)
    d_min     = flex.min(d)
    d_max     = flex.max(d)
    n         = d.size()
    fc        = f_calc.select(sel)
    fo        = mc.select(sel)
    cc        = fc.map_correlation(other = fo)
    if(sgn==1):
      fo_masked = mc_masked.select(sel)
      cc_masked = fc.map_correlation(other = fo_masked)
      if(cc_masked is not None and cc is not None):
        fmt2="%2d: %7.3f-%-7.3f %7.4f %7.4f"
        print >> log, prefix, fmt2%(i_bin, d_max, d_min, cc, cc_masked)
      else:
        fmt2="%2d: %7.3f-%-7.3f %s %s"
        print >> log, prefix, fmt2%(i_bin, d_max, d_min, "none", "none")
    else:
      print >> log, prefix, fmt1%(i_bin, d_max, d_min, cc)

def assert_same_gridding(map_1, map_2):
  assert map_1.focus()==map_2.focus()
  assert map_1.origin()==map_2.origin()
  assert map_1.all()==map_2.all()

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
      fft_map.apply_sigma_scaling()
      self.map_model = fft_map.real_map_unpadded()
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
