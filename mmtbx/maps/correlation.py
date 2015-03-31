from __future__ import division
from libtbx import adopt_init_args
from scitbx.array_family import flex
import mmtbx.map_tools
from cctbx import maptbx
from cctbx import miller

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
      mc = mmtbx.map_tools.electron_density_map(
        fmodel=self.fmodel).map_coefficients(
          map_type         = map_type,
          isotropize       = True,
          fill_missing     = False)
      crystal_gridding = self.fmodel.f_obs().crystal_gridding(
        d_min              = self.fmodel.f_obs().d_min(),
        resolution_factor  = resolution_factor)
      fft_map = miller.fft_map(
        crystal_gridding     = crystal_gridding,
        fourier_coefficients = mc)
      self.map_data = fft_map.real_map_unpadded()
    # get model map
    if(self.fmodel is not None):
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
      self.weights    = self.fmodel.xray_structure.atomic_weights()
    else:
      self.sites_cart = self.xray_structure.sites_cart()
      self.sites_frac = self.xray_structure.sites_frac()
      self.weights    = self.xray_structure.atomic_weights()

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
