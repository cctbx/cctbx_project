from __future__ import division
from libtbx import adopt_init_args
from scitbx.array_family import flex
import mmtbx.map_tools
from cctbx import maptbx
from cctbx import miller

class from_map_and_xray_structure_or_fmodel(object):

  def __init__(self,
        xray_structure=None,
        fmodel=None,
        map_data=None,
        d_min=None):
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
    # get map_data defined
    if(self.fmodel is not None):
      mc = mmtbx.map_tools.electron_density_map(
        fmodel=self.fmodel).map_coefficients(
          map_type         = "2mFo-DFc",
          isotropize       = True,
          fill_missing     = False)
      crystal_gridding = self.fmodel.f_obs().crystal_gridding(
        d_min              = self.fmodel.f_obs().d_min(),
        resolution_factor  = 1./3)
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
      self.map_model = fft_map.real_map_unpadded()
    if(self.fmodel is not None):
      self.sites_cart = self.fmodel.xray_structure.sites_cart()
      self.sites_frac = self.fmodel.xray_structure.sites_frac()
      self.weights    = self.fmodel.xray_structure.atomic_weights()
    else:
      self.sites_cart = self.xray_structure.sites_cart()
      self.sites_frac = self.xray_structure.sites_frac()
      self.weights    = self.xray_structure.atomic_weights()

  def cc(self, selections=None, selection=None, atom_radius=2.0):
    assert [selections, selection].count(None) == 1
    def compute(sites_cart):
      sel = maptbx.grid_indices_around_sites(
        unit_cell  = self.xray_structure.unit_cell(),
        fft_n_real = self.map_data.focus(),
        fft_m_real = self.map_data.all(),
        sites_cart = sites_cart,
        site_radii = flex.double(sites_cart.size(),atom_radius))
      return flex.linear_correlation(
        x=self.map_data.select(sel).as_1d(),
        y=self.map_model.select(sel).as_1d()).coefficient()
    if(selections is not None):
      result = []
      for s in selections:
        result.append(compute(sites_cart=self.sites_cart.select(s)))
      return result
    else:
      return compute(sites_cart=self.sites_cart.select(selection))

  def map_value(self, selections=None, selection=None):
    assert [selections, selection].count(None) == 1
    def compute(sites_frac, weights):
      result = 0
      for sf, w in zip(sites_frac, weights):
        result += self.map_data.eight_point_interpolation(sf)*w
      return result/flex.sum(weights)
    if(selections is not None):
      result = []
      for s in selections:
        result.append(compute(sites_frac=self.sites_frac.select(s),
          weights = self.weights.select(s)))
      return result
    else:
      return compute(sites_frac=self.sites_frac.select(selection),
        weights = self.weights.select(selection))
