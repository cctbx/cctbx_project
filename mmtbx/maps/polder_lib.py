from __future__ import division
#import time
#import os, sys
import sys
import mmtbx.f_model
import mmtbx.utils
import mmtbx.masks
import iotbx.pdb
from mmtbx import map_tools
from iotbx import ccp4_map
#from iotbx import file_reader
#from iotbx import phil
#from iotbx import reflection_file_utils
#from iotbx import crystal_symmetry_from_any
#from cStringIO import StringIO
from cctbx import maptbx
from cctbx import miller
from cctbx.array_family import flex
from libtbx import group_args
#from libtbx.utils import Sorry
#from libtbx.utils import multi_out
#from libtbx.math_utils import ifloor, iceil


master_params_str = """
polder {
  resolution_factor = 0.25
    .type = float
    .short_caption = Resolution factor
    .help = Used to determine the grid step = resolution_factor * high resolution
  sphere_radius = 5
    .type = float
    .short_caption = Solvent exclusion radius
    .help = Radius of sphere around atoms where solvent mask is reset to zero
}
"""

def master_params():
  return iotbx.phil.parse(master_params_str, process_includes = False)

class compute_polder_map():
  def __init__(self,
               f_obs,
               r_free_flags,
               xray_structure,
               pdb_hierarchy,
               params,
               selection_bool):
    self.f_obs = f_obs
    self.r_free_flags = r_free_flags
    self.xray_structure = xray_structure
    self.pdb_hierarchy = pdb_hierarchy
    self.params = params
    self.selection_bool = selection_bool
    #
    self.resolution_factor = self.params.resolution_factor
    self.sphere_radius = self.params.sphere_radius
    #

  def validate(self):
    assert not None in [self.f_obs, self.xray_structure, self.pdb_hierarchy,
      self.params, self.selection_bool]

  def run(self):
    # when extracting cartesian coordinates, xray_structure needs to be in P1:
    self.sites_cart_ligand = self.xray_structure.select(
      self.selection_bool).expand_to_p1(sites_mod_positive = True).sites_cart()
    # xray structure object without ligand/selection
    self.xray_structure_noligand = self.xray_structure.select(~self.selection_bool)
    self.crystal_gridding = self.f_obs.crystal_gridding(
      d_min             = self.f_obs.d_min(),
      symmetry_flags    = maptbx.use_space_group_symmetry,
      resolution_factor = self.resolution_factor)
    n_real = self.crystal_gridding.n_real()
    # mask using all atoms
    self.mask_data_all = self.mask_from_xrs_unpadded(
      xray_structure = self.xray_structure,
      n_real         = n_real)
    # mask if ligand is not in model
    self.mask_data_omit = self.mask_from_xrs_unpadded(
      xray_structure = self.xray_structure_noligand,
      n_real         = n_real)
    # polder mask
    self.mask_data_polder = self.modify_mask(
      mask_data     = self.mask_data_all.deep_copy(),
      sites_cart    = self.sites_cart_ligand)
    # compute fmodel and map coeffs for input, biased, polder, omit case
    self.fmodel_input = mmtbx.f_model.manager(
     f_obs          = self.f_obs,
     r_free_flags   = self.r_free_flags,
     xray_structure = self.xray_structure)
    self.fmodel_input.update_all_scales()
    self.fmodel_biased, self.mc_biased = self.get_fmodel_and_map_coefficients(
        xray_structure = self.xray_structure_noligand,
        mask_data      = self.mask_data_all)
    self.fmodel_polder, self.mc_polder = self.get_fmodel_and_map_coefficients(
      xray_structure = self.xray_structure_noligand,
      mask_data      = self.mask_data_polder)
    self.fmodel_omit, self.mc_omit = self.get_fmodel_and_map_coefficients(
      xray_structure = self.xray_structure_noligand,
      mask_data      = self.mask_data_omit)
    self.validate_polder_map()

  def modify_mask(self, mask_data, sites_cart):
    sel = maptbx.grid_indices_around_sites(
      unit_cell  = self.f_obs.crystal_symmetry().unit_cell(),
      fft_n_real = mask_data.focus(),
      fft_m_real = mask_data.all(),
      sites_cart = sites_cart,
      site_radii = flex.double(sites_cart.size(), self.sphere_radius))
    mask = mask_data.as_1d()
    mask.set_selected(sel, 0)
    mask.reshape(mask_data.accessor())
    return mask

  def mask_from_xrs_unpadded(self, xray_structure, n_real):
    mask_params = mmtbx.masks.mask_master_params.extract()
    mask = mmtbx.masks.mask_from_xray_structure(
      xray_structure           = xray_structure,
      p1                       = True,
      shrink_truncation_radius = mask_params.shrink_truncation_radius,
      solvent_radius           = mask_params.solvent_radius,
      for_structure_factors    = True,
      n_real                   = n_real).mask_data
    maptbx.unpad_in_place(map = mask)
    return mask

  def get_fmodel_and_map_coefficients(self, xray_structure, mask_data):
    f_calc = self.f_obs.structure_factors_from_scatterers(
      xray_structure = xray_structure).f_calc()
    mask = self.f_obs.structure_factors_from_map(
      map            = mask_data,
      use_scale      = True,
      anomalous_flag = False,
      use_sg         = False)
    # To check: is it really use_sg = false?
    fmodel = mmtbx.f_model.manager(
      f_obs        = self.f_obs,
      r_free_flags = self.r_free_flags,
      f_calc       = f_calc,
      f_mask       = mask)
    fmodel.update_all_scales()
    mc_fofc = map_tools.electron_density_map(fmodel = fmodel).map_coefficients(
      map_type     = "mFo-DFc",
      isotropize   = True,
      fill_missing = False)
    return fmodel, mc_fofc

  def get_polder_diff_map(self, f_obs, r_free_flags, f_calc, f_mask, xrs_selected):
    fmodel = mmtbx.f_model.manager(
      f_obs        = f_obs,
      r_free_flags = r_free_flags,
      f_calc       = f_calc,
      f_mask       = f_mask)
    fmodel.update_all_scales(remove_outliers=False)
    mc_diff = map_tools.electron_density_map(
      fmodel = fmodel).map_coefficients(
        map_type         = "mFo-DFc",
        isotropize       = True,
        fill_missing     = False)
    fft_map = miller.fft_map(
      crystal_gridding     = self.crystal_gridding,
      fourier_coefficients = mc_diff)
    fft_map.apply_sigma_scaling()
    map_data = fft_map.real_map_unpadded()
    return mmtbx.utils.extract_box_around_model_and_map(
      xray_structure = xrs_selected,
      map_data       = map_data,
      box_cushion    = 2.1)

  def get_results(self):
    return group_args(
      fmodel_input     = self.fmodel_input,
      fmodel_biased    = self.fmodel_biased,
      fmodel_omit      = self.fmodel_omit,
      fmodel_polder    = self.fmodel_polder,
      mask_data_all    = self.mask_data_all,
      mask_data_omit   = self.mask_data_omit,
      mask_data_polder = self.mask_data_polder,
      mc_polder        = self.mc_polder,
      mc_biased        = self.mc_biased,
      mc_omit          = self.mc_omit,
      box_1            = self.box_1,
      box_2            = self.box_2,
      box_3            = self.box_3)

  def validate_polder_map(self):
  # Significance check
    fmodel = mmtbx.f_model.manager(
     f_obs          = self.f_obs,
     r_free_flags   = self.r_free_flags,
     xray_structure = self.xray_structure)
    fmodel.update_all_scales(
      remove_outliers = False,
      fast            = True)
    f_obs_1 = abs(fmodel.f_model())
    fmodel.update_xray_structure(
      xray_structure      = self.xray_structure_noligand,
      update_f_calc       = True,
      update_f_mask       = True,
      force_update_f_mask = True)
  ## PVA: do we need it? fmodel.update_all_scales(remove_outliers=False)
    f_obs_2 = abs(fmodel.f_model())
    pdb_hierarchy_selected = self.pdb_hierarchy.select(self.selection_bool)
    xrs_selected = pdb_hierarchy_selected.extract_xray_structure(
      crystal_symmetry = self.f_obs.crystal_symmetry())
    f_calc = fmodel.f_obs().structure_factors_from_scatterers(
      xray_structure = self.xray_structure_noligand).f_calc()
    f_mask = fmodel.f_obs().structure_factors_from_map(
      map            = self.mask_data_polder,
      use_scale      = True,
      anomalous_flag = False,
      use_sg         = False)
    self.box_1 = self.get_polder_diff_map(
      f_obs = f_obs_1,
      r_free_flags = fmodel.r_free_flags(),
      f_calc = f_calc,
      f_mask = f_mask,
      xrs_selected = xrs_selected)
    self.box_2 = self.get_polder_diff_map(
      f_obs = f_obs_2,
      r_free_flags = fmodel.r_free_flags(),
      f_calc = f_calc,
      f_mask = f_mask,
      xrs_selected = xrs_selected)
    self.box_3 = self.get_polder_diff_map(
      f_obs = fmodel.f_obs(),
      r_free_flags = fmodel.r_free_flags(),
      f_calc = f_calc,
      f_mask = f_mask,
      xrs_selected = xrs_selected)

#-----------------------------------------------

#if (__name__ == "__main__"):
#  run(args=sys.argv[1:])

