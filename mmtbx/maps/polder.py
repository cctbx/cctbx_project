from __future__ import absolute_import, division, print_function
import mmtbx.f_model
import mmtbx.utils
import mmtbx.masks
from mmtbx import map_tools
from iotbx import phil
from cctbx import maptbx
from cctbx import miller
from cctbx.array_family import flex
from libtbx import group_args
from libtbx.utils import Sorry
from libtbx.math_utils import ifloor, iceil
from six.moves import zip
from six.moves import range


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
  box_buffer = None
    .type = float
    .short_caption = Selection box buffer
    .help = Buffer around selection box: Increase the box for resetting the mask \
     by a buffer.
  compute_box = False
    .type = bool
    .short_caption = Use box
    .help = Reset mask within a box (parallel to unit cell axes) defined by an \
     atom selection
}
"""

def master_params():
  return phil.parse(master_params_str, process_includes = False)

# =============================================================================

class compute_polder_map():
  def __init__(self,
               f_obs,
               r_free_flags,
               model,
               params,
               selection_bool):
    self.f_obs = f_obs
    self.r_free_flags = r_free_flags
    self.xray_structure = model.get_xray_structure()
    self.pdb_hierarchy = model.get_hierarchy()
    self.params = params
    self.selection_bool = selection_bool
    #
    self.resolution_factor = self.params.resolution_factor
    self.sphere_radius = self.params.sphere_radius
    self.validation_results = None
    self.fmodel_biased, self.mc_biased = None, None

  # ---------------------------------------------------------------------------

  def validate(self):
    assert not None in [self.f_obs, self.xray_structure, self.pdb_hierarchy,
      self.params, self.selection_bool]

  # ---------------------------------------------------------------------------

  def run(self):
    # When extracting cartesian coordinates, xray_structure needs to be in P1:
    sites_cart_ligand_expanded = self.xray_structure.select(
      self.selection_bool).expand_to_p1(sites_mod_positive = True).sites_cart()
    sites_frac_ligand_expanded = self.xray_structure.select(
      self.selection_bool).expand_to_p1(sites_mod_positive = False).sites_frac()
    # xray_structure object without ligand/selection
    if (self.params.compute_box):
      self.xray_structure_noligand = self.xray_structure
    else:
      self.xray_structure_noligand = self.xray_structure.select(~self.selection_bool)
    self.crystal_gridding = self.f_obs.crystal_gridding(
      d_min             = self.f_obs.d_min(),
      symmetry_flags    = maptbx.use_space_group_symmetry,
      resolution_factor = self.resolution_factor)
    n_real = self.crystal_gridding.n_real()
    # Mask using all atoms
    self.mask_data_all = self.mask_from_xrs_unpadded(
      xray_structure = self.xray_structure,
      n_real         = n_real)
    # Mask if ligand is not in model
    self.mask_data_omit = self.mask_from_xrs_unpadded(
      xray_structure = self.xray_structure_noligand,
      n_real         = n_real)
    # Polder mask
    if (self.params.compute_box):
      # TODO: check if mask_omit = mask_all
      self.mask_data_polder = self.modify_mask_box(
        mask_data  = self.mask_data_all.deep_copy(),
        sites_frac =  sites_frac_ligand_expanded)
    else:
      self.mask_data_polder = self.modify_mask(
        mask_data     = self.mask_data_all.deep_copy(),
        sites_cart    = sites_cart_ligand_expanded)
    # Compute fmodel and map coeffs for input, biased, polder, omit case
    # Input model
    self.fmodel_input = mmtbx.f_model.manager(
     f_obs          = self.f_obs,
     r_free_flags   = self.r_free_flags,
     xray_structure = self.xray_structure)
    self.fmodel_input.update_all_scales()
    # Biased map
    if (not self.params.compute_box):
      self.fmodel_biased, self.mc_biased = self.get_fmodel_and_map_coefficients(
          xray_structure = self.xray_structure_noligand,
          mask_data      = self.mask_data_all)
    # Polder map
    self.fmodel_polder, self.mc_polder = self.get_fmodel_and_map_coefficients(
      xray_structure = self.xray_structure_noligand,
      mask_data      = self.mask_data_polder)
    # OMIT map
    self.fmodel_omit, self.mc_omit = self.get_fmodel_and_map_coefficients(
      xray_structure = self.xray_structure_noligand,
      mask_data      = self.mask_data_omit)
    # Validation only applies if selection present in model:
    if not(self.params.compute_box):
      self.validate_polder_map()

  # ---------------------------------------------------------------------------

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

  # ---------------------------------------------------------------------------

  def modify_mask_box(self, mask_data, sites_frac):
    box_buffer = self.params.box_buffer
    # Number of selected atoms
    n_selected = self.selection_bool.count(True)
    na = mask_data.all()
    n_selected_p1 = sites_frac.size()
    n_boxes = int(n_selected_p1/n_selected)
    box_list = [[] for i in range(n_boxes)]
    for n_box in range(n_boxes):
      for i in range(n_selected):
        box_list[n_box].append(sites_frac[n_box + n_boxes*i])
    na = self.mask_data_all.all()
    k = 0
    for box in box_list:
      k+=1
      x_min = min(frac[0] for frac in box)
      y_min = min(frac[1] for frac in box)
      z_min = min(frac[2] for frac in box)
      x_max = max(frac[0] for frac in box)
      y_max = max(frac[1] for frac in box)
      z_max = max(frac[2] for frac in box)
      frac_min = [x_min, y_min, z_min]
      frac_max = [x_max, y_max, z_max]

      cs = self.xray_structure.crystal_symmetry()

      # Add buffer to box if indicated.
      if (box_buffer is not None):
        cushion = flex.double(cs.unit_cell().fractionalize((box_buffer,)*3))
        frac_min = list(flex.double(frac_min) - cushion)
        frac_max = list(flex.double(frac_max) + cushion)

      gridding_first = [ifloor(f * n) for f,n in zip(frac_min, na)]
      gridding_last  = [iceil(f * n) for f,n in zip(frac_max, na)]

      for j in range(3):
        if (gridding_last[j] - gridding_first[j] >= na[j]):
          raise Sorry("The box is too big. Decrease box_buffer or use a " +
                      "different selection")

      maptbx.set_box(
        value         = 0,
        map_data_to   = mask_data,
        start         = gridding_first,
        end           = gridding_last)
    return mask_data

  # ---------------------------------------------------------------------------

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

  # ---------------------------------------------------------------------------

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

  # ---------------------------------------------------------------------------

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

  # ---------------------------------------------------------------------------

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
      validation_results = self.validation_results)

  # ---------------------------------------------------------------------------

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
    box_1 = self.get_polder_diff_map(
      f_obs = f_obs_1,
      r_free_flags = fmodel.r_free_flags(),
      f_calc = f_calc,
      f_mask = f_mask,
      xrs_selected = xrs_selected)
    box_2 = self.get_polder_diff_map(
      f_obs = f_obs_2,
      r_free_flags = fmodel.r_free_flags(),
      f_calc = f_calc,
      f_mask = f_mask,
      xrs_selected = xrs_selected)
    box_3 = self.get_polder_diff_map(
      f_obs = fmodel.f_obs(),
      r_free_flags = fmodel.r_free_flags(),
      f_calc = f_calc,
      f_mask = f_mask,
      xrs_selected = xrs_selected)

    sites_cart_box = box_1.xray_structure_box.sites_cart()
    sel = maptbx.grid_indices_around_sites(
      unit_cell  = box_1.xray_structure_box.unit_cell(),
      fft_n_real = box_1.map_box.focus(),
      fft_m_real = box_1.map_box.all(),
      sites_cart = sites_cart_box,
      site_radii = flex.double(sites_cart_box.size(), 2.0))
    b1 = box_1.map_box.select(sel).as_1d()
    b2 = box_2.map_box.select(sel).as_1d()
    b3 = box_3.map_box.select(sel).as_1d()
    # Map 1: calculated Fobs with ligand
    # Map 2: calculated Fobs without ligand
    # Map 3: real Fobs data
    cc12 = flex.linear_correlation(x=b1,y=b2).coefficient()
    cc13 = flex.linear_correlation(x=b1,y=b3).coefficient()
    cc23 = flex.linear_correlation(x=b2,y=b3).coefficient()
    #### D-function
    b1 = maptbx.volume_scale_1d(map=b1, n_bins=10000).map_data()
    b2 = maptbx.volume_scale_1d(map=b2, n_bins=10000).map_data()
    b3 = maptbx.volume_scale_1d(map=b3, n_bins=10000).map_data()
    cc12_peak = flex.linear_correlation(x=b1,y=b2).coefficient()
    cc13_peak = flex.linear_correlation(x=b1,y=b3).coefficient()
    cc23_peak = flex.linear_correlation(x=b2,y=b3).coefficient()
    #### Peak CC:
    cutoffs = flex.double(
      [i/10. for i in range(1,10)]+[i/100 for i in range(91,100)])
    d12 = maptbx.discrepancy_function(map_1=b1, map_2=b2, cutoffs=cutoffs)
    d13 = maptbx.discrepancy_function(map_1=b1, map_2=b3, cutoffs=cutoffs)
    d23 = maptbx.discrepancy_function(map_1=b2, map_2=b3, cutoffs=cutoffs)
    pdb_hierarchy_selected.adopt_xray_structure(box_1.xray_structure_box)
    self.validation_results = group_args(
      box_1 = box_1,
      box_2 = box_2,
      box_3 = box_3,
      cc12  = cc12,
      cc13  = cc13,
      cc23  = cc23,
      cc12_peak = cc12_peak,
      cc13_peak = cc13_peak,
      cc23_peak = cc23_peak,
      d12 = d12,
      d13 = d13,
      d23 = d23,
      cutoffs = cutoffs,
      ph_selected = pdb_hierarchy_selected
      )

# =============================================================================

#if (__name__ == "__main__"):
#  run(args=sys.argv[1:])
