from __future__ import absolute_import, division, print_function
import mmtbx.f_model
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
  altloc_scale = False
    .type = bool
    .short_caption = apply scale for altlocs
    .help = Apply special scaling procedure for alternate conformations \
     (only one residue, no mix)
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
               selection_string):
    self.f_obs = f_obs
    self.r_free_flags = r_free_flags
    self.xray_structure = model.get_xray_structure()
    self.pdb_hierarchy = model.get_hierarchy()
    self.model = model
    self.params = params
    self.selection_string = selection_string
    self.cs = self.xray_structure.crystal_symmetry()
    #
    self.resolution_factor = self.params.resolution_factor
    self.sphere_radius = self.params.sphere_radius
    # New scaling method for alternative conformations
    self.apply_scale_for_altloc = False

  # ---------------------------------------------------------------------------

  def validate(self):
    assert not None in [self.f_obs, self.xray_structure, self.pdb_hierarchy,
      self.params, self.selection_string]

  # ---------------------------------------------------------------------------

  def run(self):
    selection_bool = self.pdb_hierarchy.atom_selection_cache().selection(
      string = self.selection_string)
    ph_selected = self.pdb_hierarchy.select(selection_bool)
    altloc_indices = ph_selected.altloc_indices()
    if (len(altloc_indices) == 2 and '' not in list(altloc_indices) and
      not self.params.compute_box and self.params.altloc_scale):
      self.apply_scale_for_altloc = True
    else:
      self.apply_scale_for_altloc = False

    computed_results = self.compute(selection_bool = selection_bool)
    if not self.apply_scale_for_altloc:
      self.computed_results = computed_results
    else:
      # to be changed with modified code but use this to make it run
      self.computed_results = self.compute_if_altloc(
                                altloc_indices      = altloc_indices,
                                computed_results_AB = computed_results)

  # ---------------------------------------------------------------------------

  def compute_if_altloc(self, altloc_indices, computed_results_AB):
    #
    def combine(mA, mB, sc):
      #m = (mA+mB*sc)/2
      m = maptbx.combine_and_maximize_maps(map_data_1=mA, map_data_2=mB*sc, n_real=mB.all())
      return m/m.sample_standard_deviation()
    def scale(x, y):
      x = flex.abs(x)
      y = flex.abs(y)
      d = flex.sum(y*y)
      return flex.sum(x*y)/d
    def get_map(mc):
      fft_map = mc.fft_map(resolution_factor=0.25)
      #fft_map.apply_volume_scaling()
      fft_map.apply_sigma_scaling()
      m = fft_map.real_map_unpadded()
      return m
    # Not sure if this definition of sel strings can be prone to errors...
    sel_A_string = self.selection_string + ' altloc %s' % list(altloc_indices)[0]
    sel_B_string = self.selection_string + ' altloc %s' % list(altloc_indices)[1]
    sel_AB_string = self.selection_string
    sel_A = self.model.selection(string=sel_A_string)
    sel_B = self.model.selection(string=sel_B_string)
    sel_AB = self.model.selection(string=sel_AB_string)
    sf_A = self.model.get_sites_frac().select(sel_A)
    sf_B = self.model.get_sites_frac().select(sel_B)
    sf_AB = self.model.get_sites_frac().select(sel_AB)
    #
    computed_results_A = self.compute(selection_bool = sel_A)
    computed_results_B = self.compute(selection_bool = sel_B)
    mc_A_polder = computed_results_A.mc_polder
    mc_B_polder = computed_results_B.mc_polder
    mc_AB_polder = computed_results_AB.mc_polder
    m_A  = get_map(mc_A_polder)
    m_B  = get_map(mc_B_polder)
    m_AB = get_map(mc_AB_polder)
    mv_A = flex.double()
    mv_B = flex.double()
    for sf_A_, sf_B_ in zip(sf_A, sf_B):
      mv_A.append(m_A.tricubic_interpolation(sf_A_))
      mv_B.append(m_B.tricubic_interpolation(sf_B_))
    sc = scale(mv_A, mv_B)
    m = combine(m_A, m_B, sc)
    ###
    mv_A = flex.double()
    mv_B = flex.double()
    for sf_ in sf_AB:
      mv_A.append(m.tricubic_interpolation(sf_))
      mv_B.append(m_AB.tricubic_interpolation(sf_))
    sc = scale(mv_A, mv_B)
    m = combine(m, m_AB, sc)
    #
    mc = maptbx.map_to_map_coefficients(m=m, cs=self.cs, d_min=mc_A_polder.d_min())
#    mtz_dataset = mc.as_mtz_dataset(column_root_label='polder')
#    mtz_object = mtz_dataset.mtz_object()
#    mtz_object.write(file_name = "%s.mtz"%'polder_cl')
    #
    computed_results = computed_results_AB
    computed_results.mc_polder = mc

    return computed_results

  # ---------------------------------------------------------------------------

  def compute(self, selection_bool):
    # When extracting cartesian coordinates, xray_structure needs to be in P1:
    sites_cart_ligand_expanded = self.xray_structure.select(
      selection_bool).expand_to_p1(sites_mod_positive = True).sites_cart()
    sites_frac_ligand_expanded = self.xray_structure.select(
      selection_bool).expand_to_p1(sites_mod_positive = False).sites_frac()
    # xray_structure object without ligand/selection
    if (self.params.compute_box):
      xray_structure_noligand = self.xray_structure
    else:
      xray_structure_noligand = self.xray_structure.select(~selection_bool)
    self.crystal_gridding = self.f_obs.crystal_gridding(
      d_min             = self.f_obs.d_min(),
      symmetry_flags    = maptbx.use_space_group_symmetry,
      resolution_factor = self.resolution_factor)
    n_real = self.crystal_gridding.n_real()
    # Mask using all atoms
    mask_data_all = self.mask_from_xrs_unpadded(
      xray_structure = self.xray_structure,
      n_real         = n_real)
    # Mask if ligand is not in model
    mask_data_omit = self.mask_from_xrs_unpadded(
      xray_structure = xray_structure_noligand,
      n_real         = n_real)
    # Polder mask
    if (self.params.compute_box):
      # TODO: check if mask_omit = mask_all
      mask_data_polder = self.modify_mask_box(
        mask_data  = mask_data_all.deep_copy(),
        sites_frac =  sites_frac_ligand_expanded,
        selection_bool = selection_bool)
    else:
      mask_data_polder = self.modify_mask(
        mask_data     = mask_data_all.deep_copy(),
        sites_cart    = sites_cart_ligand_expanded)
    # Compute fmodel and map coeffs for input, biased, polder, omit case
    # Input model
    fmodel_input = mmtbx.f_model.manager(
     f_obs          = self.f_obs,
     r_free_flags   = self.r_free_flags,
     xray_structure = self.xray_structure)
    fmodel_input.update_all_scales()
    # Biased map
    if self.params.compute_box:
      fmodel_biased, mc_biased = None, None
    else:
      fmodel_biased, mc_biased = self.get_fmodel_and_map_coefficients(
          xray_structure = xray_structure_noligand,
          mask_data      = mask_data_all)
    # Polder map
    fmodel_polder, mc_polder = self.get_fmodel_and_map_coefficients(
      xray_structure = xray_structure_noligand,
      mask_data      = mask_data_polder)
    # OMIT map
    fmodel_omit, mc_omit = self.get_fmodel_and_map_coefficients(
      xray_structure = xray_structure_noligand,
      mask_data      = mask_data_omit)
    # Validation only applies if selection present in model:
    if (self.params.compute_box or self.apply_scale_for_altloc):
      validation_results = None
    else:
      validation_results = self.validate_polder_map(
        selection_bool = selection_bool,
        xray_structure_noligand = xray_structure_noligand,
        mask_data_polder = mask_data_polder)

    computed_results = group_args(
      fmodel_input     = fmodel_input,
      fmodel_biased    = fmodel_biased,
      fmodel_omit      = fmodel_omit,
      fmodel_polder    = fmodel_polder,
      mask_data_all    = mask_data_all,
      mask_data_omit   = mask_data_omit,
      mask_data_polder = mask_data_polder,
      mc_biased        = mc_biased,
      mc_polder        = mc_polder,
      mc_omit          = mc_omit,
      validation_results = validation_results)

    return computed_results

  # ---------------------------------------------------------------------------

  def get_results(self):
    return self.computed_results

  # ---------------------------------------------------------------------------

  def modify_mask(self, mask_data, sites_cart):
    sel = maptbx.grid_indices_around_sites(
      unit_cell  = self.cs.unit_cell(),
      fft_n_real = mask_data.focus(),
      fft_m_real = mask_data.all(),
      sites_cart = sites_cart,
      site_radii = flex.double(sites_cart.size(), self.sphere_radius))
    mask = mask_data.as_1d()
    mask.set_selected(sel, 0)
    mask.reshape(mask_data.accessor())
    return mask

  # ---------------------------------------------------------------------------

  def modify_mask_box(self, mask_data, sites_frac, selection_bool):
    box_buffer = self.params.box_buffer
    # Number of selected atoms
    n_selected = selection_bool.count(True)
    na = mask_data.all()
    n_selected_p1 = sites_frac.size()
    n_boxes = int(n_selected_p1/n_selected)
    box_list = [[] for i in range(n_boxes)]
    for n_box in range(n_boxes):
      for i in range(n_selected):
        box_list[n_box].append(sites_frac[n_box + n_boxes*i])
    #na = self.mask_data_all.all()
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

      #cs = self.xray_structure.crystal_symmetry()

      # Add buffer to box if indicated.
      if (box_buffer is not None):
        cushion = flex.double(self.cs.unit_cell().fractionalize((box_buffer,)*3))
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

  def get_polder_diff_map(self,
                          f_obs,
                          r_free_flags,
                          f_calc,
                          f_mask,
                          model_selected,
                          box_cushion):
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
    mm= fft_map.as_map_manager()

    from iotbx.map_model_manager import map_model_manager
    inputs=map_model_manager(
      model=model_selected.deep_copy(),
      map_manager=mm,) # no need to allow ignore_symmetry_conflicts
    return inputs

  # ---------------------------------------------------------------------------

  def validate_polder_map(self,
                          selection_bool,
                          xray_structure_noligand,
                          mask_data_polder,
                          box_cushion = 2.1):
    '''
    The parameter box_cushion is hardcoded to be 2.1
    The value is related to the site_radii used for CC calculation (box_cushion - 0.1)
    Ideally the site_radii are calculated according to resolution, atom type and B factor for each atom
    However, for the purpose of polder map validation, it is a reasonable approximation
    to use 2.0.
    If this value is changed, it will affect the values of the CCs and therefore also the
    output messages (see mmtbx/programs/polder.py --> result_message)
    So modify this value with caution.
    '''
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
      xray_structure      = xray_structure_noligand,
      update_f_calc       = True,
      update_f_mask       = True,
      force_update_f_mask = True)
  ## PVA: do we need it? fmodel.update_all_scales(remove_outliers=False)
    f_obs_2 = abs(fmodel.f_model())
    model_selected = self.model.select(selection_bool)
    pdb_hierarchy_selected = self.pdb_hierarchy.select(selection_bool)
    xrs_selected = pdb_hierarchy_selected.extract_xray_structure(
      crystal_symmetry = self.cs)
    f_calc = fmodel.f_obs().structure_factors_from_scatterers(
      xray_structure = xray_structure_noligand).f_calc()
    f_mask = fmodel.f_obs().structure_factors_from_map(
      map            = mask_data_polder,
      use_scale      = True,
      anomalous_flag = False,
      use_sg         = False)
    box_1 = self.get_polder_diff_map(
      f_obs = f_obs_1,
      r_free_flags = fmodel.r_free_flags(),
      f_calc = f_calc,
      f_mask = f_mask,
      model_selected = model_selected,
      box_cushion = box_cushion)
    box_2 = self.get_polder_diff_map(
      f_obs = f_obs_2,
      r_free_flags = fmodel.r_free_flags(),
      f_calc = f_calc,
      f_mask = f_mask,
      model_selected = model_selected,
      box_cushion = box_cushion)
    box_3 = self.get_polder_diff_map(
      f_obs = fmodel.f_obs(),
      r_free_flags = fmodel.r_free_flags(),
      f_calc = f_calc,
      f_mask = f_mask,
      model_selected = model_selected,
      box_cushion = box_cushion)

    sites_cart_box = box_1.model().get_xray_structure().sites_cart()
    sel = maptbx.grid_indices_around_sites(
      unit_cell  = box_1.model().get_xray_structure().unit_cell(),
      fft_n_real = box_1.map_manager().map_data().focus(),
      fft_m_real = box_1.map_manager().map_data().all(),
      sites_cart = sites_cart_box,
      site_radii = flex.double(sites_cart_box.size(), box_cushion-0.1))
    b1 = box_1.map_manager().map_data().select(sel).as_1d()
    b2 = box_2.map_manager().map_data().select(sel).as_1d()
    b3 = box_3.map_manager().map_data().select(sel).as_1d()
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
    pdb_hierarchy_selected.adopt_xray_structure(box_1.model().get_xray_structure())
    return group_args(
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
