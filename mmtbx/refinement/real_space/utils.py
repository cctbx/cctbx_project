from __future__ import division
from cctbx import miller
from cctbx import maptbx

class target_map(object):
  """
  Data structure for real-space refinement. Holds map and its reciprocal-space
  equivalent, and its masked version that can be updated as model changes.
  """
  def __init__(self, map_data, xray_structure, d_min, atom_radius):
    self.d_min = d_min
    self.map_data = map_data
    self.atom_radius = atom_radius
    self.f_map_diff = None # XXX rudimentary left-over, remove later
    self.crystal_gridding = maptbx.crystal_gridding( #XXX Likewise, remove later
      unit_cell             = xray_structure.unit_cell(),
      pre_determined_n_real = map_data.all(),
        space_group_info    = xray_structure.space_group_info())
    self.complete_set = miller.build_set(
      crystal_symmetry = xray_structure.crystal_symmetry(),
      anomalous_flag   = False,
      d_min            = d_min)
    self.miller_array = self.map_to_sf(map_data = self.map_data)
    self.miller_array_masked = self.update_miller_array_masked(
      xray_structure = xray_structure)

  def update_miller_array_masked(self, xray_structure):
    if(xray_structure.crystal_symmetry().space_group().type().number()!=1):
      return None
    mask_params = mmtbx.masks.mask_master_params.extract()
    mask_params.solvent_radius = 0
    mask_params.shrink_truncation_radius = 0
    mask_params.ignore_zero_occupancy_atoms = False
    mask_params.ignore_hydrogens = False
    mmtbx_masks_asu_mask_obj = mmtbx.masks.asu_mask(
      xray_structure = xray_structure,
      n_real         = self.map_data.all(),
      mask_params    = mask_params,
      atom_radius    = self.atom_radius)
    mask = mmtbx_masks_asu_mask_obj.mask_data_whole_uc()
    self.miller_array_masked = self.map_to_sf(map_data = self.map_data*mask)

  def map_to_sf(self, map_data):
    return self.complete_set.structure_factors_from_map(
      map            = map_data,
      use_scale      = True,
      anomalous_flag = False,
      use_sg         = True)
