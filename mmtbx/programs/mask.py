"""Given PDB file calculate bulk-solvent mask"""
from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME phenix.mask
from libtbx.program_template import ProgramTemplate
import iotbx.pdb
from cctbx import maptbx
from mmtbx.maps import fem
from mmtbx import masks

master_phil_str = '''
resolution=1.0
  .type=float
resolution_factor=1./4
  .type=float
'''

# ------------------------------------------------------------------------------

class Program(ProgramTemplate):
  description = '''
phenix.mask: Given PDB file calculate bulk-solvent mask

How to run:
  phenix.mask model.pdb
'''
  datatypes = ['model', 'phil']
  master_phil_str = master_phil_str

# ------------------------------------------------------------------------------

  def validate(self):
    self.data_manager.has_models(raise_sorry=True)

# ------------------------------------------------------------------------------

  def run(self):
    model = self.data_manager.get_model()
    xrs = model.get_xray_structure()
    crystal_gridding = maptbx.crystal_gridding(
      unit_cell         = xrs.unit_cell(),
      d_min             = self.params.resolution,
      resolution_factor = self.params.resolution_factor,
      symmetry_flags    = maptbx.use_space_group_symmetry,
      space_group_info  = xrs.space_group().info())
    mp = masks.mask_master_params.extract()
    mp.n_real = crystal_gridding.n_real()
    mp.step = None
    mmtbx_masks_asu_mask_obj = masks.asu_mask(
      xray_structure = xrs.expand_to_p1(sites_mod_positive=True),
      mask_params    = mp)
    bulk_solvent_mask = mmtbx_masks_asu_mask_obj.mask_data_whole_uc()
    #
    fem.ccp4_map(
      cg        = crystal_gridding,
      file_name = "mask.ccp4",
      map_data  = bulk_solvent_mask)
