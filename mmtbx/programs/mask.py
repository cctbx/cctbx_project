"""Given PDB file calculate bulk-solvent mask"""
from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME phenix.mask
from libtbx.program_template import ProgramTemplate
import iotbx.pdb
from cctbx import maptbx
from mmtbx.maps import fem
from mmtbx import masks

master_phil_str = '''
n_real = None
  .type=ints(3)
solvent_radius = None
  .type = float
invert = False
  .type = bool
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
    assert self.params.n_real is not None
    assert [type(_) is int for _ in self.params.n_real].count(True)==3

# ------------------------------------------------------------------------------

  def run(self):
    model = self.data_manager.get_model()
    xrs = model.get_xray_structure()
    crystal_gridding = maptbx.crystal_gridding(
      unit_cell             = xrs.unit_cell(),
      pre_determined_n_real = self.params.n_real,
      symmetry_flags        = maptbx.use_space_group_symmetry,
      space_group_info      = xrs.space_group().info()
      )
    mp = masks.mask_master_params.extract()
    mp.n_real = crystal_gridding.n_real()
    mp.step = None
    if self.params.solvent_radius is not None:
      mp.solvent_radius = self.params.solvent_radius
    mmtbx_masks_asu_mask_obj = masks.asu_mask(
      xray_structure = xrs.expand_to_p1(sites_mod_positive=True),
      mask_params    = mp)
    bulk_solvent_mask = mmtbx_masks_asu_mask_obj.mask_data_whole_uc()
    if self.params.invert:
      s = bulk_solvent_mask < 0.5
      bulk_solvent_mask.set_selected(s, 1)
      bulk_solvent_mask.set_selected(~s, 0)
    #
    fem.ccp4_map(
      cg        = crystal_gridding,
      file_name = "mask.mrc",
      map_data  = bulk_solvent_mask)
