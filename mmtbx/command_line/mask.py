from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME phenix.mask

import sys
import iotbx.pdb
from libtbx.utils import Sorry
from cctbx import maptbx
import mmtbx.utils
from mmtbx.maps import fem

legend = """phenix.mask: Given PDB file calculate bulk-solvent mask

How to run:
  phenix.mask model.pdb
"""

master_params_str = """
resolution=1.0
  .type=float
resolution_factor=1./4
  .type=float
"""

def master_params():
  return iotbx.phil.parse(master_params_str)

def run(args, log=sys.stdout):
  print("-"*79, file=log)
  print(legend, file=log)
  print("-"*79, file=log)
  inputs = mmtbx.utils.process_command_line_args(args = args,
    master_params = master_params())
  file_names = inputs.pdb_file_names
  if(len(file_names) != 1): raise Sorry("A PDB file is expected.")
  xrs = iotbx.pdb.input(file_name = file_names[0]).xray_structure_simple()
  params = inputs.params.extract()
  #
  crystal_gridding = maptbx.crystal_gridding(
    unit_cell         = xrs.unit_cell(),
    d_min             = params.resolution,
    resolution_factor = params.resolution_factor,
    symmetry_flags    = maptbx.use_space_group_symmetry,
    space_group_info  = xrs.space_group().info())
  mmtbx_masks_asu_mask_obj = mmtbx.masks.asu_mask(
    xray_structure = xrs.expand_to_p1(sites_mod_positive=True),
    n_real         = crystal_gridding.n_real())
  bulk_solvent_mask = mmtbx_masks_asu_mask_obj.mask_data_whole_uc()
  #
  fem.ccp4_map(cg=crystal_gridding, file_name="mask.ccp4",
    map_data=bulk_solvent_mask)

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
