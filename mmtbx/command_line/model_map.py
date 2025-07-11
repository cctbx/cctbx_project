"""Generate model map from coordinates"""
from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME phenix.model_map

import sys
import iotbx.pdb
from libtbx.utils import Sorry
from cctbx import maptbx
from mmtbx.maps import fem
import mmtbx.real_space
import mmtbx.utils

legend = """phenix.model_map: Given PDB file calculate model map

How to run:
  phenix.model_map model.pdb
"""

master_params_str = """
grid_step=0.3
  .type=float
output_file_name_prefix = None
  .type=str
scattering_table = *n_gaussian wk1995 it1992 electron neutron
  .type = choice(multi=False)
  .help = Choices of scattering table for structure factors calculations
"""

def master_params():
  return iotbx.phil.parse(master_params_str)

def run(args, log=sys.stdout):
  print("-"*79, file=log)
  print(legend, file=log)
  print("-"*79, file=log)
  inputs = mmtbx.utils.process_command_line_args(args = args,
    master_params = master_params())
  inputs.params.show(prefix="  ", out=log)
  print(file=log)
  file_names = inputs.pdb_file_names
  if(len(file_names) != 1): raise Sorry("A PDB file is expected.")
  pdb_inp = iotbx.pdb.input(file_name = file_names[0])
  awl = list(pdb_inp.atoms_with_labels())
  xrs = pdb_inp.xray_structure_simple().expand_to_p1(sites_mod_positive=True)
  # Check for B=0
  bs = xrs.extract_u_iso_or_u_equiv()
  sel_zero = bs<1.e-3
  n_zeros = sel_zero.count(True)
  if(n_zeros>0):
    print("Atoms with B=0:")
    for i_seq in sel_zero.iselection():
      print(awl[i_seq].format_atom_record())
    raise Sorry("Input model contains %d atoms with B=0"%n_zeros)
  #
  params = inputs.params.extract()
  mmtbx.utils.setup_scattering_dictionaries(
    scattering_table = params.scattering_table,
    xray_structure   = xrs,
    d_min            = 0.5)
  #
  crystal_gridding = maptbx.crystal_gridding(
    unit_cell        = xrs.unit_cell(),
    space_group_info = xrs.space_group_info(),
    symmetry_flags   = maptbx.use_space_group_symmetry,
    step             = params.grid_step)
  m = mmtbx.real_space.sampled_model_density(
    xray_structure = xrs,
    n_real         = crystal_gridding.n_real())
  map_data = m.data()
  #
  prefix = "model_map"
  if(params.output_file_name_prefix is not None):
    prefix = params.output_file_name_prefix
  #
  m.write_as_xplor_map(file_name = "%s.xplor"%prefix)
  fem.ccp4_map(cg=crystal_gridding, file_name="%s.ccp4"%prefix,
    map_data=map_data)

if (__name__ == "__main__"):
  run(args=sys.argv[1:])

