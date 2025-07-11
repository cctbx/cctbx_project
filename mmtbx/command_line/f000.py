"""Compute F000 using atomic model and bulk-solvent contribution"""
from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME phenix.f000

import sys
import iotbx.pdb
from libtbx.utils import Sorry
import mmtbx.utils

legend = """phenix.f000: Given PDB file estimate F(0,0,0)

How to run:
  phenix.f000 model.pdb
  phenix.f000 model.pdb mean_solvent_density=0.44
"""

master_params_str = """
mean_solvent_density=0.35
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
  if(len(file_names) != 1): raise Sorry("One PDB file is expected.")
  xrs = iotbx.pdb.input(file_name = file_names[0]).xray_structure_simple()
  msd = inputs.params.extract().mean_solvent_density
  f_000 = mmtbx.utils.f_000(xray_structure=xrs,
    mean_solvent_density=msd)
  f_000_str = str("%-13.3f"%f_000.f_000).strip()
  msd_str = str("%-6.3f"%msd).strip()
  sf_str = str("%-6.3f"%f_000.solvent_fraction).strip()
  msg = """
Estimate of F(0,0,0)=%s given mean bulk-solvent density %s and fraction %s
"""
  print(msg%(f_000_str, msd, sf_str), file=log)

if (__name__ == "__main__"):
  run(args=sys.argv[1:])

