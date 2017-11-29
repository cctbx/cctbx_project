from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME mmtbx.fix_rotamer_outliers

import sys
import iotbx.pdb
import mmtbx.monomer_library.server
import iotbx.phil
import iotbx.pdb
from libtbx.utils import Sorry
from mmtbx.utils import fix_rotamer_outliers
from mmtbx.rotamer.rotamer_eval import RotamerEval
import mmtbx.model


master_phil_str = """
show_all_params = False
  .type = bool
  .style = hidden
radius = 5.
  .type = float
  .help = radius for determining surrounding atoms around the residue
  .style = hidden
output_prefix = rotamer_fixed
  .type = str
file_name = None
  .type = path
  .optional = True
  .style = hidden
"""

master_phil = iotbx.phil.parse(master_phil_str, process_includes=True)

def show_usage():
  help_msg = """\
phenix.fix_rotamer_outliers: tool for fixing rotamer outliers.
    For every rotamer outlier it will choose the rotamer with minimum
    clashes with surrounding atoms.

Usage examples:
  mmtbx.fix_rotamer_outliers model.pdb
  mmtbx.fix_rotamer_outliers model.pdb
  mmtbx.fix_rotamer_outliers model.pdb ligands.cif output_prefix=rot_fixed

Full scope of parameters:
  """
  print help_msg
  master_phil.show()

def run(args, params=None, out=sys.stdout, log=sys.stderr):
  if ( ((len(args) == 0) and (params is None)) or
       ((len(args) > 0) and ((args[0] == "-h") or (args[0] == "--help"))) ):
    show_usage()
    return

  if (params is None):
    pcl = iotbx.phil.process_command_line_with_files(
      args=args,
      master_phil_string=master_phil_str,
      pdb_file_def="file_name")
    work_params = pcl.work.extract()
  # or use parameters defined by GUI
  else:
    work_params = params
  pdb_inp = iotbx.pdb.input(file_name=work_params.file_name)

  pdb_int_params = mmtbx.model.manager.get_default_pdb_interpretation_params()
  pdb_int_params.pdb_interpretation.clash_guard.nonbonded_distance_threshold=None
  model = mmtbx.model.manager(
      model_input = pdb_inp,
      pdb_interpretation_params = pdb_int_params,
      build_grm = True)

  fixed_pdb_h = model.get_hierarchy().deep_copy()
  fixed_pdb_h.reset_atom_i_seqs()

  fixed_pdb_h = fix_rotamer_outliers(
      model = model,
      radius=work_params.radius,
      )

  res_pdb_str = model.model_as_pdb()
  with open("%s.pdb" % work_params.output_prefix, "w") as f:
    f.write(res_pdb_str)

if (__name__ == "__main__"):
  run(sys.argv[1:])
