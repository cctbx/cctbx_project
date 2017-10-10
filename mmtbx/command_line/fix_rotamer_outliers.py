from __future__ import division
from __future__ import print_function
# LIBTBX_SET_DISPATCHER_NAME mmtbx.fix_rotamer_outliers

import sys
import iotbx.pdb
import mmtbx.monomer_library.server
from mmtbx.command_line.geometry_minimization import get_geometry_restraints_manager
import iotbx.phil
from libtbx.utils import Sorry
from mmtbx.utils import fix_rotamer_outliers
from mmtbx.rotamer.rotamer_eval import RotamerEval


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
  .multiple = True
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
  print(help_msg)
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
  pdb_files = work_params.file_name

  from mmtbx.monomer_library.pdb_interpretation import grand_master_phil_str
  params = iotbx.phil.parse(
      input_string=grand_master_phil_str, process_includes=True).extract()
  params.pdb_interpretation.clash_guard.nonbonded_distance_threshold=None

  mon_lib_srv = mmtbx.monomer_library.server.server()
  ppf_srv = mmtbx.utils.process_pdb_file_srv(
      crystal_symmetry          = None,
      pdb_interpretation_params = params.pdb_interpretation,
      stop_for_unknowns         = False,
      log                       = log,
      cif_objects               = None, # need to figure out how to get them
      mon_lib_srv               = mon_lib_srv,
      ener_lib                  = None,
      use_neutron_distances     = False)
  processed_pdb_file, pdb_inp = ppf_srv.process_pdb_files(
      pdb_file_names = pdb_files)
  cs = ppf_srv.crystal_symmetry
  xrs = processed_pdb_file.xray_structure(show_summary = True)

  if(xrs is None):
    raise Sorry("Cannot extract xray_structure.")

  grm = get_geometry_restraints_manager(
      processed_pdb_file=processed_pdb_file,
      xray_structure=xrs,
      log=log)

  pdb_h = processed_pdb_file.all_chain_proxies.pdb_hierarchy
  fixed_pdb_h = pdb_h.deep_copy()
  fixed_pdb_h.reset_atom_i_seqs()

  rotamer_manager = RotamerEval(mon_lib_srv=mon_lib_srv)
  fixed_pdb_h = fix_rotamer_outliers(
      pdb_hierarchy=fixed_pdb_h,
      grm=grm.geometry,
      xrs=xrs,
      radius=work_params.radius,
      mon_lib_srv=mon_lib_srv,
      rotamer_manager=rotamer_manager,
      asc=None)

  iotbx.pdb.write_whole_pdb_file(
      file_name="%s.pdb" % work_params.output_prefix,
      processed_pdb_file=processed_pdb_file,
      pdb_hierarchy=fixed_pdb_h,
      crystal_symmetry=cs)

if (__name__ == "__main__"):
  run(sys.argv[1:])
