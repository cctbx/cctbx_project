# LIBTBX_SET_DISPATCHER_NAME mmtbx.simple_build_alt_confs

from __future__ import absolute_import, division, print_function
from libtbx.str_utils import make_header
from libtbx.utils import multi_out
import os.path
import sys

def master_phil():
  from mmtbx.command_line import generate_master_phil_with_inputs
  return generate_master_phil_with_inputs(
    phil_string="""
selection = None
  .type = atom_selection
nproc = Auto
  .type = int
output {
  file_name = alternates.pdb
    .type = path
  verbose = False
    .type = bool
  debug = False
    .type = bool
}
include scope mmtbx.building.alternate_conformations.single_residue.master_phil_str
""")

def run(args, out=None):
  if (out is None) : out = sys.stdout
  from mmtbx.building.alternate_conformations import single_residue
  import mmtbx.command_line
  cmdline = mmtbx.command_line.load_model_and_data(
    args=args,
    master_phil=master_phil(),
    process_pdb_file=True,
    create_fmodel=True,
    out=out,
    usage_string="""\
mmtbx.build_alt_confs_simple model.pdb data.mtz [options]

Simple tool for building alternate conformations by real-space refinement into
difference density.  Not intended for production use - use the program
mmtbx.build_alternate_conformations if you want refinement and post-processing.
""")
  params = cmdline.params
  validate_params(params)
  log = multi_out()
  log.register("stdout", out)
  log_file_name = os.path.splitext(params.output.file_name)[0] + ".log"
  logfile = open(log_file_name, "w")
  log.register("logfile", logfile)
  pdb_hierarchy, n_alternates = single_residue.build_cycle(
    pdb_hierarchy = cmdline.pdb_hierarchy,
    fmodel = cmdline.fmodel,
    geometry_restraints_manager = cmdline.geometry,
    params = params,
    cif_objects=cmdline.cif_objects,
    selection=params.selection,
    nproc=params.nproc,
    verbose=params.output.verbose,
    debug=params.output.debug,
    out=log)
  # TODO real-space refinement of multi-conformer model
  f = open(params.output.file_name, "w")
  f.write(pdb_hierarchy.as_pdb_string(
    crystal_symmetry=cmdline.fmodel.xray_structure))
  f.close()
  make_header("Building complete", out=out)
  print("", file=log)
  print("Wrote %s" % params.output.file_name, file=log)
  print("You MUST refine this model before using it!", file=log)

def validate_params(params):
  if (params.output.file_name is None):
    raise Sorry("Please specify an output file name.")
  return True

if (__name__ == "__main__"):
  run(sys.argv[1:])
