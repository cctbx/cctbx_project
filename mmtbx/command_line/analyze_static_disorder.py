# LIBTBX_SET_DISPATCHER_NAME mmtbx.analyze_static_disorder

from __future__ import absolute_import, division, print_function
from libtbx.str_utils import make_header, make_sub_header
from libtbx import easy_pickle
import os.path
import sys

def master_phil():
  from mmtbx.command_line import generate_master_phil_with_inputs
  return generate_master_phil_with_inputs(
    enable_stop_for_unknowns=False,
    enable_pdb_interpretation_params=True,
    phil_string="""
ignore_inconsistent_occupancy = True
  .type = bool
pickle = False
  .type = bool
verbose = False
  .type = bool
""")

def run(args, out=sys.stdout):
  from mmtbx.disorder import analyze_model
  import mmtbx.validation.molprobity
  import mmtbx.command_line
  cmdline = mmtbx.command_line.load_model_and_data(
    args=args,
    master_phil=master_phil(),
    require_data=False,
    create_fmodel=True,
    process_pdb_file=True,
    usage_string="mmtbx.analyze_static_disorder model.pdb",
    out=out)
  hierarchy = cmdline.pdb_hierarchy
  params = cmdline.params
  validation = mmtbx.validation.molprobity.molprobity(
    pdb_hierarchy=hierarchy,
    xray_structure=cmdline.xray_structure,
    fmodel=cmdline.fmodel,
    crystal_symmetry=cmdline.crystal_symmetry,
    geometry_restraints_manager=cmdline.geometry,
    header_info=None,
    keep_hydrogens=False,
    outliers_only=False,
    nuclear=False)
  segments = []
  make_header("Analyzing model", out=out)
  if (params.ignore_inconsistent_occupancy):
    print("Discontinuous occupancies will be ignored.", file=out)
  process = analyze_model.process_pdb_hierarchy(
    pdb_hierarchy=hierarchy,
    validation=validation,
    ignore_inconsistent_occupancy=params.ignore_inconsistent_occupancy,
    log=out)
  make_sub_header("MolProbity validation", out=out)
  validation.show_summary(out=out)
  make_sub_header("Disorder analysis", out=out)
  if (process.n_disordered == 0):
    print("No alternate conformations found.", file=out)
  else :
    process.show(out=out, verbose=params.verbose)
  if (params.pickle):
    file_name = os.path.basename(
      os.path.splitext(params.input.pdb.file_name[0])[0]) + ".pkl"
    easy_pickle.dump(file_name, process)
  return process

if (__name__ == "__main__"):
  run(sys.argv[1:])
