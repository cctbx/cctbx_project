# LIBTBX_SET_DISPATCHER_NAME phenix.reciprocal_space_arrays

import sys, os
import mmtbx.utils
from iotbx import reflection_file_utils
from cStringIO import StringIO
import mmtbx.f_model
from libtbx.utils import Sorry
import iotbx.phil

msg="""\

phenix.reciprocal_space_arrays:
  tool to compute various arrays such as Fcalc, Fmask, Fmodel, Fbulk, and more.

Usage:
  phenix.reciprocal_space_arrays mode.pdb data.hkl [options]

Output: MTZ file with data arrays.
"""

master_params_str="""\
f_obs_label = None
  .type = str
r_free_flags_label = None
  .type = str
remove_f_obs_outliers = True
  .type = bool
bulk_solvent_and_scaling = True
  .type = bool
hendrickson_lattman_coefficients_label = None
  .type = str
output_file_name = None
  .type = str
"""


def extract_experimental_phases(experimental_phases, f_obs):
  if(experimental_phases is not None):
    if(not f_obs.anomalous_flag()):
      if(experimental_phases.anomalous_flag()):
        experimental_phases = experimental_phases.average_bijvoet_mates()
    elif(not experimental_phases.anomalous_flag()):
      experimental_phases = experimental_phases.generate_bijvoet_mates()
    return experimental_phases.map_to_asu().matching_set(other = f_obs,
      data_substitute=(0,0,0,0))
  else: return None

def run(args):
  if(len(args)==0 or (len(args)==1 and args[0].lower() in
     ["-h","--h","-help","--help","h","help"])):
       print msg
       return
  #
  print "Default params:"
  parsed = iotbx.phil.parse(master_params_str)
  parsed.show(prefix="  ")
  processed_args = mmtbx.utils.process_command_line_args(args = args,
    log = sys.stdout, master_params = parsed)
  params = processed_args.params.extract()
  reflection_files = processed_args.reflection_files
  if(len(reflection_files) == 0):
    raise Sorry("No reflection file found.")
  crystal_symmetry = processed_args.crystal_symmetry
  if(crystal_symmetry is None):
    raise Sorry("No crystal symmetry found.")
  if(len(processed_args.pdb_file_names) == 0):
    raise Sorry("No PDB file found.")
  pdb_file_names = processed_args.pdb_file_names
  #
  rfs = reflection_file_utils.reflection_file_server(
    crystal_symmetry = crystal_symmetry,
    force_symmetry   = True,
    reflection_files = reflection_files,
    err              = StringIO())
  parameters = mmtbx.utils.data_and_flags_master_params().extract()
  parameters.labels = params.f_obs_label
  parameters.r_free_flags.label = params.r_free_flags_label
  determine_data_and_flags_result = mmtbx.utils.determine_data_and_flags(
    reflection_file_server = rfs,
    parameters             = parameters,
    keep_going             = True,
    log                    = StringIO())
  f_obs = determine_data_and_flags_result.f_obs
  print "Input data:"
  print "  Iobs or Fobs:", f_obs.info().labels
  r_free_flags = determine_data_and_flags_result.r_free_flags
  print "  Free-R flags:", r_free_flags.info().labels
  #
  parameters = mmtbx.utils.experimental_phases_params.extract()
  parameters.labels = params.hendrickson_lattman_coefficients_label
  experimental_phases_result = mmtbx.utils.determine_experimental_phases(
    reflection_file_server = rfs,
    parameters             = parameters,
    log                    = StringIO(),
    parameter_scope        = "",
    working_point_group    = None,
    symmetry_safety_check  = True,
    ignore_all_zeros       = True)
  if(experimental_phases_result is not None):
    print "  HL coefficients:", experimental_phases_result.info().labels
  experimental_phases = extract_experimental_phases(
    experimental_phases = experimental_phases_result, f_obs = f_obs)
  #
  if(r_free_flags is None):
    r_free_flags=f_obs.array(data=flex.bool(f_obs.data().size(), False))
  #
  mmtbx_pdb_file = mmtbx.utils.pdb_file(
    pdb_file_names   = pdb_file_names,
    crystal_symmetry = crystal_symmetry,
    use_elbow        = False,
    log              = sys.stdout)
  if(len(mmtbx_pdb_file.pdb_inp.xray_structures_simple())>1): #XXX support multi-models
    raise Sorry("Multiple model file not supported in this tool.")
  # XXX Twining not supported
  xray_structure = mmtbx_pdb_file.pdb_inp.xray_structure_simple()
  print "Input model:"
  print "  number of atoms:", xray_structure.scatterers().size()
  fmodel = mmtbx.f_model.manager(
    xray_structure = xray_structure,
    r_free_flags   = r_free_flags,
    f_obs          = f_obs,
    abcd           = experimental_phases)
  if(params.remove_f_obs_outliers): fmodel = fmodel.remove_outliers()
  if(params.bulk_solvent_and_scaling): fmodel.update_solvent_and_scale()
  print "Overall statistics:"
  fmodel.info().show_all()
  #
  print "Output data:"
  if(params.output_file_name is not None):
    output_file_name = params.output_file_name
  else:
    pdb_file_bn = os.path.basename(pdb_file_names[0])
    hkl_file_bn = os.path.basename(reflection_files[0].file_name())
    try: pdb_file_prefix = pdb_file_bn[:pdb_file_bn.index(".")]
    except ValueError: pdb_file_prefix = pdb_file_bn
    try:
      hkl_file_prefix = hkl_file_bn[:hkl_file_bn.index(".")]
    except ValueError: hkl_file_prefix = hkl_file_bn
    out_file_name = "%s_%s.mtz"%(pdb_file_prefix, hkl_file_prefix)
  print "  file name:", out_file_name
  print "  to see the contnt of %s:"%out_file_name
  print "    phenix.mtz.dump %s"%out_file_name
  out = open(out_file_name,"w")
  fmodel.export(out = out)
  out.close()
  print "All done."

if(__name__ == "__main__"):
  run(sys.argv[1:])
