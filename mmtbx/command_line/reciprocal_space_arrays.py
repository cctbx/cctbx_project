# LIBTBX_SET_DISPATCHER_NAME phenix.reciprocal_space_arrays

import mmtbx.utils
import mmtbx.f_model
from iotbx import reflection_file_utils
from iotbx import file_reader
import iotbx.phil
from libtbx import runtime_utils
from libtbx.utils import Sorry
from cStringIO import StringIO
import sys, os

msg="""\

phenix.reciprocal_space_arrays:
compute various arrays such as Fcalc, Fmask, Fmodel, Fbulk, and more.

Inputs:
  - File with reflection data (Fobs or Iobs), R-free flags, and optionally HL
    coefficients. It can be in most of known formats and spread across
    multiple files;
  - label(s) selecting which reflection data arrays should be used (in case
    there are multiple choices in input file, there is no need to provide labels
    otherwise);
  - PDB file with input model.

Usage examples:
  1. phenix.reciprocal_space_arrays model.pdb data.hkl f_obs_label="IOBS"
  2. phenix.reciprocal_space_arrays model.pdb data.hkl r_free_flags_label="FREE"

Output:
  MTZ file with data arrays. Depending on the input data, the following arrays
  may be present:

  - FOBS         : data from input reflection file
  - SIGFOBS      : corresponding sigmas
  - R_FREE_FLAGS : R-free flags (0 - work, 1 - test)
  - FMODEL       : total model structure factor. See phenix.fmodel for details
  - PHIFMODEL    : corresponding phases
  - FCALC        : Fcalc from atomic model
  - PHIFCALC     : corresponding phases
  - FMASK        : Fmask from bulk-solvent mask. See phenix.fmodel for details
  - PHIFMASK     : corresponding phases
  - FBULK        : bulk-solvent contribution. See phenix.fmodel for details
  - PHIFBULK     : corresponding phases
  - FB_CART      : overall anisotropic scale factor
  - FOM          : figures of merit
  - ALPHA        : ML parameter used in m&D calculation for 2mFo-DFc maps
  - BETA         : ML parameter used in m&D calculation for 2mFo-DFc maps
  - HLA          : HL coefficients from from input reflection file
  - HLB          : HL coefficients from from input reflection file
  - HLC          : HL coefficients from from input reflection file
  - HLD          : HL coefficients from from input reflection file
  - HLmodelA     : HL coefficients from the model (C=D=0)
  - HLmodelB     : HL coefficients from the model (C=D=0)
  - HLmodelC     : HL coefficients from the model (C=D=0)
  - HLmodelD     : HL coefficients from the model (C=D=0)
  - HLcombA      : combined HL: model + input
  - HLcombB      : combined HL: model + input
  - HLcombC      : combined HL: model + input
  - HLcombD      : combined HL: model + input
  - RESOLUTION   : resolution per reflection
"""

master_params_str="""\
hkl_file = None
  .type = path
  .short_caption = Experimental data
  .style = bold file_type:hkl input_file process_hkl child:fobs:f_obs_label \
    child:rfree:r_free_flags_label child:space_group:space_group \
    child:unit_cell:unit_cell \
    child:hl_coeffs:hendrickson_lattman_coefficients_label
pdb_file = None
  .type = path
  .short_caption = PDB file
  .style = bold file_type:pdb input_file
f_obs_label = None
  .type = str
  .short_caption = Data labels
  .input_size = 160
  .style = bold renderer:draw_fobs_label_widget
r_free_flags_label = None
  .type = str
  .short_caption = R-free flags
  .input_size = 160
  .style = bold renderer:draw_rfree_label_widget
remove_f_obs_outliers = True
  .type = bool
  .short_caption = Remove F-obs outliers
bulk_solvent_and_scaling = True
  .type = bool
  .short_caption = Bulk solvent correction and scaling
hendrickson_lattman_coefficients_label = None
  .type = str
  .short_caption = Hendrickson-Lattman coefficients
  .input_size = 160
  .style = renderer:draw_hl_label_widget
output_file_name = None
  .type = path
  .style = bold new_file
space_group = None
  .type = space_group
unit_cell = None
  .type = unit_cell
include scope libtbx.phil.interface.tracking_params
"""

def defaults(log):
  print >> log, "Default params::\n"
  parsed = iotbx.phil.parse(master_params_str, process_includes=True)
  parsed.show(prefix="  ", out=log)
  print >> log
  return parsed

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

def run(args, log = sys.stdout):
  if(len(args)==0):
    print >> log, msg
    defaults(log=log)
    return
  #
  parsed = defaults(log=log)
  processed_args = mmtbx.utils.process_command_line_args(args = args,
    log = sys.stdout, master_params = parsed)
  params = processed_args.params.extract()
  reflection_files = processed_args.reflection_files
  if(len(reflection_files) == 0):
    if (params.hkl_file is None) :
      raise Sorry("No reflection file found.")
    else :
      hkl_in = file_reader.any_file(params.hkl_file, force_type="hkl")
      hkl_in.assert_file_type("hkl")
      reflection_files = [ hkl_in.file_object ]
  crystal_symmetry = processed_args.crystal_symmetry
  if(crystal_symmetry is None):
    if (params.space_group is not None) and (params.unit_cell is not None) :
      from cctbx import crystal
      crystal_symmetry = crystal.symmetry(
        space_group_info=params.space_group,
        unit_cell=params.unit_cell)
    else :
      raise Sorry("No crystal symmetry found.")
  if(len(processed_args.pdb_file_names) == 0):
    if (params.pdb_file is None) :
      raise Sorry("No PDB file found.")
    else :
      pdb_file_names = [ params.pdb_file ]
  else :
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
  xray_structure = mmtbx_pdb_file.pdb_inp.xray_structure_simple(
    crystal_symmetry=crystal_symmetry,
    weak_symmetry=True)
  if (not xray_structure.unit_cell().is_similar_to(f_obs.unit_cell())) :
    raise Sorry("The unit cells in the PDB and reflections files are not "+
      "isomorphous.")
  print "Input model:"
  print "  number of atoms:", xray_structure.scatterers().size()
  fmodel = mmtbx.f_model.manager(
    xray_structure = xray_structure,
    r_free_flags   = r_free_flags,
    f_obs          = f_obs,
    abcd           = experimental_phases)
  fmodel.update_all_scales(
    remove_outliers = params.remove_f_obs_outliers,
    bulk_solvent_and_scaling = params.bulk_solvent_and_scaling)
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
    output_file_name = "%s_%s.mtz"%(pdb_file_prefix, hkl_file_prefix)
  print "  file name:", output_file_name
  print "  to see the contnt of %s:"%output_file_name
  print "    phenix.mtz.dump %s"%output_file_name
  out = open(output_file_name,"w")
  fmodel.export(out = out)
  out.close()
  print "All done."
  return output_file_name

def validate_params (params) :
  if (params.hkl_file is None) :
    raise Sorry("No reflections file provided.")
  elif (params.pdb_file is None) :
    raise Sorry("No PDB file provided.")
  elif (params.output_file_name is None) :
    raise Sorry("No output file name provided.")
  elif (params.f_obs_label is None) :
    raise Sorry("No data label selected.")
  elif (params.r_free_flags_label is None) :
    raise Sorry("No R-free flags label selected.")
  elif (params.space_group is None) or (params.unit_cell is None) :
    raise Sorry("Missing or incomplete symmetry information.")
  return True

def finish_job (result) :
  output_files, stats = [], []
  if (result is not None) and (os.path.isfile(result)) :
    output_files.append((result, "Data and phases"))
  return output_files, stats

class launcher (runtime_utils.target_with_save_result) :
  def run (self) :
    return run(args=self.args, log=sys.stdout)

if(__name__ == "__main__"):
  run(sys.argv[1:])
