from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME phenix.polder
import time
import os, sys
#import mmtbx.f_model
import mmtbx.utils
#import mmtbx.masks
import iotbx.pdb
#-----
import mmtbx.maps.polder_lib
#-----
#from mmtbx import map_tools
#from iotbx import ccp4_map
from iotbx import file_reader
#from iotbx import phil
from iotbx import reflection_file_utils
from iotbx import crystal_symmetry_from_any
from cStringIO import StringIO
#from cctbx import maptbx
#from cctbx import miller
#from cctbx.array_family import flex
#from libtbx import group_args
from libtbx.utils import Sorry
from libtbx.utils import multi_out
#from libtbx.math_utils import ifloor, iceil

legend = """\

Computes omit maps by excluding the bulk solvent in the area around a
selection. An example of application are ligand omit maps. Polder omit maps
can be helpful if the ligand density is weak and obscured by bulk solvent
in conventional omit maps (where the ligand is deleted from the model).

Inputs:
  - File with reflection data (Fobs or Iobs) and R-free flags. It can
    be in most of known formats and spread across multiple files;
  - label(s) selecting which reflection data arrays should be used (in case
    there are multiple choices in input file; otherwise there is no need to
    provide labels);
  - Model file;
  - Atom selection (such as ligand)

Usage examples:
  1. phenix.polder model.cif data.mtz selection="chain A and resseq 1"
  2. phenix.polder model.pdb data.hkl data_labels="FP" selection="chain A"
  3. phenix.polder a.hkl b.hkl model.pdb selection="resseq 435"

Output:
  MTZ file with map coefficients for:
  Polder map:
  - mFo-DFc_polder    : polder difference map coefficients
  - PHImFo-DFc_polder : corresponding phases
  Omit map:
  For this map, the OMIT selection is deleted from the model and bulk solvent
  enters the area:
  - mFo-DFc_omit      : omit difference map coefficients
  - PHImFo-DFc_omit   : corresponding phases

Optional output:
  CCP4 files with mask data:
  - mask_all.ccp4    : mask of original model
  - mask_omit.ccp4   : mask when ligand is omitted
  - mask_polder.ccp4 : mask obtained by polder procedure
"""

#include scope libtbx.phil.interface.tracking_params
master_params_str = """
model_file_name = None
  .type = path
  .short_caption = Model file
  .multiple = False
  .help = Model file name
  .style = file_type:pdb bold input_file
solvent_exclusion_mask_selection = None
  .type = str
  .short_caption = Omit selection
  .help = Atoms around which bulk solvent mask is set to zero
  .input_size = 400
reflection_file_name = None
  .type = path
  .short_caption = Data file
  .help = File with experimental data (most of formats: CNS, SHELX, MTZ, etc).
  .style = file_type:hkl bold input_file process_hkl child:fobs:data_labels \
           child:rfree:r_free_flags_labels child:d_min:high_resolution \
           child:d_max:low_resolution
data_labels = None
  .type = str
  .short_caption = Data labels
  .help = Labels for experimental data.
  .style = renderer:draw_fobs_label_widget parent:file_name:reflection_file_name
r_free_flags_labels = None
  .type = str
  .short_caption = Rfree labels
  .help = Labels for free reflections.
  .style = renderer:draw_rfree_label_widget parent:file_name:reflection_file_name
high_resolution = None
  .type = float
  .short_caption = High resolution
  .help = High resolution limit
low_resolution = None
  .type = float
  .short_caption = Low resolution
  .help = Low resolution limit
scattering_table = *n_gaussian wk1995 it1992 neutron electron
  .type = choice
  .short_caption = Scattering table
  .help = Scattering table for structure factors calculations
box_buffer = None
  .type = float
  .short_caption = Buffer around selection box
  .help = Buffer around selection box: Increase the box for resetting the mask \
   by a buffer.
compute_box = False
  .type = bool
  .short_caption = Reset mask within a box defined by atom selection
  .help = Reset mask within a box (parallel to unit cell axes) defined by an \
   atom selection
gui
  .help = "GUI-specific parameter required for output directory"
{
  output_dir = None
  .type = path
  .style = output_dir
}
include scope mmtbx.maps.polder_lib.master_params
"""

def master_params():
  return iotbx.phil.parse(master_params_str, process_includes=True)

def print_legend_and_usage(log):
  print >> log, "-"*79
  print >> log, "                               phenix.polder"
  print >> log, "-"*79
  print >> log, legend
  print >> log, "-"*79
  print >> log, master_params().show()

# validation for GUI
def validate_params(params):
  if (params.solvent_exclusion_mask_selection is None):
    raise Sorry("No selection for mask calculation found.")
  if (params.polder.sphere_radius < 3):
    raise Sorry("Sphere radius out of range: must be larger than 3 A")
  if (params.model_file_name is None):
    raise Sorry("Model file should be given")
  if (params.reflection_file_name is None):
    raise Sorry("Data file should be given")
  # check if file type is OK
  file_reader.any_file(
    file_name = params.model_file_name).check_file_type(expected_type = 'pdb')
  file_reader.any_file(
    file_name = params.reflection_file_name).check_file_type(
      expected_type = 'hkl')
  if (params.data_labels is None):
    raise Sorry("Data labels should be given")
  if (params.polder.resolution_factor < 0.0):
    raise Sorry(
      "Please use a positive value for the resolution gridding factor.")
  return True

def prepare_f_obs_and_flags(f_obs, r_free_flags):
  sel = f_obs.data()>0
  f_obs = f_obs.select(sel)
  r_free_flags = r_free_flags.select(sel)
  #
  merged = f_obs.as_non_anomalous_array().merge_equivalents()
  f_obs = merged.array().set_observation_type(f_obs)
  #
  merged = r_free_flags.as_non_anomalous_array().merge_equivalents()
  r_free_flags = merged.array().set_observation_type(r_free_flags)
  #
  f_obs, r_free_flags = f_obs.common_sets(r_free_flags)
  return f_obs, r_free_flags

def obtain_cs_if_gui_input(model_file_name, reflection_file_name):
  crystal_symmetries = []
  for f in [str(model_file_name), str(reflection_file_name)]:
    cs = crystal_symmetry_from_any.extract_from(f)
    if(cs is not None): crystal_symmetries.append(cs)
  #
  if(len(crystal_symmetries) == 1): crystal_symmetry = crystal_symmetries[0]
  elif(len(crystal_symmetries) == 0):
    raise Sorry("No crystal symmetry found.")
  else:
    if(not crystal_symmetries[0].is_similar_symmetry(crystal_symmetries[1])):
      raise Sorry("Crystal symmetry mismatch between different files.")
    crystal_symmetry = crystal_symmetries[0]
  return crystal_symmetry

def cmd_run(args, validated = False):
  # processing command-line (out of the polder class)
  log = multi_out()
  log.register("stdout", sys.stdout)
  log_file_name = "polder.log"
  logfile = open(log_file_name, "w")
  log.register("logfile", logfile)
  if len(args) == 0:
    print_legend_and_usage(log)
    return
  print >> log, "phenix.polder is running..."
  print >> log, "input parameters:\n", args
  inputs = mmtbx.utils.process_command_line_args(
    args                             = args,
    master_params                    = master_params(),
    suppress_symmetry_related_errors = True)
  params = inputs.params.extract()
  #
  # Check model file
  if (len(inputs.pdb_file_names) == 0 and (params.model_file_name is None)):
    raise Sorry("No model file found.")
  elif (len(inputs.pdb_file_names) == 1):
    params.model_file_name = inputs.pdb_file_names[0]
  else:
    raise Sorry("Only one model file should be given")
  #
  # Check reflection file(s)
  reflection_files = inputs.reflection_files
  if (len(reflection_files) == 0):
    if (params.reflection_file_name is None):
      raise Sorry("No reflection file found.")
    else:
      hkl_in = file_reader.any_file(params.reflection_file_name,
        force_type="hkl")
      hkl_in.assert_file_type("hkl")
      reflection_files = [ hkl_in.file_object ]
  #
  # Get crystal symmetry
  crystal_symmetry = None
  crystal_symmetry = inputs.crystal_symmetry
  if (crystal_symmetry is None):
    crystal_symmetry = obtain_cs_if_gui_input(
      model_file_name      = params.model_file_name,
      reflection_file_name = params.reflection_file_name)
  print >> log, "Working crystal symmetry after inspecting all inputs:"
  crystal_symmetry.show_summary(f=log, prefix="  ")
  #
  # Get data labels
  f_obs, r_free_flags = None, None
  rfs = reflection_file_utils.reflection_file_server(
    crystal_symmetry = crystal_symmetry,
    force_symmetry   = True,
    reflection_files = reflection_files,
    err              = StringIO())
  parameters = mmtbx.utils.data_and_flags_master_params().extract()
  if (params.data_labels is not None):
    parameters.labels = params.data_labels
  if (params.r_free_flags_labels is not None):
    parameters.r_free_flags.label = params.r_free_flags_labels
  determined_data_and_flags = mmtbx.utils.determine_data_and_flags(
    reflection_file_server = rfs,
    parameters             = parameters,
    keep_going             = True,
    working_point_group = crystal_symmetry.space_group().build_derived_point_group(),
    log                    = StringIO(),
    symmetry_safety_check  = True)
  f_obs = determined_data_and_flags.f_obs
  if (params.data_labels is None):
    params.data_labels = f_obs.info().label_string()
  if (params.reflection_file_name is None):
    params.reflection_file_name = parameters.file_name
  r_free_flags = determined_data_and_flags.r_free_flags
  assert f_obs is not None
  print >> log,  "Input data:"
  print >> log, "  Iobs or Fobs:", f_obs.info().labels
  if (r_free_flags is not None):
    print >> log, "  Free-R flags:", r_free_flags.info().labels
    params.r_free_flags_labels = r_free_flags.info().label_string()
  else:
    print >> log, "  Free-R flags: Not present"
  model_basename = os.path.basename(params.model_file_name.split(".")[0])
  if (len(model_basename) > 0 and
    params.polder.output_file_name_prefix is None):
    params.polder.output_file_name_prefix = model_basename
  #print params.output_file_name_prefix
  new_params =  master_params().format(python_object = params)
  print >> log, "*"*79
  new_params.show()
  if (not validated):
    validate_params(params)
  pdb_input = iotbx.pdb.input(file_name = params.model_file_name)
  pdb_hierarchy = pdb_input.construct_hierarchy()
  xray_structure = pdb_hierarchy.extract_xray_structure(
    crystal_symmetry = crystal_symmetry)
  # DON'T USE:
  # xray_structure = pdb_input.xray_structure_simple()
  # because the atom order might be wrong
  mmtbx.utils.setup_scattering_dictionaries(
    scattering_table = params.scattering_table,
    xray_structure   = xray_structure,
    d_min            = f_obs.d_min())
  f_obs = f_obs.resolution_filter(
    d_min = params.high_resolution,
    d_max = params.low_resolution)
  if (r_free_flags is not None):
    r_free_flags = r_free_flags.resolution_filter(
      d_min = params.high_resolution,
      d_max = params.low_resolution)
  #
  # If data are anomalous
  if (f_obs.anomalous_flag()):
    f_obs, r_free_flags = prepare_f_obs_and_flags(
      f_obs        = f_obs,
      r_free_flags = r_free_flags)
  #
  # Check if selection syntax makes sense
  print >> log, "*"*79
  print >> log, "selecting atoms..."
  print >> log, "Selection string:", params.solvent_exclusion_mask_selection
  selection_bool = pdb_hierarchy.atom_selection_cache().selection(
    string = params.solvent_exclusion_mask_selection)
  n_selected = selection_bool.count(True)
  if(n_selected == 0):
    raise Sorry("No atoms where selected. Check selection syntax again.")
  print >> log, "Number of atoms selected:", n_selected
  pdb_hierarchy_selected = pdb_hierarchy.select(selection_bool)
  ligand_str = pdb_hierarchy_selected.as_pdb_string()
  print >> log, "Atoms selected:\n", ligand_str
  print >> log, "*"*79
  #
  # Calculate polder map
  polder_object = mmtbx.maps.polder_lib.compute_polder_map(
    f_obs             = f_obs,
    r_free_flags      = r_free_flags,
    xray_structure    = xray_structure,
    pdb_hierarchy     = pdb_hierarchy,
    params            = params.polder,
    selection_bool    = selection_bool)
  polder_object.validate()
  polder_object.run()
  polder_object.print_rfactors()
  polder_object.write_files()
  polder_object.validate_polder_map()
  #
  print >> log, '*' * 79
  print >> log, "Finished."
  return True

# =============================================================================
# GUI-specific class for running command
from libtbx import runtime_utils
class launcher (runtime_utils.target_with_save_result) :
  def run (self) :
    import os
    from wxGUI2 import utils
    utils.safe_makedirs(self.output_dir)
    os.chdir(self.output_dir)
    result = cmd_run(args=self.args, validated=True, out=sys.stdout)
    return result

# =============================================================================
#if (__name__ == "__main__"):
#  assert cmd_run(args=sys.argv[1:]) is None # assert here is intentional
if __name__ == "__main__":
  t0 = time.time()
  cmd_run(args = sys.argv[1:])
  print "Time:", round(time.time()-t0, 2)
  print 'polder_new'
