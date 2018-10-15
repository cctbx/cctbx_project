from __future__ import division

#from __future__ import division, print_function

from iotbx.cli_parser import run_program
from mmtbx.programs import polder

if __name__ == '__main__':
  run_program(program_class=polder.Program)

# =============================================================================
# old code - maybe necessary until GUI is updated

# LIBTBX_SET_DISPATCHER_NAME phenix.polder
import time
import os, sys
import mmtbx.utils
import iotbx.pdb
#from mmtbx import map_tools
from iotbx import ccp4_map
from iotbx import file_reader
#from iotbx import phil
from iotbx import reflection_file_utils
from iotbx import crystal_symmetry_from_any
from cStringIO import StringIO
from cctbx import maptbx
#from cctbx import miller
from cctbx.array_family import flex
from libtbx import group_args
from libtbx.utils import Sorry
from libtbx.utils import multi_out
#-----
import mmtbx.maps.polder
#-----

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
output_file_name_prefix = None
  .type = str
  .short_caption = Output prefix
  .help = Prefix for output filename
mask_output = False
  .type = bool
  .short_caption = Output masks
  .help = Additional output: ccp4 maps containing the solvent mask for inital \
   model (mask_all.ccp4), when ligand is omitted (mask_omit.ccp4) and the mask \
   used for polder (mask_polder.ccp4)
debug = False
  .type = bool
  .expert_level = 3
  .short_caption = Output biased map
  .help = Additional output: biased omit map (ligand used for mask calculation \
   but omitted from model)
gui
  .help = "GUI-specific parameter required for output directory"
{
  output_dir = None
  .type = path
  .style = output_dir
}
include scope mmtbx.maps.polder.master_params
include scope libtbx.phil.interface.tracking_params
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
  if (params.polder.box_buffer is not None and (params.polder.box_buffer < 0 or
      params.polder.box_buffer > 5)):
    raise Sorry("Box buffer out of range: must be between 0 and 5")
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

def format_print_rfactors(log, r_work, r_free):
  print >> log, "R_work = %6.4f R_free = %6.4f" % (r_work, r_free)
  print >> log, "*"*79

def print_rfactors(debug, results, log):
  fmodel_input  = results.fmodel_input
  fmodel_biased = results.fmodel_biased
  fmodel_omit   = results.fmodel_omit
  fmodel_polder = results.fmodel_polder
  print >> log, "R factors for unmodified input model and data:"
  format_print_rfactors(log, fmodel_input.r_work(), fmodel_input.r_free())
  if (debug):
    print >> log, "R factor when ligand is used for mask calculation (biased map):"
    format_print_rfactors(log, fmodel_biased.r_work(), fmodel_biased.r_free())
  print >> log, "R factor for polder map"
  format_print_rfactors(log, fmodel_polder.r_work(), fmodel_polder.r_free())
  print >> log, "R factor for OMIT map (ligand is excluded for mask calculation):"
  format_print_rfactors(log, fmodel_omit.r_work(), fmodel_omit.r_free())

def result_message(cc12, cc13, cc23):
  if (cc13 < 0.7 or
      (cc23 > cc12 and cc23 > cc13) or (cc13 < cc12 and cc13 < cc23)):
    msg = """The polder map is very likely to show bulk-solvent or noise."""
  elif (cc13 >= 0.8):
    msg = 'The polder map is likely to show the ligand.'
  elif (cc13 >= 0.7 and cc13 < 0.8):
    if (cc23 < 0.7*cc13):
      msg = """The polder map is more likely to show ligand than bulk solvent.
It is recommended to carefully inspect the maps to confirm."""
    else:
      msg = """The polder map is more likely to show bulk-solvent or noise
instead of the ligand. But it is recommended to inspect the maps to confirm."""
  return msg

def print_validation(log, results, debug, pdb_hierarchy_selected):
  box_1 = results.box_1
  box_2 = results.box_2
  box_3 = results.box_3
  sites_cart_box = box_1.xray_structure_box.sites_cart()
  sel = maptbx.grid_indices_around_sites(
    unit_cell  = box_1.xray_structure_box.unit_cell(),
    fft_n_real = box_1.map_box.focus(),
    fft_m_real = box_1.map_box.all(),
    sites_cart = sites_cart_box,
    site_radii = flex.double(sites_cart_box.size(), 2.0))
  b1 = box_1.map_box.select(sel).as_1d()
  b2 = box_2.map_box.select(sel).as_1d()
  b3 = box_3.map_box.select(sel).as_1d()
  print >> log, "Map 1: calculated Fobs with ligand"
  print >> log, "Map 2: calculated Fobs without ligand"
  print >> log, "Map 3: real Fobs data"
  cc12 = flex.linear_correlation(x=b1,y=b2).coefficient()
  cc13 = flex.linear_correlation(x=b1,y=b3).coefficient()
  cc23 = flex.linear_correlation(x=b2,y=b3).coefficient()
  print >> log, "CC(1,2): %6.4f" % cc12
  print >> log, "CC(1,3): %6.4f" % cc13
  print >> log, "CC(2,3): %6.4f" % cc23
  #### D-function
  b1 = maptbx.volume_scale_1d(map=b1, n_bins=10000).map_data()
  b2 = maptbx.volume_scale_1d(map=b2, n_bins=10000).map_data()
  b3 = maptbx.volume_scale_1d(map=b3, n_bins=10000).map_data()
  print >> log, "Peak CC:"
  print >> log, "CC(1,2): %6.4f"%flex.linear_correlation(x=b1,y=b2).coefficient()
  print >> log, "CC(1,3): %6.4f"%flex.linear_correlation(x=b1,y=b3).coefficient()
  print >> log, "CC(2,3): %6.4f"%flex.linear_correlation(x=b2,y=b3).coefficient()
  cutoffs = flex.double(
    [i/10. for i in range(1,10)]+[i/100 for i in range(91,100)])
  d12 = maptbx.discrepancy_function(map_1=b1, map_2=b2, cutoffs=cutoffs)
  d13 = maptbx.discrepancy_function(map_1=b1, map_2=b3, cutoffs=cutoffs)
  d23 = maptbx.discrepancy_function(map_1=b2, map_2=b3, cutoffs=cutoffs)
  print >> log, "q    D(1,2) D(1,3) D(2,3)"
  for c,d12_,d13_,d23_ in zip(cutoffs,d12,d13,d23):
    print >> log, "%4.2f %6.4f %6.4f %6.4f"%(c,d12_,d13_,d23_)
  ###
  if(debug):
    #box_1.write_ccp4_map(file_name="box_1_polder.ccp4")
    #box_2.write_ccp4_map(file_name="box_2_polder.ccp4")
    #box_3.write_ccp4_map(file_name="box_3_polder.ccp4")
    write_map_box(
      box      = box_1,
      filename = "box_1_polder.ccp4")
    write_map_box(
      box      = box_2,
      filename = "box_2_polder.ccp4")
    write_map_box(
      box      = box_3,
      filename = "box_3_polder.ccp4")
    pdb_hierarchy_selected.adopt_xray_structure(box_1.xray_structure_box)
    pdb_hierarchy_selected.write_pdb_file(file_name="box_polder.pdb",
      crystal_symmetry=box_1.box_crystal_symmetry)
  #
  print >> log, '*' * 79
  message = result_message(cc12 = cc12, cc13 = cc13, cc23 = cc23)
  print >> log, message
  return message

def write_map_box(box, filename):
    ccp4_map.write_ccp4_map(
      file_name   = filename,
      unit_cell   = box.xray_structure_box.unit_cell(),
      space_group = box.xray_structure_box.space_group(),
      map_data    = box.map_box,
      labels      = flex.std_string([""]))


def write_files(results, mask_output, debug, f_obs, prefix, log):
  # write mask files (if specified)
  if (mask_output):
    masks = [results.mask_data_all, results.mask_data_omit, results.mask_data_polder]
    filenames = ["all", "omit", "polder"]
    for mask_data, filename in zip(masks, filenames):
      ccp4_map.write_ccp4_map(
        file_name   = "mask_" + filename + ".ccp4",
        unit_cell   = f_obs.unit_cell(),
        space_group = f_obs.space_group(),
        map_data    = mask_data,
        labels      = flex.std_string([""]))
  mtz_dataset = results.mc_polder.as_mtz_dataset(
    column_root_label = "mFo-DFc_polder")
  # add map coeffs for biased map if debug=True
  if (debug):
    mtz_dataset.add_miller_array(
      miller_array      = results.mc_biased,
      column_root_label = "mFo-DFc_bias_omit")
  mtz_dataset.add_miller_array(
    miller_array      = results.mc_omit,
    column_root_label = "mFo-DFc_omit")
  mtz_object = mtz_dataset.mtz_object()
  polder_file_name = "polder_map_coeffs.mtz"
  if (prefix is not None):
    polder_file_name = prefix + "_" + polder_file_name
  mtz_object.write(file_name = polder_file_name)
  print >> log, 'File %s was written.' % polder_file_name

def get_inputs(args, log, master_params, validated):
  inputs = mmtbx.utils.process_command_line_args(
    args                             = args,
    master_params                    = master_params,
    suppress_symmetry_related_errors = True)
  params = inputs.params.extract()
  print params.model_file_name
  # Check model file
  if (len(inputs.pdb_file_names) == 0 and (params.model_file_name is None)):
    raise Sorry("No model file found.")
  elif (len(inputs.pdb_file_names) == 1):
    params.model_file_name = inputs.pdb_file_names[0]
  elif (len(inputs.pdb_file_names) > 1):
  #else:
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
    params.output_file_name_prefix is None):
    params.output_file_name_prefix = model_basename
  new_params =  master_params.format(python_object = params)
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
  return group_args(
    f_obs             = f_obs,
    r_free_flags      = r_free_flags,
    xray_structure    = xray_structure,
    pdb_hierarchy     = pdb_hierarchy,
    params            = params)

def run(args, validated = False, out=sys.stdout):
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
  #
  inputs = get_inputs(
    args          = args,
    log           = log,
    master_params = master_params(),
    validated     = validated)
  params = inputs.params
  #
  # Check if selection syntax makes sense
  print >> log, "*"*79
  print >> log, "selecting atoms..."
  print >> log, "Selection string:", params.solvent_exclusion_mask_selection
  selection_bool = inputs.pdb_hierarchy.atom_selection_cache().selection(
    string = params.solvent_exclusion_mask_selection)
  n_selected = selection_bool.count(True)
  n_selected_all = inputs.pdb_hierarchy.atom_selection_cache().selection(
    string = 'all').count(True)
  if(n_selected == 0):
    raise Sorry("No atoms where selected. Check selection syntax again.")
  if (n_selected/n_selected_all > 0.5):
    raise Sorry("""More than half of total number of atoms selected. Omit
      selection should be smaller, such as one ligand or a few residues.""")
  print >> log, "Number of atoms selected:", n_selected
  pdb_hierarchy_selected = inputs.pdb_hierarchy.select(selection_bool)
  ligand_str = pdb_hierarchy_selected.as_pdb_string()
  print >> log, "Atoms selected:\n", ligand_str
  print >> log, "*"*79
  #
  # Calculate polder map
  polder_object = mmtbx.maps.polder.compute_polder_map(
    f_obs             = inputs.f_obs,
    r_free_flags      = inputs.r_free_flags,
    xray_structure    = inputs.xray_structure,
    pdb_hierarchy     = inputs.pdb_hierarchy,
    params            = inputs.params.polder,
    selection_bool    = selection_bool)
  polder_object.validate()
  polder_object.run()
  results = polder_object.get_results()
  print_rfactors(
    debug   = params.debug,
    results = results,
    log     = log)
  write_files(
    results     = results,
    mask_output = params.mask_output,
    debug       = params.debug,
    f_obs       = inputs.f_obs,
    prefix      = params.output_file_name_prefix,
    log         = log)
  message = None
  if (not params.polder.compute_box):
    message = print_validation(
      log     = log,
      pdb_hierarchy_selected = pdb_hierarchy_selected,
      results = results,
      debug   = params.debug)
  #
  print >> log, '*' * 79
  print >> log, "Finished."
  # results object not returned because it contains maps
  return group_args(message=message)

# =============================================================================
# GUI-specific class for running command
from libtbx import runtime_utils
class launcher (runtime_utils.target_with_save_result) :
  def run (self) :
    import os
    from wxGUI2 import utils
    utils.safe_makedirs(self.output_dir)
    os.chdir(self.output_dir)
    result = run(args=self.args, validated=True, out=sys.stdout)
    return result

# =============================================================================
#if (__name__ == "__main__"):
#  assert cmd_run(args=sys.argv[1:]) is None # assert here is intentional
#if __name__ == "__main__":
#  t0 = time.time()
#  run(args = sys.argv[1:])
#  print "Time:", round(time.time()-t0, 2)
