from __future__ import absolute_import, division, print_function
from six.moves import zip
try:
  from phenix.program_template import ProgramTemplate
except ImportError:
  from libtbx.program_template import ProgramTemplate
import os
from libtbx.utils import null_out, Sorry
import mmtbx.maps.polder
from iotbx import crystal_symmetry_from_any
import mmtbx.utils
from iotbx import mrcfile
from libtbx import group_args
from cctbx.array_family import flex

master_phil_str = '''
include scope libtbx.phil.interface.tracking_params
include scope mmtbx.maps.polder.master_params
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
'''

# =============================================================================

class Program(ProgramTemplate):

  description = '''

Computes omit maps by excluding the bulk solvent in the area around a
selection. An example of application are ligand omit maps. Polder omit maps
can be helpful if the ligand density is weak and obscured by bulk solvent
in conventional omit maps (where the ligand is deleted from the model).

Inputs:
  - Reflection file: It can be in most of known formats and data can be
    spread across multiple files (Fobs in one file, Rfree in another)
  - Model file
  - Omit selection (such as for a ligand)
  - optional: label(s) for data arrays to be used

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
  enters the area.
  - mFo-DFc_omit      : omit difference map coefficients
  - PHImFo-DFc_omit   : corresponding phases

Optional output:
  CCP4 files with mask data:
  - mask_all.ccp4    : mask of original model
  - mask_omit.ccp4   : mask when ligand is omitted
  - mask_polder.ccp4 : mask obtained by polder procedure

'''

  datatypes = ['model', 'phil', 'miller_array']

  master_phil_str = master_phil_str
  known_article_ids = ['phenix.polder']

  # ---------------------------------------------------------------------------

  def validate(self):
    print('Validating inputs:\n', file=self.logger)
    self.data_manager.has_models(
      raise_sorry = True,
      expected_n  = 1,
      exact_count = True)
    self.data_manager.has_miller_arrays(raise_sorry=True)
    if (len(self.data_manager.get_miller_array_names()) > 2):
      raise Sorry('Dont input more than 2 reflection files.')
    if (self.params.reflection_file_name is None):
      self.params.reflection_file_name = self.data_manager.get_default_miller_array_name()
    if (self.params.model_file_name is None):
      self.params.model_file_name = self.data_manager.get_default_model_name()

    if (self.params.solvent_exclusion_mask_selection is None):
      raise Sorry('''Selection for atoms to be omitted is required.

  Try something like

    solvent_exclusion_mask_selection=LIG
  ''')
    if (self.params.polder.sphere_radius < 3):
      raise Sorry("Sphere radius out of range: must be larger than 3 A")
    if (self.params.polder.box_buffer is not None and
      (self.params.polder.box_buffer < 0 or self.params.polder.box_buffer > 5)):
      raise Sorry("Box buffer out of range: must be between 0 and 5")

    if (self.params.polder.resolution_factor < 0.0):
      raise Sorry('Use a positive value for the resolution gridding factor.')

  # ---------------------------------------------------------------------------

  def get_crystal_symmetry(self):
    crystal_symmetries = []
    files = self.data_manager.get_model_names() + \
      self.data_manager.get_miller_array_names()
    for f in files:
      cs = crystal_symmetry_from_any.extract_from(f)
      if (cs is not None):
        crystal_symmetries.append(cs)

    if (len(crystal_symmetries) == 1): crystal_symmetry = crystal_symmetries[0]
    elif (len(crystal_symmetries) == 0):
     raise Sorry("No crystal symmetry found.")
    else:
      if (not crystal_symmetries[0].is_similar_symmetry(crystal_symmetries[1])):
        raise Sorry("Crystal symmetry mismatch between different files.")
        # TODO what if 3 symmetries?
      crystal_symmetry = crystal_symmetries[0]
    return crystal_symmetry

  # ---------------------------------------------------------------------------

  def prepare_f_obs_and_flags_if_anomalous(self, f_obs, r_free_flags):
    sel = f_obs.data()>0
    f_obs = f_obs.select(sel)
    merged = f_obs.as_non_anomalous_array().merge_equivalents()
    f_obs = merged.array().set_observation_type(f_obs)
    if r_free_flags:
      r_free_flags = r_free_flags.select(sel)
      merged_free = r_free_flags.as_non_anomalous_array().merge_equivalents()
      r_free_flags = merged_free.array().set_observation_type(r_free_flags)
      f_obs, r_free_flags = f_obs.common_sets(r_free_flags)
    return f_obs, r_free_flags

  # ---------------------------------------------------------------------------

  def get_fobs_rfree(self, crystal_symmetry):
    f_obs, r_free_flags = None, None

    rfs = self.data_manager.get_reflection_file_server(
      filenames = self.data_manager.get_miller_array_names(),
      crystal_symmetry = crystal_symmetry,
      logger=null_out())

    parameters = mmtbx.utils.data_and_flags_master_params().extract()
    if (self.params.data_labels is not None):
      parameters.labels = self.params.data_labels
    if (self.params.r_free_flags_labels is not None):
      parameters.r_free_flags.label = self.params.r_free_flags_labels
    determined_data_and_flags = mmtbx.utils.determine_data_and_flags(
      reflection_file_server = rfs,
      parameters             = parameters,
      keep_going             = True,
      working_point_group = crystal_symmetry.space_group().build_derived_point_group(),
      log                    = null_out(),
      symmetry_safety_check  = True)

    f_obs = determined_data_and_flags.f_obs
    r_free_flags = determined_data_and_flags.r_free_flags
    assert (f_obs is not None)
    if (self.params.data_labels is None):
      self.params.data_labels = f_obs.info().label_string()
    if (r_free_flags is not None):
      self.params.r_free_flags_labels = r_free_flags.info().label_string()

    return f_obs, r_free_flags

  # ---------------------------------------------------------------------------

  def prepare_fobs_rfree(self, f_obs, r_free_flags):
    f_obs = f_obs.resolution_filter(
      d_min = self.params.high_resolution,
      d_max = self.params.low_resolution)
    if (r_free_flags is not None):
      r_free_flags = r_free_flags.resolution_filter(
        d_min = self.params.high_resolution,
        d_max = self.params.low_resolution)

    if (f_obs.anomalous_flag()):
      f_obs, r_free_flags = self.prepare_f_obs_and_flags_if_anomalous(
        f_obs        = f_obs,
        r_free_flags = r_free_flags)

    return f_obs, r_free_flags

  # ---------------------------------------------------------------------------

  def check_selection(self, pdb_hierarchy):
   print("*"*79, file=self.logger)
   print('Selecting atoms...\n', file=self.logger)
   print('Selection string:', self.params.solvent_exclusion_mask_selection)

   selection_bool = pdb_hierarchy.atom_selection_cache().selection(
     string = self.params.solvent_exclusion_mask_selection)
   n_selected = selection_bool.count(True)
   n_selected_all = pdb_hierarchy.atom_selection_cache().selection(
     string = 'all').count(True)
   if(n_selected == 0):
     raise Sorry("No atoms where selected. Check selection syntax again.")
   if (n_selected/n_selected_all > 0.5):
     raise Sorry("""More than half of total number of atoms selected. Omit
       selection should be smaller, such as one ligand or a few residues.""")

   print('Number of atoms selected:', n_selected, file=self.logger)
   pdb_hierarchy_selected = pdb_hierarchy.select(selection_bool)
   ligand_str = pdb_hierarchy_selected.as_pdb_string()
   print(ligand_str, file=self.logger)
   print("*"*79, file=self.logger)

   return selection_bool

  # ---------------------------------------------------------------------------

  def broadcast_rfactors(self, r_work, r_free):
    print('R_work = %6.4f R_free = %6.4f' % (r_work, r_free))
    print ('*'*79, file=self.logger)

  # ---------------------------------------------------------------------------

  def print_rfactors(self, results):
    print ('*'*79, file=self.logger)
    fmodel_input  = results.fmodel_input
    fmodel_biased = results.fmodel_biased
    fmodel_omit   = results.fmodel_omit
    fmodel_polder = results.fmodel_polder
    print('R factors for unmodified input model and data:', file=self.logger)
    self.broadcast_rfactors(fmodel_input.r_work(), fmodel_input.r_free())
    if (self.params.debug):
      print('R factor when ligand is used for mask calculation (biased map):',
        file=self.logger)
      self.broadcast_rfactors(fmodel_biased.r_work(), fmodel_biased.r_free())
    print('R factor for polder map', file=self.logger)
    self.broadcast_rfactors(fmodel_polder.r_work(), fmodel_polder.r_free())
    print('R factor for OMIT map (ligand is excluded for mask calculation):',
      file=self.logger)
    self.broadcast_rfactors(fmodel_omit.r_work(), fmodel_omit.r_free())

  # ---------------------------------------------------------------------------

  def write_files(self, results, f_obs):
    if (self.params.mask_output):
      masks = [results.mask_data_all, results.mask_data_omit, results.mask_data_polder]
      filenames = ["all", "omit", "polder"]
      for mask_data, filename in zip(masks, filenames):
        mrcfile.write_ccp4_map(
          file_name   = "mask_" + filename + ".ccp4",
          unit_cell   = f_obs.unit_cell(),
          space_group = f_obs.space_group(),
          map_data    = mask_data,
          labels      = flex.std_string([""]))
    mtz_dataset = results.mc_polder.as_mtz_dataset(
      column_root_label = "mFo-DFc_polder")
    mtz_dataset.add_miller_array(
      miller_array      = results.mc_omit,
      column_root_label = "mFo-DFc_omit")
    if (self.params.debug):
      mtz_dataset.add_miller_array(
        miller_array      = results.mc_biased,
        column_root_label = "mFo-DFc_bias_omit")
    mtz_object = mtz_dataset.mtz_object()
    polder_file_name = "polder_map_coeffs.mtz"
    if (self.params.output_file_name_prefix is not None):
      polder_file_name = self.params.output_file_name_prefix + "_" + polder_file_name
    mtz_object.write(file_name = polder_file_name)
    print('File %s was written.' % polder_file_name, file=self.logger)

  # ---------------------------------------------------------------------------
  def print_validation(self, results):
    vr = results.validation_results
    print('Map 1: calculated Fobs with ligand')
    print('Map 2: calculated Fobs without ligand')
    print('Map 3: real Fobs data')
    print('CC(1,2): %6.4f' % vr.cc12)
    print('CC(1,3): %6.4f' % vr.cc13)
    print('CC(2,3): %6.4f' % vr.cc23)
    print('Peak CC:')
    print('CC(1,2): %6.4f' % vr.cc12_peak)
    print('CC(1,3): %6.4f' % vr.cc13_peak)
    print('CC(2,3): %6.4f' % vr.cc23_peak)
    print('q    D(1,2) D(1,3) D(2,3)')
    for c,d12_,d13_,d23_ in zip(vr.cutoffs,vr.d12,vr.d13,vr.d23):
      print('%4.2f %6.4f %6.4f %6.4f'%(c,d12_,d13_,d23_))
    ###
    if(self.params.debug):
      #box_1.write_ccp4_map(file_name="box_1_polder.ccp4")
      self.write_map_box(
        box      = vr.box_1,
        filename = "box_1_polder.ccp4")
      self.write_map_box(
        box      = vr.box_2,
        filename = "box_2_polder.ccp4")
      self.write_map_box(
        box      = vr.box_3,
        filename = "box_3_polder.ccp4")
      vr.ph_selected.write_pdb_file(file_name="box_polder.pdb",
        crystal_symmetry=vr.box_1.box_crystal_symmetry)
    #
    print ('*'*79, file=self.logger)
    message = self.result_message(cc12 = vr.cc12, cc13 = vr.cc13, cc23 = vr.cc23)
    print(message, file=self.logger)
    return message

  # ---------------------------------------------------------------------------

  def write_map_box(self, box, filename):
      mrcfile.write_ccp4_map(
        file_name   = filename,
        unit_cell   = box.xray_structure_box.unit_cell(),
        space_group = box.xray_structure_box.space_group(),
        map_data    = box.map_box,
        labels      = flex.std_string([""]))

  # ---------------------------------------------------------------------------

  def result_message(self, cc12, cc13, cc23):
    if (cc13 < 0.7 or
        (cc23 > cc12 and cc23 > cc13) or (cc13 < cc12 and cc13 < cc23)):
      msg = """The polder map is very likely to show bulk-solvent or noise."""
    elif (cc13 >= 0.8):
      msg = 'The polder map is likely to show the omitted atoms.'
    elif (cc13 >= 0.7 and cc13 < 0.8):
      if (cc23 < 0.7*cc13):
        msg = """The polder map is more likely to show the omitted atoms than bulk solvent.
  It is recommended to carefully inspect the maps to confirm."""
      else:
        msg = """The polder map is more likely to show bulk-solvent or noise
  instead of the omitted atoms. But it is recommended to inspect the maps to confirm."""
    return msg

  # ---------------------------------------------------------------------------

  def run(self):

    print('Using model file:', self.params.model_file_name)
    print('Using reflection file(s):', self.data_manager.get_miller_array_names())

    cs = self.get_crystal_symmetry()

    model = self.data_manager.get_model()
    ph = model.get_hierarchy()
    xrs = model.get_xray_structure()
    selection_bool = self.check_selection(pdb_hierarchy = ph)

    f_obs, r_free_flags = self.get_fobs_rfree(crystal_symmetry = cs)
    print('Input data...', file=self.logger)
    print('  Reflection data:', f_obs.info().labels, file=self.logger)
    if (r_free_flags is not None):
      print('  Free-R flags:', r_free_flags.info().labels, file=self.logger)
    else:
      print('  Free-R flags: not present or not found', file=self.logger)
    print('\nWorking crystal symmetry after inspecting all inputs:', file=self.logger)
    cs.show_summary(f=self.logger)

    f_obs, r_free_flags = self.prepare_fobs_rfree(
        f_obs        = f_obs,
        r_free_flags = r_free_flags)

    model_basename = os.path.basename(self.params.model_file_name.split(".")[0])
    if (len(model_basename) > 0 and self.params.output_file_name_prefix is None):
      self.params.output_file_name_prefix = model_basename

    mmtbx.utils.setup_scattering_dictionaries(
      scattering_table = self.params.scattering_table,
      xray_structure   = xrs,
      d_min            = f_obs.d_min())

    polder_object = mmtbx.maps.polder.compute_polder_map(
      f_obs          = f_obs,
      r_free_flags   = r_free_flags,
      model          = model,
      params         = self.params.polder,
      selection_bool = selection_bool)
    polder_object.validate()
    polder_object.run()
    results = polder_object.get_results()

    self.print_rfactors(results = results)
    self.write_files(
      results = results,
      f_obs   = f_obs)
    self.message = None
    if (not self.params.polder.compute_box):
      self.message = self.print_validation(results = results)

    print ('*'*79, file=self.logger)
    print ('Finished', file=self.logger)
    # results object not returned because it contains maps


  def get_results(self):
    return group_args(message=self.message)
