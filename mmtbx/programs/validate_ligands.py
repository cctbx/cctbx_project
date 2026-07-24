"""Validate ligands in a model"""
from __future__ import absolute_import, division, print_function
import os
import traceback
try:
  from phenix.program_template import ProgramTemplate
except ImportError:
  from libtbx.program_template import ProgramTemplate
import mmtbx.validation.ligands
from mmtbx.validation import validate_ligands
from mmtbx.hydrogens import place_and_optimize_hydrogens
from mmtbx.hydrogens import reduce_hydrogen
import iotbx.pdb
from libtbx.utils import null_out, Sorry
from libtbx.str_utils import make_sub_header
from libtbx import group_args


master_phil_str = """
include scope mmtbx.validation.validate_ligands.master_params
include scope mmtbx.probe.Helpers.probe_phil_parameters
scattering_table = *n_gaussian wk1995 it1992 neutron electron
  .type = choice
  .short_caption = Scattering table
  .help = Scattering table for structure factors calculations
run_reduce2 = True
  .type = bool
save_reduce2_model = False
  .type = bool
save_map_coeffs = False
  .type = bool
  .short_caption = Save map coefficients (MTZ)
  .help = "Write 2mFo-DFc and mFo-DFc map coefficients (from fmodel) to <basename>_map_coeffs.mtz. Requires reflection data."
verbose = False
  .type = bool
gui
  .help = GUI-specific parameters
{
  output_dir = None
    .type = path
    .style = output_dir

  data_column_label = None
    .type = str
    .style = noauto renderer:draw_any_label_widget
    .input_size = 200

  free_column_label = None
    .type = str
    .style = noauto renderer:draw_any_label_widget
    .input_size = 200
}
"""

# =============================================================================

class Program(ProgramTemplate):

  description = '''
mmtbx.development.validate_ligands model.pdb data.mtz
mmtbx.development.validate_ligands model.pdb

Print out basic statistics for residue(s) with the given code(s), including
RSCC.
'''

  datatypes = ['model', 'phil', 'restraint', 'miller_array', 'real_map']

  master_phil_str = master_phil_str
  data_manager_options = ['model_skip_expand_with_mtrix',
                          'model_skip_ss_annotations']

  # ---------------------------------------------------------------------------

  def validate(self):
    print('Validating inputs...\n', file=self.logger)

    # allow only one model
    self.data_manager.has_models(
      raise_sorry = True,
      expected_n  = 1,
      exact_count = True)

  # ---------------------------------------------------------------------------

  def add_hydrogens(self, model):
    '''
    Place H atoms with reduce2
    '''
    make_sub_header('Placing hydrogen atoms with reduce2', out=self.logger)
    try:
      self.working_model = place_and_optimize_hydrogens(
        model           = model,
        keep_existing_H = False,
        probe_phil      = self.params.probe,
        stop_for_unknowns = False,
        raise_on_missing = False,
        log             = self.logger)
      self.working_model.unset_riding_h_manager()
    except Exception:
      msg = traceback.format_exc()
      print('Reduce2 failed.\n' + msg, file=self.logger)
      raise Sorry('Hydrogen placement (reduce2) failed; '
                  'see log for details.')

  # ---------------------------------------------------------------------------

  @staticmethod
  def _remove_element_x(model):
    '''Drop pseudo-atoms with element "X"; they choke pdb_interpretation.'''
    return model.select(~model.selection('element X'))

  # ---------------------------------------------------------------------------

  def check_ligands(self, model):
    make_sub_header('Check if input model has ligands', out=self.logger)
    get_class = iotbx.pdb.common_residue_names_get_class
    exclude = ["common_amino_acid", "modified_amino_acid", "common_rna_dna",
               "modified_rna_dna", "ccp4_mon_lib_rna_dna", "common_water",
                "common_element"]
    self.has_ligands = False
    model = self._remove_element_x(model)
    for chain in model.chains():
      for rg in chain.residue_groups():
        for ag in rg.atom_groups():
          if (get_class(name=ag.resname) in exclude): continue
          print('Found ligand: ', ag.resname, file=self.logger)
          self.has_ligands = True
          mlq, cif_object = reduce_hydrogen.mon_lib_query(
                              residue     = ag,
                              mon_lib_srv = model.get_mon_lib_srv(),
                              raise_sorry = False)

  # ---------------------------------------------------------------------------

  def run(self):
    has_miller = False
    has_map = False
    fmodel = None
    map_manager = None
    self.ligand_manager = None
    self.model_fn_reduce2 = None
    model_fn = self.data_manager.get_default_model_name()
    self._original_model_fn = model_fn
    data_fn = self.data_manager.get_default_miller_array_name()
    map_fn = self.data_manager.get_default_real_map_name()

    print('Using model file:', model_fn, file=self.logger)
    if data_fn is not None:
      print('Using reflection file:', data_fn, file=self.logger)
      has_miller = True
    if map_fn is not None:
      print('Using map file', map_fn, file=self.logger)
      has_map = True

    # get model object from input file
    m = self.data_manager.get_model()
    m.set_log(log = null_out())
    if self.data_manager.has_restraints():
      m.set_stop_for_unknowns(False)
      #m.set_log(log = null_out())
      m.process(make_restraints=False)

    # stop if multi-model file
    if(len(m.get_hierarchy().models())>1):
      raise Sorry('Multi-model files currently not supported.')

    self.check_ligands(model = m)
    if not self.has_ligands:
      print('No ligands found. Exiting.', file=self.logger)
      return

    if self.params.validate_ligands.ligand_code:
      print('\nFocusing on the following ligands only:', file=self.logger)
      for lc in self.params.validate_ligands.ligand_code:
        print('\t', lc, file=self.logger)

    # get rid of element X as it will choke pdb_interpretation
    if ' X' in m.get_hierarchy().atoms().extract_element():
      print('\nFound atoms with element "X" in model. Removing...',
        file=self.logger)
      m = self._remove_element_x(m)

    self.working_model = None

    if self.params.run_reduce2:
      self.add_hydrogens(model = m)
    else:
      self.working_model = m

    if self.working_model is None:
      raise Sorry('Could not create model object.')

    # Register the working model (with reduce2 H placed and element-X atoms
    # removed) so downstream managers use it instead of the original input
    # model. This does not change the DataManager default model.
    # In-memory DataManager key for the working model; never written to disk.
    _model_fn = '%s_validate_ligands_working.pdb' % os.path.splitext(
      os.path.basename(model_fn))[0]
    self.data_manager.add_model(_model_fn, self.working_model)

    if has_map:
      mmm = self.data_manager.get_map_model_manager(
        model_file = _model_fn)
      map_manager = mmm.map_manager()
      self.working_model = mmm.model()
      # Map inputs are cryo-EM; electron scattering is physically correct here
      # regardless of the phil default (scattering_table only applies to the
      # fmodel/reflection path).
      self.working_model.setup_scattering_dictionaries(scattering_table='electron')
      # keep the registered model in sync with the boxed map model
      self.data_manager.add_model(_model_fn, self.working_model)

    ro = self.working_model.get_restraint_objects()
    if ro is None: ro=[]
    ro_no_duplicates = []
    if ro:
      seen = set()
      for name, _ro in ro:
        if name not in seen:
          ro_no_duplicates.append((name, _ro))
          seen.add(name)

    #  print(_ro[1]['comp_list']['_chem_comp.three_letter_code'][0])

    # get fmodel object if reflection data were provided
    if has_miller:
      make_sub_header(' Creating fmodel object ', out=self.logger)
      fmodel_params = self.data_manager.get_fmodel_params()
      fmodel_params.xray_data.r_free_flags.required = False
      fmodel_params.xray_data.r_free_flags.ignore_r_free_flags = True
      self.data_manager.set_fmodel_params(fmodel_params)
      fmodel = self.data_manager.get_fmodel(
        scattering_table = self.params.scattering_table,
        model_filename   = _model_fn)
      print('\n', file = self.logger)
      fmodel.update_all_scales()
      fmodel.show(log=self.logger, show_header=False)
      print ("r_work=%6.4f r_free=%6.4f"%(fmodel.r_work(), fmodel.r_free()),
        file=self.logger)

    #self.data_manager.write_real_map_file(map_manager,filename="my_map.map")
    #self.data_manager.write_model_file(self.working_model,filename="my_model.pdb")

    self.working_model.set_restraint_objects(ro_no_duplicates)
    self.working_model.set_stop_for_unknowns(False)
    pi = self.working_model.get_current_pdb_interpretation_params()
    pi.pdb_interpretation.allow_polymer_cross_special_position = True
    try:
      self.working_model.process(
        make_restraints=True,
        pdb_interpretation_params = pi)
    except Exception as e:
      print(e, file=self.logger)
      print('Could not process model to create restraints.', file=self.logger)
      return


    # split(".")[0] strips any remaining extension after splitext, e.g. the
    # ".ent" left over from a compound extension like "pdb1avd.ent.gz".
    basename = os.path.splitext(os.path.basename(model_fn))[0].split(".")[0]
    self.model_fn_reduce2 = "%s_newH.cif" % basename
    if self.params.save_reduce2_model:
      self.data_manager.set_overwrite(True)
      self.data_manager.write_model_file(self.working_model,filename=self.model_fn_reduce2, format='cif')

    if self.params.save_map_coeffs:
      if fmodel is not None:
        map_coeffs_fn = "%s_map_coeffs.mtz" % basename
        mtz_object = validate_ligands.map_coefficients_as_mtz_object(fmodel)
        self.data_manager.set_overwrite(True)
        self.data_manager.write_miller_array_file(
          mtz_object, filename=map_coeffs_fn)
        print('Wrote map coefficients:', map_coeffs_fn, file=self.logger)
      else:
        print('save_map_coeffs requested but no reflection data were provided; '
              'skipping map output.', file=self.logger)

    ligand_manager = validate_ligands.manager(
      model = self.working_model,
      fmodel = fmodel,
      map_manager = map_manager,
      params = self.params.validate_ligands,
      log   = self.logger)
    ligand_manager.run()
    ligand_manager.show_ligand_counts()
    ligand_manager.show_fragmentation()
    ligand_manager.show_sites_within()
    ligand_manager.show_table(out=self.logger)

    self.ligand_manager = ligand_manager

  # ---------------------------------------------------------------------------

  def get_results(self):
    if (self.params.run_reduce2
        and self.params.save_reduce2_model
        and self.model_fn_reduce2 is not None):
      model_to_open = os.path.abspath(self.model_fn_reduce2)
    else:
      model_to_open = getattr(self, '_original_model_fn', None)
    ligand_results = None
    if self.ligand_manager is not None:
      ligand_results = [lr.as_picklable_snapshot()
                        for lr in self.ligand_manager]
    return group_args(
      working_model_fn = model_to_open,
      ligand_manager   = self.ligand_manager,
      ligand_results   = ligand_results)
