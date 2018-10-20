from __future__ import division, print_function
try:
  from phenix.program_template import ProgramTemplate
except ImportError:
  from libtbx.program_template import ProgramTemplate
import mmtbx.validation.ligands
from mmtbx.validation import validate_ligands
from libtbx.utils import null_out, Sorry
from iotbx import crystal_symmetry_from_any
from libtbx.str_utils import make_sub_header


master_phil_str = '''
ligand_code = None
  .type = str
  .multiple = True
reference_structure = None
  .type = path
only_segid = None
  .type = str
verbose = False
  .type = bool
'''

# =============================================================================

class Program(ProgramTemplate):

  description = '''
mmtbx.validate_ligands model.pdb data.mtz LIGAND_CODE [...]

Print out basic statistics for residue(s) with the given code(s), including
electron density values/CC.
'''

  datatypes = ['model', 'phil', 'miller_array']

  master_phil_str = master_phil_str

  # ---------------------------------------------------------------------------

  def validate(self):
    print('Validating inputs...\n', file=self.logger)
    self.data_manager.has_models(
      raise_sorry = True,
      expected_n  = 1,
      exact_count = True)
#    if (self.params.ligand_code is None or self.params.ligand_code[0] is None):
#      raise Sorry("Ligand code required!")

  # ---------------------------------------------------------------------------

  def get_crystal_symmetry(self):
    crystal_symmetries = []
    files = [self.data_manager.get_default_model_name(),
      self.data_manager.get_default_miller_array_name() ]
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
      crystal_symmetry = crystal_symmetries[0]
    return crystal_symmetry

  # ---------------------------------------------------------------------------

  def get_fobs_rfree(self, crystal_symmetry):
    rfs = self.data_manager.get_reflection_file_server(
      filenames = [self.data_manager.get_default_miller_array_name()],
      crystal_symmetry = crystal_symmetry,
      logger=null_out())
    parameters = mmtbx.utils.data_and_flags_master_params().extract()
    determined_data_and_flags = mmtbx.utils.determine_data_and_flags(
      reflection_file_server = rfs,
      parameters             = parameters,
      keep_going             = True,
      working_point_group = crystal_symmetry.space_group().build_derived_point_group(),
      log                    = null_out(),
      symmetry_safety_check  = True)
    f_obs = determined_data_and_flags.f_obs
    r_free_flags = determined_data_and_flags.r_free_flags
    return f_obs, r_free_flags

  # ---------------------------------------------------------------------------

  def print_results(self, results, ph):
    print('\nThe following ligands were found in the input model:')
    for i, ligand_result in results.items():
      print(ligand_result.resname, ligand_result.id_str)
      #isel = ligand_result.isel
      #for rg in ph.select(isel).residue_groups():
      #  rn = ",".join(rg.unique_resnames())
      #  print(rn, rg.id_str())

    print('\nOccupancies')
    for i, ligand_result in results.items():
      occ = ligand_result.get_occupancies()
      print('Ligand: %s mean: %s' %(ligand_result.resname, occ.mean))

  # ---------------------------------------------------------------------------

  def run(self):

    print('Using model file:', self.data_manager.get_default_model_name())
    print('Using reflection file:', self.data_manager.get_default_miller_array_name())

    cs = self.get_crystal_symmetry()
    model = self.data_manager.get_model()
    ph = model.get_hierarchy()
    xrs = model.get_xray_structure()
    f_obs, r_free_flags = self.get_fobs_rfree(crystal_symmetry = cs)

    print('\nInput data...', file=self.logger)
    print('  Reflection data:', f_obs.info().labels, file=self.logger)
    if (r_free_flags is not None):
      print('  Free-R flags:', r_free_flags.info().labels, file=self.logger)
    else:
      print('  Free-R flags: not present or not found', file=self.logger)
    print('\nWorking crystal symmetry after inspecting all inputs:', file=self.logger)
    cs.show_summary(f=self.logger)

    fmodel = mmtbx.f_model.manager(
     f_obs          = f_obs,
     r_free_flags   = r_free_flags,
     xray_structure = xrs)
    #fmodel.update_all_scales()

    # This is the new class, currently a stub but will be developed
    # winter 2018/spring 2019 by DL and NWM
    results = validate_ligands.manager(model = model)
    results.run()
#    validation_obj = validate_ligands.validate_ligands(model = model)
#    validation_obj.validate_inputs()
#    validation_obj.run()

    self.print_results(
      results = results,
      ph      = ph)

    # TODO
    # DL: Eventually, delete "old" call below, but leave it for now to keep the
    # funcitonality alive, just in case
    if (not(self.params.ligand_code is None or self.params.ligand_code[0] is None)):
      make_sub_header("Validating ligands", out=self.logger)
      for ligand_code in self.params.ligand_code :
        validations = mmtbx.validation.ligands.validate_ligands(
          pdb_hierarchy       = ph,
          fmodel              = fmodel,
          ligand_code         = ligand_code,
          reference_structure = self.params.reference_structure,
          only_segid          = self.params.only_segid)
        if (validations is None) :
          raise Sorry("No ligands named '%s' found." % ligand_code)
        mmtbx.validation.ligands.show_validation_results(validations=validations,
          out     = self.logger,
          verbose = self.params.verbose)
