from __future__ import absolute_import, division, print_function
try:
  from phenix.program_template import ProgramTemplate
except ImportError:
  from libtbx.program_template import ProgramTemplate
import mmtbx.validation.ligands
from mmtbx.validation import validate_ligands
from libtbx.utils import null_out, Sorry
from iotbx import crystal_symmetry_from_any
from libtbx.str_utils import make_sub_header
from iotbx import extract_xtal_data


master_phil_str = """
include scope mmtbx.validation.validate_ligands.master_params
ligand_code = None
  .type = str
  .multiple = True
reference_structure = None
  .type = path
only_segid = None
  .type = str

verbose = False
  .type = bool
update_scales = True
  .type = bool

"""
# TODO update_scales if for development only, delete for production!

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
    files = [self.data_manager.get_default_model_name()]
    if self.data_manager.get_default_miller_array_name() is not None:
      files.append(self.data_manager.get_default_miller_array_name())
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
    parameters = extract_xtal_data.data_and_flags_master_params().extract()
    determined_data_and_flags = extract_xtal_data.run(
      reflection_file_server = rfs,
      parameters             = parameters,
      keep_going             = True,
      working_point_group = crystal_symmetry.space_group().build_derived_point_group())
    f_obs = determined_data_and_flags.f_obs
    r_free_flags = determined_data_and_flags.r_free_flags
    return f_obs, r_free_flags

  # ---------------------------------------------------------------------------

  def run(self):

    model_fn = self.data_manager.get_default_model_name()
    print('Using model file:', model_fn, file=self.logger)
    print('Using reflection file:',
      self.data_manager.get_default_miller_array_name(), file=self.logger)

    cs = self.get_crystal_symmetry()
    model = self.data_manager.get_model()
    #grm = model.get_restraints_manager()
    ph = model.get_hierarchy()
    xrs = model.get_xray_structure()

    fmodel = None
    if self.data_manager.get_default_miller_array_name():
      f_obs, r_free_flags = self.get_fobs_rfree(crystal_symmetry = cs)
      print('\nInput data...', file=self.logger)
      print('  Reflection data:', f_obs.info().labels, file=self.logger)
      if (r_free_flags is not None):
        print('  Free-R flags:', r_free_flags.info().labels, file=self.logger)
      else:
        print('  Free-R flags: not present or not found', file=self.logger)
      fmodel = mmtbx.f_model.manager(
       f_obs          = f_obs,
       r_free_flags   = r_free_flags,
       xray_structure = xrs)
      print('\n', file = self.logger)
      fmodel.show(log=self.logger, show_header=False)
      # TODO: delete this keyword for production
      #if self.params.update_scales:
      fmodel.update_all_scales()
      fmodel.show(log=self.logger, show_header=False)

    print('\nWorking crystal symmetry after inspecting all inputs:', file=self.logger)
    cs.show_summary(f=self.logger)

    # This is the new class, currently a stub but will be developed
    # spring 2019 by DL and NWM
    #t0 = time.time()
    # TODO: Decide if H should be placed here or in the class
    # if readyset is used, filename is needed
    # if readyset can be run as class, filename could be avoided
    ligand_manager = validate_ligands.manager(
      model = model,
#      model_fn = model_fn,
      fmodel = fmodel,
      params = self.params.validate_ligands,
      log   = self.logger)
    ligand_manager.run()
    ligand_manager.show_ligand_counts()
    ligand_manager.show_ligand_occupancies()
    ligand_manager.show_adps()
    ligand_manager.show_ccs()
    ligand_manager.show_nonbonded_overlaps()
    #print('time running manager: ', time.time()-t0)

    # TODO
    # DL: Eventually, delete "old" call below, but leave it for now to keep the
    # funcitonality alive, just in case
    if self.params.ligand_code and self.data_manager.get_default_miller_array_name() is not None:
      if (not(self.params.ligand_code is None or self.params.ligand_code[0] is None)):
        make_sub_header("Validating ligands", out=self.logger)
        for ligand_code in self.params.ligand_code :
          validations = mmtbx.validation.ligands.validate_ligands(
            pdb_hierarchy       = ph,
            fmodel              = fmodel,
            ligand_code         = ligand_code,
            reference_structure = self.params.reference_structure,
            only_segid          = self.params.only_segid)
          if (validations is None):
            raise Sorry("No ligands named '%s' found." % ligand_code)
          mmtbx.validation.ligands.show_validation_results(validations=validations,
            out     = self.logger,
            verbose = self.params.verbose)
