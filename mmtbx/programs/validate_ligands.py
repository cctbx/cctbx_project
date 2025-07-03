"""Validate ligands in a model"""
from __future__ import absolute_import, division, print_function
import os
try:
  from phenix.program_template import ProgramTemplate
except ImportError:
  from libtbx.program_template import ProgramTemplate
import mmtbx.validation.ligands
from mmtbx.validation import validate_ligands
from libtbx.utils import null_out, Sorry
from iotbx import crystal_symmetry_from_any
from libtbx.str_utils import make_sub_header
from libtbx import group_args


master_phil_str = """
include scope mmtbx.validation.validate_ligands.master_params
ligand_code = None
  .type = str
  .multiple = True
scattering_table = *n_gaussian wk1995 it1992 neutron electron
  .type = choice
  .short_caption = Scattering table
  .help = Scattering table for structure factors calculations
run_reduce2 = True
  .type = bool
verbose = False
  .type = bool
"""

# old params from Nat
#reference_structure = None
#  .type = path
#only_segid = None
#  .type = str

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

  # ---------------------------------------------------------------------------

  def get_crystal_symmetry(self, model_fn, data_fn):
    crystal_symmetries = []
    for f in [model_fn, data_fn]:
      if f is None: continue
      cs = crystal_symmetry_from_any.extract_from(f)
      if (cs is not None):
        crystal_symmetries.append(cs)
    if (len(crystal_symmetries) == 0):
     raise Sorry("No crystal symmetry found.")
    elif (len(crystal_symmetries) > 1):
      if (not crystal_symmetries[0].is_similar_symmetry(crystal_symmetries[1])):
        raise Sorry("Crystal symmetry mismatch between different files.")
    crystal_symmetry = crystal_symmetries[0]
    return crystal_symmetry

  # ---------------------------------------------------------------------------

  def add_hydrogens(self):
    '''
    Place H atoms with reduce2
    '''
    make_sub_header(' Placing H with reduce2 ', out=self.logger)
    model_reduce2 = None
    model_fn = self.data_manager.get_default_model_name()
    basename = os.path.splitext(os.path.basename(model_fn))[0]
    model_fn_reduce2 = "%s_newH.cif" % basename
    from iotbx.cli_parser import run_program
    from mmtbx.programs import reduce2 as reduce2
    args=["overwrite=True",
          "%s" % model_fn,
          #"use_neutron_distances=True",
          "output.filename=%s" % model_fn_reduce2]
    print("mmtbx.reduce2 %s" %(" ".join(args)), file=self.logger)
    try:
      result = run_program(program_class=reduce2.Program,args=args,
       logger = null_out())
      model_reduce2 = result.model
    except Exception as e:
      msg = traceback.format_exc()
      self.success   = False
      self.write_log(step = 'Reduce2', msg  = msg)
      print('Reduce2 failed.\n' + msg, file=self.logger)
    self.data_manager.add_model(model_fn_reduce2, model_reduce2)
    self.working_model_fn = model_fn_reduce2
    self.working_model = model_reduce2


  # ---------------------------------------------------------------------------

  def run(self):
    has_data = False
    fmodel = None
    model_fn = self.data_manager.get_default_model_name()
    data_fn = self.data_manager.get_default_miller_array_name()
    print('Using model file:', model_fn, file=self.logger)
    if data_fn is not None:
      print('Using reflection file:', data_fn, file=self.logger)
      has_data = True

    if self.params.run_reduce2:
      self.add_hydrogens()
    else:
      self.working_model_fn = model_fn
      m = self.data_manager.get_model()
      m.process(make_restraints=True)
      self.working_model = m

    # get fmodel object
    if has_data:
      make_sub_header(' Creating fmodel object ', out=self.logger)
      fmodel = self.data_manager.get_fmodel(
        scattering_table = self.params.scattering_table,
        model_filename   = self.working_model_fn)
      print('\n', file = self.logger)
      fmodel.update_all_scales()
      fmodel.show(log=self.logger, show_header=False)

#
#    cs = self.get_crystal_symmetry(model_fn, data_fn)
#    model = self.data_manager.get_model()
#    print('\nWorking crystal symmetry after inspecting all inputs:', file=self.logger)
#    cs.show_summary(f=self.logger)

    #t0 = time.time()
    ligand_manager = validate_ligands.manager(
      model = self.working_model,
      fmodel = fmodel,
      params = self.params.validate_ligands,
      log   = self.logger)
    ligand_manager.run()
    ligand_manager.show_ligand_counts()
    ligand_manager.show_ligand_occupancies()
    ligand_manager.show_adps()
    ligand_manager.show_ccs()
    ligand_manager.show_nonbonded_overlaps()

    self.ligand_manager = ligand_manager
    #print('time running manager: ', time.time()-t0)

    # old class from Nat
    #if self.params.ligand_code and self.data_manager.get_default_miller_array_name() is not None:
    #  if (not(self.params.ligand_code is None or self.params.ligand_code[0] is None)):
    #    make_sub_header("Validating ligands", out=self.logger)
    #    for ligand_code in self.params.ligand_code :
    #      validations = mmtbx.validation.ligands.validate_ligands(
    #        pdb_hierarchy       = ph,
    #        fmodel              = fmodel,
    #        ligand_code         = ligand_code,
    #        reference_structure = self.params.reference_structure,
    #        only_segid          = self.params.only_segid)
    #      if (validations is None):
    #        raise Sorry("No ligands named '%s' found." % ligand_code)
    #      mmtbx.validation.ligands.show_validation_results(validations=validations,
    #        out     = self.logger,
    #        verbose = self.params.verbose)

  def get_results(self):
    return group_args(
      ligand_manager     = self.ligand_manager)
