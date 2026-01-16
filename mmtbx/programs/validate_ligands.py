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
from mmtbx.hydrogens import reduce_hydrogen
import iotbx.pdb
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
  data_manager_options = ['model_skip_expand_with_mtrix',
                          'model_skip_ss_annotations']

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

  def add_hydrogens(self, model_fn):
    '''
    Place H atoms with reduce2
    '''
    make_sub_header(' Placing H with reduce2 ', out=self.logger)
    model_reduce2 = None
    basename = os.path.splitext(os.path.basename(model_fn))[0]
    model_fn_reduce2 = "%s_newH.cif" % basename.split(".")[0]
    from iotbx.cli_parser import run_program
    from mmtbx.programs import reduce2 as reduce2
    args=["overwrite=True",
          "%s" % model_fn,
          "ignore_missing_restraints=True",
          #"use_neutron_distances=True",
          "output.filename=%s" % model_fn_reduce2]
    print("mmtbx.reduce2 %s" %(" ".join(args)), file=self.logger)
    try:
      result = run_program(program_class=reduce2.Program,args=args,
       logger = null_out())
      #model_reduce2 = self.data_manager.get_model(model_fn_reduce2)
      model_reduce2 = result.model
      model_reduce2.unset_riding_h_manager()
    except Exception as e:
      msg = traceback.format_exc()
      print('Reduce2 failed.\n' + msg, file=self.logger)
      return

    self.data_manager.add_model(model_fn_reduce2, model_reduce2)
    self.working_model_fn = model_fn_reduce2
    self.working_model = model_reduce2

  # ---------------------------------------------------------------------------

  def check_ligands(self, model):
    make_sub_header('Check if input model has ligands', out=self.logger)
    get_class = iotbx.pdb.common_residue_names_get_class
    exclude = ["common_amino_acid", "modified_amino_acid", "common_rna_dna",
               "modified_rna_dna", "ccp4_mon_lib_rna_dna", "common_water",
                "common_element"]
    self.has_ligands = False
    model = model.select(~model.selection('element X'))
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
          if cif_object:
            self.additional_ro.append(('auto_%s' % ag.resname, cif_object))
          elif mlq is None:
            print('\tNo restraints available for %s. No H atoms will be added.'
              % ag.resname, file=self.logger)

  # ---------------------------------------------------------------------------

  def run(self):
    has_data = False
    fmodel = None
    self.additional_ro = []
    self.working_model_fn = None
    self.ligand_manager = None
    model_fn = self.data_manager.get_default_model_name()
    data_fn = self.data_manager.get_default_miller_array_name()
    print('Using model file:', model_fn, file=self.logger)
    if data_fn is not None:
      print('Using reflection file:', data_fn, file=self.logger)
      has_data = True

    m = self.data_manager.get_model()
    #
    if(len(m.get_hierarchy().models())>1):
      raise Sorry('Multi-model files currently not supported.')
    #
    self.check_ligands(model = m)
    if not self.has_ligands:
      print('No ligands found. Exiting.', file=self.logger)
      return
    #
    if ' X' in m.get_hierarchy().atoms().extract_element():
      m = m.select(~m.selection('element X'))
      basename = os.path.splitext(os.path.basename(model_fn))[0]
      model_fn = "%s_noX.cif" % basename
      self.data_manager.write_model_file(
        model_str = m.model_as_mmcif(),
        filename  = model_fn,
        overwrite = True)

    self.working_model = None
    if self.params.run_reduce2:
      self.add_hydrogens(model_fn = model_fn)
    else:
      self.working_model_fn = model_fn
      m = self.data_manager.get_model()
      self.working_model = m

    if self.working_model is None:
      raise Sorry('Could not create model object')
    self.working_model.set_log(log = null_out())
    if self.additional_ro:
      ro = self.working_model.get_restraint_objects()
      if ro is None: ro=[]
      ro.extend(self.additional_ro)
      self.working_model.set_restraint_objects(ro)
    self.working_model.set_stop_for_unknowns(False)
    pi = self.working_model.get_current_pdb_interpretation_params()
    pi.pdb_interpretation.allow_polymer_cross_special_position = True
    try:
      self.working_model.process(
        make_restraints=True,
        pdb_interpretation_params = pi)
    except Exception as e:
      print(e, file=self.logger)
      print('Could not process model to create restraints.')
      return

    # get fmodel object
    if has_data:
      make_sub_header(' Creating fmodel object ', out=self.logger)
      fmodel_params = self.data_manager.get_fmodel_params()
      fmodel_params.xray_data.r_free_flags.required = False
      fmodel_params.xray_data.r_free_flags.ignore_r_free_flags = True
      self.data_manager.set_fmodel_params(fmodel_params)
      fmodel = self.data_manager.get_fmodel(
        scattering_table = self.params.scattering_table,
        model_filename   = self.working_model_fn)
      print('\n', file = self.logger)
      fmodel.update_all_scales()
      fmodel.show(log=self.logger, show_header=False)
      print ("r_work=%6.4f r_free=%6.4f"%(fmodel.r_work(), fmodel.r_free()),
        file=self.logger)

    #t0 = time.time()
    ligand_manager = validate_ligands.manager(
      model = self.working_model,
      fmodel = fmodel,
      params = self.params.validate_ligands,
      log   = self.logger)
    ligand_manager.run()
    ligand_manager.show_ligand_counts()
    ligand_manager.show_sites_within()
    ligand_manager.show_table(out=self.logger)
    # ligand_manager.show_ligand_occupancies()
    # ligand_manager.show_adps()
    # ligand_manager.show_ccs()
    # ligand_manager.show_nonbonded_overlaps()

    self.ligand_manager = ligand_manager
    #print('time running manager: ', time.time()-t0)

  # ---------------------------------------------------------------------------

  def get_results(self):
    return group_args(
      working_model_fn = self.working_model_fn,
      ligand_manager   = self.ligand_manager)
