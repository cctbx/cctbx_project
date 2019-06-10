from __future__ import division, print_function
try:
  from phenix.program_template import ProgramTemplate
except ImportError:
  from libtbx.program_template import ProgramTemplate
import libtbx.load_env
#import time
import cctbx.geometry_restraints.nonbonded_overlaps as nbo
import cctbx.geometry_restraints.process_nonbonded_proxies as pnp
from elbow.command_line.ready_set import model_interface as ready_set_model_interface


master_phil_str = """

  keep_temp = False
    .type = bool
    .help = '''Keep temporary readyset folder'''

  add_hydrogen = True
    .type = bool
    .help = '''Add H atoms to input model.'''

  nuclear = False
    .type = bool
    .help = '''Use nuclear hydrogen positions'''
"""

# =============================================================================

class Program(ProgramTemplate):

  description = '''
phenix.nonbonded_overlaps file.pdb [params.eff] [options ...]

Options:

  keep_hydrogens=True       keep input hydrogen files (otherwise regenerate)
  add_hydrogen=True         Place H atoms with ready_set
  nuclear=False             use nuclear x-H distances and vdW radii

Example:

>>> mmtbx.nonbonded_overlaps xxxx.pdb keep_hydrogens=True

>>> mmtbx.nonbonded_overlaps xxxx.pdb add_hydrogen=False

>>> mmtbx.nonbonded_overlaps xxxx.pdb yyyyy.cif
'''

  datatypes = ['model', 'phil', 'restraint']

  master_phil_str = master_phil_str


  def validate(self):
    print('Validating inputs:\n', file=self.logger)
    self.data_manager.has_models(
      raise_sorry = True,
      expected_n  = 1,
      exact_count = True)


  def run(self):
    pdb_fn = self.data_manager.get_default_model_name()
    print('Using model file:', pdb_fn, file=self.logger)
    restraint_objects = list()
    if self.data_manager.has_restraints():
      print('Using restraints files', self.data_manager.get_restraint_names(),
        file=self.logger)
      for filename in self.data_manager.get_restraint_names():
        restraint_objects.append((filename, self.data_manager.get_restraint(filename)))

    model = self.data_manager.get_model()

    pi_params = model.get_default_pdb_interpretation_params()
    pi_params.pdb_interpretation.clash_guard.nonbonded_distance_threshold = None

    # TODO: below are some non-defaults for pdb_interpretation,
    # but they cannot be changed
    # Do we need the non-defaults?
    # assume_hydrogens_all_missing --> default is False, but becomes True if H are present
    # hard_minimum_nonbonded_distance --> default is 0.001 but is 0 in NBO
    # substitute_non_crystallographic_unit_cell_if_necessary --> not necessary

    # add H atoms with readyset
    if self.params.add_hydrogen:
      params=["add_h_to_water=False",
              "optimise_final_geometry_of_hydrogens=False"]
      assert (libtbx.env.has_module(name="reduce"))
      assert (libtbx.env.has_module(name="elbow"))
      readyset_model = ready_set_model_interface(
            model  = model,
            params = params,
            keep_temp = self.params.keep_temp)
    else:
      readyset_model = model

    readyset_model.set_pdb_interpretation_params(pi_params)
    readyset_model.set_restraint_objects(restraint_objects)
    readyset_model.get_restraints_manager()

    # TODO: do we need macro_mol_sel, do we care?
    # If we use model.select(), we don't need it.
    proxies = readyset_model.all_chain_proxies
    cache = proxies.pdb_hierarchy.atom_selection_cache()
    macro_mol_sel = proxies.selection(
      cache  = cache,
      string = 'protein or dna or rna')

    #t0 = time.time()
    nb_overlaps = nbo.info(
      model = readyset_model,
      macro_molecule_selection=macro_mol_sel)
    #t1 = time.time()

    nb_overlaps.show(
      log=self.logger,
      nbo_type='all',
      normalized_nbo=True)

    #t2 = time.time()
    processed_nbps = pnp.manager(model = readyset_model)
    clashes = processed_nbps.get_clashes()
    #t3 = time.time()
    clashes.show(log=self.logger)

    hbonds = processed_nbps.get_hbonds()
    hbonds.show(log=self.logger)

    #print("OLD time: %8.3f"%(t1-t0))
    #print("NEW time: %8.3f"%(t3-t2))
