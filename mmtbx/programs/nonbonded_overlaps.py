from __future__ import division, print_function
try:
  from phenix.program_template import ProgramTemplate
except ImportError:
  from libtbx.program_template import ProgramTemplate
import cctbx.geometry_restraints.nonbonded_overlaps as nbo
import cctbx.geometry_restraints.process_nonbonded_proxies as pnp
import mmtbx.validation.clashscore as mvc
from libtbx.utils import null_out
import mmtbx.model
import iotbx.pdb
#import time

master_phil_str = """
  verbose = True
    .type = bool

  keep_hydrogens = True
    .type = bool
    .help = '''Keep hydrogens in input file
    (if there are no hydrogens in input file they will be added)'''

  skip_hydrogen_test = False
    .type = bool
    .help = '''Ignore hydrogen considerations, check NBO on PDB file as-is '''

  nuclear = False
    .type = bool
    .help = '''Use nuclear hydrogen positions'''

  show_overlap_type = *all sym macro selection
  .type = choice(multi=False)
  .help = '''When using cctbx method, this parameter allows selecting to show
  all clashes 'all', clashes dues to symmetry operation 'sym' or clashes in
  the macro molecule 'macro'.'''

  show_normalized_nbo = False
    .type = bool
    .help = When True, will show non-bonded overlaps per 1000 atoms

  show_non_binary_overlap_values = True
    .type = bool
    .help = use a function
"""

# =============================================================================

class Program(ProgramTemplate):

  description = '''
phenix.nonbonded_overlaps file.pdb [params.eff] [options ...]

Options:

  keep_hydrogens=True       keep input hydrogen files (otherwise regenerate)
  skip_hydrogen_test=False  Ignore hydrogen considerations,
                            check NBO on model file as-is
  nuclear=False             use nuclear x-H distances and vdW radii
  verbose=True              verbose text output
  show_overlap_type=all     what type of overlaps to show (all, sym or macro)
  show_normalized_nbo=False Show non-bonded overlaps per 1000 atoms


Example:

>>> mmtbx.nonbonded_overlaps xxxx.pdb keep_hydrogens=True

>>> mmtbx.nonbonded_overlaps xxxx.pdb verbose=false

>>> mmtbx.nonbonded_overlaps xxxx.pdb yyyyy.cif
'''

  datatypes = ['model', 'phil', 'restraint']

  master_phil_str = master_phil_str

  # ---------------------------------------------------------------------------

  def validate(self):
    print('Validating inputs:\n', file=self.logger)
    self.data_manager.has_models(
      raise_sorry = True,
      expected_n  = 1,
      exact_count = True)

  # ---------------------------------------------------------------------------

  def run(self):
    pdb_fn = self.data_manager.get_default_model_name()
    print('Using model file:', pdb_fn, file=self.logger)
    if self.data_manager.has_restraints():
      print('Using restraints files', self.data_manager.get_restraint_names(),
        file=self.logger)

    model = self.data_manager.get_model()

    pi_params = model.get_default_pdb_interpretation_params()
    pi_params.pdb_interpretation.clash_guard.nonbonded_distance_threshold = None

    # TODO: below are some non-defaults for pdb_interpretation,
    # but they cannot be changed
    # Do we need the non-defaults?
    # assume_hydrogens_all_missing --> default is False, but becomes True if H are present
    # hard_minimum_nonbonded_distance --> default is 0.001 but is 0 in NBO
    # substitute_non_crystallographic_unit_cell_if_necessary --> not necessary

    # add H atoms with reduce
    # TODO: This should be replaced with readyset in the future
    if not self.params.skip_hydrogen_test:
      pdb_with_h, h_were_added = mvc.check_and_add_hydrogen(
          file_name      = pdb_fn,
          model_number   = 0,
          nuclear        = self.params.nuclear,
          verbose        = self.params.verbose,
          keep_hydrogens = self.params.keep_hydrogens,
          allow_multiple_models = False,
          log            = self.logger)
      if h_were_added:
        pdb_fn = pdb_fn.replace('.pdb','_with_h.pdb')
        open(pdb_fn,'w').write(pdb_with_h)
        pdb_inp = iotbx.pdb.input(file_name=pdb_fn)

        model = mmtbx.model.manager(
          model_input = pdb_inp,
          pdb_interpretation_params = pi_params,
          stop_for_unknowns = False,
          log         = null_out())
        if self.data_manager.has_restraints():
          restraint_objects = list()
          for filename in self.data_manager.get_restraint_names():
            restraint_objects.append((filename, self.data_manager.get_restraint(filename)))
          model.set_restraint_objects(restraint_objects)

    model.set_pdb_interpretation_params(pi_params)
    model.get_restraints_manager()

    # TODO: do we need macro_mol_sel, do we care?
    # If we use model.select(), we don't need it.
    proxies = model.all_chain_proxies
    cache = proxies.pdb_hierarchy.atom_selection_cache()
    macro_mol_sel = proxies.selection(
      cache  = cache,
      string = 'protein or dna or rna')

    #t0 = time.time()
    nb_overlaps = nbo.info(
      model = model,
      macro_molecule_selection=macro_mol_sel)
    #t1 = time.time()


    if self.params.verbose:
      nb_overlaps.show(
        log=self.logger,
        nbo_type=self.params.show_overlap_type,
        normalized_nbo=self.params.show_normalized_nbo)

    #t2 = time.time()
    processed_nbps = pnp.manager(model = model)
    clashes = processed_nbps.get_clashes()
    #t3 = time.time()
    clashes.show()

    #print("OLD time: %8.3f"%(t1-t0))
    #print("NEW time: %8.3f"%(t3-t2))

