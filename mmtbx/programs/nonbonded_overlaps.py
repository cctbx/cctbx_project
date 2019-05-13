from __future__ import division, print_function
try:
  from phenix.program_template import ProgramTemplate
except ImportError:
  from libtbx.program_template import ProgramTemplate
import cctbx.geometry_restraints.nonbonded_overlaps as nbo
import mmtbx.validation.clashscore as mvc
from libtbx.utils import null_out
import mmtbx.model

master_phil_str = """
  model = None
    .type = path
    .optional = False
    .help = '''input PDB file'''

  cif = None
    .type = path
    .optional = True
    .help = '''Optional Crystallographic Information File (CIF)'''

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

  time_limit = 120
    .type = int
    .help = '''Time limit (sec) for Reduce optimization'''

  show_overlap_type = *all sym macro selection
  .type = choice(multi=False)
  .help = '''When using cctbx method, this parameter allows selecting to show
  all clashes 'all', clashes dues to symmetry operation 'sym' or clashes in
  the macro molecule 'macro'.'''

  substitute_non_crystallographic_unit_cell_if_necessary = False
    .type = bool
    .help = '''\
    Will replace the crystallographic unit cell when the model
    crystallographic information is bad'''

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

  model=input_file          input PDB file
  cif=input_file            input CIF file for additional model information
  keep_hydrogens=True       keep input hydrogen files (otherwise regenerate)
  skip_hydrogen_test=False  Ignore hydrogen considerations,
                            check NBO on PDB file as-is
  nuclear=False             use nuclear x-H distances and vdW radii
  verbose=True              verbose text output
  time_limit=120            Time limit (sec) for Reduce optimization
  show_overlap_type=all     what type of overlaps to show (all,sym or macro)
  show_normalized_nbo=False Show non-bonded overlaps per 1000 atoms
  substitute_non_crystallographic_unit_cell_if_necessary=false
                            fix CRYST1 records if needed

Example:

>>> mmtbx.nonbonded_overlaps xxxx.pdb keep_hydrogens=True

>>> mmtbx.nonbonded_overlaps xxxx.pdb verbose=false
'''

  # TODO: read in cif files
  datatypes = ['model', 'phil']

  master_phil_str = master_phil_str

  # ---------------------------------------------------------------------------

  def validate(self):
    print('Validating inputs:\n', file=self.logger)
    self.data_manager.has_models(
      raise_sorry = True,
      expected_n  = 1,
      exact_count = True)
    if (self.params.model is None):
      self.params.model = self.data_manager.get_default_model_name()

  # ---------------------------------------------------------------------------

  def run(self):

    print('Using model file:', self.params.model, file=self.logger)
    model = self.data_manager.get_model()

    pi_params = model.get_default_pdb_interpretation_params()

    # TODO: below are some non-defaults for pdb_interpretation,
    # but they cannot be changed once model is created.
    # Do we need the non-defaults?
    #print(pi_params.pdb_interpretation.assume_hydrogens_all_missing) default is False, but becomes True if H are present
    #print(pi_params.pdb_interpretation.hard_minimum_nonbonded_distance)
    pi_params.pdb_interpretation.clash_guard.nonbonded_distance_threshold = None
    #print(pi_params.pdb_interpretation.substitute_non_crystallographic_unit_cell_if_necessary)


    # TODO: This should be replaced with readyset in the future
    # add H atoms with reduce
    pdb_fn = self.params.model
    if not self.params.skip_hydrogen_test:
      pdb_with_h, h_were_added = mvc.check_and_add_hydrogen(
          file_name      = pdb_fn,
          model_number   = 0,
          nuclear        = self.params.nuclear,
          verbose        = self.params.verbose,
          time_limit     = self.params.time_limit,
          keep_hydrogens = self.params.keep_hydrogens,
          allow_multiple_models = False,
          log            = self.logger)
      if h_were_added:
        pdb_fn = pdb_fn.replace('.pdb','_with_h.pdb')
        open(pdb_fn,'w').write(pdb_with_h)
        pdb_inp = iotbx.pdb.input(file_name=pdb_fn)
        model = mmtbx.model.manager(
          model_input = pdb_inp,
          #restraint_objects = cif_objects,
          pdb_interpretation_params = pi_params,
          log         = null_out())

    model.set_pdb_interpretation_params(pi_params)
    geometry = model.get_restraints_manager().geometry

    xrs = model.get_xray_structure()
    sites_cart = model.get_sites_cart()
    site_labels = xrs.scatterers().extract_labels()
    hd_sel = model.get_hd_selection()

    # TODO: do we need macro_mol_sel, do we care?
    proxies = model.all_chain_proxies
    cache = proxies.pdb_hierarchy.atom_selection_cache()
    macro_mol_sel = proxies.selection(
      cache  = cache,
      string = 'protein or dna or rna')

    # TODO replace input parameters with model object
    nb_overlaps = nbo.info(
      geometry_restraints_manager=geometry,
      macro_molecule_selection=macro_mol_sel,
      sites_cart=sites_cart,
      site_labels=site_labels,
      hd_sel=hd_sel)

    if self.params.verbose:
      nb_overlaps.show(
        log=self.logger,
        nbo_type=self.params.show_overlap_type,
        normalized_nbo=self.params.show_normalized_nbo)

