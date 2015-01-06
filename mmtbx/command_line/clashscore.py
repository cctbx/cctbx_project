# LIBTBX_SET_DISPATCHER_NAME phenix.clashscore

from __future__ import division
import mmtbx.monomer_library.pdb_interpretation as pdb_inter
import cctbx.geometry_restraints.clash_score as clash_score
import mmtbx.validation.clashscore
from libtbx.utils import null_out
from libtbx.utils import Usage
import iotbx.phil
import sys

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
    .help = '''Keep hydrogens in input file'''

  nuclear = False
    .type = bool
    .help = '''Use nuclear hydrogen positions'''

  time_limit = 120
    .type = int
    .help = '''Time limit (sec) for Reduce optimization'''

  b_factor_cutoff = None
    .type = int
    .help = '''B factor cutoff for use with MolProbity'''

  method = *molprobity cctbx
    .type = choice(multi=False)
    .help = '''Select molprobity or ccbtx nonbonded clashscore method'''

  hard_minimum_nonbonded_distance = 0.0
    .type = float
    .help = '''Clash guard minimum distance'''

  nonbonded_distance_threshold = None
    .type = float

  substitute_non_crystallographic_unit_cell_if_necessary = True
    .type = bool

  assume_hydrogens_all_missing = True
    .type = bool
    .help = '''geometry restraints manager assume that model does not
    contains hydrogens.'''
"""

usage_string = """\
phenix.clashscore file.pdb [params.eff] [options ...]

Options:

  model=input_file          input PDB file
  method=molprobity         Non-bonded clashscore method: molprobity or cctbx
Option for molprobity clashscore:
  keep_hydrogens=True       keep input hydrogen files (otherwise regenerate)
  nuclear=False             use nuclear x-H distances and vdW radii
  verbose=True              verbose text output
  b_factor_cutoff=40        B factor cutoff for clash analysis
Option for cctbx clashscore:
  cif=input_file            optional cif file
  assume_hydrogens_all_missing=True
                            geometry restraints manager assume that model
                            does not contains hydrogens.

Example:

  phenix.clashscore model=1ubq.pdb keep_hydrogens=True

"""

def run (args, out=sys.stdout, quiet=None) :
  """
  Calculates nonbonded clashscore using MolProbity (PROBE) or CCTBX
  """
  cmdline = iotbx.phil.process_command_line_with_files(
    args=args,
    master_phil_string=master_phil_str,
    pdb_file_def="model",
    cif_file_def="cif",
    usage_string=usage_string)
  params = cmdline.work.extract()

  if (params.model is None) :
    raise Usage(usage_string)
  if params.method == 'molprobity':
    pdb_in = cmdline.get_file(params.model, force_type="pdb")
    hierarchy = pdb_in.file_object.hierarchy
    result = mmtbx.validation.clashscore.clashscore(
      pdb_hierarchy=hierarchy,
      keep_hydrogens=params.keep_hydrogens,
      nuclear=params.nuclear,
      out=out,
      verbose=params.verbose and not quiet,
      b_factor_cutoff=params.b_factor_cutoff)
    if params.verbose:
      result.show_old_output(out=out)
  else:
    assert params.method == 'cctbx'
    pdb_file_name = [x for x in args if x.endswith('.pdb')]
    cif_file_name = [x for x in args if x.endswith('.cif')]
    assert pdb_file_name
    files = [pdb_file_name[0]]
    if cif_file_name:
        files.append(cif_file_name[0])

    pdb_processed_file = pdb_inter.run(
      args=files,
      assume_hydrogens_all_missing=params.assume_hydrogens_all_missing,
      hard_minimum_nonbonded_distance=params.hard_minimum_nonbonded_distance,
      nonbonded_distance_threshold=params.nonbonded_distance_threshold,
      substitute_non_crystallographic_unit_cell_if_necessary
      =params.substitute_non_crystallographic_unit_cell_if_necessary,
      log=null_out())

    grm = pdb_processed_file.geometry_restraints_manager()
    xrs = pdb_processed_file.xray_structure()
    sites_cart = xrs.sites_cart()
    site_labels = xrs.scatterers().extract_labels()
    hd_sel = xrs.hd_selection()

    clash_score.info(
      geometry_restraints_manager=grm,
      sites_cart=sites_cart,
      site_labels=site_labels,
      hd_sel=hd_sel).show(log=out)

if (__name__ == "__main__") :
  run(sys.argv[1:])
