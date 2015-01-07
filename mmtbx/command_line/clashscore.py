# LIBTBX_SET_DISPATCHER_NAME phenix.clashscore

from __future__ import division
from cctbx.geometry_restraints.clash_score import check_and_add_hydrogen
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
"""

usage_string = """\
phenix.clashscore file.pdb [params.eff] [options ...]

Options:

  model=input_file          input PDB file
  method=molprobity         Non-bonded clashscore method: molprobity or cctbx
  keep_hydrogens=True       keep input hydrogen files (otherwise regenerate)
  nuclear=False             use nuclear x-H distances and vdW radii
  verbose=True              verbose text output
  b_factor_cutoff=40        B factor cutoff for clash analysis
  cif=input_file            optional cif file (if needed for cctbx clashscore)

Example:

  phenix.clashscore model=1ubq.pdb keep_hydrogens=True

"""

def run (args, out=sys.stdout, quiet=None) :
  """
  Calculates nonbonded clashscore using MolProbity (PROBE) or CCTBX

  Returns:
    When verbose=True the function print detailed results to log
    When verbose=False it will print:
      for molprobity (float): clashscore
      for cctbx (list of floats):
        [simple_cctbx_clashscore (without symmetry and solvent),
         symmetry_cctbx_clashscore,
         solvent_cctbx_clashscore,
         total_cctbx_clashscore]
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
      print >> out,round(result.get_clashscore(),2)
  else:
    assert params.method == 'cctbx'
    pdb_file_name = [x for x in args if x.endswith('.pdb')]
    cif_file_name = [x for x in args if x.endswith('.cif')]
    assert pdb_file_name
    pdb_file_name = pdb_file_name[0]
    pdb_with_h, h_were_added = check_and_add_hydrogen(
        file_name=pdb_file_name,
        model_number=0,
        nuclear=params.nuclear,
        verbose=params.verbose,
        time_limit=params.time_limit,
        keep_hydrogens=params.keep_hydrogens,
        log=out)
    if h_were_added:
      pdb_file_name = pdb_file_name.replace('.pdb','_with_h.pdb')
      open(pdb_file_name,'w').write(pdb_with_h)
    files = [pdb_file_name]
    if cif_file_name:
        files.append(cif_file_name[0])

    pdb_processed_file = pdb_inter.run(
      args=files,
      assume_hydrogens_all_missing=False,
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

    cctbx_clashscore = clash_score.info(
      geometry_restraints_manager=grm,
      sites_cart=sites_cart,
      site_labels=site_labels,
      hd_sel=hd_sel)
    if params.verbose:
      cctbx_clashscore.show(log=out)
    else:
      all = cctbx_clashscore.result.nb_clashscore_all_clashes
      simple = cctbx_clashscore.result.nb_clashscore_simple
      sym = cctbx_clashscore.result.nb_clashscore_due_to_sym_op
      solv = cctbx_clashscore.result.nb_clashscore_solvent_solvent
      out_list = [simple,sym,solv,all]
      print >> out,map(lambda x: round(x,2),out_list)

if (__name__ == "__main__") :
  run(sys.argv[1:])
