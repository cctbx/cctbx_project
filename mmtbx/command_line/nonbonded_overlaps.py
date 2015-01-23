# LIBTBX_SET_DISPATCHER_NAME mmtbx.nonbonded_overlaps

from __future__ import division
import mmtbx.monomer_library.pdb_interpretation as pdb_inter
import cctbx.geometry_restraints.nonbonded_overlaps as nbo
import mmtbx.validation.clashscore as mvc
from libtbx.utils import null_out
from libtbx.utils import Sorry
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

  show_overlap_type = *all sym macro selection
  .type = choice(multi=False)
  .help = '''When using cctbx method, this parameter allows selecting to show
  all clashes 'all', clashes dues to symmetry operation 'sym' or clashes in
  the macro molecule 'macro'.'''

  substitute_non_crystallographic_unit_cell_if_necessary = False
    .type = bool
    .help = '''\
    Will replace the crystallographic unit cell when the model
    crystallographic information is bad
    '''
"""

usage_string = """\
phenix.clashscore file.pdb [params.eff] [options ...]

Options:

  model=input_file          input PDB file
  cif=input_file            input CIF file for additional model information
  keep_hydrogens=True       keep input hydrogen files (otherwise regenerate)
  nuclear=False             use nuclear x-H distances and vdW radii
  verbose=True              verbose text output
  time_limit=120            Time limit (sec) for Reduce optimization
  show_overlap_type=all     what type of overlaps to show
  substitute_non_crystallographic_unit_cell_if_necessary=false
                            fix CRYST1 records if needed

Example:

  phenix.clashscore xxxx.pdb keep_hydrogens=True
"""

def run (args, out=sys.stdout) :
  """
  Calculates number of non-bonded atoms overlaps in a model

  Returns:
    When verbose=True the function print detailed results to log
    When verbose=False it will print:
        nb_overlaps_macro_molecule,
        nb_overlaps_due_to_sym_op,
        nb_overlaps_all
  """
  if not args:
    print >> out,usage_string
    return None
  cmdline = iotbx.phil.process_command_line_with_files(
    args=args,
    master_phil_string=master_phil_str,
    pdb_file_def="model",
    cif_file_def="cif",
    usage_string=usage_string)
  params = cmdline.work.extract()
  if (params.model is None) :
    raise Usage(usage_string)

  pdb_file_name = [x for x in args if x.endswith('.pdb')]
  cif_file_name = [x for x in args if x.endswith('.cif')]
  assert pdb_file_name
  pdb_file_name = pdb_file_name[0]
  pdb_with_h, h_were_added = mvc.check_and_add_hydrogen(
      file_name=pdb_file_name,
      model_number=0,
      nuclear=params.nuclear,
      verbose=params.verbose,
      time_limit=params.time_limit,
      keep_hydrogens=params.keep_hydrogens,
      allow_multiple_models=False,
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
    hard_minimum_nonbonded_distance=0.0,
    nonbonded_distance_threshold=None,
    substitute_non_crystallographic_unit_cell_if_necessary=
    params.substitute_non_crystallographic_unit_cell_if_necessary,
    log=null_out()    )
  # test that CRYST1 records are ok
  sps = pdb_processed_file.all_chain_proxies.special_position_settings
  if not sps:
    msg = 'None valid CRSYT1 records.\n'
    msg += 'Consider running phenix.clashscore with the option:\n'
    msg += 'substitute_non_crystallographic_unit_cell_if_necessary=true'
    raise Sorry(msg)
  grm = pdb_processed_file.geometry_restraints_manager()
  xrs = pdb_processed_file.xray_structure()
  sites_cart = xrs.sites_cart()
  site_labels = xrs.scatterers().extract_labels()
  hd_sel = xrs.hd_selection()
  macro_mol_sel = nbo.get_macro_mol_sel(pdb_processed_file)
  nb_overlaps = nbo.info(
    geometry_restraints_manager=grm,
    macro_molecule_selection=macro_mol_sel,
    sites_cart=sites_cart,
    site_labels=site_labels,
    hd_sel=hd_sel)
  if params.verbose:
    nb_overlaps.show(log=out,nbo_type=params.show_overlap_type)
  else:
    all = nb_overlaps.result.nb_overlaps_all
    macro_molecule = nb_overlaps.result.nb_overlaps_macro_molecule
    sym = nb_overlaps.result.nb_overlaps_due_to_sym_op
    out_list = map(lambda x: str(round(x,2)),[macro_molecule,sym,all])
    print >> out,', '.join(out_list)

if (__name__ == "__main__") :
  run(sys.argv[1:])

