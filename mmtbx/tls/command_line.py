from iotbx import pdb
from cctbx.array_family import flex
import sys
from mmtbx_tls_ext import tlso
import os
from cctbx import adptbx
from iotbx.option_parser import iotbx_option_parser
import iotbx.phil
from mmtbx import utils
from libtbx.utils import Sorry
import mmtbx.tls.tools


master_params = iotbx.phil.parse("""\
selection = None
  .type = str
  .multiple = True
  .help = Selection for TLS groups. Extracted from PDB file header or provided \
          below.
combine_tls = None
  .type = bool
  .help = Combine ADP contribution from TLS matrices (from PDB file header) \
          with ADP in ATOM records.
extract_tls = None
  .type = bool
  .help = Fit TLS matrices to B-factors of selected sets of atoms.
output_file_name = None
  .type = str
  .help = Prefix for output files.
enforce_positive_definite_TL = True
  .type = bool
  .help = Make sure the output T and L matrices are positive definite.
""")

banner = \
"""Author:
  Pavel Afonine

Pupose:
  - Fit TLS matrices to B-factors of selected sets of atoms. Example:
      phenix.tls model.pdb extract_tls=true selection='chain A' selection='chain B'
  - Combine B-factors computed from TLS matrices (extracted from PDB file header)
    with the B-factors from ATOM records of selected sets of atoms. Example:
      phenix.tls model.pdb selection='chain A' selection='chain B' combine_tls=true
  - The output PDB file can be of two types: contain total ADP in ATOM records
    (B_tls + B_residual) and B_residual in ATOM records.

Questions and problems:
  phenixbb@phenix-online.org"""

def get_xrs_helper(mmtbx_pdb_file, log, silent):
  mmtbx_pdb_file.set_ppf()
  xsfppf = utils.xray_structures_from_processed_pdb_file(
    processed_pdb_file = mmtbx_pdb_file.processed_pdb_file,
    scattering_table   = "wk1995",
    d_min              = None)
  xray_structure = xsfppf.xray_structure_all
  if(not silent):
    xray_structure.show_summary(f = log)
    print >> log
  return xray_structure

def run(args, command_name = "phenix.tls"):
  if(len(args) == 0): args = ["--help"]
  usage_fmt = "%s pdb_file [parameters: file or command line string]"
  des_fmt = "Example: %s model.pdb fit_tls_to.selection='%s' fit_tls_to.selection='%s'"
  command_line = (iotbx_option_parser(
    usage = usage_fmt % command_name,
    description = banner)
    .option("--show_defaults",
      action="store_true",
      help="Do not output to the screen (except errors).")
    .option("--silent",
      action="store_true",
      help="Suppress output to the screen.")
    .option("--use_elbow",
      action="store_false",
      help="Use eLBOW if there are unknown ligands in PDB file (available in PHENIX only).")
    ).process(args=args)
  #
  log = sys.stdout
  if(not command_line.options.silent):
    utils.print_header("TLS tools", out = log)
  if(command_line.options.show_defaults):
    master_params.show(out = log)
    print >> log
    return
  if(not command_line.options.silent):
    print >> log, banner
  #
  processed_args = utils.process_command_line_args(args = command_line.args,
    master_params = master_params, log = log)
  reflection_files = processed_args.reflection_files
  if(processed_args.crystal_symmetry is None):
    raise Sorry("No crystal symmetry found.")
  if(len(processed_args.pdb_file_names) == 0):
    raise Sorry("No PDB file found.")
  params = processed_args.params
  if(not command_line.options.silent):
    utils.print_header("Input parameters", out = log)
    params.show(out = log)
  params = params.extract()
  #
  if(processed_args.crystal_symmetry.unit_cell() is None or
     processed_args.crystal_symmetry.space_group() is None):
    raise Sorry("No CRYST1 record found.")
  mmtbx_pdb_file = utils.pdb_file(
    pdb_file_names   = processed_args.pdb_file_names,
    cif_objects      = processed_args.cif_objects,
    crystal_symmetry = processed_args.crystal_symmetry,
    use_elbow        = command_line.options.use_elbow,
    log              = log)
  #
  if(not command_line.options.silent):
    utils.print_header("TLS groups from PDB file header", out = log)
  pdb_inp_tls = mmtbx.tls.tools.tls_from_pdb_inp(
    remark_3_records = mmtbx_pdb_file.pdb_inp.extract_remark_iii_records(3),
    pdb_hierarchy = mmtbx_pdb_file.pdb_inp.construct_hierarchy())
  #
  tls_groups = []
  if(pdb_inp_tls.tls_present):
    if(pdb_inp_tls.error_string is not None):
      raise Sorry(pdb_inp_tls.error_string)
    mmtbx_pdb_file.set_ppf()
    xray_structure = get_xrs_helper(mmtbx_pdb_file = mmtbx_pdb_file, log = log,
      silent = command_line.options.silent)
    pdb_tls = mmtbx.tls.tools.extract_tls_from_pdb(
      pdb_inp_tls       = pdb_inp_tls,
      all_chain_proxies = mmtbx_pdb_file.processed_pdb_file.all_chain_proxies,
      xray_structure    = xray_structure)
    tls_groups = pdb_tls.pdb_inp_tls.tls_params
  #
  tls_selections_strings = []
  #
  if(len(tls_groups) == 0 and not command_line.options.silent):
    print >> log, "No TLS groups found in PDB file header."
  else:
    for i_seq, tls_group in enumerate(tls_groups):
      tls_selections_strings.append(tls_group.selection_string)
      if(not command_line.options.silent):
        print >> log, "TLS group %d: %s" % (i_seq+1, tls_group.selection_string)
        mmtbx.tls.tools.show_tls_one_group(tlso = tls_group, out = log)
        print >> log
  #
  if(len(tls_selections_strings) > 0 and len(params.selection) > 0):
    raise Sorry("Two TLS selection sources found: PDB file header and parameters.")
  if(len(params.selection) > 0):
    tls_selections_strings = params.selection
    xray_structure = get_xrs_helper(mmtbx_pdb_file = mmtbx_pdb_file, log = log,
      silent = command_line.options.silent)
  if([params.combine_tls, params.extract_tls].count(True) > 1):
    raise Sorry("Cannot simultaneously pereform: combine_tls and extract_tls")
  if([params.combine_tls, params.extract_tls].count(True) > 0):
    if(len(tls_selections_strings)==0):
      raise Sorry("No TLS selections found.")
  #
  if(len(tls_selections_strings)):
    if(not command_line.options.silent):
      utils.print_header("TLS groups selections", out = log)
    selections = utils.get_atom_selections(
      all_chain_proxies = mmtbx_pdb_file.processed_pdb_file.all_chain_proxies,
      selection_strings = tls_selections_strings,
      xray_structure    = xray_structure)
    if(not command_line.options.silent):
      print >> log, "Number of TLS groups: ", len(selections)
      print >> log, "Number of atoms: %d" % xray_structure.scatterers().size()
    n_atoms_in_tls = 0
    for sel_a in selections:
      n_atoms_in_tls += sel_a.size()
    if(not command_line.options.silent):
      print >> log, "Number of atoms in TLS groups: %d" % n_atoms_in_tls
      print >> log
    assert len(tls_selections_strings) == len(selections)
    if(not command_line.options.silent):
      for sel_a, sel_s in zip(selections,tls_selections_strings):
        print >> log, "Selection string:\n%s" % sel_s
        print >> log, "selects %d atoms." % sel_a.size()
        print >> log
      print >> log, "Ready-to-use in phenix.refine:\n"
      for sel_a, sel_s in zip(selections,tls_selections_strings):
        print >> log, sel_s
  #
  ofn = params.output_file_name
  if(ofn is None):
    ofn = os.path.splitext(os.path.basename(processed_args.pdb_file_names[0]))[0]
    if(len(processed_args.pdb_file_names) > 1):
      ofn = ofn+"_el_al"
    if(params.combine_tls):
      ofn = ofn+"_combine_tls.pdb"
    elif(params.extract_tls):
      ofn = ofn+"_extract_tls.pdb"
    else: ofn = None
  if(ofn is not None):
    ofo = open(ofn, "w")
  #
  if(params.extract_tls):
    utils.print_header(
      "Fit TLS matrices to B-factors of selected sets of atoms", out = log)
    tlsos = mmtbx.tls.tools.generate_tlsos(
      selections     = selections,
      xray_structure = xray_structure,
      value          = 0.0)
    for rt,rl,rs in [[1,0,1],[1,1,1],[0,1,1],
                     [1,0,0],[0,1,0],[0,0,1],[1,1,1],
                     [0,0,1]]*10:
      tlsos = mmtbx.tls.tools.tls_from_uanisos(
        xray_structure               = xray_structure,
        selections                   = selections,
        tlsos_initial                = tlsos,
        number_of_macro_cycles       = 10,
        max_iterations               = 100,
        refine_T                     = rt,
        refine_L                     = rl,
        refine_S                     = rs,
        enforce_positive_definite_TL = params.enforce_positive_definite_TL,
        verbose                      = -1,
        out                          = log)
      mmtbx.tls.tools.show_tls(tlsos = tlsos, out = log)
    u_cart_from_tls = mmtbx.tls.tools.u_cart_from_tls(
      sites_cart = xray_structure.sites_cart(),
      selections = selections,
      tlsos      = tlsos)
    unit_cell = xray_structure.unit_cell()
    for i_seq, sc in enumerate(xray_structure.scatterers()):
      if(u_cart_from_tls[i_seq] != (0,0,0,0,0,0)):
        u_star_tls = adptbx.u_cart_as_u_star(unit_cell,
          tuple(u_cart_from_tls[i_seq]))
        sc.u_star = tuple(flex.double(sc.u_star) - flex.double(u_star_tls))
    for sel in selections:
      xray_structure.convert_to_isotropic(selection = sel)
    mmtbx.tls.tools.remark_3_tls(tlsos = tlsos,
      selection_strings = tls_selections_strings, out = ofo)
  #
  if(params.combine_tls):
    utils.print_header("Combine B_tls with B_residual", out = log)
    mmtbx.tls.tools.combine_tls_and_u_local(xray_structure = xray_structure,
      tls_selections = selections, tls_groups = tls_groups)
    print >> log, "All done."
  #
  if(ofn is not None):
    utils.print_header("Write output PDB file %s"%ofn, out = log)
    utils.write_pdb_file(
      xray_structure = xray_structure,
      pdb_hierarchy  = mmtbx_pdb_file.processed_pdb_file.all_chain_proxies.pdb_hierarchy,
      out            = ofo)
    ofo.close()
    print >> log, "All done."
