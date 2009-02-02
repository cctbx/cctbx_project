import sys, os, time
from cctbx.array_family import flex
from iotbx import pdb
from cctbx import adptbx
from iotbx.option_parser import iotbx_option_parser
from libtbx.utils import Sorry
from iotbx import reflection_file_utils
from mmtbx import monomer_library
import mmtbx.monomer_library.server
import mmtbx.monomer_library.pdb_interpretation
import mmtbx.f_model
from iotbx.pdb import crystal_symmetry_from_pdb
from iotbx import crystal_symmetry_from_any
from iotbx import reflection_file_utils
from iotbx import reflection_file_reader
from libtbx.str_utils import format_value
import iotbx
from mmtbx import utils
from iotbx import pdb
from libtbx import easy_run
from cStringIO import StringIO
from mmtbx import model_statistics
from iotbx.pdb import extract_rfactors_resolutions_sigma
import iotbx.pdb.remark_3_interpretation
import mmtbx.bulk_solvent.bulk_solvent_and_scaling as bss
from iotbx.pdb import combine_unique_pdb_files


def get_program_name(file_lines):
  result = None
  for line in file_lines:
    line = line.strip()
    result = iotbx.pdb.remark_3_interpretation.get_program(st = line)
    if(result is not None): return result
  if(result is not None):
    result = "_".join(result.split())
  return result

def show_geometry(processed_pdb_file, scattering_table, pdb_inp,
                  model_selections, show_geometry_statistics):
  xray_structures = processed_pdb_file.xray_structure()
  if(show_geometry_statistics):
    from phenix.command_line.ramalyze import ramalyze
    from phenix.command_line.rotalyze import rotalyze
    from phenix.command_line.cbetadev import cbetadev
    sctr_keys = \
      xray_structures.scattering_type_registry().type_count_dict().keys()
    has_hd = "H" in sctr_keys or "D" in sctr_keys
    geometry = processed_pdb_file.geometry_restraints_manager(
      show_energies                = False,
      plain_pairs_radius           = 5.0,
      edits                        = None,
      assume_hydrogens_all_missing = not has_hd)
    restraints_manager_all = mmtbx.restraints.manager(
      geometry      = geometry,
      normalization = True)
  hierarchy = pdb_inp.construct_hierarchy()
  models = hierarchy.models()
  print "  Number of models:", len(models)
  for i_seq, model_selection in enumerate(model_selections):
    print "  Model #%s:"%str("%d"%(i_seq+1)).strip()
    hierarchy_i_seq = pdb.hierarchy.root()
    hierarchy_i_seq.append_model(models[i_seq].detached_copy())
    #
    overall_counts_i_seq = hierarchy_i_seq.overall_counts()
    #
    print "    Number of residues in alternative conformations:", \
      overall_counts_i_seq.n_alt_conf_pure + \
      overall_counts_i_seq.n_alt_conf_proper + \
      overall_counts_i_seq.n_alt_conf_improper
    #
    rc = overall_counts_i_seq.resname_classes
    print "    Residue content:"
    len_max = 0
    for k in rc.keys():
      k=k.replace("common_","")
      if(len(k)>len_max): len_max = len(k)
    fmt = "      %-"+str(len_max)+"s : %d"
    for k,v in zip(rc.keys(), rc.values()):
      print fmt%(k.replace("common_",""), v)
    #
    xray_structure = xray_structures.select(model_selection)
    hd_sel = xray_structure.hd_selection()
    if(hd_sel.count(True) > 0 and scattering_table != "neutron"):
      show_xray_structure_statistics(xray_structure =
        xray_structure.select(~hd_sel))
    else:
      show_xray_structure_statistics(xray_structure = xray_structure)
    if(show_geometry_statistics):
      # exclude hydrogens
      if(hd_sel.count(True) > 0 and scattering_table != "neutron"):
        xray_structure = xray_structure.select(~hd_sel)
        model_selection = model_selection.select(~hd_sel)
        geometry = restraints_manager_all.geometry.select(selection = ~hd_sel)
      model_selection_as_bool = flex.bool(xray_structures.scatterers().size(),
        model_selection)
      geometry = restraints_manager_all.geometry.select(selection =
        model_selection_as_bool)
      restraints_manager = mmtbx.restraints.manager(
        geometry      = geometry,
        normalization = True)
      restraints_manager.geometry.pair_proxies(sites_cart =
        xray_structure.sites_cart())
      result = model_statistics.geometry(
        sites_cart         = xray_structure.sites_cart(),
        hd_selection       = hd_sel,
        ignore_hd          = False,
        restraints_manager = restraints_manager)
      print "    Stereochemistry statistics (mean, max, count):"
      print "      bonds            : %8.4f %8.4f %d" % (result.b_mean, result.b_max, result.b_number)
      print "      angles           : %8.4f %8.4f %d" % (result.a_mean, result.a_max, result.a_number)
      print "      dihedrals        : %8.4f %8.4f %d" % (result.d_mean, result.d_max, result.d_number)
      print "      chirality        : %8.4f %8.4f %d" % (result.c_mean, result.c_max, result.c_number)
      print "      planarity        : %8.4f %8.4f %d" % (result.p_mean, result.p_max, result.p_number)
      print "      non-bonded (min) : %8.4f" % (result.n_min)
      need_ramachandran = False
      rc = overall_counts_i_seq.resname_classes
      n_residues = 0
      for k in rc.keys():
        if(k.count('amino_acid')):
          need_ramachandran = True
          n_residues = int(rc[k])
          break
      if(need_ramachandran):
        ramalyze_obj = ramalyze()
        output, output_list = ramalyze_obj.analyze_pdb(hierarchy =
          hierarchy_i_seq, outliers_only = False)
        outl = ramalyze_obj.get_outliers_count_and_fraction()
        allo = ramalyze_obj.get_allowed_count_and_fraction()
        favo = ramalyze_obj.get_favored_count_and_fraction()
        gene = ramalyze_obj.get_general_count_and_fraction()
        glyc = ramalyze_obj.get_gly_count_and_fraction()
        prol = ramalyze_obj.get_pro_count_and_fraction()
        prep = ramalyze_obj.get_prepro_count_and_fraction()
        assert n_residues != 0
        print "      Ramachandran plot, number of:"
        print "        outliers : %-5d (%-5.2f %s)"%(outl[0],outl[1]*100.,"%")
        print "        general  : %-5d (%-5.2f %s)"%(gene[0],gene[1]*100.,"%")
        print "        allowed  : %-5d (%-5.2f %s)"%(allo[0],allo[1]*100.,"%")
        print "        favored  : %-5d (%-5.2f %s)"%(favo[0],favo[1]*100.,"%")
        print "        glycine  : %-5d (%-5.2f %s)"%(glyc[0],glyc[1]*100.,"%")
        print "        proline  : %-5d (%-5.2f %s)"%(prol[0],prol[1]*100.,"%")
        print "        prepro   : %-5d (%-5.2f %s)"%(prep[0],prep[1]*100.,"%")

def show_xray_structure_statistics(xray_structure):
  b_isos = xray_structure.extract_u_iso_or_u_equiv()
  n_aniso = xray_structure.use_u_aniso().count(True)
  n_not_positive_definite = xray_structure.is_positive_definite_u().count(False)
  b_mean = format_value("%-6.1f",adptbx.u_as_b(flex.mean(b_isos))).strip()
  b_min = format_value("%-6.1f",adptbx.u_as_b(flex.min(b_isos))).strip()
  b_max = format_value("%-6.1f",adptbx.u_as_b(flex.max(b_isos))).strip()
  n_atoms = format_value("%-8d",xray_structure.scatterers().size()).strip()
  n_npd = format_value("%-8s",n_not_positive_definite).strip()
  occ = xray_structure.scatterers().extract_occupancies()
  o_mean = format_value("%-6.2f",flex.mean(occ)).strip()
  o_min = format_value("%-6.2f",flex.min(occ)).strip()
  o_max = format_value("%-6.2f",flex.max(occ)).strip()
  atom_counts = xray_structure.scattering_types_counts_and_occupancy_sums()
  atom_counts_strs = []
  for ac in atom_counts:
    atom_counts_strs.append("%s:%s:%s"%(ac.scattering_type,str(ac.count),
      str("%10.2f"%ac.occupancy_sum).strip()))
  atom_counts_str = " ".join(atom_counts_strs)
  print "    Atoms:"
  print "      atom_number_(type:count:occ_sum) : %s (%s)"%(n_atoms,atom_counts_str)
  print "      ADP_(min,max,mean)               : %s %s %s"%(b_min,b_max,b_mean)
  print "      occupancies_(min,max,mean)       : %s %s %s"%(o_min,o_max,o_mean)
  print "      number_of_anisotropic            : "+format_value("%-7s",n_aniso)
  print "      number_of_non_positive_definite  : %s"%n_npd

def reflection_file_server(crystal_symmetry, reflection_files):
  return reflection_file_utils.reflection_file_server(
    crystal_symmetry=crystal_symmetry,
    force_symmetry=True,
    reflection_files=reflection_files,
    err=StringIO())

def show_data(fmodel, n_outl, test_flag_value, f_obs_labels):
  info = fmodel.info()
  flags_pc = \
   fmodel.r_free_flags.data().count(True)*1./fmodel.r_free_flags.data().size()
  print "  Data:"
  try: f_obs_labels = f_obs_labels[:f_obs_labels.index(",")]
  except ValueError: pass
  result = " \n    ".join([
    "data_label              : %s"%f_obs_labels,
    "high_resolution         : "+format_value("%-5.2f",info.d_min),
    "low_resolution          : "+format_value("%-6.2f",info.d_max),
    "completeness_in_range   : "+format_value("%-6.2f",info.completeness_in_range),
    "completeness(d_min-inf) : "+format_value("%-6.2f",info.completeness_d_min_inf),
    "completeness(6A-inf)    : "+format_value("%-6.2f",info.completeness_6_inf),
    "wilson_b                : "+format_value("%-6.1f",fmodel.wilson_b()),
    "number_of_reflections   : "+format_value("%-8d",  info.number_of_reflections),
    "test_set_size           : "+format_value("%-8.4f",flags_pc),
    "test_flag_value         : "+format_value("%-d",   test_flag_value),
    "number_of_Fobs_outliers : "+format_value("%-8d",  n_outl),
    "anomalous_flag          : "+format_value("%-6s",  fmodel.f_obs.anomalous_flag())])
  print "   ", result

def show_model_vs_data(fmodel):
  d_max, d_min = fmodel.f_obs.d_max_min()
  flags_pc = fmodel.r_free_flags.data().count(True)*100./\
    fmodel.r_free_flags.data().size()
  if(flags_pc == 0): r_free = None
  else: r_free = fmodel.r_free()
  k_sol = format_value("%-5.2f",fmodel.k_sol())
  b_sol = format_value("%-7.2f",fmodel.b_sol())
  b_cart = " ".join([("%8.2f"%v).strip() for v in fmodel.b_cart()])
  print "  Model_vs_Data:"
  result = " \n    ".join([
    "r_work(re-computed)                : "+format_value("%-6.4f",fmodel.r_work()),
    "r_free(re-computed)                : "+format_value("%-6.4f",r_free),
    "bulk_solvent_(k_sol,b_sol)         : %s %s"%(k_sol,b_sol),
    "overall_anisotropic_scale_(b_cart) : "+format_value("%-s",b_cart)])
  print "   ", result

def run(args,
        command_name             = "mmtbx.model_vs_data",
        show_geometry_statistics = True,
        model_size_max_atoms     = 80000,
        data_size_max_reflections= 1000000,
        unit_cell_max_dimension  = 800.):
  if(len(args) == 0): args = ["--help"]
  command_line = (iotbx_option_parser(
    usage="%s reflection_file pdb_file [options]" % command_name,
    description='Example: %s data.mtz model.pdb'%command_name)
    .option(None, "--f_obs_label",
      action="store",
      default=None,
      type="string",
      help="Label for F-obs (or I-obs).")
    .option(None, "--r_free_flags_label",
      action="store",
      default=None,
      type="string",
      help="Label for free R flags.")
    .option(None, "--twin_law",
      action="store",
      default=None,
      type="string",
      help="Provide twin law operator if twinned data used (for example: h,-h-k,-l).")
    .option(None, "--scattering_table",
      action="store",
      default="n_gaussian",
      type="string",
      help="Choice for scattering table: n_gaussian (default) or wk1995 or it1992 or neutron.")
    .option("--ignore_giant_models_and_datasets",
      action="store_true",
      help="Ignore too big models and data files to avoid potential memory problems.")
    ).process(args=args)
  if(command_line.options.scattering_table not in ["n_gaussian","wk1995",
     "it1992","neutron"]):
    raise Sorry("Incorrect scattering_table.")
  #
  processed_args = utils.process_command_line_args(args = args,
    log = sys.stdout)
  reflection_files = processed_args.reflection_files
  if(len(reflection_files) == 0):
    raise Sorry("No reflection file found.")
  crystal_symmetry = processed_args.crystal_symmetry
  if(crystal_symmetry is None):
    raise Sorry("No crystal symmetry found.")
  if(len(processed_args.pdb_file_names) == 0):
    raise Sorry("No PDB file found.")
  pdb_file_names = processed_args.pdb_file_names
  #
  rfs = reflection_file_server(
    crystal_symmetry = crystal_symmetry,
    reflection_files = reflection_files)
  parameters = utils.data_and_flags.extract()
  if(command_line.options.f_obs_label is not None):
    parameters.labels = command_line.options.f_obs_label
  if(command_line.options.r_free_flags_label is not None):
    parameters.r_free_flags.label = command_line.options.r_free_flags_label
  determine_data_and_flags_result = utils.determine_data_and_flags(
    reflection_file_server  = rfs,
    parameters              = parameters,
    data_parameter_scope    = "refinement.input.xray_data",
    flags_parameter_scope   = "refinement.input.xray_data.r_free_flags",
    data_description        = "X-ray data",
    keep_going              = True,
    log                     = StringIO())
  f_obs = determine_data_and_flags_result.f_obs
  number_of_reflections = f_obs.indices().size()
  if(command_line.options.ignore_giant_models_and_datasets and
     number_of_reflections > data_size_max_reflections):
    raise Sorry("Too many reflections: %d"%number_of_reflections)
  #
  max_unit_cell_dimension = max(f_obs.unit_cell().parameters()[:3])
  if(command_line.options.ignore_giant_models_and_datasets and
     max_unit_cell_dimension > unit_cell_max_dimension):
    raise Sorry("Too large unit cell (max dimension): %s"%
      str(max_unit_cell_dimension))
  #
  r_free_flags = determine_data_and_flags_result.r_free_flags
  test_flag_value = determine_data_and_flags_result.test_flag_value
  if(r_free_flags is None):
    r_free_flags=f_obs.array(data=flex.bool(f_obs.data().size(), False))
    test_flag_value=None
  #
  mmtbx_pdb_file = mmtbx.utils.pdb_file(
    pdb_file_names = pdb_file_names,
    cif_objects    = processed_args.cif_objects,
    cryst1         = pdb.format_cryst1_record(crystal_symmetry = crystal_symmetry),
    use_elbow      = show_geometry_statistics,
    log            = sys.stdout)
  processed_pdb_file = mmtbx_pdb_file.processed_pdb_file
  pdb_raw_records = mmtbx_pdb_file.pdb_raw_records
  pdb_inp = mmtbx_pdb_file.pdb_inp
  #
  xsfppf = mmtbx.utils.xray_structures_from_processed_pdb_file(
    processed_pdb_file = processed_pdb_file,
    scattering_table   = command_line.options.scattering_table,
    d_min              = f_obs.d_min())
  xray_structures = xsfppf.xray_structures
  model_selections = xsfppf.model_selections
  #
  print "  Unit cell:       ", f_obs.unit_cell()
  print "  Space group:     ", f_obs.crystal_symmetry().space_group_info().\
    symbol_and_number(), "number of symmetry operators:",\
    f_obs.crystal_symmetry().space_group_info().type().group().order_z()
  print "  Unit cell volume: %-15.4f" % f_obs.unit_cell().volume()
  #
  show_geometry(
    processed_pdb_file       = processed_pdb_file,
    scattering_table         = command_line.options.scattering_table,
    pdb_inp                  = pdb_inp,
    model_selections         = model_selections,
    show_geometry_statistics = show_geometry_statistics)
  #
  bss_params = bss.master_params.extract()
  bss_params.k_sol_max = 0.8
  bss_params.k_sol_min = 0.0
  bss_params.b_sol_max = 500.0
  bss_params.b_sol_min = 0.0
  bss_params.k_sol_grid_search_max = 0.6
  bss_params.k_sol_grid_search_min = 0.0
  bss_params.b_sol_grid_search_max = 80.0
  bss_params.b_sol_grid_search_min = 20.0
  bss_params.k_sol_step = 0.3
  bss_params.b_sol_step = 30.0
  fmodel = utils.fmodel_simple(xray_structures = xray_structures,
                               f_obs           = f_obs,
                               r_free_flags    = r_free_flags,
                               twin_law        = command_line.options.twin_law,
                               bss_params      = bss_params)
  n_outl = f_obs.data().size() - fmodel.f_obs.data().size()
  show_data(fmodel          = fmodel,
            n_outl          = n_outl,
            test_flag_value = test_flag_value,
            f_obs_labels    = f_obs.info().label_string())
  show_model_vs_data(fmodel)
  # Extract information from PDB file header and output (if any)
  published_results = extract_rfactors_resolutions_sigma.extract(
    file_name = pdb_file_names[0])
  if(published_results is not None):
    published_results_r_work = published_results.r_work
    published_results_r_free = published_results.r_free
    published_results_high   = published_results.high
    published_results_low    = published_results.low
    published_results_sigma  = published_results.sigma
    program_name             = get_program_name(file_lines = pdb_raw_records)
    published_results_result =  [published_results_r_work,
                          published_results_r_free,
                          published_results_high  ,
                          published_results_low   ,
                          published_results_sigma ,
                          program_name            ,
                          pdb_inp.extract_header_year()]
    if(len(published_results_result) != published_results_result.count(None)):
      print "  Information extracted from PDB file header:"
      print "    program_name    : ", program_name
      print "    year            : ", published_results_result[6]
      print "    r_work          : ", published_results.r_work
      print "    r_free          : ", published_results.r_free
      print "    high_resolution : ", published_results.high
      print "    low_resolution  : ", published_results.low
      print "    sigma_cutoff    : ", published_results.sigma
