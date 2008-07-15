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
from libtbx import smart_open
from iotbx import reflection_file_utils
from iotbx import reflection_file_reader
from libtbx.str_utils import format_value
import iotbx
from cStringIO import StringIO
from mmtbx import utils
from iotbx import pdb
from libtbx import easy_run
from cStringIO import StringIO
from mmtbx import model_statistics


def get_fmodel_object(xray_structure, f_obs, r_free_flags):
  sel = f_obs.d_spacings().data() > 0.25
  f_obs = f_obs.select(sel)
  r_free_flags = r_free_flags.select(sel)
  n_outl = sel.count(False)
  fmodel = mmtbx.f_model.manager(
    xray_structure = xray_structure,
    r_free_flags   = r_free_flags,
    target_name    = "ml",
    f_obs          = f_obs)
  sel = fmodel.outlier_selection()
  fmodel.update_xray_structure(update_f_calc = True, update_f_mask = True)
  fmodel.update_solvent_and_scale(verbose = -1)
  fmodel = fmodel.select(selection = sel)
  n_outl += sel.count(False)
  return fmodel, n_outl

def show_geometry(processed_pdb_file, scattering_table):
  xray_structure = processed_pdb_file.xray_structure()
  sctr_keys = \
    xray_structure.scattering_type_registry().type_count_dict().keys()
  has_hd = "H" in sctr_keys or "D" in sctr_keys
  hd_sel = xray_structure.hd_selection()
  geometry = processed_pdb_file.geometry_restraints_manager(
    show_energies = False,
    plain_pairs_radius = 5.0,
    edits = None,
    assume_hydrogens_all_missing = not has_hd)
  restraints_manager = mmtbx.restraints.manager(
    geometry = geometry,
    normalization = True)
  if(hd_sel.count(True) > 0 and scattering_table != "neutron"):
    xray_structure = xray_structure.select(~hd_sel)
    geometry = restraints_manager.geometry.select(selection = ~hd_sel)
    restraints_manager = mmtbx.restraints.manager(
      geometry      = geometry,
      normalization = True)
    restraints_manager.geometry.pair_proxies(sites_cart =
      xray_structure.sites_cart())
  result = model_statistics.geometry(
    sites_cart         = xray_structure.sites_cart(),
    restraints_manager = restraints_manager)
  print "  Stereochemistry statistics (mean, max, count):"
  print "    bonds            : %6.2f %6.2f %d" % (result.b_mean, result.b_max, result.b_number)
  print "    angles           : %6.2f %6.2f %d" % (result.a_mean, result.a_max, result.a_number)
  print "    dihedrals        : %6.2f %6.2f %d" % (result.d_mean, result.d_max, result.d_number)
  print "    chirality        : %6.2f %6.2f %d" % (result.c_mean, result.c_max, result.c_number)
  print "    planarity        : %6.2f %6.2f %d" % (result.p_mean, result.p_max, result.p_number)
  print "    non-bonded (min) : %6.2f" % (result.n_min)

def get_processed_pdb_file(pdb_file_name, cryst1, show_geometry_statistics):
  pdb_raw_records = smart_open.for_reading(
    file_name=pdb_file_name).read().splitlines()
  pdb_raw_records.append(cryst1)
  #############################################################################
  cif_objects = None
  if(show_geometry_statistics):
    t = time.ctime().split() # to make it safe to remove files
    time_stamp = "_"+t[4]+"_"+t[1].upper()+"_"+t[2]+"_"+t[3][:-3].replace(":","h")
    prefix = os.path.basename(pdb_file_name)
    from elbow.command_line import builder
    from elbow.command_line import join_cif_files
    del sys.argv[1:]
    molecules = builder.run(pdb_file_name,
                            do_all=True,
                            output="%s" % (prefix+time_stamp),
                            no_output=True,
                            silent=True,
                            )
    cif_file = prefix+time_stamp+".cif"
    if molecules:
      if type(molecules)!=type([]): molecules = [molecules]
      cif_lines = []
      for molecule in molecules:
        #print molecule.DisplayBrief()
        lines = molecule.WriteCIF2String()
        cif_lines.append(lines)
      rc = join_cif_files.run(cif_lines, cif_file, no_file_access=True)
      if rc:
        f=file(cif_file, "wb")
        f.write(rc)
        f.close()
    cif_objects = None
    if(os.path.isfile(cif_file)):
      cif_objects = []
      cif_objects.append((cif_file,
        mmtbx.monomer_library.server.read_cif(file_name = cif_file)))
  #############################################################################
  processed_pdb_files_srv = utils.process_pdb_file_srv(cif_objects =
    cif_objects, log = StringIO())
  processed_pdb_file, pdb_inp = \
    processed_pdb_files_srv.process_pdb_files(raw_records = pdb_raw_records)
  return processed_pdb_file

def get_xray_structures(processed_pdb_file, scattering_table, d_min):
  if 0: # XXX remove later once tested
    print list(processed_pdb_file.all_chain_proxies.pdb_inp.model_ids())
    print list(processed_pdb_file.all_chain_proxies.pdb_inp.model_indices())
  xray_structure = processed_pdb_file.xray_structure()
  if(xray_structure is None or xray_structure.scatterers().size()==0):
    raise Sorry("Cannot extract xray_structure.")
  if(scattering_table != "neutron"):
    xray_structure.scattering_type_registry(
      table = scattering_table,
      d_min = d_min,
      types_without_a_scattering_contribution=["?"])
  else:
    xray_structure.scattering_type_registry(
      types_without_a_scattering_contribution=["?"])
    xray_structure.switch_to_neutron_scattering_dictionary()
  model_indices = processed_pdb_file.all_chain_proxies.pdb_inp.model_indices()
  if(len(model_indices)>1):
     result = []
     model_indices_padded = flex.size_t([0])
     model_indices_padded.extend(model_indices)
     ranges = []
     for i, v in enumerate(model_indices_padded):
       try: ranges.append([model_indices_padded[i], model_indices_padded[i+1]])
       except IndexError: pass
     for ran in ranges:
       sel = flex.size_t(range(ran[0],ran[1]))
       result.append(xray_structure.select(sel))
     return result
  else:
    return [xray_structure]

def reflection_file_server(crystal_symmetry, reflection_files):
  return reflection_file_utils.reflection_file_server(
    crystal_symmetry=crystal_symmetry,
    force_symmetry=True,
    reflection_files=reflection_files,
    err=StringIO())

def show_model(xray_structure, serial):
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
  print "  Model #%s:"%str("%d"%serial).strip()
  result = " \n    ".join([
    "atom_number_(type:count:occ_sum) : %s (%s)"%(n_atoms,atom_counts_str),
    "ADP_(min,max,mean)               : %s %s %s"%(b_min,b_max,b_mean),
    "occupancies_(min,max,mean)       : %s %s %s"%(o_min,o_max,o_mean),
    "number_of_anisotropic            : "+format_value("%-7s",n_aniso),
    "number_of_non_positive_definite  : %s"%n_npd])
  print "   ", result

def show_data(fmodel, n_outl, test_flag_value, f_obs_labels):
  info = fmodel.info()
  flags_pc = \
   fmodel.r_free_flags.data().count(True)*100./fmodel.r_free_flags.data().size()
  print "  Data:"
  try: f_obs_labels = f_obs_labels[:f_obs_labels.index(",")]
  except ValueError: pass
  result = " \n    ".join([
    "data_label              : %s"%f_obs_labels,
    "high_ressolution        : "+format_value("%-5.2f",info.d_min),
    "low_resolution          : "+format_value("%-6.2f",info.d_max),
    "completeness_in_range   : "+format_value("%-6.2f",info.completeness_in_range),
    "wilson_b                : "+format_value("%-6.1f",fmodel.wilson_b()),
    "number_of_reflections   : "+format_value("%-8d",  info.number_of_reflections),
    "test_set_size(%)        : "+format_value("%-8.1f",flags_pc),
    "test_flag_value         : "+format_value("%-d",   test_flag_value),
    "is_twinned              : "+format_value("%-5s",  fmodel.twin_test()),
    "number_of_Fobs_outliers : "+format_value("%-8d",  n_outl),
    "anomalous_flag          : "+format_value("%-6s",  fmodel.f_obs.anomalous_flag())])
  print "   ", result

def show_model_vs_data(fmodel):
  d_max, d_min = fmodel.f_obs.d_max_min()
  flags_pc = \
   fmodel.r_free_flags.data().count(True)*100./fmodel.r_free_flags.data().size()
  k_sol = format_value("%-5.2f",fmodel.k_sol())
  b_sol = format_value("%-7.2f",fmodel.b_sol())
  b_cart = " ".join([("%8.2f"%v).strip() for v in fmodel.b_cart()])
  print "  Model_vs_Data:"
  result = " \n    ".join([
    "r_work                             : "+format_value("%-6.4f",fmodel.r_work()),
    "r_free                             : "+format_value("%-6.4f",fmodel.r_free()),
    "bulk_solvent_(k_sol,b_sol)         : %s %s"%(k_sol,b_sol),
    "overall_anisotropic_scale_(b_cart) : "+format_value("%-s",b_cart)])
  print "   ", result

def run(args, command_name = "phenix.model_vs_data",
        show_geometry_statistics = True):
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
    .option(None, "--scattering_table",
      action="store",
      default="n_gaussian",
      type="string",
      help="Choice for scattering table: n_gaussian (default) or wk1995 or it1992 or neutron.")
    ).process(args=args)
  if(command_line.options.scattering_table not in ["n_gaussian","wk1995",
     "it1992","neutron"]):
    raise Sorry("Incorrect scattering_table.")
  crystal_symmetry = None
  crystal_symmetry_data = None
  crystal_symmetry_model = None
  hkl_file_name = None
  pdb_file_name = None
  reflection_file = None
  for arg in command_line.args:
    arg_is_processed = False
    if(not os.path.isfile(arg)):
      raise Sorry("The command line argument %s is not a file."%arg)
    else:
      if(pdb.is_pdb_file(file_name=arg)):
        pdb_file_name = arg
        arg_is_processed = True
        try:
          crystal_symmetry_model = crystal_symmetry_from_pdb.extract_from(
            file_name=arg)
        except RuntimeError, e:
          if(str(e) == "No CRYST1 record."): pass
      if(not arg_is_processed):
        reflection_file = reflection_file_reader.any_reflection_file(
          file_name=arg, ensure_read_access=False)
        if(reflection_file.file_type() is not None):
          hkl_file_name = arg
          arg_is_processed = True
          crystal_symmetry_data = crystal_symmetry_from_any.extract_from(arg)
      if(not arg_is_processed):
        raise Sorry(
          "The command line argument %s is not a valid PDB or data file."%arg)
  print "Model and data files: %s %s"%(
    format_value("%5s",os.path.basename(pdb_file_name)),
    format_value("%5s",os.path.basename(hkl_file_name)))
  if([crystal_symmetry_model, crystal_symmetry_data].count(None)==0):
    if(not crystal_symmetry_model.is_similar_symmetry(crystal_symmetry_data)):
      raise Sorry("Crystal symmetry mismatch between data and PDB files.")
  if(crystal_symmetry_model is not None):
    crystal_symmetry = crystal_symmetry_model
  if(crystal_symmetry_data is not None):
    crystal_symmetry = crystal_symmetry_data
  if(crystal_symmetry is None): raise Sorry("Crystal symmetry is not defined.")
  if(pdb_file_name is None): raise Sorry("No PDB file given.")
  if(reflection_file is None): raise Sorry("No reflection file given.")
  rfs = reflection_file_server(
    crystal_symmetry = crystal_symmetry,
    reflection_files = [reflection_file])
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
  r_free_flags = determine_data_and_flags_result.r_free_flags
  test_flag_value = determine_data_and_flags_result.test_flag_value
  if(r_free_flags is None):
    r_free_flags=f_obs.array(data=flex.bool(f_obs.data().size(), False))
    test_flag_value=1
  processed_pdb_file = get_processed_pdb_file(
    pdb_file_name = pdb_file_name,
    cryst1 = pdb.format_cryst1_record(crystal_symmetry = crystal_symmetry),
    show_geometry_statistics = show_geometry_statistics)
  xray_structures = get_xray_structures(
    processed_pdb_file       = processed_pdb_file,
    scattering_table         = command_line.options.scattering_table,
    d_min                    = f_obs.d_min())
  if(not xray_structures[0].crystal_symmetry().is_similar_symmetry(
     f_obs.crystal_symmetry())):
    raise Sorry("Inconsistent crystal symmetry.")
  #
  print "  unit cell:       ", f_obs.unit_cell()
  print "  space group:     ", f_obs.crystal_symmetry().space_group_info().\
    symbol_and_number()
  if(len(xray_structures) == 1):
    fmodel, n_outl = get_fmodel_object(xray_structure = xray_structures[0],
                                       f_obs          = f_obs,
                                       r_free_flags   = r_free_flags)
    show_data(fmodel          = fmodel,
              n_outl          = n_outl,
              test_flag_value = test_flag_value,
              f_obs_labels    = f_obs.info().label_string())
    show_model(xray_structure = fmodel.xray_structure, serial = 0)
    if(show_geometry_statistics):
      show_geometry(processed_pdb_file = processed_pdb_file,
                    scattering_table   = command_line.options.scattering_table)
    show_model_vs_data(fmodel)
  else:
    f_model_data = None
    for i_seq, xray_structure in enumerate(xray_structures):
      fmodel, n_outl = get_fmodel_object(xray_structure = xray_structure,
                                         f_obs          = f_obs,
                                         r_free_flags   = r_free_flags)
      if(i_seq == 0):
        show_data(fmodel          = fmodel,
                  n_outl          = n_outl,
                  test_flag_value = test_flag_value,
                  f_obs_labels    = f_obs.info().label_string())
        f_model_data = fmodel.f_model_scaled_with_k1().data()
      else:
        f_model_data += fmodel.f_model_scaled_with_k1().data()
      show_model(xray_structure = fmodel.xray_structure, serial = i_seq)
    fmodel_average = fmodel.f_obs.array(data = f_model_data)
    fmodel_result = mmtbx.f_model.manager(
     r_free_flags = fmodel.r_free_flags,
     target_name  = "ml",
     f_obs        = fmodel.f_obs,
     f_mask       = fmodel.f_mask(),
     f_calc       = fmodel_average)
    print
    import mmtbx.bulk_solvent.bulk_solvent_and_scaling as bss
    params = bss.master_params.extract()
    params.bulk_solvent=False
    fmodel_result.update_solvent_and_scale(params = params, verbose = -1)
    if(show_geometry_statistics):
      show_geometry(processed_pdb_file = processed_pdb_file,
                    scattering_table   = command_line.options.scattering_table)
    show_model_vs_data(fmodel_result)
