# LIBTBX_SET_DISPATCHER_NAME phenix.model_vs_data

import sys, os
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

def get_xray_structures(pdb_file_name, scattering_table, cryst1, d_min):
  pdb_raw_records = smart_open.for_reading(
    file_name=pdb_file_name).read().splitlines()
  mon_lib_srv = monomer_library.server.server()
  ener_lib = monomer_library.server.ener_lib()
  pdb_raw_records.append(cryst1)
  processed_pdb_file = monomer_library.pdb_interpretation.process(
    mon_lib_srv = mon_lib_srv,
    ener_lib    = ener_lib,
    raw_records = pdb_raw_records)
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
  print "  Model #%s:"%str("%d"%serial).strip()
  result = " \n    ".join([
    "number_of_atoms                 : %s"%n_atoms,
    "ADP_(min,max,mean)              : %s %s %s"%(b_min,b_max,b_mean),
    "occupancies_(min,max,mean)      : %s %s %s"%(o_min,o_max,o_mean),
    "number_of_anisotropic           : "+format_value("%-7s",n_aniso),
    "number_of_non_positive_definite : %s"%n_npd])
  print "   ", result

def show_data(fmodel, n_outl, test_flag_value):
  d_max, d_min = fmodel.f_obs.d_max_min()
  flags_pc = \
   fmodel.r_free_flags.data().count(True)*100./fmodel.r_free_flags.data().size()
  print "  Data:"
  result = " \n    ".join([
    "high_ressolution      : "+format_value("%-5.2f",d_min),
    "low_resolution        : "+format_value("%-6.2f",d_max),
    "wilson_b              : "+format_value("%-6.1f",fmodel.wilson_b()),
    "number_of_reflections : "+format_value("%-8d",fmodel.f_obs.data().size()),
    "test_set_size(%)      : "+format_value("%-8.1f",flags_pc),
    "test_flag_value       : "+format_value("%-d",test_flag_value),
    "is_twinned            : "+format_value("%-5s",fmodel.twin_test()),
    "number_of_Fobs_outl   : "+format_value("%-8d",n_outl),
    "anomalous_flag        : "+format_value("%-6s",fmodel.f_obs.anomalous_flag())])
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

def run(args, command_name = "phenix.model_vs_data"):
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
          arg_is_processed = True
          crystal_symmetry_data = crystal_symmetry_from_any.extract_from(arg)
      if(not arg_is_processed):
        raise Sorry(
          "The command line argument %s is not a valid PDB or data file."%arg)
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
  xray_structures = get_xray_structures(
    pdb_file_name = pdb_file_name,
    scattering_table = command_line.options.scattering_table,
    d_min  = f_obs.d_min(),
    cryst1 = pdb.format_cryst1_record(crystal_symmetry = crystal_symmetry))
  if(not xray_structures[0].crystal_symmetry().is_similar_symmetry(
     f_obs.crystal_symmetry())):
    raise Sorry("Inconsistent crystal symmetry.")
  #
  print "Model file: "+format_value("%5s",os.path.basename(pdb_file_name))
  print "  unit cell:       ", f_obs.unit_cell()
  print "  space group:     ",f_obs.crystal_symmetry().space_group_info().\
    symbol_and_number()
  if(len(xray_structures) == 1):
    fmodel, n_outl = get_fmodel_object(xray_structure = xray_structures[0],
                                       f_obs          = f_obs,
                                       r_free_flags   = r_free_flags)
    show_data(fmodel = fmodel, n_outl = n_outl, test_flag_value=test_flag_value)
    show_model(xray_structure = fmodel.xray_structure, serial = 0)
    show_model_vs_data(fmodel)
  else:
    f_model_data = None
    for i_seq, xray_structure in enumerate(xray_structures):
      fmodel, n_outl = get_fmodel_object(xray_structure = xray_structure,
                                         f_obs          = f_obs,
                                         r_free_flags   = r_free_flags)
      if(i_seq == 0):
        show_data(fmodel = fmodel, n_outl = n_outl, test_flag_value=test_flag_value)
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
    show_model_vs_data(fmodel_result)

if(__name__ == "__main__"):
   run(sys.argv[1:])
